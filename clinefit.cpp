#include <phypp.hpp>

// Structure to define a line group to be fitted simultaneously
struct line_t {
    line_t(vec1d lam, vec1d ra) : lambda(lam), ratio(ra) {
        ratio /= ratio[0];
    }

    vec1d lambda; // wavelengths of the lines
    vec1d ratio;  // flux ratios of the lines relative to the first
};

// Local functions, defined at the end of the file
void print_help(const std::map<std::string,line_t>& db);
void print_available_lines(const std::map<std::string,line_t>& db);

int main(int argc, char* argv[]) {
    // Build the line data base (you can add your own there!)
    std::map<std::string,line_t> linedb = {
        {"o2",     line_t({0.3727}, {1.0})},
        {"o3",     line_t({0.5007, 0.4959}, {1.0, 0.3})},
        {"hbeta",  line_t({0.4861}, {1.0})},
        {"halpha", line_t({0.6563}, {1.0})},
        {"n2",     line_t({0.6584}, {1.0})}
    };

    if (argc < 2) {
        print_help(linedb);
        return 0;
    }

    double z0 = dnan;
    double dz = 0.01;
    double width = 0.0007;
    bool subtract_continuum = true;
    uint_t continuum_width = 100;
    bool fit_background = false;
    double spatial_smooth = 0.0; // radius of the smoothing kernel
    double minsnr = 3.0;
    std::string expmap;
    double minexp = 1;
    bool velocity = false;
    bool verbose = false;
    double oh_threshold = dnan;
    bool allow_absorption = false;
    std::string tline;

    // Read command line arguments
    read_args(argc-1, argv+1, arg_list(z0, dz, name(tline, "line"), width,
        subtract_continuum, continuum_width, minsnr, fit_background, expmap, minexp,
        velocity, spatial_smooth, verbose, oh_threshold, allow_absorption
    ));

    // Check validity of input
    bool bad = false;
    if (!is_finite(z0)) {
        error("please provide the fiducial redshift z0=...");
        bad = true;
    }
    if (tline.empty()) {
        error("please provide the name of an emission line to fit with line=...");
        note("available lines:");
        print_available_lines(linedb);
        bad = true;
    } else if (linedb.find(tline) == linedb.end()) {
        error("unknown line '", tline, "'");
        note("available lines:");
        print_available_lines(linedb);
        bad = true;
    }

    if (bad) return 1;

    // Get line in database
    const line_t& line = linedb.find(tline)->second;

    // Read cube
    if (verbose) note("read input cube...");
    std::string infile = argv[1];
    vec3d cflx, cerr;
    fits::input_image fimg(infile);
    fimg.reach_hdu(1);
    fimg.read(cflx);
    fimg.reach_hdu(2);
    fimg.read(cerr);

    // Build wavelength axis
    uint_t nlam = cflx.dims[0];
    double cdelt = 1, crpix = 1, crval = 1;
    if (!fimg.read_keyword("CDELT3", cdelt) ||
        !fimg.read_keyword("CRPIX3", crpix) ||
        !fimg.read_keyword("CRVAL3", crval)) {
        error("could not read WCS information for wavelength axis");
        return 1;
    }

    vec1d lam = crval + cdelt*(findgen(nlam) + (1 - crpix));

    // Read 2D astrometry
    std::string ctype1, ctype2;
    std::string cunit1, cunit2;
    double crpix1, crpix2, crval1, crval2, cdelt1, cdelt2;
    double cd11, cd12, cd21, cd22;
    if (!fimg.read_keyword("CRPIX1", crpix1) || !fimg.read_keyword("CRPIX2", crpix2) ||
        !fimg.read_keyword("CRVAL1", crval1) || !fimg.read_keyword("CRVAL2", crval2) ||
        !fimg.read_keyword("CDELT1", cdelt1) || !fimg.read_keyword("CDELT2", cdelt2) ||
        !fimg.read_keyword("CTYPE1", ctype1) || !fimg.read_keyword("CTYPE2", ctype2) ||
        !fimg.read_keyword("CUNIT1", cunit1) || !fimg.read_keyword("CUNIT2", cunit2) ||
        !fimg.read_keyword("CTYPE1", ctype1) || !fimg.read_keyword("CTYPE2", ctype2) ||
        !fimg.read_keyword("CD1_1",  cd11)   || !fimg.read_keyword("CD1_2",  cd12) ||
        !fimg.read_keyword("CD2_1",  cd21)   || !fimg.read_keyword("CD2_2",  cd22)) {
        error("could not read WCS information for spatial axes");
        return 1;
    }

    // Subtract continuum (optional)
    // NB: we do it before reducing the wavelength range, in order to use the largest possible
    // part of the spectrum to properly subtract the continuum. Avoids border effects close to
    // the edges of the wavelength domain.
    if (subtract_continuum) {
        if (verbose) note("estimate and subtract continuum...");

        for (uint_t y : range(cflx.dims[1]))
        for (uint_t x : range(cflx.dims[2])) {
            vec1d tmp = cflx(_,y,x);
            for (uint_t l : range(cflx.dims[0])) {
                uint_t l0 = max(0, int_t(l)-int_t(continuum_width/2));
                uint_t l1 = min(cflx.dims[0]-1, int_t(l)+int_t(continuum_width/2));
                cflx(l,y,x) -= median(tmp[l0-_-l1]);
            }
        }
    }

    // Build a simple 1D spectrum of the cube to identify covered regions
    vec1d cf1d = partial_median(1, reform(cflx, cflx.dims[0], cflx.dims[1]*cflx.dims[2]));

    // Select a wavelength domain centered on the line(s)
    vec1u idl = where(lam > min(line.lambda)*(1.0+z0-2*dz)
        && lam < max(line.lambda)*(1.0+z0+2*dz)
        && is_finite(cf1d));

    if (idl.empty()) {
        error("the line '", tline, "' is not covered by the provided cube at z=", z0, " +/- ", dz);
        idl = where(is_finite(cf1d));
        note("the cube covers ", min(lam[idl]), " to ", max(lam[idl]));
        note("your redshift search for this line requires ", min(line.lambda)*(1.0+z0-2*dz),
            " to ", max(line.lambda)*(1.0+z0+2*dz));
        return 1;
    }

    cflx = cflx(idl,_,_);
    cerr = cerr(idl,_,_);
    lam = lam[idl];

    // Read exposure map if provided and flag out baddly covered pixels (optional)
    if (!expmap.empty() && is_finite(minexp)) {
        if (verbose) note("taking into account exposure map...");

        vec2d expo;
        fits::input_image(expmap).read(expo);

        phypp_check(expo.dims[0] == cflx.dims[1] &&
                    expo.dims[1] == cflx.dims[2],
            "incompatible dimensions for cube and exposure map");

        vec1u idb = where(expo < minexp);
        reform(cflx(_,_,_), cflx.dims[0], cflx.dims[1]*cflx.dims[2])(_,idb) = dnan;
        reform(cerr(_,_,_), cflx.dims[0], cflx.dims[1]*cflx.dims[2])(_,idb) = dnan;
    }

    // Flag out wavelengths that are too noisy because of OH lines (optional)
    if (is_finite(oh_threshold)) {
        if (verbose) note("flag out OH lines regions...");

        vec1d ce1d = partial_median(1, reform(cerr, cflx.dims[0], cflx.dims[1]*cflx.dims[2]));
        double clear_err = median(ce1d);
        double rms_err = 1.48*mad(ce1d);

        idl = where((ce1d - clear_err) < oh_threshold*rms_err);
        if (verbose) {
            note("flagging ", round(100 - 100*idl.size()/float(ce1d.size())),
                "% of the selected wavelength range");
        }

        cflx = cflx(idl,_,_);
        cerr = cerr(idl,_,_);
        lam = lam[idl];
    }

    if (verbose) {
        note("fitting ", cflx.dims[0], " spectral elements between ",
           lam[0], " and ", lam.back(), " um...");
    }

    // Apply spatial smoothing (optional)
    if (is_finite(spatial_smooth) && spatial_smooth > 0) {
        if (verbose) note("perform spatial convolution...");

        for (uint_t l : range(cflx.dims[0])) {
            // Find good pixels remaining
            vec2d tmp = cflx(l,_,_)/cerr(l,_,_);
            vec1u idg = where(is_finite(tmp));
            if (idg.empty()) continue;

            // Smooth kernel dimension (must be an odd number)
            uint_t npix = 10*spatial_smooth;
            if (npix % 2 == 0) npix += 1;
            vec2d kernel = gaussian_profile({{npix, npix}}, spatial_smooth);
            kernel /= total(kernel);

            tmp = cflx(l,_,_);

            // Put invalid pixels to zero before convolving
            vec1u idb = where(!is_finite(tmp));
            tmp[idb] = 0;

            tmp = convolve2d(tmp, kernel);

            // Bring them back as invalid afterwards
            tmp[idb] = dnan;

            cflx(l,_,_) = tmp;

            // Decrease the uncertainty accordingly
            cerr(l,_,_) *= sqrt(total(sqr(kernel)));
        }
    }

    // Define redshift grid so as to have two redshift samples per wavelength element
    double dlam = lam[1]-lam[0];
    uint_t nz = 2*(2*dz)/(dlam/line.lambda[0]);

    vec1d zs = rgen(z0-dz, z0+dz, nz);

    // Perform a redshift search for each pixel of the 2D map
    if (verbose) note("redshift search...");

    vec2d chi2(cflx.dims[1], cflx.dims[2]); chi2[_] = dinf;
    vec2d z(cflx.dims[1], cflx.dims[2]); z[_] = dnan;
    vec2d flux(cflx.dims[1], cflx.dims[2]); flux[_] = dnan;
    vec2d flux_err(cflx.dims[1], cflx.dims[2]); flux_err[_] = dnan;
    for (uint_t y : range(cflx.dims[1]))
    for (uint_t x : range(cflx.dims[2])) {
        vec1d flx = cflx(_,y,x);
        vec1d err = cerr(_,y,x);

        for (uint_t iz : range(zs)) {
            vec1d model(lam.dims);
            for (uint_t il : range(line.lambda)) {
                model += line.ratio[il]*exp(
                    -sqr(line.lambda[il]*(1.0+zs[iz]) - lam)/(2.0*sqr(width))
                )/(sqrt(2.0*dpi)*width*1e4);
            }

            linfit_result res;
            if (fit_background) {
                res = linfit(flx, err, model, 1.0);
            } else {
                res = linfit(flx, err, model);
            }

            if (res.chi2 < chi2(y,x) && (res.params[0] > 0 || allow_absorption)) {
                chi2(y,x) = res.chi2;
                z(y,x) = zs[iz];
                flux(y,x) = res.params[0];
                flux_err(y,x) = res.errors[0];
            }
        }
    }

    // Truncate the fit to a minimum SNR (optional)
    vec2d snr = flux/flux_err;
    if (is_finite(minsnr)) {
        if (verbose) note("truncate low SNR...");

        vec1u ids = where(snr < minsnr);
        z[ids] = dnan;
    }

    // Convert redshift into velocity based on the reference redshift (optional)
    if (velocity) {
        if (verbose) note("compute velocity map...");
        z = 3e5*(z-z0);
    }

    // Write the result
    if (verbose) note("write to disk...");
    fits::output_image oimg(file::remove_extension(file::get_basename(argv[1]))+"_lfit_"+tline+".fits");
    auto write_wcs = [&]() {
        oimg.write_keyword("CRPIX1", crpix1);
        oimg.write_keyword("CRPIX2", crpix2);
        oimg.write_keyword("CRVAL1", crval1);
        oimg.write_keyword("CRVAL2", crval2);
        oimg.write_keyword("CDELT1", cdelt1);
        oimg.write_keyword("CDELT2", cdelt2);
        oimg.write_keyword("CTYPE1", ctype1);
        oimg.write_keyword("CTYPE2", ctype2);
        oimg.write_keyword("CUNIT1", cunit1);
        oimg.write_keyword("CUNIT2", cunit2);
        oimg.write_keyword("CD1_1",  cd11);
        oimg.write_keyword("CD1_2",  cd12);
        oimg.write_keyword("CD2_1",  cd21);
        oimg.write_keyword("CD2_2",  cd22);
    };

    oimg.write(vec2d(0,0)); // empty primary extension, KMOS convention

    oimg.reach_hdu(1);
    oimg.write(z);
    write_wcs();

    oimg.reach_hdu(2);
    oimg.write(flux);
    write_wcs();

    oimg.reach_hdu(3);
    oimg.write(flux_err);
    write_wcs();

    oimg.reach_hdu(4);
    oimg.write(chi2);
    write_wcs();

    oimg.reach_hdu(5);
    oimg.write(snr);
    write_wcs();
    return 0;
}

void print_help(const std::map<std::string,line_t>& db) {
    using namespace format;

    print("clinefit v1.0");
    print("usage: clinefit <kmos_cube.fits> z0=... line=... [options]");
    print("");
    print("Main parameters:");
    paragraph("'kmos_cube.fits' should be a cube created by the KMOS pipeline, with 3 "
        "extensions: the first is empty (KMOS convention), the second contains the flux, "
        "while the third contains the uncertainty. Within this cube, the program will "
        "search and fit for the line given in 'line' (see the list of available lines "
        "and their code name below) around a redshift of 'z0'. It will do the fit pixel "
        "by pixel, so you end up with a series of 2D maps saved in the "
        "'*_lfit_<line>.fits' file, where '<line>' is replaced by the code name of the "
        "line that was fitted. This file contains multiple extensions: the first is empty "
        "(KMOS convention), the second contains the redshift (or velocity) map, the third "
        "contains the line flux, the fourth is the uncertainty, the fifth is the chi2, and "
        "lastly the sixth extension contains the S/N (ratio of flux to uncertainty).");
    print("Available lines:");
    print_available_lines(db);
    print("\nAvailable options (in order of importance):");
    bullet("dz=...", "Must be a number. Defines the redshift search window around the "
        "fiducial value 'z0'. In other words, the program will fit the lines within "
        "'z0-dz' and 'z0+dz'. Default is dz=0.01. You want this value to be large enough "
        "to fully encompass the range of redshifts of your object(s) and the width of the "
        "line(s), but not too large so as to avoid fitting other nearby lines or "
        "additional noise.");
    bullet("width=...", "Must be a number. Defines the width of the line in microns. This "
        "value is kept fixed in the fitting process to reduce complexity and avoid "
        "degeneracies in low S/N situations. Default is 0.0007 um. It is up to you to "
        "find out the adequate width of your lines.");
    bullet("velocity", "Set this flag if you want the program to generate a velocity map "
        "rather than a redshift map.");
    bullet("spatial_smooth=...", "Must be a number. It defines the size in pixels of the "
        "Gaussian profile that will be used for spatial convolution, in order to increase "
        "the S/N. The width of this Gaussian should ideally be similar to the width of "
        "the Point Spread Function of KMOS, which depends on the seeing at which your "
        "data was observed. Default value is 0 pixels, which does not perform spatial "
        "smoothing. Be careful that enabling this option will increase the 'effective' "
        "PSF of your image, which will result in more blurry velocity/redshift/flux maps. "
        "Only enable this option if your data are too noisy.");
    bullet("subtract_continuum", "Set this flag to zero if you do not want the program "
        "to estimate the continuum emission of your target(s) from the spectrum. The "
        "continuum is estimated from the average flux over multiple large spectral "
        "windows; large enough to wash out any line emission (the size of this window can "
        "be controlled with the parameter 'continuum_width'). By default this step is "
        "enabled, and it is recommended to leave it on unless you know that the continuum "
        "emission of your target(s) is very weak or weakly varying and will not perturb "
        "the line identification.");
    bullet("continuum_width=...", "Must be a integer. It defines the number of spectral "
        "elements that will be averaged to compute the continuum emission. Default is 100 "
        "pixels.");
    bullet("fit_background", "Set this flag to fit a constant wavelength-independent "
        "background level to the measured spectrum. Enabling this option is recommended "
        "only if you decided to disable the continuum subtraction, since it will "
        "introduce additional uncertainty to the fit (probably not much though).");
    bullet("expmap=...", "Must be a 2D FITS file containing the exposure map of the IFU "
        "in units of exposure time or number of exposures (if exposure time per exposure "
        "is constant). This map will be used with the 'minexp' option to limit the "
        "analyzed region of the IFU. Default is not to use any exposure map, but it is "
        "recommended to use one if you have it, and flag out the pixels with poor coverage.");
    bullet("minexp=...", "Must be a number, and is only useful if 'expmap' is also "
        "provided. It defines the minimum value in 'expmap' that will be considered in "
        "the analysis. Pixels below this value will be flagged out and discarded. This "
        "can be used to discard regions with poor coverage, for which the uncertainty "
        "or the background level can be hard to estimate properly. The default is to use "
        "all valid pixels. It is recommended to use this option and flag out the regions "
        "with only a handful of exposures, where strong outliers can be found (cosmic "
        "rays, detector hot pixels, etc).");
    bullet("minsnr=...", "Must be a number. This value defines the minimum S/N ratio in "
        "line flux; pixels below this threshold will not appear in the redshift (or "
        "velocity) maps. They will be kept in the flux, error, chi2 and S/N maps though. "
        "Default is 3.");
    bullet("oh_threshold", "Must be a number. Give this parameter a value if you want to "
        "exclude from the fit the wavelength regions with strong OH residuals (as "
        "estiamted from the uncertainty spectrum). The value is such that pixels with an "
        "uncertainty 'E' are flagged out if 'E' is larger than the median uncertainty "
        "by 'oh_threshold' sigma (sigma being the standard deviation of the uncertainties). "
        "A value of 3 or more should be adequate.");
    bullet("allow_absorption", "Enable this flag if you also want to fit absorption "
        "features, i.e., lines of negative fluxes. By default only positive fluxes are "
        "considered.");
    bullet("verbose", "Set this flag to print the progress of the detection process in "
        "the terminal. Can be useful if something goes wrong, or just to understand what "
        "is going on.");
}

void print_available_lines(const std::map<std::string,line_t>& db) {
    for (auto& l : db) {
        if (l.second.lambda.size() == 1) {
            print("  - ", l.first, ", lambda=", l.second.lambda[0]);
        } else {
            print("  - ", l.first, ", lambda=", l.second.lambda, ", ratios=", l.second.ratio);
        }
    }
}