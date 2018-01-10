#include <phypp.hpp>

// Structure to define a line group to be fitted simultaneously
struct line_t {
    line_t() = default;
    line_t(std::string n, vec1d lam, vec1d ra) : name(n), lambda(lam), ratio(ra) {
        ratio /= ratio[0];
    }

    std::string name; // identifier of the line
    vec1d lambda;     // wavelengths of the lines
    vec1d ratio;      // flux ratios of the lines relative to the first
};

// Local functions (defined at the end of the file)
vec2u grow_within(vec2u map, vec2b mask);
void print_help(const std::map<std::string,line_t>& db);
void print_available_lines(const std::map<std::string,line_t>& db);

int phypp_main(int argc, char* argv[]) {
    // Build the line data base (you can add your own there!)
    // NB: the line grouping and theoretical ratios are not used here
    // the whole data base will be "flattended" later on. This format
    // is just used for consistency with clinefit and slinefit.
    std::map<std::string,line_t> linedb = {
        // name            name        wavelength        flux ratio
        {"lyalpha", line_t("lyalpha", {0.12157},        {1.0})},
        {"c4",      line_t("c4",      {0.15495},        {1.0})},
        {"c3",      line_t("c3",      {0.19087},        {1.0})},
        {"mg2",     line_t("mg2",     {0.2799},         {1.0})},
        {"o2",      line_t("o2",      {0.3727},         {1.0})},
        {"ne3",     line_t("ne3",     {0.3869},         {1.0})},
        {"o3",      line_t("o3",      {0.5007, 0.4959}, {1.0, 0.3})},
        {"hdelta",  line_t("hdelta",  {0.4103},         {1.0})},
        {"hgamma",  line_t("hgamma",  {0.4342},         {1.0})},
        {"hbeta",   line_t("hbeta",   {0.4861},         {1.0})},
        {"halpha",  line_t("halpha",  {0.6563},         {1.0})},
        {"n2",      line_t("n2",      {0.6584},         {1.0})},
        {"s2",      line_t("s2",      {0.6718, 0.6733}, {1.0, 0.75})},
        {"palpha",  line_t("palpha",  {1.875},          {1.0})},
        {"o1",      line_t("o1",      {145.525},        {1.0})},
        {"c2",      line_t("c2",      {157.741},        {1.0})},
        {"c1_370",  line_t("c1_370",  {370.42},         {1.0})},
        {"c1_609",  line_t("c1_609",  {609.14},         {1.0})},
        {"co98",    line_t("co98",    {298.12},         {1.0})},
        {"co87",    line_t("co87",    {325.23},         {1.0})},
        {"co76",    line_t("co76",    {371.65},         {1.0})},
        {"co65",    line_t("co65",    {433.57},         {1.0})},
        {"co54",    line_t("co54",    {520.23},         {1.0})},
        {"co43",    line_t("co43",    {650.25},         {1.0})},
        {"co32",    line_t("co32",    {866.96},         {1.0})},
        {"co21",    line_t("co21",    {1300.40},        {1.0})},
        {"co10",    line_t("co10",    {2600.75},        {1.0})}
    };

    if (argc < 2) {
        print_help(linedb);
        return 0;
    }

    std::string expmap;
    double minexp = 1;
    double spatial_smooth = 1.2; // radius of the smoothing kernel
    uint_t spectral_bin = 1; // number of spectral pixels to sum
    bool save_cubes = false;
    double snr_det = 5.0; // SNR threshold for detections
    double snr_source = 3.0; // lower SNR threshold for the extents of the source
    bool verbose = false;
    double error_scale = 1.0;
    uint_t lambda_pad = 5;
    vec1d zhint;
    double maxdv = 200.0;
    double maxdpos = 5;
    double maxdist = 4;
    double qflag2_snr_threshold = 5.0;
    bool disable_zsearch = false;
    uint_t minqflag = 0;
    bool single_source = false;
    std::string outdir;
    bool ascii = false;
    vec1s lines;
    std::string semethod = "stddevneg";

    read_args(argc-1, argv+1, arg_list(expmap, minexp, spatial_smooth, save_cubes,
        snr_det, snr_source, verbose, spectral_bin, zhint, maxdv, maxdpos, single_source,
        qflag2_snr_threshold, minqflag, error_scale, outdir, ascii, disable_zsearch,
        lambda_pad, lines, name(semethod, "emethod"), maxdist
    ));

    if (!outdir.empty()) {
        outdir = file::directorize(outdir);
        file::mkdir(outdir);
    }

    if (zhint.size() != 0 && zhint.size() != 2) {
        error("'zhint' must contain two values: the minimum and maximum allowed redshift");
        return 1;
    }

    if (!lines.empty()) {
        bool bad = false;
        std::map<std::string, line_t> tdb;
        for (auto& l : lines) {
            auto iter = linedb.find(l);
            if (iter == linedb.end()) {
                error("unknown line '", l, "'");
                bad = true;
            } else {
                tdb.insert(*iter);
            }
        }

        if (bad) {
            print_available_lines(linedb);
            return 1;
        }

        std::swap(linedb, tdb);
    }

    // Find out which method to use to compute the uncertainties
    enum class error_method {
        pipeline,
        mad,
        stddev,
        madneg,
        stddevneg
    } emethod = error_method::stddevneg;

    vec<1,std::pair<std::string,error_method>> semethod_dict = {
        {"pipeline",  error_method::pipeline},
        {"mad",       error_method::mad},
        {"stddev",    error_method::stddev},
        {"madneg",    error_method::madneg},
        {"stddevneg", error_method::stddevneg},
    };

    bool found_method = false;
    for (auto& d : semethod_dict) {
        if (d.first == semethod) {
            emethod = d.second;
            found_method = true;
            break;
        }
    }

    if (!found_method) {
        error("unknown method '", semethod, "' to compute uncertainties");

        auto pair_first = vectorize_lambda([](const std::pair<std::string,error_method>& p) {
            return p.first;
        });

        note("must be one of: ", collapse(pair_first(semethod_dict), ", "));

        return 1;
    }

    // Read cube
    std::string infile = argv[1];
    vec3f cube;
    fits::input_image fimg(infile);
    fimg.reach_hdu(1);
    fimg.read(cube);

    vec3f perror;
    if (emethod == error_method::pipeline) {
        fimg.reach_hdu(2);
        fimg.read(perror);
    }

    // Read wavelength/frequency axis WCS
    uint_t nlam = cube.dims[0];
    double cdelt = 1, crpix = 1, crval = 1;
    vec1s missing;
    if (!fimg.read_keyword("CDELT3", cdelt)) missing.push_back("CDELT3");
    if (!fimg.read_keyword("CRPIX3", crpix)) missing.push_back("CRPIX3");
    if (!fimg.read_keyword("CRVAL3", crval)) missing.push_back("CRVAL3");
    if (!missing.empty()) {
        error("could not read WCS information for wavelength axis");
        note("missing keyword", missing.size() > 1 ? "s " : " ", collapse(missing, ", "));
        return 1;
    }

    bool frequency = false;
    std::string ctype;
    if (fimg.read_keyword("CTYPE3", ctype)) {
        if (ctype == "FREQ") {
            // Cube in frequency units, assuming in Hz
            frequency = true;
        } else if (ctype == "WAVE") {
            // Cube in wavelength units, assuming in micron
            frequency = false;
        }
    }

    // Read 2D astrometry
    std::string ctype1, ctype2;
    std::string cunit1, cunit2;
    double crpix1, crpix2, crval1, crval2, cdelt1, cdelt2;
    double cd11, cd12, cd21, cd22;

    if (!fimg.read_keyword("CRPIX1", crpix1)) missing.push_back("CRPIX1");
    if (!fimg.read_keyword("CRPIX2", crpix2)) missing.push_back("CRPIX2");
    if (!fimg.read_keyword("CRVAL1", crval1)) missing.push_back("CRVAL1");
    if (!fimg.read_keyword("CRVAL2", crval2)) missing.push_back("CRVAL2");
    if (!fimg.read_keyword("CDELT1", cdelt1)) missing.push_back("CDELT1");
    if (!fimg.read_keyword("CDELT2", cdelt2)) missing.push_back("CDELT2");
    if (!fimg.read_keyword("CTYPE1", ctype1)) missing.push_back("CTYPE1");
    if (!fimg.read_keyword("CTYPE2", ctype2)) missing.push_back("CTYPE2");
    if (!fimg.read_keyword("CUNIT1", cunit1)) missing.push_back("CUNIT1");
    if (!fimg.read_keyword("CUNIT2", cunit2)) missing.push_back("CUNIT2");
    if (!fimg.read_keyword("CD1_1",  cd11))   missing.push_back("CD1_1");
    if (!fimg.read_keyword("CD1_2",  cd12))   missing.push_back("CD1_2");
    if (!fimg.read_keyword("CD2_1",  cd21))   missing.push_back("CD2_1");
    if (!fimg.read_keyword("CD2_2",  cd22))   missing.push_back("CD2_2");
    if (!missing.empty()) {
        error("could not read WCS information for spatial axes");
        note("missing keyword", missing.size() > 1 ? "s " : " ", collapse(missing, ","));
        return 1;
    }

    // Flag out the pixels at the border of the spectrum
    if (lambda_pad != 0) {
        if (2*lambda_pad >= cube.dims[0]-1) {
            error("the input cube is smaller than the requested lambda pad ("
                "the cube contains ", cube.dims[0], " elements)");
            return 1;
        }

        // Build a simple 1D spectrum to identify valid regions
        vec1d tmp_spec = partial_median(1, partial_median(2, cube));

        // Locate first and last valid spectral elements
        uint_t lambda_min = 0, lambda_max = cube.dims[0]-1;
        for (uint_t i : range(tmp_spec)) {
            if (is_finite(tmp_spec[i])) {
                lambda_min = i;
                break;
            }
        }
        for (uint_t i : range(tmp_spec)) {
            if (is_finite(tmp_spec[cube.dims[0]-1-i])) {
                lambda_max = cube.dims[0]-1-i;
                break;
            }
        }

        // Apply padding
        lambda_min = (lambda_min+lambda_pad > cube.dims[0]-1 ?
            cube.dims[0]-1 : lambda_min + lambda_pad);
        lambda_max = (lambda_max < lambda_pad ?
            0 : lambda_max - lambda_pad);

        if (lambda_max <= lambda_min) {
            error("the cube does not contain any valid spectral slice");
            return 1;
        }

        // Flag bad spectral elements
        cube(_-lambda_min,_,_) = fnan;
        cube(lambda_max-_,_,_) = fnan;
    }

    // Bin spectral pixels and start building exposure map
    vec3f exposure;
    if (spectral_bin > 1) {
        if (verbose) note("rebinning spectral axis...");

        // Rebin SNR array
        nlam = floor(nlam/double(spectral_bin));

        vec3f ocube = cube;
        cube.resize(nlam, cube.dims[1], cube.dims[2]);
        cube[_] = 0;

        vec3f operror;
        if (emethod == error_method::pipeline) {
            operror = perror;
            perror.resize(cube.dims);
            perror[_] = 0;
        }

        // Keep track of how many pixels were used in the sum
        exposure.resize(cube.dims);

        // Compute the sum of the pixels in each bin
        for (uint_t l : range(nlam)) {
            for (uint_t i : range(spectral_bin)) {
                vec2f tmp = ocube(l*spectral_bin + i,_,_);
                vec1u idg = where(is_finite(tmp));
                cube(l,_,_)[idg] += tmp[idg];
                exposure(l,_,_)[idg] += 1;

                if (emethod == error_method::pipeline) {
                    tmp = operror(l*spectral_bin + i,_,_);
                    perror(l,_,_)[idg] += sqr(tmp[idg]);
                }
            }
        }

        cube /= exposure;
        if (emethod == error_method::pipeline) {
            perror = sqrt(perror)/exposure;
        }

        exposure = 1/sqrt(exposure);

        // Update wavelength array
        crpix = (crpix - 0.5)/spectral_bin + 0.5;
        cdelt *= spectral_bin;
    } else {
        exposure = replicate(1.0f, cube.dims);
    }

    // Build wavelength axis
    vec1d  lam = crval + cdelt*(findgen(nlam) + (1 - crpix));
    double dspec = abs(cdelt);
    double dlam = dspec;
    if (frequency) {
        double lam0 = 3e14/mean(lam);
        lam = 3e14/lam;
        dlam = dlam*sqr(lam0)/3e14;
    }

    // Read exposure map if provided
    vec1u badexp;
    if (!expmap.empty()) {
        if (verbose) note("taking into account exposure map...");

        vec2d tmp;
        fits::input_image(expmap).read(tmp);

        phypp_check(tmp.dims[0] == cube.dims[1] &&
                    tmp.dims[1] == cube.dims[2],
            "incompatible dimensions for cube and exposure map");

        if (is_finite(minexp)) {
            badexp = where(tmp < minexp);
            tmp[badexp] = 0;
        }

        for (uint_t l : range(nlam)) {
            exposure(l,_,_) *= sqrt(tmp/max(tmp));
        }
    }

    // For each spectral slice...
    vec3f err;
    if (emethod == error_method::pipeline) {
        err = std::move(perror);
    } else {
        err = replicate(fnan, cube.dims);
    }

    if (verbose) note("filtering wavelength slices...");

    vec3f cube_ns, err_ns;
    if (spatial_smooth > 0) {
        cube_ns.resize(cube.dims);
        err_ns.resize(cube.dims);
    }

    for (uint_t l : range(nlam)) {
        // 1) Flag baddly covered areas with exposure map (optional)
        if (!badexp.empty()) {
            cube(l,_,_)[badexp] = dnan;
        }

        // Find good pixels remaining
        vec2d tmp = cube(l,_,_)/exposure(l,_,_);
        vec1u idg = where(is_finite(tmp));
        if (!idg.empty()) {
            double img_rms = 0.0;

            switch (emethod) {
                // Use errors computed by the pipeline, nothing to do
                case error_method::pipeline : break;

                // 2) Compute pixel fluctuations of exposure-renormalized fluxes
                // and estimate error cube by this value and the exposure

                // Fluctuations from standard deviation
                // -> will be more sensitive toward outliers and the actual source in
                // the map, will tend to overestimate the flux fluctuations
                case error_method::stddev :
                    img_rms = stddev(tmp[idg]);
                    break;

                // Fluctuations from median absolute deviation
                // -> not enough pixels in the IFU to use this one reliably, it tends
                // to underestimate the actual flux fluctuation in some situation
                case error_method::mad :
                    img_rms = 1.48*mad(tmp[idg]);
                    break;

                // Fluctuations from median of negative pixels in median subtracted map
                case error_method::madneg :
                    tmp -= median(tmp[idg]);
                    idg = where(is_finite(tmp) && tmp < 0);
                    img_rms = -1.48*median(tmp[idg]);
                    break;

                // Fluctuations from RMS of negative pixels in median subtracted map
                // -> more sensitive to outliers, but only from noise, not the sources
                case error_method::stddevneg :
                    tmp -= median(tmp[idg]);
                    idg = where(is_finite(tmp) && tmp < 0);
                    img_rms = rms(tmp[idg]);
                    break;
            }

            if (emethod != error_method::pipeline) {
                err(l,_,_) = (img_rms*error_scale)*exposure(l,_,_);
            } else {
                err(l,_,_) *= error_scale;
            }

            if (spatial_smooth > 0) {
                cube_ns(l,_,_) = cube(l,_,_);
                err_ns(l,_,_) = err(l,_,_);

                // 3) Apply spatial smoothing
                // Smooth kernel dimension (must be an odd number)
                uint_t npix = 10*spatial_smooth;
                if (npix % 2 == 0) npix += 1;
                vec2d kernel = gaussian_profile({{npix, npix}}, spatial_smooth);
                kernel /= total(kernel);

                tmp = cube(l,_,_);

                // Put invalid pixels to zero before convolving
                vec1u idb = where(!is_finite(tmp));
                tmp[idb] = 0;

                tmp = convolve2d(tmp, kernel);

                // Bring them back as invalid afterwards
                tmp[idb] = dnan;

                cube(l,_,_) = tmp;

                // 4) Decrease the uncertainty accordingly
                tmp = err(l,_,_);
                idb = where(!is_finite(tmp));
                tmp[idb] = 0;
                tmp = sqrt(convolve2d(sqr(tmp), sqr(kernel)));
                tmp[idb] = dnan;
                err(l,_,_) = tmp;
            }
        } else {
            cube(l,_,_) = fnan;
            err(l,_,_) = fnan;

            if (spatial_smooth > 0) {
                cube_ns(l,_,_) = fnan;
                err_ns(l,_,_) = fnan;
            }
        }
    }

    // Free some unused cubes
    { vec3f eraser = std::move(exposure); }

    // If asked, save the processed cubes to disk
    std::string ofilebase = outdir+file::remove_extension(file::get_basename(infile));
    if (save_cubes) {
        fits::output_image oimg(ofilebase+"_filt.fits");

        auto write_wcs = [&]() {
            oimg.write_header(fimg.read_header());
            if (spectral_bin > 1) {
                oimg.write_keyword("CDELT3", cdelt);
                oimg.write_keyword("CRPIX3", crpix);
                oimg.write_keyword("CD3_3", cdelt);
            }
        };

        // Empty primary array (KMOS convention)
        oimg.write(vec2d(0,0));

        // Flux
        oimg.reach_hdu(1);
        oimg.write(cube);
        write_wcs();

        // Uncertainty
        oimg.reach_hdu(2);
        oimg.write(err);
        write_wcs();

        // S/N
        oimg.reach_hdu(3);
        oimg.write(cube/err);
        write_wcs();
    }

    // Now we are ready to try detecting stuff above the noise
    // Build a detection list and segmentation map at each wavelength slice
    // independently

    vec3u seg(cube.dims);
    uint_t nsrc_tot = 0;

    vec2d cx = generate_img({{cube.dims[1], cube.dims[2]}}, [](double y, double x) { return x; });
    vec2d cy = generate_img({{cube.dims[1], cube.dims[2]}}, [](double y, double x) { return y; });

    vec1u id;
    vec1d lambda;
    vec1u lpix;
    vec1d x, y;
    vec1u npix;
    vec1d flux;
    vec1d flux_err;

    vec2u cseg(cube.dims[1], cube.dims[2]);
    double hsize = cube.dims[1]/2;
    vec2b is_central;
    if (is_finite(maxdist)) {
        is_central = sqr(cx - hsize) + sqr(cy - hsize) < sqr(maxdist);
    } else {
        is_central = replicate(true, cseg.dims);
    }

    if (verbose) note("finding and segmenting detections...");
    for (uint_t l : range(nlam)) {
        // Flag pixels above the S/N threshold
        vec2d snr = cube(l,_,_)/err(l,_,_);
        vec2u det = vec2u(snr > snr_det && is_central);

        cseg += det;

        // Segment them into individual sources
        segment_output so;
        segment_params sp;
        sp.first_id = nsrc_tot + 1;
        det = segment(det, so, sp);
        nsrc_tot += so.id.size();

        if (!so.id.empty()) {
            // We have a new source(s)!
            // Make them grow a little toward lower S/N since they must be real
            if (snr_source < snr_det) {
                det = grow_within(det, snr > snr_source);
            }

            // Save that into the segmentation map
            seg(l,_,_) = det;

            // For each source, save its identifiers and some information
            for (uint_t i : range(so.id)) {
                vec1u idd = where(det == so.id[i]);

                id.push_back(so.id[i]);
                npix.push_back(idd.size());

                lpix.push_back(l+1);
                lambda.push_back(lam[l]);

                // If we applied spatial smoothing, the errors are correlated
                // and the standard error propagation will under estimate the
                // true uncertainty. So we rather measure the flux on the
                // non-smoothed image
                if (spatial_smooth > 0) {
                    flux.push_back(total(cube_ns(l,_,_)[idd])*dspec);
                    flux_err.push_back(sqrt(total(sqr(err_ns(l,_,_)[idd])))*dspec);
                } else {
                    flux.push_back(total(cube(l,_,_)[idd])*dspec);
                    flux_err.push_back(sqrt(total(sqr(err(l,_,_)[idd])))*dspec);
                }

                x.push_back(total((snr*cx)[idd])/total(snr[idd]));
                y.push_back(total((snr*cy)[idd])/total(snr[idd]));
            }
        }
    }

    if (verbose) note("found ", nsrc_tot, " source", (nsrc_tot > 1 ? "s" : ""), " in the cube");

    fits::output_image fseg(ofilebase+"_seg.fits");

    fseg.write(vec2d(0,0)); // empty primary

    // Write segmentation cube
    fseg.reach_hdu(1);
    fseg.write(seg);
    fseg.write_header(fimg.read_header());
    if (spectral_bin > 1) {
        fseg.write_keyword("CDELT3", cdelt);
        fseg.write_keyword("CRPIX3", crpix);
        fseg.write_keyword("CD3_3", cdelt);
    }

    // Write collapsed segmentation
    fseg.reach_hdu(2);
    fseg.write(cseg);
    fseg.write_keyword("CRPIX1", crpix1);
    fseg.write_keyword("CRPIX2", crpix2);
    fseg.write_keyword("CRVAL1", crval1);
    fseg.write_keyword("CRVAL2", crval2);
    fseg.write_keyword("CDELT1", cdelt1);
    fseg.write_keyword("CDELT2", cdelt2);
    fseg.write_keyword("CTYPE1", ctype1);
    fseg.write_keyword("CTYPE2", ctype2);
    fseg.write_keyword("CUNIT1", cunit1);
    fseg.write_keyword("CUNIT2", cunit2);
    fseg.write_keyword("CD1_1",  cd11);
    fseg.write_keyword("CD1_2",  cd12);
    fseg.write_keyword("CD2_1",  cd21);
    fseg.write_keyword("CD2_2",  cd22);
    fseg.close();

    vec1d ra, dec;
    astro::xy2ad(astro::wcs(fimg.read_header()), x+1, y+1, ra, dec);

    std::string lambda_name = (frequency ? "frequency" : "lambda");
    std::string lpix_name = (frequency ? "fpix" : "lpix");
    std::string lambda_unit = (frequency ? "[Hz]" : "[um]");
    std::string flux_unit = (frequency ? "flux [Jy Hz]" : "flux [erg/s/cm2]");

    // Adjust flux unit
    if (frequency) {
        // Northing to do
        // Assuming cube is given in Jy
    } else {
        // Convert dspec from micron to Angstroms
        // Assuming cube is given in erg/s/cm2/A
        flux     *= 1e4;
        flux_err *= 1e4;
    }

    if (ascii) {
        vec1s hdr = {"ID", "x", "y", "RA [deg]", "Dec [deg]", "Npix",
            lambda_name+" "+lambda_unit, lambda_name+" [pix]", flux_unit, "error", "SNR"};

        vec1s slambda = (frequency ? strna_sci(3e14/lambda) : strna(lambda));

        ascii::write_table_hdr(ofilebase+"_cat.cat", 18, hdr,
            id, x, y, ra, dec, npix, slambda, lpix, strna_sci(flux), strna_sci(flux_err), flux/flux_err
        );
    } else {
        fits::write_table(ofilebase+"_cat.fits", ftable(
            id, x, y, ra, dec, npix, name(lambda, lambda_name), name(lpix, lpix_name),
            flux, flux_err
        ));
    }

    // If no detection, there nothing more we can do
    if (x.empty() || disable_zsearch) return 0;

    // Combine multiple wavelengths into a single source, based on distance of emission
    // centroid
    vec1d xg = x, yg = y, group = id, ngroup = replicate(1, x.size());
    std::vector<vec1u> gsids;
    for (uint_t i : range(x)) {
        gsids.push_back(vec1u{i});
    }

    while (true) {
        // First compute the nearest neighbors between groups, within 'maxdpos'
        vec1d dg = replicate(dinf, xg.size());
        for (uint_t i : range(xg))
        for (uint_t j : range(i+1, xg.size())) {
            double d = sqr(xg[j]-xg[i]) + sqr(yg[j]-yg[i]);
            if (single_source || (d < sqr(maxdpos) && d < dg[i] && d < dg[j])) {
                group[j] = group[i];
                dg[j] = dg[i] = d;
            }
        }

        vec1u uid = unique_ids(group);

        // If no group was merged, then we have converged
        if (uid.size() == group.size()) break;

        // Some groups were merged, compute their barycenter, group size and member ids
        uint_t nn = uid.size();
        vec1d nxg(nn), nyg(nn), ng = uindgen(nn), nng(nn);
        std::vector<vec1u> ngi(nn);
        for (uint_t iu : range(nn)) {
            vec1u idd = where(group == group[uid[iu]]);
            nng[iu] = total(ngroup[idd]);
            nxg[iu] = total(ngroup[idd]*xg[idd])/nng[iu];
            nyg[iu] = total(ngroup[idd]*yg[idd])/nng[iu];
            for (uint_t i : range(idd)) {
                append<0>(ngi[iu], gsids[idd[i]]);
            }
        }

        // Swap with previous group variables and try again
        std::swap(nxg, xg); std::swap(nyg, yg);
        std::swap(ng, group); std::swap(nng, ngroup);
        std::swap(ngi, gsids);
    }

    if (verbose) {
        note("merged ", x.size(), " individual spectral detections into ", xg.size(),
            " independent spatial sources");
    }

    // Identify possible lines & redshift
    vec1u gqf1(xg.size());
    vec1u gqf2(xg.size());
    vec1d gfm(xg.size());
    vec1d gdv(xg.size());
    vec1d gz = replicate(dnan, xg.size());
    vec1d gze = replicate(dnan, xg.size());
    vec1s gl(xg.size());

    // Flatten the line database
    vec1d ldb;
    vec1s ndb;
    for (auto& l : linedb)
    for (uint_t ill : range(l.second.lambda)) {
        ldb.push_back(l.second.lambda[ill]);

        if (l.second.lambda.size() == 1) {
            ndb.push_back(l.first);
        } else {
            ndb.push_back(l.first+"-"+strn(ill+1));
        }
    }

    // Sort it by wavelength
    {
        vec1u ids = sort(ldb);
        if (frequency) ids = reverse(ids);
        ldb = ldb[ids]; ndb = ndb[ids];
    }

    // Now explore each group
    for (uint_t ig : range(xg.size())) {
        vec1d lg(gsids[ig].size());
        vec1d snrg(gsids[ig].size());

        for (uint_t j : range(gsids[ig])) {
            lg[j] = lambda[gsids[ig][j]];
            snrg[j] = (flux/flux_err)[gsids[ig][j]];
        }

        vec1d zs;
        vec1u zfrom;
        vec1u zline;
        vec1d fmatch;
        vec1d dvmax;
        vec1d dzs;
        vec1u qflag1;
        vec1u qflag2;
        vec1s lids;

        if (verbose) {
            note("group ", strn(ig+1), ": ", lg.size(), " detection", (lg.size() > 1 ? "s" : ""));
        }

        // For each line in the group, first blindly try the line alone
        // against the line database to list all the possible redshifts
        std::vector<vec1u> found_ill;
        std::vector<vec1u> found_izl;
        vec1u              found_il;
        for (uint_t il : range(lg)) {
            found_ill.push_back(vec1u{});
            found_izl.push_back(vec1u{});

            for (uint_t ill : range(ldb)) {
                double ztry = lg[il]/ldb[ill] - 1.0;

                // Discard negative redshits...
                if (ztry < 0) continue;

                // If the user provided a redshift hint, we can only list
                // the possibilities within the available range
                if (!zhint.empty() && (ztry < zhint[0] || ztry > zhint[1])) continue;

                zs.push_back(ztry);
                zfrom.push_back(il);
                zline.push_back(ill);
                fmatch.push_back(1.0/lg.size());
                dvmax.push_back(0.0);
                dzs.push_back(0.5*dlam/lg[il]);
                lids.push_back(ndb[ill]);

                uint_t tqflag2 = snrg[il] > qflag2_snr_threshold ? 1 : 0;
                qflag1.push_back(1);
                qflag2.push_back(tqflag2);

                found_ill.back().push_back(ill);
                found_izl.back().push_back(zs.size()-1);
            }

            if (found_ill.back().empty()) {
                found_ill.pop_back();
                found_izl.pop_back();
            } else {
                found_il.push_back(il);
            }
        }

        // Now, if we have more than one line, try to see if we can
        // strengthen the redshifts by including additional detected lines
        // at the same redshift, within the uncertainties
        if (lg.size() > 1) {
            vec1d ozs = zs, odzs = dzs;
            for (uint_t iz : range(zs)) {
                // Pick only the other lines that correspond to the same redshift
                vec1d dz = ozs - zs[iz];
                vec1u idz = where(zfrom >= zfrom[iz] &&
                    abs(dz) <= sqrt(sqr(maxdv/3e5) + sqr(odzs) + sqr(dzs[iz]))
                );

                // Make sure we only combine lines that come from unique
                // spectral features. Indeed, one spectral feature may be
                // associated to several lines, within the allowed uncertainty,
                // but we want it to only count once regardless.
                vec1u idus = idz[unique_ids(zfrom[idz])];

                if (idus.size() > 2) {
                    vec1d tz = zs[idus];
                    vec1d tze = dzs[idus];
                    vec1d tw = snrg[zfrom[idus]];

                    // Compute average redshift of the group, only counting
                    // each spectral feature once
                    double totw = total(tw);
                    zs.push_back(total(tz*tw)/totw);
                    dvmax.push_back((max(dz[idz]) - min(dz[idz]))*3e5);
                    fmatch.push_back(idus.size()/double(lg.size()));
                    dzs.push_back(sqrt(total(sqr(tze*tw)))/totw);
                    qflag2.push_back(total(qflag2[idus]));

                    // Build line list
                    vec1u ids = idus[sort(snrg[zfrom[idus]])];
                    vec1s tmpl(idus.size());
                    for (uint_t i : range(ids)) {
                        vec1u idl = idz[where(zfrom[idz] == zfrom[ids[i]])];
                        tmpl[i] = strn(zfrom[ids[i]]+1)+":"+collapse(lids[idl], "/");
                    }

                    lids.push_back(collapse(tmpl,","));

                    // Now compute the number of independent lines for the quality flag.
                    // To do so, we only count the lines that are unambiguously
                    // observed, i.e., not counting those that are systematically
                    // blended with the same lines. Note that the same line can be
                    // observed in multiple spectral elements if the 'maxdv' is large
                    // enough, so we have to be careful.

                    // First count the number of unique lines
                    vec1u idul = ids[unique_ids(zline[ids])];

                    vec1u blend_group(idul.size());
                    uint_t groupid = 0;
                    uint_t nline = 0;
                    for (uint_t iu1 : range(idul)) {
                        // For each line not already grouped
                        if (blend_group[iu1] != 0) continue;
                        // ... pick up all its observed spectral elements
                        vec1u zfrom1 = zfrom[idz[where(zline[idz] == zline[idul[iu1]])]];

                        // Give it a group
                        ++groupid;
                        blend_group[iu1] = groupid;
                        // fully_blended_group.push_back(false);
                        bool fully_blended_group = false;

                        for (uint_t iu2 : range(idul)) {
                            // For each other line
                            if (iu2 == iu1) continue;
                            // ... pick up all its observed spectral elements
                            vec1u zfrom2 = zfrom[idz[where(zline[idz] == zline[idul[iu2]])]];

                            // See if the two lines overlap
                            vec1u id1, id2;
                            match(zfrom1, zfrom2, id1, id2);

                            // If no overlap, go to the next line
                            if (id1.empty()) continue;

                            // If the first line is fully included inside the second
                            if (id1.size() == zfrom1.size()) {
                                if (id2.size() == zfrom2.size()) {
                                    // If the second is also included in the first, these
                                    // two lines are fully blended so group them together
                                    blend_group[iu2] = blend_group[iu1];
                                } else {
                                    // Else, forget about this line and its group
                                    fully_blended_group = true;
                                }
                            }
                        }

                        if (!fully_blended_group) ++nline;
                    }

                    qflag1.push_back(nline);
                }
            }
        }

        if (!zs.empty()) {
            lids[_-(zfrom.size()-1)] = strna(zfrom+1)+":"+lids[_-(zfrom.size()-1)];
        }

        // Sort all the solutions by decreasing quality flag and increasing dvmax
        vec1u ids = uindgen(zs.size());
        inplace_sort(ids, [&](uint_t i1, uint_t i2) {
            if (qflag1[i1] > qflag1[i2]) return true;
            if (qflag1[i1] < qflag1[i2]) return false;
            if (qflag2[i1] > qflag2[i2]) return true;
            if (qflag2[i1] < qflag2[i2]) return false;
            if (fmatch[i1] > fmatch[i2]) return true;
            if (fmatch[i1] < fmatch[i2]) return false;
            if (dvmax[i1] < dvmax[i2]) return true;
            if (dvmax[i1] > dvmax[i2]) return false;
            return zs[i1] < zs[i2];
        });

        // Apply qflag cut
        ids = ids[where(qflag1[ids] >= minqflag)];

        zs = zs[ids]; dzs = dzs[ids]; fmatch = fmatch[ids]; dvmax = dvmax[ids];
        qflag1 = qflag1[ids]; qflag2 = qflag2[ids]; lids = lids[ids];

        if (zs.empty()) {
            if (verbose) {
                note("group ", ig+1, ": no reliable redshift");
            }
        } else {
            gqf1[ig] = qflag1[0];
            gqf2[ig] = qflag2[0];
            gfm[ig] = fmatch[0];
            gdv[ig] = dvmax[0];
            gz[ig] = zs[0];
            gze[ig] = dzs[0];
            gl[ig] = lids[0];

            if (verbose) {
                note("group ", ig+1, ": most likely redshift is ", zs[0],
                    " +/- ", dzs[0], " (nline: ", qflag1[0], ", ngline: ", qflag2[0], ", "
                    "match: ", round(100*fmatch[0]), "%, maxdv: ", round(dvmax[0]), " km/s)");
                note("from lines: ", lids[0]);
            }
        }

        // Save group
        if (!zs.empty()) {
            if (ascii) {
                vec1s hdr = {"qflag1", "qflag2", "fmatch", "maxdv", "z", "zerr", "lines"};

                // Manual handling of column widths
                vec1s sq1 = strna(qflag1);
                vec1s sq2 = strna(qflag2);
                vec1s sfm = strna(round(100.0*fmatch)/100.0);
                vec1s sdv = strna(round(dvmax));
                vec1s sz = strna(zs);
                vec1s sdz = strna(dzs);
                vec1s sl = lids;

                sq1 = align_right(sq1, max(length(sq1))+8);
                hdr[0] = align_right(hdr[0], sq1[0].size());
                sq2 = align_right(sq2, max(length(sq2))+8);
                hdr[1] = align_right(hdr[1], sq2[0].size());
                sfm = align_right(sfm, max(length(sfm))+8);
                hdr[2] = align_right(hdr[2], sfm[0].size());
                sdv = align_right(sdv, max(length(sdv))+8);
                hdr[3] = align_right(hdr[3], sdv[0].size());
                sz = align_right(sz, max(length(sz))+4);
                hdr[4] = align_right(hdr[4], sz[0].size());
                sdz = align_right(sdz, max(length(sdz))+4);
                hdr[5] = align_right(hdr[5], sdz[0].size());
                sl  = "  "+align_left(sl, max(length(sl)));
                hdr[6] = "  "+align_left(hdr[6], sl[0].size()-2);

                ascii::write_table_hdr(ofilebase+"_gcat_"+strn(ig+1)+".cat", 0,
                    hdr, sq1, sq2, sfm, sdv, sz, sdz, sl
                );
            } else {
                fits::write_table(ofilebase+"_gcat_"+strn(ig+1)+".fits",
                    "qflag1", qflag1, "qflag2", qflag2, "fmatch", fmatch,
                    "maxdv", dvmax, "z", zs, "zerr", dzs, "lines", lids
                );
            }
        }
    }

    if (ascii) {
        vec1s hdr = {"qflag1", "qflag2", "fmatch", "maxdv", "z", "zerr", "lines"};

        // Manual handling of column widths
        vec1s sq1 = strna(gqf1);
        vec1s sq2 = strna(gqf2);
        vec1s sfm = strna(round(100.0*gfm)/100.0);
        vec1s sdv = strna(round(gdv));
        vec1s sz = strna(gz);
        vec1s sdz = strna(gze);
        vec1s sl = gl;

        sq1 = align_right(sq1, max(length(sq1))+8);
        hdr[0] = align_right(hdr[0], sq1[0].size());
        sq2 = align_right(sq2, max(length(sq2))+8);
        hdr[1] = align_right(hdr[1], sq2[0].size());
        sfm = align_right(sfm, max(length(sfm))+8);
        hdr[2] = align_right(hdr[2], sfm[0].size());
        sdv = align_right(sdv, max(length(sdv))+8);
        hdr[3] = align_right(hdr[3], sdv[0].size());
        sz = align_right(sz, max(length(sz))+4);
        hdr[4] = align_right(hdr[4], sz[0].size());
        sdz = align_right(sdz, max(length(sdz))+4);
        hdr[5] = align_right(hdr[5], sdz[0].size());
        sl  = "  "+align_left(sl, max(length(sl)));
        hdr[6] = "  "+align_left(hdr[6], sl[0].size()-2);

        ascii::write_table_hdr(ofilebase+"_gcat.cat", 0,
            hdr, sq1, sq2, sfm, sdv, sz, sdz, sl
        );
    } else {
        fits::write_table(ofilebase+"_gcat.fits",
            "qflag1", gqf1, "qflag2", gqf2, "fmatch", gfm,
            "maxdv", gdv, "z", gz, "zerr", gze, "lines", gl
        );
    }

    return 0;
}

// Function to grow a segmentation map within an allowed mask. The segmentation map
// should contain integer values above 0 indicating the IDs of different segments. Each
// segment has its own dedicated ID, which must simply be unique in this map. IDs can
// start from any integer value, and can be disjoint. The mask must be set to 'true' in
// the regions where growth is allowed. When two segments compete for growth over the
// same pixels of the mask, the biggest segment (in terms of number of pixels) will be
// prefered.
vec2u grow_within(vec2u map, vec2b mask) {
    phypp_check(map.dims == mask.dims, "incompatible dimensions between map and mask "
        "(", map.dims, " vs. ", mask.dims, ")");

    struct obj_state {
        uint_t id;
        uint_t npix = 0;
        std::vector<uint_t> oy, ox;
    };

    std::vector<obj_state> states;

    // Initialize states: identify segments and their boundaries where growth is allowed
    std::vector<uint_t> toy, tox;

    for (uint_t y : range(map.dims[0]))
    for (uint_t x : range(map.dims[1])) {
        if (map.safe(y,x) == 0) continue;

        // Found a guy
        states.push_back(obj_state());
        auto& state = states.back();
        state.id = map.safe(y,x);

        vec2u tmap = map;
        vec2b tmask = mask;

        toy.clear(); tox.clear();

        auto process_point = [&toy,&tox,&state,&tmap,&tmask](uint_t ty, uint_t tx) {
            tmap.safe(ty,tx) = 0;
            ++state.npix;

            auto check_add = [&toy,&tox,&state,&tmap,&tmask](uint_t tty, uint_t ttx) {
                if (tmap.safe(tty,ttx) == state.id) {
                    toy.push_back(tty);
                    tox.push_back(ttx);
                } else if (tmask.safe(tty,ttx)) {
                    tmask.safe(tty,ttx) = false;
                    state.oy.push_back(tty);
                    state.ox.push_back(ttx);
                }
            };

            if (ty != 0)              check_add(ty-1,tx);
            if (ty != tmap.dims[0]-1) check_add(ty+1,tx);
            if (tx != 0)              check_add(ty,tx-1);
            if (tx != tmap.dims[1]-1) check_add(ty,tx+1);
        };

        process_point(y, x);

        while (!tox.empty()) {
            uint_t ty = toy.back(); toy.pop_back();
            uint_t tx = tox.back(); tox.pop_back();
            process_point(ty, tx);
        }
    }

    // Priority given to the biggest
    std::sort(states.begin(), states.end(), [](const obj_state& s1, const obj_state& s2) {
        return s1.npix > s2.npix;
    });

    // Grow each component by one pixel at a time, in the pre-established order
    bool starved = false;
    while (!starved) {
        starved = true;
        for (uint_t i : range(states)) {
            auto& state = states[i];

            if (state.ox.empty()) continue;
            starved = false;

            auto process_point = [&state,&map,&mask](uint_t ty, uint_t tx) {
                map.safe(ty,tx) = state.id;
                ++state.npix;

                auto check_add = [&state,&map,&mask](uint_t tty, uint_t ttx) {
                    if (map.safe(tty,ttx) == 0 && mask.safe(tty,ttx)) {
                        mask.safe(tty,ttx) = false;
                        state.oy.push_back(tty);
                        state.ox.push_back(ttx);
                    }
                };

                if (ty != 0)             check_add(ty-1,tx);
                if (ty != map.dims[0]-1) check_add(ty+1,tx);
                if (tx != 0)             check_add(ty,tx-1);
                if (tx != map.dims[1]-1) check_add(ty,tx+1);
            };

            toy = state.oy;   tox = state.ox;
            state.oy.clear(); state.ox.clear();

            while (!tox.empty()) {
                uint_t ty = toy.back(); toy.pop_back();
                uint_t tx = tox.back(); tox.pop_back();
                process_point(ty, tx);
            }
        }
    }

    return map;
}

void print_help(const std::map<std::string,line_t>& db) {
    using namespace format;

    print("cdetect v1.0");
    print("usage: cdetect <kmos_cube.fits> [options]");
    print("");
    print("Main parameter and program workflow:");
    paragraph("'kmos_cube.fits' should be a cube created by the KMOS pipeline, with at "
        "least 2 extensions: the first is empty (KMOS convention), and the second "
        "contains the flux. If a third extension is present and contains the uncertainty, "
        "it is only used if 'emethod=pipeline' (see below). Else, the uncertainty is "
        "derived from the flux cube itself. Then, the program filters the cube, binning "
        "spectrally and smoothing spatially if requested, to detect high S/N features in "
        "the IFU. It generates a segmentation cube containing the position and extents of "
        "the detections (*_seg.fits), and a catalog listing these detections (*_cat.fits). "
        "In a second pass, it combines these independent spectral detections into "
        "spatially coherent sources (based on their distance on the sky) and tries to "
        "match the detected spectral features with known emission lines to determine the "
        "redshift of each source. The most likely redshift of each source is saved in a "
        "catalog (*_gcat.fits), and all the possible redshifts of one particular source "
        "'X' is saved in a dedicated catalog (*_gcat_X.fits). To each redshift solution "
        "is associated a number of quality flags and values. 'qflag1' is the number of "
        "unique emission lines that are used to determine the redshift (blended lines, "
        "such as Halpha and [NII], count as one). Solutions with 'qflag1=1' are not "
        "robust unless you have strong external constraints from the broad band "
        "photometry. 'qflag2' is the number of matched spectral feature with S/N greater "
        "than 5. NB: one line can be matched to multiple spectral features if the line is "
        "spectrally extended. 'fmatch' is the fraction of the detected spectral features "
        "that were matched to known emission lines. 'maxdv' is the largest "
        "velocity offset between two matched lines in km/s. An uncertainty on the derived "
        "redshift is computed based solely on 1) the spectral resolution of the binned "
        "spectrum and 2) the number of matched lines. Lastly, the list of matched lines "
        "is also given, with the following format: 'match1,match2,...', where 'matchX' "
        "lists the lines that are matched to one spectral feature, and the format is "
        "'id:line1/line2/line3/...' where 'id' is the ID of the spectral feature and "
        "'lineY' is the code name of the line that was matched. The 'matchX' are sorted "
        "by decreasing S/N.");
    print("Available lines for matching:");
    print_available_lines(db);
    paragraph("\nNote: you can add your own lines either by modifying the source code of the "
        "program, or directly into the command line arguments. In the 'lines=[...]' "
        "parameter, you can indeed create a new line with the synthax 'name:lambda'. In this "
        "case, 'name' can be whatever you want (should not contain spaces), 'lambda' must be "
        "the rest-frame wavelength of the line in microns. For example, to add Lyman alpha: "
        "lines=[lyalpha:0.12157].");
    print("Available options for source detection (in order of importance):");
    bullet("expmap=...", "Must be a 2D FITS file containing the exposure map of the IFU "
        "in units of exposure time or number of exposures (if exposure time per exposure "
        "is constant). This map will be used to estimate more finely the uncertainties. "
        "In addition, it can be combined with the 'minexp' option to limit the analyzed "
        "region of the IFU. Default is not to use any exposure map, but it is recommended "
        "to use one if you have it.");
    bullet("minexp=...", "Must be a number, and is only useful if 'expmap' is also "
        "provided. It defines the minimum value in 'expmap' that will be considered in "
        "the analysis. Pixels below this value will be flagged out and discarded. This "
        "can be used to discard regions with poor coverage, for which the uncertainty "
        "can be hard to estimate properly. The default is to use all valid pixels. It is "
        "recommended to use this option and flag out the regions with only a handful of "
        "exposures, where strong outliers can be found (cosmic rays, detector hot pixels, "
        "etc).");
    bullet("maxdist=...", "Must be a number. It defines the maximum distance from the IFU "
        "center for a detection (in pixels). All detections further than this distance will "
        "be discarded. This helps removing spurious detections, assuming the IFU was well "
        "centered on the target and assuming that the line emission is not substantially "
        "shifted from the center of mass of the target. Default is 3 pixels.");
    bullet("emethod=...", "Must be a string. It defines the method used to estimate the "
        "uncertainty on each pixel in the cube. Possible values are the following. "
        "'pipeline': use the uncertainty estimated by the KMOS pipeline, which is "
        "provided with the cube. 'stddev': compute the uncertainty of each wavelength "
        "slice from the standard deviation of the pixels (renormalized for exposure). "
        "'mad': same as stddev, but using the median absolute deviation. "
        "'stddevneg': same as stddev, but only using the negative pixels in the map. "
        "'madneg': same as stddevneg but using the median absolute deviation. Generally "
        "speaking, the standard deviation is more sensitive to outliers than the median "
        "absolute deviation, so it will tend to give more conservative (if not "
        "clearly overestimated) uncertainties. It will also tend to be biased high by the "
        "presence of genuine sources in the map. For this reason, there are alternative "
        "version of both stddev and mad that only use the negative pixels; i.e., are free "
        "from contamination by sources. Except for the 'pipeline' estimate, all the "
        "methods assume that your source(s) only occupy a small fraction of the space in "
        "the IFU.");
    bullet("spectral_bin=...", "Must be an integer number. It defines the number of "
        "spectral pixels that should be combined for the detection. The larger the value, "
        "the more pixels will be used for the detection, hence the higher the S/N. "
        "However, if an emission line is very narrow and only present in one or a few "
        "spectral elements, using too large binning values will dilute the line's "
        "signal, and you will actually loose S/N. The optimal value therefore depends on "
        "the expected line width. Default is no binning, but you should experiment.");
    bullet("spatial_smooth=...", "Must be a number. It defines the size of the Gaussian "
        "profile that will be used for spatial convolution, in order to increase the S/N. "
        "The width of this Gaussian should be at least the width of the Point Spread "
        "Function of KMOS, which depends on the seeing at which your data was observed. "
        "It can be chosen larger than that if you expect extended emission. As for "
        "spectral binning, the optimal value of this parameter depends on your targets. "
        "Set this value to zero to disable the convolution. Default value is 1.2 pixels.");
    bullet("error_scale=...", "Must be a number. It can be used to scale up or down the "
        "estimated uncertainties globally. It usually happens that uncertainties are "
        "underestimated, and the program produces many false positives. In this case "
        "you can raise this value, typically to 1.3. Default is 1, so no rescaling.");
    bullet("snr_det=...", "Must be a number. It defines the minimum S/N ratio to consider "
        "for a detection. Default is 5.");
    bullet("snr_source=...", "Must be a number. It defines the minimum S/N ratio to consider "
        "to define the spatial extents of a detection. Default is 3.");
    bullet("lambda_pad", "Must be an integer. It defines the number of wavelength element "
        "that are ignored both at the beginning and end of the spectrum. Default is 5 "
        "elements. This is used to flag out invalid and poorly covered spectral regions "
        "which often trigger spurious detections.");
    bullet("save_cubes", "Set this flag to write to disk the filtered flux and uncertainty "
        "cubes that are used to perform the detection. They will be saved in *_filt.fits.");
    print("\nAvailable options for redshift identification (in order of importance):");
    bullet("disable_zsearch", "Set this flag to skip the redshift identification step "
        "and only list the spectral detections. You would typically set this flag on if "
        "you are looking for continuum emission.");
    bullet("maxdpos=...", "Must be a number. It defines the maximum distance on the sky (in "
        "pixels) between two spectral sources to be considered as part of the same "
        "object. Default is 5 pixels, which is fairly generous. You may want to decrease "
        "this value if you know that your sources are very compact and that there should "
        "not be spatial offsets between the sources and their line emission (which could "
        "happen, e.g., if there is a cloud of outflowing gas). Using a too large value "
        "will lead to more spurious associations of noise fluctuations.");
    bullet("single_source", "Set this flag if you want all the spectral detections to be "
        "associated to a unique source (this is equivalent to choosing 'maxdpos' to a "
        "value larger than the size of the IFU).");
    bullet("zhint=[...]", "Must be a vector of two numbers. It defines the allowed redshift "
        "range when searching for line emission. The first number is the lowest redshift "
        "and the second number is the highest redshift. Example: 'zhint=[1,2]' will only "
        "search for lines between z=1 and z=2. Default is to consider all possible "
        "positive redshifts. This can be used to narrow down the possibilities, for "
        "example if you already have strong constraints from the broad band photometry or "
        "low resolution spectroscopy.");
    bullet("maxdv=...", "Must be a number. It defines the maximum allowed velocity offset "
        "between two lines that belong to the same source. Default is 200 km/s. This is, "
        "if you wish, the tolerence threshold on the line identification. Increasing this "
        "number will allow new redshift solutions where lines can be strongly offset from "
        "one another, likely because of outflows or multiple components with different "
        "velocities. It will also increase the chance of spurious detection.");
    bullet("minqflag=...", "Must be an integer. It defines the minimum allowed value of "
        "'qflag1'; redshift solutions with a lower quality flag will be discarded as "
        "unreliable. Default is to keep all solutions.");
    print("\nAvailable generic options:");
    bullet("outdir=...", "Name of the directory into which the output files should be created. "
        "Default is the current directory.");
    bullet("ascii", "Set this flag to save the output catalogs in ASCII format rather "
        "than FITS.");
    bullet("verbose", "Set this flag to print the progress of the detection process in "
        "the terminal. Can be useful if something goes wrong, or just to understand what "
        "is going on.");
}

void print_available_lines(const std::map<std::string,line_t>& db) {
    for (auto& l : db) {
        if (l.second.lambda.size() == 1) {
            print("  - ", l.first, ", lambda=", l.second.lambda[0]);
        } else {
            for (uint_t il : range(l.second.lambda)) {
                print("  - ", l.first+"-"+strn(il+1), ", lambda=", l.second.lambda[il]);
            }
        }
    }
}
