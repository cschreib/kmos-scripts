#include <phypp.hpp>

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print("usage: cube2slit <cube.fits> [out=... slit_width=... verbose]");
        return 0;
    }

    double slit_width = 0.7; // [arcsec]
    bool verbose = false;
    uint_t hdu = 1;
    uint_t error_hdu = 2;
    double center_x = dnan;
    double center_y = dnan;
    std::string out_file;

    read_args(argc-1, argv+1, arg_list(
        slit_width, verbose, name(out_file, "out"), hdu, error_hdu, center_x, center_y
    ));

    std::string infile = argv[1];
    if (out_file.empty()) {
        out_file = file::remove_extension(file::get_basename(infile))+"_slit.fits";
    }

    fits::input_image fimg(infile);
    fimg.reach_hdu(hdu);

    // Read wavelength/frequency axis WCS
    uint_t nlam = 0;
    std::string cunit;
    double cdelt = 1, crpix = 1, crval = 1;
    vec1s missing;
    fimg.read_keyword("CUNIT3", cunit); // optional
    if (!fimg.read_keyword("NAXIS3", nlam))  missing.push_back("NAXIS3");
    if (!fimg.read_keyword("CDELT3", cdelt)) {
        if (!fimg.read_keyword("CD3_3", cdelt)) {
            missing.push_back("CDELT3 or CD3_3");
        }
    }
    if (!fimg.read_keyword("CRPIX3", crpix)) missing.push_back("CRPIX3");
    if (!fimg.read_keyword("CRVAL3", crval)) missing.push_back("CRVAL3");
    if (!missing.empty()) {
        error("could not read WCS information for wavelength axis (axis 3)");
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

    if (verbose) {
        double l0 = crval + cdelt*(1 - crpix);
        double l1 = crval + cdelt*(nlam - crpix);
        if (l0 > l1) std::swap(l0, l1);

        std::string unit = cunit.empty() ? "[unknown unit?]" : cunit;
        note("cube is covering ", l0, " to ", l1, " ", unit);
    }

    // Find pixel scale
    double aspix = 0;
    if (!fimg.read_keyword("CDELT1", aspix)) {
        if (!fimg.read_keyword("CD1_1", aspix)) {
            error("could not read WCS information for sky axis (axis 1)");
            note("missing keyword CDELT1 or CD1_1");
            return 1;
        }
    }

    aspix = abs(aspix*3600);
    if (verbose) note("cube has a pixel scale of ", aspix, "\"");

    uint_t navg = ceil(slit_width/aspix)/2;
    if (verbose) note("using a slit of width ", 2*navg + 1, " pixels");

    fits::output_image foimg(out_file);
    foimg.write_empty();

    vec3d flx3d, err3d;
    fimg.reach_hdu(hdu);
    fimg.read(flx3d);
    fimg.reach_hdu(error_hdu);
    fimg.read(err3d);

    uint_t nslit = flx3d.dims[1];
    if (2*navg+1 > flx3d.dims[2]) {
        error("slit would be larger than data cube (", flx3d.dims[2], " pixels, or ",
            flx3d.dims[2]*aspix, "\")");
        note("please choose a smaller value for 'slit_width'");
        return 1;
    }

    uint_t i0, i1;
    if (!is_finite(center_x)) {
        i0 = flx3d.dims[2]/2 - navg;
        i1 = flx3d.dims[2]/2 + navg;
    } else {
        double cp = round(center_x);
        if (cp-navg < 0 || cp+navg > flx3d.dims[2]-1) {
            error("slit position reaches the edge of the data cube");
            return 1;
        }

        i0 = uint_t(cp) - navg;
        i1 = uint_t(cp) + navg;
    }

    vec2d flx2d(nslit, flx3d.dims[0]);
    vec2d err2d(nslit, flx3d.dims[0]);
    for (uint_t l : range(flx3d.dims[0]))
    for (uint_t p : range(flx3d.dims[1])) {
        vec1d lflx = flx3d(l,p,i0-_-i1);
        vec1d lerr = err3d(l,p,i0-_-i1);
        vec1u idl = where(is_finite(lflx) && is_finite(lerr) && lerr > 0);
        auto pp = optimal_mean(lflx[idl], lerr[idl]);
        flx2d(p,l) = pp.first;
        err2d(p,l) = pp.second;
    }

    foimg.reach_hdu(1);
    foimg.write(flx2d);

    auto write_wcs = [&]() {
        if (!cunit.empty()) foimg.write_keyword("CUNIT1", cunit);
        if (!ctype.empty()) foimg.write_keyword("CTYPE1", ctype);
        foimg.write_keyword("CDELT1", cdelt);
        foimg.write_keyword("CRPIX1", crpix);
        foimg.write_keyword("CRVAL1", crval);

        foimg.write_keyword("CTYPE2", "SLIT");
        foimg.write_keyword("CUNIT2", "arcsec");
        foimg.write_keyword("CDELT2", aspix);
        foimg.write_keyword("CRPIX2", flx2d.dims[0]/2 + 1);
        foimg.write_keyword("CRVAL2", 0.0);
    };

    write_wcs();

    foimg.reach_hdu(2);
    foimg.write(err2d);

    return 0;
}
