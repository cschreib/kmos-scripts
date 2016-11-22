#include <phypp.hpp>
#include <phypp/math/mpfit.hpp>

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print("usage: make_profile <cube.fits> [verbose]");
        return 0;
    }

    bool verbose = false;
    uint_t wave_padding = 10;
    uint_t edge_padding = 2;
    uint_t hdu = npos;
    uint_t err_hdu = npos;
    bool no_error = false;

    bool do_sigma_clip = true;          // enable/disable sigma clipping of outliers
    double sigma_clip_threshold = 5.0;  // significance threshold for rejecting outliers

    read_args(argc-1, argv+1, arg_list(verbose, wave_padding, hdu, err_hdu, no_error,
        do_sigma_clip, sigma_clip_threshold));

    if (hdu != npos && err_hdu == npos && !no_error) {
        error("please specify both hdu=... and err_hdu=... (unless no_error is set)");
        return 1;
    }

    if (hdu == npos && err_hdu == npos) {
        hdu = 1; err_hdu = 2;
    }

    vec1d lines = {
        // H band
        1.4564, 1.4606, 1.4666, 1.4699, 1.4784, 1.4832, 1.4864, 1.4888,
        1.4909, 1.4932, 1.4803, 1.5054, 1.5057, 1.5069, 1.5088, 1.5188,
        1.5241, 1.5287, 1.5332, 1.5395, 1.5433, 1.5505, 1.5515, 1.5541,
        1.5544, 1.5570, 1.5598, 1.5632, 1.5656, 1.5702, 1.5832, 1.5848,
        1.5869, 1.5973, 1.6031, 1.6080, 1.6129, 1.6195, 1.6236, 1.6317,
        1.6352, 1.6388, 1.6415, 1.6444, 1.6477, 1.6503, 1.6555, 1.6611,
        1.6690, 1.6693, 1.6706, 1.6709, 1.6733, 1.6841, 1.6904, 1.6955,
        1.7009, 1.7079, 1.7124, 1.7211, 1.7249, 1.7283, 1.7330, 1.7350,
        1.7358, 1.7386, 1.7428, 1.7450, 1.7505, 1.7652, 1.7810, 1.7879,
        1.7991, 1.7995, 1.8067, 1.8120, 1.8209, 1.8216,

        // H+K gap
        1.8645, 1.8809, 1.8843,

        // K band
        1.9200, 1.9207, 1.9246, 1.9250, 1.9277, 1.9283, 1.9309, 1.9350,
        1.9399, 1.9560, 1.9642, 1.9701, 1.9753, 1.9774, 1.9841, 2.0010,
        2.0034, 2.0278, 2.0342, 2.0414, 2.0500, 2.0566, 2.0731, 2.0862,
        2.0909, 2.1177, 2.1233, 2.1252, 2.1509, 2.1541, 2.1710, 2.1804,
        2.1874, 2.1957, 2.2125
    };

    vec1d bands_low = {1.8256, 1.8678, 1.8898, 1.8981, 1.9109};
    vec1d bands_up  = {1.8609, 1.8779, 1.8939, 1.9072, 1.9169};

    std::string infile = argv[1];
    std::string outfile_model = file::remove_extension(file::get_basename(infile))+"_profile.fits";
    std::string outfile_collapsed = file::remove_extension(file::get_basename(infile))+"_collapsed.fits";

    // Read cube
    vec3d flx3d, err3d;
    fits::input_image fimg(infile);
    fimg.reach_hdu(hdu);

    if (fimg.image_dims().size() != 3) {
        error("this HDU does not contain a data cube");
        return 1;
    }

    fimg.read(flx3d);

    // Read 2D astrometry
    std::string ctype1, ctype2;
    std::string cunit1, cunit2;
    double crpix1, crpix2, crval1, crval2, cdelt1, cdelt2;
    double cd11, cd12, cd21, cd22;

    vec1s missing;
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

    double cdelt = 1, crpix = 1, crval = 1;
    if (!fimg.read_keyword("CDELT3", cdelt)) missing.push_back("CDELT3");
    if (!fimg.read_keyword("CRPIX3", crpix)) missing.push_back("CRPIX3");
    if (!fimg.read_keyword("CRVAL3", crval)) missing.push_back("CRVAL3");
    if (!missing.empty()) {
        error("could not read WCS information for wavelength axis");
        note("missing keyword", missing.size() > 1 ? "s " : " ", collapse(missing, ", "));
        return 1;
    }

    uint_t nlam = flx3d.dims[0];
    vec1d lam = crval + cdelt*(dindgen(nlam) + (1 - crpix));

    vec1b flagged(lam.dims);
    for (uint_t l : range(lines)) {
        flagged = flagged || (lam >= lines[l]-0.0006 && lam <= lines[l]+0.0006);
    }
    for (uint_t l : range(bands_low)) {
        flagged = flagged || (lam >= bands_low[l]-cdelt && lam <= bands_up[l]+cdelt);
    }

    if (!no_error) {
        fimg.reach_hdu(err_hdu);
        fimg.read(err3d);
    }

    // Flag edges
    uint_t i0 = where_first(is_finite(flx3d(_,flx3d.dims[1]/2,flx3d.dims[2]/2)));
    uint_t i1 = where_last(is_finite(flx3d(_,flx3d.dims[1]/2,flx3d.dims[2]/2)));
    if (i0 == npos) {
        error("cube does not contain any valid pixel");
        return 1;
    }

    if (wave_padding != 0) {
        i0 += wave_padding-1;
        i1 -= wave_padding-1;
        flagged[_-i0] = true;
        flagged[i1-_] = true;
    }

    // Compute collapsed profile
    vec2d flx2d(flx3d.dims[1], flx3d.dims[2]);
    vec2d err2d(flx3d.dims[1], flx3d.dims[2]);
    for (uint_t y : range(flx3d.dims[1]))
    for (uint_t x : range(flx3d.dims[2])) {
        vec1d tflx = flx3d(_,y,x);
        if (no_error) {
            vec1u idl = where(is_finite(tflx) && !flagged);
            flx2d(y,x) = median(tflx[idl]);
            err2d(y,x) = 1.0;
        } else {
            vec1d terr = err3d(_,y,x);
            vec1u idl = where(is_finite(tflx) && is_finite(terr) && terr > 0 && !flagged);

            // Apply sigma clipping
            if (do_sigma_clip) {
                vec1d res = abs(tflx - median(tflx[idl]))/terr;
                double rms = 1.48*median(res[idl]);
                idl = where(is_finite(tflx) && is_finite(terr) && terr > 0 && !flagged && res < sigma_clip_threshold*rms);
            }

            auto pp = optimal_mean(tflx[idl], terr[idl]);
            flx2d(y,x) = pp.first;
            err2d(y,x) = pp.second;
        }
    }

    // Flag edges
    if (edge_padding != 0) {
        flx2d(_-(edge_padding-1),_) = dnan;
        flx2d((flx2d.dims[0]-edge_padding)-_,_) = dnan;
        flx2d(_,_-(edge_padding-1)) = dnan;
        flx2d(_,(flx2d.dims[1]-edge_padding)-_) = dnan;
    }

    // Fit profile
    struct base_t {
        vec1d x, y;
    } xx;

    xx.x = flatten(generate_img(flx2d.dims, [&](int_t,    int_t tx) { return tx - flx2d.dims[1]/2.0; }));
    xx.y = flatten(generate_img(flx2d.dims, [&](int_t ty, int_t)    { return ty - flx2d.dims[0]/2.0; }));

    vec1u idf = where(is_finite(flx2d) && is_finite(err2d) && err2d > 0);
    vec1d tflx2d = 1e17*flx2d[idf]; vec1d terr2d = 1e17*err2d[idf];
    xx.x = xx.x[idf]; xx.y = xx.y[idf];

    vec1d startp = {0.0, max(tflx2d), 0.0, 0.0, 1.8};
    auto model = [](const base_t& x, const vec1d& p) {
        return p[0] + abs(p[1])*exp(-(sqr(x.x - p[2]) + sqr(x.y - p[3]))/(2.0*sqr(p[4])));
    };

    auto res = mpfitfun(tflx2d, terr2d, xx, model, startp);

    print(res.params);
    print(res.params/res.errors);

    // Build model
    xx.x = flatten(generate_img(flx2d.dims, [&](int_t,    int_t tx) { return tx - flx2d.dims[1]/2.0; }));
    xx.y = flatten(generate_img(flx2d.dims, [&](int_t ty, int_t)    { return ty - flx2d.dims[0]/2.0; }));
    vec2d model2d = reform(model(xx, res.params), flx2d.dims);

    // Save observed profile + error + residual
    fits::output_image foimg_obs(outfile_collapsed);
    foimg_obs.write_empty();

    auto write_wcs = [&](fits::output_image& oimg) {
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

    foimg_obs.reach_hdu(1);
    foimg_obs.write(flx2d);
    write_wcs(foimg_obs);
    foimg_obs.reach_hdu(2);
    foimg_obs.write(err2d);
    write_wcs(foimg_obs);
    foimg_obs.reach_hdu(3);
    foimg_obs.write(flx2d - 1e-17*model2d);
    write_wcs(foimg_obs);

    print("RMS error: ", stddev((flx2d[idf] - 1e-17*model2d[idf])/err2d[idf]));

    // Normalize model to unit flux + remove the fitted background
    res.params[0] = 0;
    res.params[1] = 1;
    model2d = reform(model(xx, res.params), flx2d.dims);
    model2d /= 2.0*dpi*sqr(res.params[4]);

    // Save model
    fits::output_image foimg(outfile_model);
    foimg.write(model2d);
    write_wcs(foimg);

    return 0;
}
