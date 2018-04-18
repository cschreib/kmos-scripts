#include <phypp.hpp>
#include <phypp/astro/qstack.hpp>
#include <phypp/math/mpfit.hpp>

// Integral of a Gaussian profile between x0 and x1
template <typename TX0, typename TX1>
auto gauss_integral(TX0&& x0, TX1&& x1, double xc, double xw, double amp = 1.0) -> decltype(sqrt(x0)+sqrt(x1)) {
    // Average this from t = x-0.5 to x+0.5
    // exp(-sqr(t - xc)/(2.0*sqr(xw)))/(sqrt(2*dpi)*xw)

    return 0.5*(amp/(x1 - x0))*(
        erf((x1 - xc)/(sqrt(2.0)*xw)) -
        erf((x0 - xc)/(sqrt(2.0)*xw))
    );
};

// Integral of a 2D Gaussian profile between [x0,x1] and [y0,y1]
template <typename TX0, typename TX1>
auto gauss_integral_2d(TX0&& x0, TX1&& x1, TX0&& y0, TX1&& y1, double xc, double yc, double w, double amp = 1.0) -> decltype(sqrt(x0)+sqrt(x1)) {
    return amp*gauss_integral(x0,x1,xc,w,1.0)*gauss_integral(y0,y1,yc,w,1.0);
};

template <typename T>
double get_radius_raw(const vec<2,T>& img, double y0, double x0) {
    double tot = 0.0;
    for (uint_t y : range(img.dims[0]))
    for (uint_t x : range(img.dims[1])) {
        tot += img.safe(y,x)*(sqr(y - y0) + sqr(x - x0));
    }

    return sqrt(tot);
}

double get_radius(vec2d img, double y0, double x0) {
    // Compute first estimate
    img /= total(img); // image needs to be normalized first
    double rr = get_radius_raw(img, y0, x0);

    // Simulate various profiles
    vec1d td = rgen_step(0.1, max(img.dims[0], img.dims[1]), 0.3);
    vec1d md(td.dims);

    vec2d timg(img.dims);
    for (uint_t i : range(td)) {
        for (uint_t y : range(timg.dims[0]))
        for (uint_t x : range(timg.dims[1])) {
            timg.safe(y,x) = gauss_integral_2d(y-0.5, y+0.5, x-0.5, x+0.5, y0, x0, td[i]);
        }

        md[i] = get_radius_raw(timg, y0, x0);
    }

    // Correct measured value using simulation
    // The additional factor of sqrt(log(4.0)) is the conversion
    // from a Gaussian sigma into a half-light radius (=1.177)
    return interpolate(td, md, rr)*sqrt(log(4.0));
}

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print("usage: make_profile <cube.fits> [verbose]");
        return 0;
    }

    bool verbose = false;
    uint_t wave_padding = 10;
    uint_t edge_padding = 2;
    uint_t hdu = npos;
    uint_t error_hdu = npos;
    double gwidth_default = 1.2;
    double seeing = dnan;
    double seeing_hri = dnan;
    double dpix_step = 0.1;
    double max_dpix = 2;
    bool fix_gwidth = false;
    bool fix_pos = false;
    bool no_error = false;
    std::string hri, segmentation;

    bool do_sigma_clip = true;          // enable/disable sigma clipping of outliers
    double sigma_clip_threshold = 5.0;  // significance threshold for rejecting outliers

    read_args(argc-1, argv+1, arg_list(verbose, wave_padding, hdu, error_hdu, no_error,
        do_sigma_clip, sigma_clip_threshold, gwidth_default, fix_gwidth, hri, fix_pos,
        seeing, seeing_hri, segmentation, dpix_step, max_dpix, edge_padding));

    if (hdu != npos && error_hdu == npos && !no_error) {
        error("please specify both hdu=... and error_hdu=... (unless no_error is set)");
        return 1;
    }

    if (hdu == npos) {
        hdu = 1;
    }
    if (error_hdu == npos) {
        error_hdu = 2;
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
    fimg.reach_hdu(hdu);
    astro::wcs wcs(fimg.read_header());

    std::string ctype1, ctype2;
    std::string cunit1, cunit2;
    double crpix1, crpix2, crval1, crval2, cdelt1, cdelt2;
    double cd11, cd12, cd21, cd22;
    bool pos_use_cd = false;
    bool wave_use_cd = false;

    vec1s missing;
    if (!fimg.read_keyword("CRPIX1", crpix1)) missing.push_back("CRPIX1");
    if (!fimg.read_keyword("CRPIX2", crpix2)) missing.push_back("CRPIX2");
    if (!fimg.read_keyword("CRVAL1", crval1)) missing.push_back("CRVAL1");
    if (!fimg.read_keyword("CRVAL2", crval2)) missing.push_back("CRVAL2");
    if (!fimg.read_keyword("CDELT1", cdelt1) || !fimg.read_keyword("CDELT2", cdelt2)) {
        if (!fimg.read_keyword("CD1_1",  cd11))   missing.push_back("CDELT1 or CD1_1");
        if (!fimg.read_keyword("CD1_2",  cd12))   missing.push_back("CD1_2");
        if (!fimg.read_keyword("CD2_1",  cd21))   missing.push_back("CD2_1");
        if (!fimg.read_keyword("CD2_2",  cd22))   missing.push_back("CDELT2 or CD2_2");
        pos_use_cd = true;
    }
    if (!fimg.read_keyword("CTYPE1", ctype1)) missing.push_back("CTYPE1");
    if (!fimg.read_keyword("CTYPE2", ctype2)) missing.push_back("CTYPE2");
    if (!fimg.read_keyword("CUNIT1", cunit1)) missing.push_back("CUNIT1");
    if (!fimg.read_keyword("CUNIT2", cunit2)) missing.push_back("CUNIT2");
    if (!missing.empty()) {
        error("could not read WCS information for spatial axes");
        note("missing keyword", missing.size() > 1 ? "s " : " ", collapse(missing, ","));
        return 1;
    }

    double cdelt = 1, crpix = 1, crval = 1;
    if (!fimg.read_keyword("CRPIX3", crpix)) missing.push_back("CRPIX3");
    if (!fimg.read_keyword("CRVAL3", cdelt)) missing.push_back("CRVAL3");
    if (!fimg.read_keyword("CDELT3", crval)) {
        if (!fimg.read_keyword("CD3_3",  cdelt)) missing.push_back("CDELT3 or CD3_3");
        wave_use_cd = true;
    }
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
        fimg.reach_hdu(error_hdu);
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
    vec2d flx2da(flx3d.dims[1], flx3d.dims[2]);
    vec2d flx2db(flx3d.dims[1], flx3d.dims[2]);
    vec2d err2d(flx3d.dims[1], flx3d.dims[2]);
    vec2d err2da(flx3d.dims[1], flx3d.dims[2]);
    vec2d err2db(flx3d.dims[1], flx3d.dims[2]);
    vec3b sigma_clipped(flx3d.dims);
    for (uint_t y : range(flx3d.dims[1]))
    for (uint_t x : range(flx3d.dims[2])) {
        vec1d tflx = flx3d(_,y,x);
        if (no_error) {
            vec1u idl = where(is_finite(tflx) && !flagged);
            flx2d(y,x) = median(tflx[idl]);
            err2d(y,x) = 1.0;
        } else {
            vec1d terr = err3d(_,y,x);
            vec1b basesel = is_finite(tflx) && is_finite(terr) && terr > 0 && !flagged;
            vec1u idl = where(basesel);

            // Apply sigma clipping
            if (do_sigma_clip && !idl.empty()) {
                vec1d res = abs(tflx - median(tflx[idl]))/terr;
                double rms = 1.48*median(res[idl]);
                sigma_clipped(_,y,x) = res > sigma_clip_threshold*rms;
                basesel = basesel && !sigma_clipped(_,y,x);
                idl = where(basesel);
            }

            auto pp = optimal_mean(tflx[idl], terr[idl]);
            flx2d(y,x) = pp.first;
            err2d(y,x) = pp.second;

            idl = where(basesel && lam <= mean(lam));
            pp = optimal_mean(tflx[idl], terr[idl]);
            flx2da(y,x) = pp.first;
            err2da(y,x) = pp.second;

            idl = where(basesel && lam > mean(lam));
            pp = optimal_mean(tflx[idl], terr[idl]);
            flx2db(y,x) = pp.first;
            err2db(y,x) = pp.second;
        }
    }

    // Flag edges
    if (edge_padding != 0) {
        auto mask_2d = [edge_padding](vec2d& img) {
            img(_-(edge_padding-1),_) = dnan;
            img((img.dims[0]-edge_padding)-_,_) = dnan;
            img(_,_-(edge_padding-1)) = dnan;
            img(_,(img.dims[1]-edge_padding)-_) = dnan;
        };

        mask_2d(flx2d);
        mask_2d(err2d);
        if (!no_error) {
            mask_2d(flx2da);
            mask_2d(flx2db);
            mask_2d(err2da);
            mask_2d(err2db);
        }
    }

    // Build jackniffe
    if (!no_error) {
        vec2d jcrms = (flx2da - flx2db)/sqrt(sqr(err2da) + sqr(err2db));
        double rms_rescale = max(1.0, 1.48*mad(jcrms[where(is_finite(jcrms))]));
        if (verbose) {
            note("rescaling uncertainties by ", rms_rescale);
        }

        err2d *= rms_rescale;
        err3d *= rms_rescale;
    } else {
        err2d = replicate(1e-20, flx2d.dims);
        err3d = replicate(1e-19, flx3d.dims);
    }

    // Build model
    vec2d model2d;
    vec2d profile2d;
    std::function<void(fits::output_image&)> add_keywords;
    if (!hri.empty()) {
        // Extract cutout from HRI
        double ra0, dec0, ra1, dec1;
        astro::xy2ad(wcs, wcs.dims[2]/2.0 + 1, wcs.dims[1]/2.0 + 1, ra0, dec0);
        astro::xy2ad(wcs, 1, 1, ra1, dec1);
        double cutsize = angdist(ra0, dec0, ra1, dec1);
        double hri_pixsize;
        if (!astro::get_pixel_size(hri, hri_pixsize)) {
            error("could not read astrometry from HRI");
            return 1;
        }
        double cube_pixsize;
        if (!astro::get_pixel_size(wcs, cube_pixsize)) {
            error("could not read astrometry from cube");
            return 1;
        }

        double dpp = cube_pixsize/hri_pixsize;

        uint_t hsize = ceil(cutsize/hri_pixsize);
        if (verbose) {
            note("extracting ", 2*hsize+1, "x", 2*hsize+1, " cutout from HRI");
        }

        vec3d ccutout;
        vec1u cids;
        qstack_params qopts;
        qopts.save_offsets = true;
        auto qout = qstack({ra0}, {dec0}, hri, hsize, ccutout, cids, qopts);

        if (cids.empty()) {
            error("this galaxy is not covered by the HRI");
            return 1;
        }

        vec2d cutout = ccutout(0,_,_);

        // Create astrometry of cutout
        fits::header hri_hdr = fits::read_header(hri);
        fits::setkey(hri_hdr, "CRVAL1", ra0);
        fits::setkey(hri_hdr, "CRVAL2", dec0);
        fits::setkey(hri_hdr, "NAXIS1", cutout.dims[1]);
        fits::setkey(hri_hdr, "NAXIS2", cutout.dims[0]);
        fits::setkey(hri_hdr, "CRPIX1", hsize+1+qout.dx[0]);
        fits::setkey(hri_hdr, "CRPIX2", hsize+1+qout.dy[0]);
        astro::wcs hri_wcs(hri_hdr);

        // If asked, mask nearby sources from the HRI using a segmentation map
        if (!segmentation.empty()) {
            if (verbose) {
                note("reading segmentation and flagging nearby sources");
            }

            double seg_pixsize;
            if (!astro::get_pixel_size(segmentation, seg_pixsize)) {
                error("could not read astrometry from segmentation map");
                return 1;
            }

            uint_t shsize = ceil(cutsize/seg_pixsize);
            vec3d scutout;
            vec1u sids;
            auto sqout = qstack({ra0}, {dec0}, segmentation, shsize, scutout, sids, qopts);

            if (sids.empty()) {
                warning("this galaxy is not covered by the segmentation map");
            } else {
                vec2d seg = scutout(0,_,_);

                // Mask out the target from the segmentation
                double icent = seg(shsize,shsize);
                if (icent == 0) {
                    // Try the neighboring pixels
                    vec1d tdx = {-1, -1, 0, 1, 1, 1, 0, -1};
                    vec1d tdy = {0, 1, 1, 1, 0, -1, -1, -1};
                    double ddp = cube_pixsize/seg_pixsize;

                    for (uint_t i : range(tdx)) {
                        icent = seg(shsize+int_t(round(tdy[i]*ddp)), shsize+int_t(round(tdx[i]*ddp)));
                        if (icent > 0) break;
                    }
                }
                if (icent > 0) {
                    seg[where(seg == icent)] = 0;
                }

                // Regrid segmentation to HRI
                fits::header seg_hdr = fits::read_header(segmentation);
                fits::setkey(seg_hdr, "CRVAL1", ra0);
                fits::setkey(seg_hdr, "CRVAL2", dec0);
                fits::setkey(seg_hdr, "NAXIS1", seg.dims[1]);
                fits::setkey(seg_hdr, "NAXIS2", seg.dims[0]);
                fits::setkey(seg_hdr, "CRPIX1", shsize+1+sqout.dx[0]);
                fits::setkey(seg_hdr, "CRPIX2", shsize+1+sqout.dy[0]);
                astro::wcs seg_wcs(seg_hdr);

                astro::regrid_interpolate_params ropts;
                ropts.method = interpolation_method::nearest;
                seg = astro::regrid_interpolate(seg, seg_wcs, hri_wcs, ropts);

                // Mask out other sources
                vec1u idb = where(is_finite(seg) && seg > 0);
                cutout[idb] = 0;
            }
        }

        if (verbose) {
            note("SNR of HRI: ", total(cutout)/(1.48*mad(cutout)));
        }

        // Smooth it by the seeing
        if (is_finite(seeing)) {
            double conv_sigma = seeing/2.335/hri_pixsize;
            bool no_convolve = false;
            if (is_finite(seeing_hri)) {
                if (seeing_hri < seeing) {
                    conv_sigma = sqrt(sqr(conv_sigma) - sqr(seeing_hri/2.335/hri_pixsize));
                } else {
                    // Seeing of the HRI is higher than that of the cube, do not convolve.
                    no_convolve = true;
                }
            }

            if (!no_convolve) {
                if (verbose) {
                    note("convolving HRI to seeing of cube");
                }

                vec2d kernel = gaussian_profile(cutout.dims, conv_sigma);
                cutout = convolve2d(cutout, kernel);
            }
        }

        if (verbose) {
            note("fitting HRI model to cube");
        }

        if (count(is_finite(cutout)) == 0) {
            error("HRI cutout only contains invalid pixels");
            return 1;
        }

        cutout[where(!is_finite(cutout))] = 0;
        cutout /= total(cutout[where(cutout > 2*1.48*mad(cutout))]);

        vec1u idf = where(is_finite(flx2d) && is_finite(err2d) && err2d > 0);
        vec1d tflx2d = 1e17*flx2d[idf]; vec1d terr2d = 1e17*err2d[idf];

        double bchi2 = dinf;
        double bamp = dnan;
        double berr = dnan;
        double bdx = dnan;
        double bbg = dnan;
        double bbge = dnan;
        double bdy = dnan;

        if (!fix_pos) {
            double dymin = -max_dpix;
            double dymax = +max_dpix;
            uint_t ndy = ceil(dymax - dymin);
            double dxmin = -max_dpix;
            double dxmax = +max_dpix;
            uint_t ndx = ceil(dxmax - dxmin);

            if (verbose) {
                note("position grid is ", ndy, "x", ndx, " with step ", 1.0);
            }

            // Regrid to cube pixel size and astrometry
            auto pg = progress_start(ndy*ndx);
            for (double dy : rgen(dymin, dymax, ndy))
            for (double dx : rgen(dxmin, dxmax, ndx)) {
                hri_wcs.w->crpix[0] = hsize+1+qout.dx[0]+dx*dpp;
                hri_wcs.w->crpix[1] = hsize+1+qout.dy[0]+dy*dpp;
                vec2d tc = astro::regrid_drizzle(cutout, hri_wcs, wcs);
                tc[where(!is_finite(tc))] = 0;

                auto res = linfit(tflx2d, terr2d, tc[idf], 1.0);
                if (res.chi2 < bchi2 && res.params[0] > 0) {
                    bchi2 = res.chi2;
                    bamp = res.params[0];
                    berr = res.errors[0];
                    bbg = res.params[1];
                    bbge = res.errors[1];
                    bdx = dx;
                    bdy = dy;
                }

                if (verbose) progress(pg, 131);
            }

            uint_t rfd = 1;
            dymin = bdy - rfd;
            dxmin = bdx - rfd;
            dymax = bdy + rfd;
            dxmax = bdx + rfd;
            ndy = ceil((dymax - dymin)/dpix_step);
            ndx = ceil((dxmax - dxmin)/dpix_step);

            if (verbose) {
                note("position grid is ", ndy, "x", ndx, " with step ", dpix_step);
            }

            pg = progress_start(ndy*ndx);
            for (double dy : rgen(dymin, dymax, ndy))
            for (double dx : rgen(dxmin, dxmax, ndx)) {
                hri_wcs.w->crpix[0] = hsize+1+qout.dx[0]+dx*dpp;
                hri_wcs.w->crpix[1] = hsize+1+qout.dy[0]+dy*dpp;
                vec2d tc = astro::regrid_drizzle(cutout, hri_wcs, wcs);
                tc[where(!is_finite(tc))] = 0;

                auto res = linfit(tflx2d, terr2d, tc[idf], 1.0);
                if (res.chi2 < bchi2 && res.params[0] > 0) {
                    bchi2 = res.chi2;
                    bamp = res.params[0];
                    berr = res.errors[0];
                    bbg = res.params[1];
                    bbge = res.errors[1];
                    bdx = dx;
                    bdy = dy;
                }

                if (verbose) progress(pg, 131);
            }

            if (verbose) {
                note("fit parameters (SNR):");
                note("  background: ", bbg, " (", abs(bbg/bbge), ")");
                note("  flux: ", bamp, " (", abs(bamp/berr), ")");
                note("  dx: ", bdx);
                note("  dy: ", bdy);
            }
        } else {
            bdx = 0.0;
            bdy = 0.0;

            vec2d tc = astro::regrid_drizzle(cutout, hri_wcs, wcs);
            tc[where(!is_finite(tc))] = 0;

            auto res = linfit(tflx2d, terr2d, tc[idf], 1.0);

            bchi2 = res.chi2;
            bamp = res.params[0];
            berr = res.errors[0];
            bbg = res.params[1];
            bbge = res.errors[1];
        }

        hri_wcs.w->crpix[0] = hsize+1+qout.dx[0]+bdx*dpp;
        hri_wcs.w->crpix[1] = hsize+1+qout.dy[0]+bdy*dpp;
        profile2d = astro::regrid_drizzle(cutout, hri_wcs, wcs);
        profile2d[where(!is_finite(profile2d))] = 0;
        model2d = 1e-17*bamp*profile2d;
        profile2d /= total(profile2d);

        // Find peak pixel of the profile
        vec1u mid = mult_ids(profile2d, max_id(profile2d));
        // Find radius of the profile
        double rad = get_radius(profile2d, mid[0], mid[1]);

        // Add fit params to header
        add_keywords = [=](fits::output_image& oimg) {
            oimg.write_keyword("SEEING", is_finite(seeing) ? seeing : 0.0, "[arcsec]");
            oimg.write_keyword("PWIDTH", is_finite(rad) ? rad : 0.0,       "[pixels]");
            oimg.write_keyword("PPEAKX", mid[1], "[pixels]");
            oimg.write_keyword("PPEAKY", mid[0], "[pixels]");
            oimg.write_keyword("SNR", abs(bamp/berr));
        };
    } else {
        // Fit profile
        struct base_t {
            vec1d x, y;
        } xx;

        xx.x = flatten(generate_img(flx2d.dims, [&](int_t,    int_t tx) { return tx - flx2d.dims[1]/2.0; }));
        xx.y = flatten(generate_img(flx2d.dims, [&](int_t ty, int_t)    { return ty - flx2d.dims[0]/2.0; }));

        vec1u idf = where(is_finite(flx2d) && is_finite(err2d) && err2d > 0);
        vec1d tflx2d = 1e17*flx2d[idf]; vec1d terr2d = 1e17*err2d[idf];
        xx.x = xx.x[idf]; xx.y = xx.y[idf];

        vec1d startp = {0.0, max(tflx2d), 0.0, 0.0, gwidth_default};
        auto model = [](const base_t& x, const vec1d& p) {
            return p[0] + gauss_integral_2d(x.x-0.5, x.x+0.5, x.y-0.5, x.y+0.5, p[2], p[3], p[4], abs(p[1]));
        };

        auto opts = mpfit_options(startp.size());
        opts.frozen[2] = fix_pos;
        opts.frozen[3] = fix_pos;
        opts.frozen[4] = fix_gwidth;
        auto res = mpfitfun(tflx2d, terr2d, xx, model, startp, opts);

        // print(res.params);
        // print(res.params/res.errors);
        if (verbose) {
            note("fit parameters (SNR):");
            vec1s pname = {"background", "flux", "dx", "dy", "width"};
            for (uint_t k : range(res.params)) {
                note("  ", pname[k], ": ", res.params[k], " (", abs(res.params[k]/res.errors[k]), ")");
            }
        }

        // Build model
        xx.x = flatten(generate_img(flx2d.dims, [&](int_t,    int_t tx) { return tx - flx2d.dims[1]/2.0; }));
        xx.y = flatten(generate_img(flx2d.dims, [&](int_t ty, int_t)    { return ty - flx2d.dims[0]/2.0; }));
        model2d = 1e-17*reform(model(xx, res.params), flx2d.dims);

        // Add fit params to header
        add_keywords = [=](fits::output_image& oimg) {
            oimg.write_keyword("PWIDTH", res.params[4], "[pixels]");
            oimg.write_keyword("PPEAKX", res.params[2], "[pixels]");
            oimg.write_keyword("PPEAKY", res.params[3], "[pixels]");
            oimg.write_keyword("SNR", abs(res.params[1]/res.errors[1]));
        };

        // Build normalized profile
        // Normalize model to unit flux + remove the fitted background
        res.params[0] = 0;
        res.params[1] = 1;
        profile2d = reform(model(xx, res.params), flx2d.dims);
        profile2d /= 2.0*dpi*sqr(res.params[4]);
    }

    // Save observed profile + error + residual
    fits::output_image foimg_obs(outfile_collapsed);
    foimg_obs.write_empty();

    auto write_wcs = [&](fits::output_image& oimg) {
        oimg.write_keyword("CRPIX1", crpix1);
        oimg.write_keyword("CRPIX2", crpix2);
        oimg.write_keyword("CRVAL1", crval1);
        oimg.write_keyword("CRVAL2", crval2);
        oimg.write_keyword("CTYPE1", ctype1);
        oimg.write_keyword("CTYPE2", ctype2);
        oimg.write_keyword("CUNIT1", cunit1);
        oimg.write_keyword("CUNIT2", cunit2);
        if (pos_use_cd) {
            oimg.write_keyword("CD1_1",  cd11);
            oimg.write_keyword("CD1_2",  cd12);
            oimg.write_keyword("CD2_1",  cd21);
            oimg.write_keyword("CD2_2",  cd22);
        } else {
            oimg.write_keyword("CDELT1", cdelt1);
            oimg.write_keyword("CDELT2", cdelt2);
        }

        add_keywords(oimg);
    };

    foimg_obs.reach_hdu(1);
    foimg_obs.write(flx2d);
    write_wcs(foimg_obs);

    foimg_obs.reach_hdu(2);
    foimg_obs.write(err2d);
    write_wcs(foimg_obs);

    foimg_obs.reach_hdu(3);
    foimg_obs.write(flx2d - model2d);
    write_wcs(foimg_obs);

    if (verbose) {
        vec1u idf = where(is_finite(flx2d) && is_finite(err2d) && err2d > 0);
        print("RMS error: ", stddev((flx2d[idf] - model2d[idf])/err2d[idf]));
    }

    // Save model
    fits::output_image foimg(outfile_model);
    foimg.write(profile2d);
    write_wcs(foimg);

    return 0;
}
