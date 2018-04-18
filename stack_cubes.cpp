#include <phypp.hpp>

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print("usage: stack_cubes files.sof [options]");
        return 1;
    }

    bool verbose = false;
    std::string stack_align = "nearest";
    std::string wave_align = "nearest";

    uint_t edge_padding = 1;

    bool do_sigma_clip = true;          // enable/disable sigma clipping of outliers
    double sigma_clip_threshold = 5.0;  // significance threshold for rejecting outliers
    double sigma_clip_error_threshold = 10.0;  // significance threshold for rejecting outliers
    uint_t sigma_clip_width = 1;        // width (in pixels) of the wavelength bin in which to define outliers

    bool bootstrap = false;             // enable/disable bootstrapping to measure uncertainties
    uint_t nbstrap = 200;               // number of bootstrap realizations
    uint_t tseed = 42;                  // random seed for bootstrap generation

    vec1s selected_targets;
    bool save_intermediate = false;

    read_args(argc-1, argv+1, arg_list(
        verbose, stack_align, wave_align,
        do_sigma_clip, sigma_clip_threshold, sigma_clip_width, edge_padding, save_intermediate,
        nbstrap, bootstrap, name(tseed, "seed"), name(selected_targets, "targets"),
        sigma_clip_error_threshold
    ));

    if (stack_align != "nearest") {
        error("unknown value for 'stack_align': ", stack_align);
        note("please choose one of: nearest");
        return 1;
    }

    if (wave_align != "nearest") {
        error("unknown value for 'wave_align': ", wave_align);
        note("please choose one of: nearest");
        return 1;
    }

    auto seed = make_seed(tseed);

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

    vec1s files;
    ascii::read_table(argv[1], ascii::find_skip(argv[1]), files);

    if (verbose) note("scanning files to build list of targets");

    bool no_error = false;
    vec1b found_targets(selected_targets.size());
    vec1s targets;
    for (auto& filename : files) {
        fits::input_image iimg(filename);

        if (iimg.hdu_count() == 25) no_error = true;

        for (uint_t i : range(0, 24)) {
            std::string target;
            iimg.reach_hdu(0);
            if (!iimg.read_keyword("ESO OCS ARM"+to_string(i+1)+" NAME", target)) {
                warning("in '", filename, "'");
                warning("missing keyword 'ESO OCS ARM"+to_string(i+1)+" NAME', is something wrong with this file?");
                continue;
            }

            // Check that HDU is not empty
            iimg.reach_hdu(no_error ? i+1 : 2*i+1);
            if (iimg.image_dims().empty()) continue;

            // Save target name
            if (!is_any_of(target, targets) &&
                (selected_targets.empty() || is_any_of(target, selected_targets))) {
                targets.push_back(target);
                if (!selected_targets.empty()) {
                    found_targets[where_first(selected_targets == target)] = true;
                }
            }
        }
    }

    if (!selected_targets.empty()) {
        for (uint_t i : range(selected_targets)) {
            if (!found_targets[i]) {
                warning("could not find any data for target '", selected_targets[i], "'");
            }
        }
    }

    if (no_error) {
        warning("the input files do not contain uncertainties");
        warning("the output uncertainties will only come from bootstrap");
    }

    if (verbose) note("found ", targets.size(), " targets among ", files.size(), " files");

    inplace_sort(targets);

    for (auto& target : targets) {
        if (verbose) note("stacking ", target, "...");

        if (verbose) note("  defining astrometry of stacked cube");

        double cd11 = dnan, cd21 = dnan, cd12 = dnan, cd22 = dnan;
        astro::wcs final_wcs;
        vec1d dx = replicate(dnan, files.size());
        vec1d dy = replicate(dnan, files.size());
        vec1d dl = replicate(dnan, files.size());
        vec1u thdu = replicate(0, files.size());
        vec1u orig_dims;

        for (uint_t f : range(files)) {
            std::string filename = files[f];
            fits::input_image iimg(filename);
            iimg.reach_hdu(0);

            // Find HDU of this target in this file
            for (uint_t i : range(0, 24)) {
                std::string t;
                if (!iimg.read_keyword("ESO OCS ARM"+to_string(i+1)+" NAME", t)) {
                    continue;
                }

                if (t == target) {
                    thdu[f] = (no_error ? i+1 : 2*i+1);
                    break;
                }
            }

            if (thdu[f] == 0) continue;

            // Check that this is a valid cube
            iimg.reach_hdu(thdu[f]);
            vec1u dims = iimg.image_dims();
            if (dims.size() != 3) {
                warning("in '", filename, "'");
                warning("data for target ", target, " is not a valid cube (dims: ", dims, ")");
                thdu[f] = 0;
                continue;
            }

            // Check that the orientation is the same as the first cube
            double tcd11 = dnan, tcd21 = dnan, tcd12 = dnan, tcd22 = dnan;
            if (!iimg.read_keyword("CD1_1", tcd11) || !iimg.read_keyword("CD2_1", tcd21) ||
                !iimg.read_keyword("CD1_2", tcd12) || !iimg.read_keyword("CD2_2", tcd22)) {
                warning("in '", filename, "'");
                warning("missing CDx_y keywords for target ", target);
                thdu[f] = 0;
                continue;
            }

            if (!is_finite(cd11)) {
                cd11 = tcd11; cd21 = tcd21; cd12 = tcd12; cd22 = tcd22;
            } else {
                if (3600*abs(cd11 - tcd11) > 0.01 || 3600*abs(cd21 - tcd21) > 0.01) {
                    // Not the same orientation, we don't support that right now
                    thdu[f] = 0;
                    continue;
                }
            }

            // Get pointing offset
            if (orig_dims.empty()) {
                orig_dims = dims;
                final_wcs = astro::wcs(iimg.read_header());
                dx[f] = 0;
                dy[f] = 0;
                dl[f] = 0;
            } else {
                astro::wcs w(iimg.read_header());

                double tra, tdec, tlam, tx, ty, ts;
                astro::xy2ad(w,         dims[2]/2+1, dims[1]/2+1, tra, tdec);
                astro::ad2xy(final_wcs, tra,         tdec,        tx,  ty);
                astro::x2w(w,         0, 1.0,  tlam);
                astro::w2x(final_wcs, 0, tlam, ts);

                dx[f] = tx - (orig_dims[2]/2+1);
                dy[f] = ty - (orig_dims[1]/2+1);
                dl[f] = ts - 1.0;
            }
        }

        vec1u idg = where(thdu > 0);
        if (idg.empty()) {
            warning("target ", target, " has no valid data to stack");
            continue;
        }

        // Define spatial extent of stacked cube
        double yshift = abs(floor(min(dy[idg])));
        double xshift = abs(floor(min(dx[idg])));
        uint_t nfile = count(thdu > 0);
        uint_t nlam = orig_dims[0];
        uint_t ny = orig_dims[1]+ceil(max(dy[idg])+yshift);
        uint_t nx = orig_dims[2]+ceil(max(dx[idg])+xshift);
        dy += yshift;
        dx += xshift;

        // Define cubes
        vec4d sflx3d = replicate(dnan, nfile, nlam, ny, nx);
        vec4d serr3d = replicate(dnan, nfile, nlam, ny, nx);
        vec4d swht3d = replicate(0.0, nfile, nlam, ny, nx);

        if (verbose) note("  stacking ", nfile, " exposures");
        if (verbose) note("  stacked cube will be ", nlam, "x", ny, "x", nx);
        if (verbose) note("  regridding");

        // Shift individual exposures on the same grid
        {
            uint_t f = 0;
            for (uint_t tf : range(files)) {
                if (thdu[tf] == 0) continue;

                std::string filename = files[tf];
                fits::input_image iimg(filename);
                iimg.reach_hdu(thdu[tf]);

                vec3d tflx3d;
                iimg.read(tflx3d);

                // Flag out edges if requested
                if (edge_padding > 0) {
                    tflx3d(_,_-(edge_padding-1),_) = dnan;
                    tflx3d(_,(tflx3d.dims[1]-edge_padding)-_,_) = dnan;
                    tflx3d(_,_,_-(edge_padding-1)) = dnan;
                    tflx3d(_,_,(tflx3d.dims[2]-edge_padding)-_) = dnan;
                }

                // Flag out OH residuals
                auto w = astro::wcs(iimg.read_header());
                vec1d lam = astro::build_axis(w, 0, astro::axis_unit::wave_um);
                double cdelt = abs(lam[1]-lam[0]);

                vec1b flagged(lam.dims);
                for (uint_t l : range(lines)) {
                    flagged = flagged || (lam >= lines[l]-0.0006 && lam <= lines[l]+0.0006);
                }
                for (uint_t l : range(bands_low)) {
                    flagged = flagged || (lam >= bands_low[l]-cdelt && lam <= bands_up[l]+cdelt);
                }

                tflx3d(where(flagged),_,_) = dnan;

                // Stack
                if (no_error) {
                    for (uint_t l0 : range(nlam)) {
                        if (l0 < dl[tf] || l0 - dl[tf] >= nlam) continue;

                        uint_t l1 = l0 - dl[tf];
                        uint_t x0 = round(dx[tf]), x1 = x0 + tflx3d.dims[2]-1;
                        uint_t y0 = round(dy[tf]), y1 = y0 + tflx3d.dims[1]-1;

                        sflx3d(f, l0, y0-_-y1, x0-_-x1) = tflx3d(l1,_,_);
                        swht3d(f, l0, y0-_-y1, x0-_-x1) = 1.0;
                    }
                } else {
                    iimg.reach_hdu(thdu[tf]+1);
                    vec3d terr3d;
                    iimg.read(terr3d);

                    terr3d(where(flagged),_,_) = dnan;

                    for (uint_t l0 : range(nlam)) {
                        if (l0 < dl[tf] || l0 - dl[tf] >= nlam) continue;

                        uint_t l1 = l0 - dl[tf];
                        uint_t x0 = round(dx[tf]), x1 = x0 + tflx3d.dims[2]-1;
                        uint_t y0 = round(dy[tf]), y1 = y0 + tflx3d.dims[1]-1;

                        sflx3d(f, l0, y0-_-y1, x0-_-x1) = tflx3d(l1,_,_);
                        serr3d(f, l0, y0-_-y1, x0-_-x1) = terr3d(l1,_,_);
                        swht3d(f, l0, y0-_-y1, x0-_-x1) = 1/sqr(terr3d(l1,_,_));
                    }
                }

                ++f;
            }
        }

        // Apply sigma clipping (if asked)
        vec4b crej(sflx3d.dims);
        if (do_sigma_clip) {
            // Compute median uncertainty for each exposure
            vec1d mederr(nfile);
            if (!no_error) {
                if (verbose) note("  determining median errors");
                for (uint_t f : range(nfile)) {
                    vec3d tmp = serr3d(f,_,_,_);
                    mederr[f] = median(tmp[where(tmp > 0 && is_finite(tmp))]);
                }
            }

            if (verbose) note("  sigma clipping");

            uint_t d1 = (sigma_clip_width-1)/2;
            uint_t d2 = sigma_clip_width-1 - d1;

            for (uint_t l : range(nlam))
            for (uint_t y : range(ny))
            for (uint_t x : range(nx)) {
                // Define wavelength region
                uint_t l0 = (l > d1 ?      l - d1 : 0);
                uint_t l1 = (l < nlam-d2 ? l + d2 : nlam-1);

                // First compute the weighted median (which we assume is unbiased)
                vec2d med(nfile, l1-l0+1);
                vec1b erej(nfile);
                for (uint_t ll : range(med.dims[1])) {
                    // First identify outliers in the weight or error spectrum
                    vec1d twei = swht3d(_,l0+ll,y,x);
                    if (!no_error) {
                        // Normalize each error value by the median of its exposure,
                        // this accounts for varying integration times and data quality
                        vec1d terr = mederr/serr3d(_,l0+ll,y,x);
                        // Then compute the median factor
                        vec1u ttid = where(terr > 0 && is_finite(terr));
                        if (ttid.size() < 5) {
                            med(_,ll) = dnan;
                            continue;
                        }

                        double tmed = median(terr[ttid]);
                        vec1d tres = abs(terr - tmed);
                        double trms = 1.48*median(tres[ttid]);
                        // And isolate data points which deviate significantly
                        erej = tres > sigma_clip_error_threshold*trms;
                        // Put their weight to zero temporarily
                        twei[where(erej)] = 0;
                    }

                    // Then compute the median
                    med(_,ll) = weighted_median(sflx3d(_,l0+ll,y,x), twei);
                }

                // Compute the absolute residuals
                vec2d res = abs(sflx3d(_,l0-_-l1,y,x) - med);

                // Compute the RMS of these using the MAD (which we assume is unbiased)
                double rms = 1.48*median(res);

                // Select significant outliers in this wavelength element
                crej(_,l,y,x) = res(_,d1) > sigma_clip_threshold*rms || erej;
            }

            if (verbose) {
                uint_t nfinite = count(is_finite(sflx3d));
                uint_t cnt = count(crej);
                note("  ", cnt, "/", nfinite, " elements sigma clipped (",
                    round(10.0*100.0*cnt/float(nfinite))/10.0, "%)");
            }
        }

        if (save_intermediate) {
            fits::output_image oimg(target+"/cube_exposures.fits");
            oimg.write_empty();

            oimg.reach_hdu(1);
            oimg.write(sflx3d);
            oimg.reach_hdu(2);
            oimg.write(serr3d);
            oimg.reach_hdu(3);
            oimg.write(crej);
        }

        // Mask invalid data points
        if (verbose) note("  masking");
        vec1u idb = where(!is_finite(sflx3d) || !is_finite(serr3d) || !is_finite(swht3d) || crej);
        sflx3d[idb] = 0; serr3d[idb] = 0; swht3d[idb] = 0;
        if (verbose) note("  ", round(10.0*100.0*idb.size()/float(sflx3d.size()))/10.0, "% of data masked");

        // Stack them
        auto stack_spectrum = [&](const vec4d& flx4d, const vec4d& err4d, const vec4d& wht4d,
            vec3d& flx, vec3d& err, vec3d& errb) {

            // vec3d wht = partial_total(0, wht4d);
            // flx = partial_total(0, flx4d*wht4d)/wht;
            // if (!no_error) {
            //     err = sqrt(partial_total(0, sqr(err4d*wht4d)))/wht;
            // }

            flx.resize(flx4d.dims[1], flx4d.dims[2], flx4d.dims[3]);
            errb.resize(flx4d.dims[1], flx4d.dims[2], flx4d.dims[3]);
            err.resize(flx4d.dims[1], flx4d.dims[2], flx4d.dims[3]);

            for (uint_t il : range(flx4d.dims[1]))
            for (uint_t iy : range(flx4d.dims[2]))
            for (uint_t ix : range(flx4d.dims[3])) {
                double wht = 0.0;
                for (uint_t i : range(flx4d.dims[0])) {
                    flx.safe(il,iy,ix) += flx4d.safe(i,il,iy,ix)*wht4d.safe(i,il,iy,ix);
                    if (!no_error) {
                        err.safe(il,iy,ix) += sqr(err4d.safe(i,il,iy,ix)*wht4d.safe(i,il,iy,ix));
                    }
                    wht += wht4d.safe(i,il,iy,ix);
                }

                flx.safe(il,iy,ix) /= wht;
                if (!no_error) {
                   err.safe(il,iy,ix) = sqrt(err.safe(il,iy,ix))/wht;
                }

                uint_t npt = 0;
                wht = 0.0;
                for (uint_t i : range(flx4d.dims[0])) {
                    if (wht4d.safe(i,il,iy,ix) > 0) {
                        ++npt;
                        wht += wht4d.safe(i,il,iy,ix);
                        errb.safe(il,iy,ix) += sqr((flx4d.safe(i,il,iy,ix) - flx.safe(il,iy,ix))*wht4d.safe(i,il,iy,ix));
                    }
                }

                errb.safe(il,iy,ix) = sqrt(errb.safe(il,iy,ix)*(npt/(npt-1.0)))/wht;
                if (no_error) {
                    // Use bootstrap error as the "formal error" because we have none...
                    err.safe(il,iy,ix) = errb.safe(il,iy,ix);
                }
            }
        };

        if (verbose) note("  stacking");
        vec3d flx3d, err3d, errb3d;
        stack_spectrum(sflx3d, serr3d, swht3d, flx3d, err3d, errb3d);

        // Compute bootstrapping
        // if (nbstrap > 0) {
        //     if (verbose) note("  bootstrapping");
        //     auto pg = progress_start(nbstrap);
        //     vec4d bres(flx3d.dims[0], flx3d.dims[1], flx3d.dims[2], nbstrap);
        //     for (uint_t i : range(nbstrap)) {
        //         // vec1u ids = shuffle(seed, uindgen(nfile))[_-(nfile/2)]; // without replacement
        //         vec1u ids = randomi(seed, 0, nfile-1, nfile/2); // with replacement

        //         // Compute flux[half] and error[half]
        //         vec3d bflx3d(flx3d.dims);
        //         vec3d berr3d;
        //         if (!no_error) berr3d.resize(flx3d.dims);

        //         stack_spectrum(sflx3d(ids,_,_,_), serr3d(ids,_,_,_), swht3d(ids,_,_,_), bflx3d, berr3d);

        //         if (!no_error) {
        //             // Subtract flux[total]
        //             bflx3d -= flx3d;
        //             // Weight by the ratio err[total]/err[half] to account for missing elements
        //             bflx3d *= err3d/berr3d;
        //         }

        //         // Save
        //         bres(_,_,_,i) = bflx3d;

        //         if (verbose) progress(pg);
        //     }

        //     if (!no_error) {
        //         // Get uncertainty from the RMS of the resulting bootstrap
        //         errb3d = partial_rms(3, bres);
        //     } else {
        //         // Get uncertainty from the stddev of the resulting bootstrap
        //         errb3d = partial_stddev(3, bres);
        //         // Use this error as the "formal error" because we have none...
        //         err3d = errb3d;
        //     }
        // }

        // Save to file
        if (verbose) note("  saving FITS file '"+target+"/cube.fits'");
        file::mkdir(target);
        fits::output_image oimg(target+"/cube.fits");
        oimg.write_empty();

        auto write_wcs = [&]() {
            double cra, cdec, clam, clam1;
            astro::xy2ad(final_wcs, orig_dims[1]/2+1, orig_dims[2]/2+1, cra, cdec);
            astro::x2w(final_wcs, 0, 1, clam, astro::axis_unit::wave_um);
            astro::x2w(final_wcs, 0, 2, clam1, astro::axis_unit::wave_um);

            oimg.write_keyword("CTYPE1", "RA---TAN");
            oimg.write_keyword("CUNIT1", "deg");
            oimg.write_keyword("CRPIX1", xshift+orig_dims[2]/2+1);
            oimg.write_keyword("CRVAL1", cra);
            oimg.write_keyword("CTYPE2", "DEC--TAN");
            oimg.write_keyword("CUNIT2", "deg");
            oimg.write_keyword("CRPIX2", yshift+orig_dims[1]/2+1);
            oimg.write_keyword("CRVAL2", cdec);
            oimg.write_keyword("CD1_1", cd11);
            oimg.write_keyword("CD2_1", cd21);
            oimg.write_keyword("CD1_2", cd12);
            oimg.write_keyword("CD2_2", cd22);

            oimg.write_keyword("CTYPE3", "WAVE");
            oimg.write_keyword("CUNIT3", "um");
            oimg.write_keyword("CRPIX3", 1.0);
            oimg.write_keyword("CRVAL3", clam);
            oimg.write_keyword("CD1_3", 0.0);
            oimg.write_keyword("CD2_3", 0.0);
            oimg.write_keyword("CD3_1", 0.0);
            oimg.write_keyword("CD3_2", 0.0);
            oimg.write_keyword("CD3_3", clam1 - clam);
        };

        oimg.reach_hdu(1);
        oimg.write(flx3d);
        write_wcs();
        oimg.reach_hdu(2);
        oimg.write(err3d);
        write_wcs();
        oimg.reach_hdu(3);
        oimg.write(errb3d);
        write_wcs();
    }

    return 0;
}
