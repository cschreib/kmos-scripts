#include <phypp.hpp>

// Local functions (defined at the end of the file)
vec2u segment(vec2u map, uint_t first_id, uint_t& nsrc);
vec2u grow_within(vec2u map, vec2b mask);
void print_help();

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    std::string expmap;
    double minexp = 1;
    double spatial_smooth = 1.2; // radius of the smoothing kernel
    uint_t spectral_bin = 1; // number of spectral pixels to sum
    bool save_cubes = false;
    double snr_det = 5.0; // SNR threshold for detections
    double snr_source = 3.0; // lower SNR threshold for the extents of the source
    bool subtract_continuum = false;
    uint_t continuum_width = 100; // in pixels of the original scale
    bool verbose = false;
    double error_scale = 1.0;
    std::string outdir;

    read_args(argc-1, argv+1, arg_list(expmap, minexp, spatial_smooth, save_cubes,
        snr_det, snr_source, subtract_continuum, continuum_width, verbose, spectral_bin,
        error_scale, outdir));

    if (!outdir.empty()) {
        outdir = file::directorize(outdir);
        file::mkdir(outdir);
    }

    // Read cube
    std::string infile = argv[1];
    vec3d cube;
    fits::input_image fimg(infile);
    fimg.read(cube);

    // Build wavelength axis
    uint_t nlam = cube.dims[0];
    double cdelt = 1, crpix = 1, crval = 1;
    if (!fimg.read_keyword("CDELT3", cdelt) ||
        !fimg.read_keyword("CRPIX3", crpix) ||
        !fimg.read_keyword("CRVAL3", crval)) {
        error("could not read WCS information for wavelength axis");
        return 1;
    }

    vec1d lam = crval + cdelt*(findgen(nlam) + (1 - crpix));

    // Bin spectral pixels and start building exposure map
    vec3d exposure;
    if (spectral_bin > 1) {
        if (verbose) note("rebinning spectral axis...");

        // Rebin SNR array
        nlam = floor(nlam/double(spectral_bin));

        vec3d ocube = cube;
        cube.resize(nlam, cube.dims[1], cube.dims[2]);
        cube[_] = 0;

        // Keep track of how many pixels were used in the sum
        vec3u cnt(cube.dims);

        // Compute the sum of the pixels in each bin
        for (uint_t l : range(nlam)) {
            for (uint_t i : range(spectral_bin)) {
                vec2d tmp = ocube(l*spectral_bin + i,_,_);
                vec1u idg = where(is_finite(tmp));
                cube(l,_,_)[idg] += tmp[idg];
                cnt(l,_,_)[idg] += 1;
            }
        }

        cube /= cnt;
        exposure = 1/sqrt(cnt);

        // Update wavelength array
        crpix -= (spectral_bin-1)*cdelt/2.0;
        cdelt *= spectral_bin;
        lam = crval + cdelt*(findgen(nlam) + (1 - crpix));
    } else {
        exposure = replicate(1.0, cube.dims);
    }

    // Read exposure map if provided
    // vec2d exposure;
    vec1u badexp;
    // vec2d rms_renorm;
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
    vec3d err = replicate(dnan, cube.dims);

    if (verbose) note("filtering wavelength slices...");
    for (uint_t l : range(nlam)) {
        // 1) Flag baddly covered areas with exposure map (optional)
        if (!badexp.empty()) {
            cube(l,_,_)[badexp] = dnan;
        }

        // Find good pixels remaining
        vec2d tmp = cube(l,_,_)/exposure(l,_,_);
        vec1u idg = where(is_finite(tmp));
        if (!idg.empty()) {
            // 2) Compute pixel fluctuations of exposure-renormalized fluxes
            // and estimate error cube by this value and the exposure

            // Fluctuations from median absolute deviation
            // double img_rms = 1.48*mad(tmp[idg]);
            // -> not enough pixels in the IFU to use this one reliably, it tends
            // to underestiamte the actual flux fluctuation in some situation

            // Fluctuations from RMS of negative pixels in median subtracted map
            tmp -= median(tmp[idg]);
            idg = where(is_finite(tmp) && tmp < 0);
            double img_rms = rms(tmp[idg]);

            err(l,_,_) = (img_rms*error_scale)*exposure(l,_,_);

            if (spatial_smooth > 0) {
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
                err(l,_,_) *= sqrt(total(sqr(kernel)));
            }
        } else {
            cube(l,_,_) = fnan;
            err(l,_,_) = fnan;
        }
    }

    // If asked, estimate and subtract the continuum emission
    if (subtract_continuum) {
        if (verbose) note("estimate continuum...");
        // Build a rebinned cube
        uint_t cbin = floor(continuum_width/double(spectral_bin));

        if (cbin <= 1) {
            error("continuum would be estimated over only one rebinned spectral element");
            error("this would remove all the flux from the cube");
            note("please use a higher value of 'continuum_width' (was ", continuum_width, ")");
            note("or disable continuum subtraction (subtract_continuum=0)");
            return 1;
        }

        vec3d ccube;
        uint_t cnlam = floor(nlam/double(cbin));
        ccube.resize(cnlam, cube.dims[1], cube.dims[2]);
        ccube[_] = 0;

        // Keep track of weights of pixels that were used in the sum
        vec3d cnt(ccube.dims);

        // Compute the weighted sum of the pixels in each bin
        for (uint_t l : range(cnlam)) {
            for (uint_t i : range(cbin)) {
                vec2d w = 1.0/err(l*spectral_bin + i,_,_);
                vec2d tmp = w*cube(l*spectral_bin + i,_,_);
                vec1u idg = where(is_finite(tmp));
                ccube(l,_,_)[idg] += tmp[idg];
                cnt(l,_,_)[idg] += w[idg];
            }
        }

        ccube /= cnt;

        // Now interpolate it back to the original pixel grid and subtract it
        if (verbose) note("subtract continuum...");
        double ccrpix = crpix - (cbin-1)*cdelt/2.0;
        double ccdelt = cdelt*cbin;
        vec1d clam = crval + ccdelt*(findgen(cnlam) + (1 - ccrpix));
        for (uint_t y : range(cube.dims[1]))
        for (uint_t x : range(cube.dims[2])) {
            cube(_,y,x) -= interpolate(ccube(_,y,x), clam, lam);
        }
    }

    // If asked, save the processed cubes to disk
    std::string ofilebase = outdir+file::remove_extension(file::get_basename(infile));
    if (save_cubes) {
        fits::output_image oimg(ofilebase+"_filt.fits");

        // Empty primary array (KMOS convention)
        oimg.write(vec2d(0,0));

        // Flux
        oimg.reach_hdu(1);
        oimg.write(cube);
        oimg.write_header(fimg.read_header());
        if (spectral_bin > 1) {
            oimg.write_keyword("CDELT3", cdelt);
            oimg.write_keyword("CRPIX3", crpix);
            oimg.write_keyword("CD3_3", cdelt);
        }

        // Uncertainty
        oimg.reach_hdu(2);
        oimg.write(err);
        oimg.write_header(fimg.read_header());
        if (spectral_bin > 1) {
            oimg.write_keyword("CDELT3", cdelt);
            oimg.write_keyword("CRPIX3", crpix);
            oimg.write_keyword("CD3_3", cdelt);
        }

        // S/N
        oimg.reach_hdu(3);
        oimg.write(cube/err);
        oimg.write_header(fimg.read_header());
        if (spectral_bin > 1) {
            oimg.write_keyword("CDELT3", cdelt);
            oimg.write_keyword("CRPIX3", crpix);
            oimg.write_keyword("CD3_3", cdelt);
        }
    }

    // Now we are ready to try detecting stuff above the noise
    // Build a detection list and segmentation map at each wavelength slice
    // independently

    vec3u seg(cube.dims);
    uint_t nsrc_tot = 0;

    vec1u id;
    vec1d lambda;
    vec1u lpix;
    vec1u npix;
    vec1d flux;
    vec1d flux_err;

    if (verbose) note("finding and segmenting detections...");
    for (uint_t l : range(nlam)) {
        // Flag pixels above the S/N threshold
        vec2d snr = cube(l,_,_)/err(l,_,_);
        vec2u det = vec2u(snr > snr_det);

        // Segment them into individual sources
        uint_t nsrc;
        uint_t first_id = nsrc_tot + 1;
        det = segment(det, first_id, nsrc);
        nsrc_tot += nsrc;

        if (nsrc != 0) {
            // We have a new source(s)!
            // Make them grow a little toward lower S/N since they must be real
            if (snr_source < snr_det) {
                det = grow_within(det, snr > snr_source);
            }

            // Save that into the segmentation map
            seg(l,_,_) = det;

            // For each source, save its identifiers and some information
            for (uint_t i : range(nsrc)) {
                vec1u idd = where(det == i+1);

                id.push_back(first_id + i);
                npix.push_back(idd.size());

                lpix.push_back(l+1);
                lambda.push_back(lam[l]);

                flux.push_back(total(cube(l,_,_)[idd]));
                flux_err.push_back(sqrt(total(sqr(err(l,_,_)[idd])))); // just indicative
            }
        }
    }

    if (verbose) note("found ", nsrc_tot, " sources in the cube");

    fits::output_image fseg(ofilebase+"_seg.fits");
    fseg.write(seg);
    fseg.write_header(fimg.read_header());
    if (spectral_bin > 1) {
        fseg.write_keyword("CDELT3", cdelt);
        fseg.write_keyword("CRPIX3", crpix);
        fseg.write_keyword("CD3_3", cdelt);
    }

    fits::write_table(ofilebase+"_cat.fits", ftable(id, npix, lambda, lpix, flux, flux_err));

    return 0;
}


// Function to segment a binary map into multiple components.
// Does no de-blending. The returned map contains IDs between 'first_id' and
// 'first_id + nsrc' (where 'nsrc' is the number of identified segments, which is
// provided in output). Values of 0 in the input binary map are also 0 in the
// segmentation map.
vec2u segment(vec2u map, uint_t first_id, uint_t& nsrc) {
    vec2u smap(map.dims);
    nsrc = 0;
    uint_t id = first_id;

    std::vector<uint_t> oy, ox;
    for (uint_t y : range(map.dims[0]))
    for (uint_t x : range(map.dims[1])) {
        if (map.safe(y,x) == 0) continue;

        // Found a guy
        // Use an A* - like algorithm to navigate around and
        // figure out its extents
        oy.clear(); ox.clear();

        auto process_point = [&ox,&oy,&map,&smap](uint_t ty, uint_t tx, uint_t tid) {
            map.safe(ty,tx) = 0;
            smap.safe(ty,tx) = tid;

            auto check_add = [&ox,&oy,&map](uint_t tty, uint_t ttx) {
                if (map.safe(tty,ttx) != 0) {
                    oy.push_back(tty);
                    ox.push_back(ttx);
                }
            };

            if (ty != 0)             check_add(ty-1,tx);
            if (ty != map.dims[0]-1) check_add(ty+1,tx);
            if (tx != 0)             check_add(ty,tx-1);
            if (tx != map.dims[1]-1) check_add(ty,tx+1);
        };

        process_point(y, x, id);

        while (!ox.empty()) {
            uint_t ty = oy.back(); oy.pop_back();
            uint_t tx = ox.back(); ox.pop_back();
            process_point(ty, tx, id);
        }

        ++id;
    }

    nsrc = id - first_id;

    return smap;
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

void print_help() {
    using namespace format;

    print("cdetect v1.0");
    print("usage: cdetect <kmos_cube.fits> [options]");
    print("");
    print("Main parameter:");
    paragraph("'kmos_cube.fits' should be a cube created by the KMOS pipeline, with at "
        "least 2 extensions: the first is empty (KMOS convention), and the second "
        "contains the flux. Even if a third extension is present and contains the "
        "uncertainty, it is not used (the uncertainty is derived from the flux cube "
        "itself.");
    print("Available options (in order of importance):");
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
    bullet("subtract_continuum", "Set this flag if you want the program to estimate the "
        "continuum emission of your target(s) from the spectrum. The continuum will be "
        "estimated from the average flux over multiple large spectral windows; large "
        "enough to wash out any line emission (the size of this window can be controlled "
        "with the parameter 'continuum_width'). Default is not to subtract the continuum, "
        "since it may be what you are after.");
    bullet("continuum_width=...", "Must be a integer. It defines the number of spectral "
        "elements that will be averaged to compute the continuum emission. Note that "
        "this number refers to the spectral elements *before* spectral binning has "
        "occured, and therefore does not depend on the chosen value of of 'spectral_bin'. "
        "Default is 100 pixels.");
    bullet("error_scale=...", "Must be a number. It can be used to scale up or down the "
        "estimated uncertainties globally. It usually happens that uncertainties are "
        "underestimated, and the program produces many false positives. In this case "
        "you can raise this value, typically to 1.3. Default is 1, so no rescaling.");
    bullet("snr_det=...", "Must be a number. It defines the minimum S/N ratio to consider "
        "for a detection. Default is 5.");
    bullet("snr_source=...", "Must be a number. It defines the minimum S/N ratio to consider "
        "to define the spatial extents of a detection. Default is 3.");
    bullet("save_cubes", "Set this flag to write to disk the filtered flux and uncertainty "
        "cubes that are used to perform the detection. They will be saved in *_filt.fits.");
    bullet("outdir", "Name of the directory into which the output files should be created. "
        "Default is the current directory.");
    bullet("verbose", "Set this flag to print the progress of the detection process in "
        "the terminal. Can be useful if something goes wrong, or just to understand what "
        "is going on.");
}
