#include <phypp.hpp>

void print_help();

int phypp_main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    vec1s suffix;
    std::string outdir;
    bool background = false;
    std::string target;
    bool verbose = false;

    bool do_sigma_clip = true;          // enable/disable sigma clipping of outliers
    double sigma_clip_threshold = 5.0;  // significance threshold for rejecting outliers
    uint_t sigma_clip_width = 1;        // width (in pixels) of the wavelength bin in which to define outliers

    read_args(argc-2, argv+2, arg_list(background, suffix, outdir, target, verbose,
        do_sigma_clip, sigma_clip_threshold, sigma_clip_width));

    if (!outdir.empty()) {
        outdir = file::directorize(outdir);
        file::mkdir(outdir);
    }

    // Read models
    if (verbose) {
        note("reading models");
    }

    vec3d models;
    astro::wcs model_wcs;
    uint_t nmodel; {
        fits::input_image fmodel(argv[2]);

        model_wcs = astro::wcs(fmodel.read_header());

        if (fmodel.axis_count() == 2) {
            nmodel = 1;
            vec2d tmodel;
            fmodel.read(tmodel);
            models = reform(tmodel, 1, tmodel.dims);
        } else {
            fmodel.read(models);
            nmodel = models.dims[0];
        }

        // Add constant background if requested
        if (background) {
            append<0>(models, replicate(1.0, 1, models.dims[1], models.dims[2]));
            if (nmodel == 1) {
                suffix = {"", "background"};
            } else if (suffix.size() == nmodel) {
                suffix.push_back("background");
            }

            ++nmodel;
        }

        if (suffix.size() != nmodel && nmodel != 1) {
            error("please provide as many suffixes as there are models (", nmodel,
                ") in suffix=[...]");
            return 1;
        }
    }

    // Adjust suffixes to include "_"
    if (nmodel == 1 && suffix.empty()) {
        suffix = {""};
    }

    for (std::string& s : suffix) {
        if (!s.empty() && s[0] != '_') {
            s = "_"+s;
        }
    }

    // Read SOF
    vec1s files;
    ascii::read_table(argv[1], ascii::find_skip(argv[1]), files);
    files = file::get_directory(argv[1])+files;

    // Read and fit spectra
    if (verbose) {
        note("reading and fitting cubes");
    }

    vec1d lam;
    vec3d cflx, cerr;
    vec1s used_files;
    for (uint_t f : range(files)) {
        fits::input_image iimg(files[f]);

        if (iimg.hdu_count() == 25) {
            error("in '", files[f], "'");
            error("this file does not contain uncertainties, which are required to do the fit");
            return 1;
        }

        // Identify HDU of this target
        uint_t ihdu = npos;
        iimg.reach_hdu(0);
        for (uint_t i : range(0, 24)) {
            std::string ttarget;
            iimg.reach_hdu(0);
            if (!iimg.read_keyword("ESO OCS ARM"+to_string(i+1)+" NAME", ttarget)) {
                warning("in '", files[f], "'");
                warning("missing keyword 'ESO OCS ARM"+to_string(i+1)+" NAME', is something wrong with this file?");
                continue;
            }

            // Check that HDU is not empty
            iimg.reach_hdu(2*i+1);
            if (iimg.image_dims().empty()) continue;

            // Save target name
            if (ttarget == target) {
                ihdu = 2*i+1;
                break;
            }
        }

        if (ihdu == npos) continue;

        used_files.push_back(files[f]);

        // Read cube
        vec3d flx, err;
        iimg.reach_hdu(ihdu);
        iimg.read(flx);

        astro::wcs cube_wcs(iimg.read_header());
        if (lam.empty()) {
            lam = astro::build_axis(cube_wcs, 0, astro::axis_unit::wave_um);
        }

        iimg.reach_hdu(ihdu+1);
        iimg.read(err);

        // Regrid models
        vec2d gmodels(nmodel, flx.dims[1]*flx.dims[2]);
        for (uint_t m : range(nmodel)) {
            vec2d tm = models(m,_,_);
            tm = astro::regrid_drizzle(tm, model_wcs, cube_wcs);

            // Make sure each model has no bad pixel and has unit integral
            tm[where(!is_finite(tm))] = 0;
            tm /= total(tm);

            gmodels(m,_) = flatten(tm);
        }

        // Do the fit wavelength by wavelength
        vec2d flx1d(nmodel, flx.dims[0]);
        vec2d err1d(nmodel, flx.dims[0]);
        for (uint_t l : range(flx.dims[0])) {
            vec1u idg = where(is_finite(flx(l,_,_)) && is_finite(err(l,_,_)) && err(l,_,_) > 0);
            auto res = linfit_pack(flx(l,_,_)[idg], err(l,_,_)[idg], gmodels(_,idg));
            flx1d(_,l) = res.params;
            err1d(_,l) = res.errors;
        }

        // Save fits
        append<0>(cflx, reform(flx1d, 1, flx1d.dims));
        append<0>(cerr, reform(err1d, 1, err1d.dims));
    }

    if (verbose) {
        note("combining ", cflx.dims[0], " exposures");
    }

    // Inverse variance weighting
    vec3d cwei = 1/sqr(cerr);

    // Apply sigma clipping (if asked)
    vec3b crej(cflx.dims);
    if (do_sigma_clip) {
        uint_t d1 = (sigma_clip_width-1)/2;
        uint_t d2 = sigma_clip_width-1 - d1;

        for (uint_t m : range(nmodel))
        for (uint_t l : range(cflx.dims[2])) {
            // Define wavelength region
            uint_t l0 = (l > d1 ?              l - d1 : 0);
            uint_t l1 = (l < cflx.dims[2]-d2 ? l + d2 : cflx.dims[2]-1);

            // First compute the weighted median (which we assume is unbiased)
            vec2d med(cflx.dims[0], l1-l0+1);
            for (uint_t ll : range(med.dims[1])) {
                med(_,ll) = weighted_median(cflx(_,m,l0+ll), cwei(_,m,l0+ll));
            }

            // Compute the absolute residuals
            vec2d res = abs(cflx(_,m,l0-_-l1) - med);

            // Compute the RMS of these using the MAD (which we assume is unbiased)
            double rms = 1.48*median(res);

            // Select significant outliers in this wavelength element
            crej(_,m,l) = res(_,d1) > sigma_clip_threshold*rms;
        }

        if (verbose) {
            uint_t nfinite = count(is_finite(cflx));
            uint_t cnt = count(crej);
            note(cnt, "/", nfinite, " elements sigma clipped (",
                round(10.0*100.0*cnt/float(nfinite))/10.0, "%)");
        }
    }

    // Down-weight bad pixels
    {
        vec1u idb = where(!is_finite(cflx) || !is_finite(cerr) || !is_finite(cwei) || crej);
        cwei[idb] = 0; cflx[idb] = 0; cerr[idb] = 0;
    }

    // Write individual spectra to disk
    if (verbose) {
        note("writing individual spectra");
    }

    std::string filebase = outdir+target;
    for (uint_t m : range(nmodel)) {
        fits::output_image fspec(filebase+suffix[m]+"_spec_indiv.fits");
        fspec.write_empty(); // empty primary extension, KMOS convention

        auto write_wcs = [&]() {
            fspec.write_keyword("CTYPE1", "WAVE");
            fspec.write_keyword("CUNIT1", "um");
            fspec.write_keyword("CRPIX1", 1.0);
            fspec.write_keyword("CRVAL1", lam[0]);
            fspec.write_keyword("CDELT1", lam[1]-lam[0]);
        };

        fspec.reach_hdu(1);
        fspec.write(cflx(_,m,_));
        write_wcs();

        fspec.reach_hdu(2);
        fspec.write(cerr(_,m,_));
        write_wcs();

        fspec.reach_hdu(3);
        fspec.write(cwei(_,m,_));
        write_wcs();
    }

    // Stack spectra
    if (verbose) {
        note("stacking spectra");
    }
    vec2d flx1d(nmodel, cflx.dims[2]);
    vec2d err1d(nmodel, cflx.dims[2]);
    vec2d err1db(nmodel, cflx.dims[2]);
    for (uint_t m : range(nmodel)) {
        vec1d wei = partial_total(0, cwei(_,m,_));
        flx1d(m,_) = partial_total(0, cflx(_,m,_)*cwei(_,m,_))/wei;
        err1d(m,_) = sqrt(partial_total(0, sqr(cerr(_,m,_)*cwei(_,m,_))))/wei;

        vec2d tmp = (cflx(_,m,_) - replicate(flx1d(m,_), cflx.dims[0]))*(replicate(err1d(m,_), cflx.dims[0])/cerr(_,m,_));
        vec2d twei = cwei(_,m,_);
        vec1u idb = where(!is_finite(tmp));
        twei[idb] = 0; tmp[idb] = 0;
        err1db(m,_) = sqrt(partial_total(0, twei*sqr(tmp))/partial_total(0, twei));
    }

    // Flag OH residuals and low transmission
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

    vec1b flagged(lam.dims);
    for (uint_t l : range(lines)) {
        flagged = flagged || (lam >= lines[l]-0.0006 && lam <= lines[l]+0.0006);
    }
    double cdelt = lam[1]-lam[0];
    for (uint_t l : range(bands_low)) {
        flagged = flagged || (lam >= bands_low[l]-cdelt && lam <= bands_up[l]+cdelt);
    }

    vec1u idflag = where(flagged);
    for (uint_t m : range(nmodel)) {
        flx1d(m,idflag) = dnan;
        err1d(m,idflag) = dnan;
        err1db(m,idflag) = dnan;
    }

    // Write spectra to disk
    if (verbose) {
        note("writing spectra");
    }

    for (uint_t m : range(nmodel)) {
        fits::output_image fspec(filebase+suffix[m]+"_spec.fits");
        fspec.write_empty(); // empty primary extension, KMOS convention

        auto write_wcs = [&]() {
            fspec.write_keyword("CTYPE1", "WAVE");
            fspec.write_keyword("CUNIT1", "um");
            fspec.write_keyword("CRPIX1", 1.0);
            fspec.write_keyword("CRVAL1", lam[0]);
            fspec.write_keyword("CDELT1", lam[1]-lam[0]);
            for (uint_t i : range(used_files)) {
                fspec.write_keyword("EXP"+align_right(to_string(i), 3, '0'), file::get_basename(used_files[i]));
            }
        };

        fspec.reach_hdu(1);
        fspec.write(flx1d(m,_));
        write_wcs();

        fspec.reach_hdu(2);
        fspec.write(err1d(m,_));
        write_wcs();

        fspec.reach_hdu(3);
        fspec.write(max(err1db(m,_), err1d(m,_)));
        write_wcs();
    }

    return 0;
}

void print_help() {
    using namespace format;

    print("multispecfit v1.0");
    print("usage: multispecfit_sof <list.sof> <models.fits> suffix=[...] [options]");
    print("");
    print("Main parameters:");
    paragraph("WIP.");
    print("Available options:");
    bullet("background", "Set this flag if you want the program to fit for a constant background "
        "level on top of your models (default: no background model). The background spectrum will "
        "be saved with suffix \"_background\".");
    bullet("outdir", "Name of the directory into which the output files should be created. "
        "Default is the current directory.");
}
