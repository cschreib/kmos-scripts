#include <phypp.hpp>

void print_help();

int phypp_main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    vec1s suffix;
    std::string fitmask;
    std::string outdir;
    bool residuals = false;
    bool background = false;
    uint_t hdu = 1;
    uint_t error_hdu = 2;
    read_args(argc-2, argv+2, arg_list(background, suffix, fitmask, residuals, outdir, hdu, error_hdu));

    if (!outdir.empty()) {
        outdir = file::directorize(outdir);
        file::mkdir(outdir);
    }

    // Read cube
    fits::input_image fimg(argv[1]);
    vec3d flx;
    fimg.reach_hdu(hdu);
    fimg.read(flx);
    // Read uncertainty
    vec3d err;
    fimg.reach_hdu(error_hdu);
    fimg.read(err);

    // Get wavelength WCS
    fimg.reach_hdu(hdu);
    astro::wcs wcs(fimg.read_header());
    vec1d lam = astro::build_axis(wcs, 0, astro::axis_unit::wave_um);

    // Read fit mask (optional)
    vec2d mask; {
        if (!fitmask.empty()) {
            fits::read(fitmask, mask);
        } else {
            mask = vec2d{partial_count(0, is_finite(flx)) > 1};
        }
    }

    // Read models
    vec2d models;
    uint_t nmodel; {
        fits::input_image fmodel(argv[2]);
        if (fmodel.axis_count() == 2) {
            nmodel = 1;
            fits::read(argv[2], models);
            models = reform(models, 1, models.size());
        } else {
            vec3d tmp;
            fits::read(argv[2], tmp);

            nmodel = tmp.dims[0];
            models = reform(tmp, nmodel, tmp.dims[1]*tmp.dims[2]);
        }

        // Add constant background if requested
        if (background) {
            append<0>(models, replicate(1.0, 1, models.dims[1]));
            if (nmodel == 1) {
                suffix = {"", "background"};
            } else if (suffix.size() == nmodel) {
                suffix.push_back("background");
            }

            ++nmodel;
        }

        // Make sure each model has unit integral
        for (uint_t i : range(nmodel)) {
            vec1u idg = where(is_finite(models(i,_)));
            models(i,_) /= total(models(i,_)[idg]);
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

    // Do the fit wavelength by wavelength
    vec2d flx1d(flx.dims[0], nmodel);
    vec2d err1d(flx.dims[0], nmodel);
    for (uint_t l : range(flx.dims[0])) {
        vec1u idg = where(is_finite(flx(l,_,_)) && is_finite(err(l,_,_)) && err(l,_,_) > 0 && mask > 0.0);
        auto res = linfit_pack(flx(l,_,_)[idg], err(l,_,_)[idg], models(_,idg));
        flx1d(l,_) = res.params;
        err1d(l,_) = res.errors;
    }

    // Write spectra to disk
    std::string filebase = outdir+file::get_basename(erase_end(argv[1], ".fits"));
    for (uint_t i : range(nmodel)) {
        fits::output_image fspec(filebase+suffix[i]+"_spec.fits");
        fspec.write_empty(); // empty primary extension, KMOS convention

        auto write_wcs = [&]() {
            fspec.write_keyword("CTYPE1", "WAVE");
            fspec.write_keyword("CUNIT1", "um");
            fspec.write_keyword("CRPIX1", 1.0);
            fspec.write_keyword("CRVAL1", lam[0]);
            fspec.write_keyword("CDELT1", lam[1]-lam[0]);
        };

        fspec.reach_hdu(1);
        fspec.write(flx1d(_,i));
        write_wcs();

        fspec.reach_hdu(2);
        fspec.write(err1d(_,i));
        write_wcs();
    }

    // Compute residuals if asked
    if (residuals) {
        fimg.reach_hdu(1);
        fits::header hdr = fimg.read_header();

        for (uint_t i : range(nmodel)) {
            vec3d res = flx;
            for (uint_t l : range(res.dims[0])) {
                flatten(res(l,_,_)) -= flx1d(l,i)*models(i,_);
            }

            fits::output_image fspec(filebase+suffix[i]+"_residual.fits");
            fspec.write_empty(); // empty primary HDU

            fspec.reach_hdu(1);
            fspec.write(res);
            fspec.write_header(hdr);

            fspec.reach_hdu(2);
            fspec.write(err);
            fspec.write_header(hdr);
        }
    }

    return 0;
}

void print_help() {
    using namespace format;

    print("multispecfit v1.0");
    print("usage: multispecfit <kmos_cube.fits> <models.fits> suffix=[...] [options]");
    print("");
    print("Main parameters:");
    paragraph("'kmos_cube.fits' must be a cube created by the KMOS pipeline, with 3 "
        "extensions: the first is empty (KMOS convention), the second contains the flux, "
        "and the third contains the uncertainty. 'models.fits' must be a cube created by "
        "yourself (e.g., with IDL or Python) containing the models to fit. Each slice of "
        "this cube corresponds to a different model. Note that you have to include a "
        "uniform model (all pixels equal to 1) if you want to fit a constant background. "
        "Finally, 'suffix' must contain the names of each model, with the following "
        "format: [name1,name1,name3,...]. There must be as many suffixes as there are "
        "models in the 'models.fits' file.");
    print("Available options:");
    bullet("fitmask=...", "Must be a 2D FITS file containing the mask to define the fitting "
        "region. Only the pixels with a non-zero value in the mask will be used "
        "(default: all valid pixels are used).");
    bullet("residuals", "Set this flag if you want the program to generate residual cubes, "
        "removing the contribution of each component separately (default: do not "
        "generate residuals)");
    bullet("background", "Set this flag if you want the program to fit for a constant background "
        "level on top of your models (default: no background model). The background spectrum will "
        "be saved with suffix \"_background\".");
    bullet("outdir", "Name of the directory into which the output files should be created. "
        "Default is the current directory.");
    bullet("hdu", "Choose from which HDU to read the flux cube (default: HDU 1, first extension).");
    bullet("error_hdu", "Choose from which HDU to read the error cube (default: HDU 2, second "
        "extension).");
}
