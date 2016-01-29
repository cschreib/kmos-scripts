#include <phypp.hpp>
#include <phypp/mpfit.hpp>

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
    double dz = 0.3;
    double width_min = 0.0003;
    double width_max = 0.003;
    double fix_width = dnan;
    bool same_width = false;
    bool brute_force_width = false;
    bool subtract_continuum = true;
    uint_t continuum_width = 100;
    bool verbose = false;
    bool save_model = false;
    std::string outdir;
    vec1s tlines;

    // Read command line arguments
    read_args(argc-1, argv+1, arg_list(z0, dz, name(tlines, "lines"), width_min, width_max,
        subtract_continuum, continuum_width, verbose, same_width, save_model, fix_width,
        brute_force_width
    ));

    if (!outdir.empty()) {
        outdir = file::directorize(outdir);
        file::mkdir(outdir);
    }

    // Check validity of input
    bool bad = false;
    if (!is_finite(z0)) {
        error("please provide the fiducial redshift z0=...");
        bad = true;
    }
    if (tlines.empty()) {
        error("please provide the name of an emission line(s) to fit with lines=[...]");
        note("available lines:");
        print_available_lines(linedb);
        bad = true;
    } else {
        bool bad2 = false;
        for (std::string& l : tlines) {
            if (linedb.find(l) == linedb.end()) {
                error("unknown line '", l, "'");
                bad = true;
                bad2 = true;
            }
        }

        if (bad2) {
            note("available lines:");
            print_available_lines(linedb);
        }
    }

    if (bad) return 1;

    // Get lines in database
    double lambda_min = finf, lambda_max = -finf;
    vec<1,line_t> lines;
    for (std::string& l : tlines) {
        lines.push_back(linedb.find(l)->second);
        auto lmima = minmax(lines.back().lambda);
        lambda_min = min(lambda_min, lmima.first);
        lambda_max = max(lambda_max, lmima.second);
    }

    // Read spectrum
    if (verbose) note("read input spectrum...");
    std::string infile = argv[1];
    vec1d flx, err;
    fits::input_image fimg(infile);
    fimg.reach_hdu(1);
    fimg.read(flx);
    fimg.reach_hdu(2);
    fimg.read(err);

    // Renormalize uncertainties to get print,

    // Build wavelength axis
    uint_t nlam = flx.size();
    double cdelt = 1, crpix = 1, crval = 1;
    if (!fimg.read_keyword("CDELT1", cdelt) ||
        !fimg.read_keyword("CRPIX1", crpix) ||
        !fimg.read_keyword("CRVAL1", crval)) {
        error("could not read WCS information for wavelength axis");
        return 1;
    }

    vec1d lam = crval + cdelt*(findgen(nlam) + (1 - crpix));

    // Subtract continuum (optional)
    if (subtract_continuum) {
        if (verbose) note("estimate and subtract continuum...");
        vec1d tmp = flx;
        for (uint_t l : range(nlam)) {
            uint_t l0 = max(0, int_t(l)-int_t(continuum_width/2));
            uint_t l1 = min(nlam-1, int_t(l)+int_t(continuum_width/2));
            flx[l] -= median(tmp[l0-_-l1]);
        }
    }

    // Select a wavelength domain centered on the line(s)
    vec1u idl = where(lam > lambda_min*(1.0+z0-2*dz)
        && lam < lambda_max*(1.0+z0+2*dz)
        && is_finite(flx) && is_finite(err) && err > 0);

    if (idl.empty()) {
        error("the chosen lines are not covered by the provided cube at z=", z0, " +/- ", dz);
        idl = where(is_finite(flx) && is_finite(err) && err > 0);
        note("the cube covers ", min(lam[idl]), " to ", max(lam[idl]));
        note("your redshift search for this line requires a range within ",
            lambda_min*(1.0+z0-2*dz), " to ", lambda_max*(1.0+z0+2*dz));
        return 1;
    }

    flx = flx[idl];
    err = err[idl];
    lam = lam[idl];

    if (verbose) {
        note("fitting ", flx.dims[0], " spectral elements between ",
           lam[0], " and ", lam.back(), " um...");
    }

    // Define redshift grid so as to have two redshift samples per wavelength element
    double dlam = lam[1]-lam[0];
    uint_t nz = ceil(2*(2*dz)/(dlam/lines[0].lambda[0]));
    vec1d zs = rgen(z0-dz, z0+dz, nz);

    // Define width grid so as to have four width samples per wavelength element
    uint_t nwidth = ceil(4*(width_max - width_min)/dlam);
    vec1d ws = rgen(width_min, width_max, nwidth);

    // Perform a redshift search
    if (verbose) note("redshift search...");

    double chi2 = dinf;
    double z = dnan;
    vec1d pz(zs.size());
    vec1d flux(lines.size());
    vec1d flux_err(lines.size());
    vec1d width(lines.size());
    vec1d width_err(lines.size());

    vec1d best_model;

    auto pg = progress_start(zs.size());
    for (uint_t iz : range(zs)) {
        // Initialize starting conditions and specify fitting constraints
        uint_t nparam = (same_width ? 1+lines.size() : 2*lines.size());
        vec1d p(nparam);

        // Function to perform a non-linear fit (if varying the line widths)
        auto try_nlfit = [&]() {
            auto model = [&](const vec1d& l, const vec1d& tp) {
                vec1d m(l.dims);
                if (same_width) {
                    for (uint_t il : range(lines))
                    for (uint_t isl : range(lines[il].lambda)) {
                        m += (tp[il+1]*lines[il].ratio[isl]/(sqrt(2.0*dpi)*tp[0]*1e4))*exp(
                            -sqr(lines[il].lambda[isl]*(1.0+zs[iz]) - l)/(2.0*sqr(tp[0]))
                        );
                    }
                } else {
                    for (uint_t il : range(lines))
                    for (uint_t isl : range(lines[il].lambda)) {
                        m += (tp[2*il]*lines[il].ratio[isl]/(sqrt(2.0*dpi)*tp[2*il+1]*1e4))*exp(
                            -sqr(lines[il].lambda[isl]*(1.0+zs[iz]) - l)/(2.0*sqr(tp[2*il+1]))
                        );
                    }
                }

                return m;
            };

            mpfit_result res = mpfitfun(flx, err, lam, model, p);

            bool better = res.success && res.chi2 < chi2;
            if (better) {
                chi2 = res.chi2;
                z = zs[iz];
                if (same_width) {
                    flux         = res.params[uindgen(lines.size())+1];
                    flux_err     = res.errors[uindgen(lines.size())+1];
                    width[_]     = res.params[0];
                    width_err[_] = res.errors[0];
                } else {
                    flux      = res.params[2*uindgen(lines.size())+0];
                    flux_err  = res.errors[2*uindgen(lines.size())+0];
                    width     = res.params[2*uindgen(lines.size())+1];
                    width_err = res.errors[2*uindgen(lines.size())+1];
                }

                if (save_model) {
                    best_model = model(lam, res.params);
                }
            }

            if (res.success) {
                pz[iz] = res.chi2;
            }

            return better;
        };

        // Function to perform a linear fit (fixed line widths)
        auto try_lfit = [&]() {
            vec2d m(lines.size(), lam.size());
            if (same_width) {
                for (uint_t il : range(lines))
                for (uint_t isl : range(lines[il].lambda)) {
                    m(il,_) += (lines[il].ratio[isl]/(sqrt(2.0*dpi)*p[0]*1e4))*exp(
                        -sqr(lines[il].lambda[isl]*(1.0+zs[iz]) - lam)/(2.0*sqr(p[0]))
                    );
                }
            } else {
                for (uint_t il : range(lines))
                for (uint_t isl : range(lines[il].lambda)) {
                    m(il,_) += (lines[il].ratio[isl]/(sqrt(2.0*dpi)*p[2*il+1]*1e4))*exp(
                        -sqr(lines[il].lambda[isl]*(1.0+zs[iz]) - lam)/(2.0*sqr(p[2*il+1]))
                    );
                }
            }

            linfit_result res = linfit_pack(flx, err, m);

            bool better = res.success && res.chi2 < chi2;
            if (better) {
                chi2 = res.chi2;
                z = zs[iz];

                flux     = res.params;
                flux_err = res.errors;

                if (same_width) {
                    p[uindgen(lines.size())+1] = flux;
                    width[_]     = p[0];
                    width_err[_] = 0.0;
                } else {
                    p[2*uindgen(lines.size())] = flux;
                    width        = p[2*uindgen(lines.size())+1];
                    width_err[_] = 0.0;
                }

                if (save_model) {
                    best_model = vec1d(lam.dims);
                    for (uint_t il : range(lines)) {
                        best_model += flux[il]*m(il,_);
                    }
                }
            }

            if (res.success) {
                pz[iz] = res.chi2;
            }

            return better;
        };

        // Estimate starting parameters and/or fix the parameter values
        // and do the fit
        if (brute_force_width && !is_finite(fix_width)) {
            if (same_width) {
                for (uint_t iw : range(ws)) {
                    p[0] = ws[iw];

                    for (uint_t il : range(lines)) {
                        p[il+1] = interpolate(flx, lam, lines[il].lambda[0]*(1.0 + zs[iz]))*
                            sqrt(2.0*dpi)*p[0]*1e4;
                    }

                    try_lfit();
                }
            } else {
                for (uint_t il : range(lines)) {
                    p[2*il+1] = 0.5*(width_max - width_min);
                    p[2*il+0] = interpolate(flx, lam, lines[il].lambda[0]*(1.0 + zs[iz]))*
                        sqrt(2.0*dpi)*p[2*il+1]*1e4;
                }

                for (uint_t il : range(lines)) {
                    for (uint_t iw : range(ws)) {
                        p[2*il+1] = ws[iw];
                        p[2*il+0] = interpolate(flx, lam, lines[il].lambda[0]*(1.0 + zs[iz]))*
                            sqrt(2.0*dpi)*p[2*il+1]*1e4;

                        try_lfit();
                    }

                    p[2*il+1] = width[il];
                }
            }
        } else {
            if (same_width) {
                if (is_finite(fix_width)) {
                    p[0] = fix_width;
                } else {
                    p[0] = 0.5*(width_max - width_min);
                }

                for (uint_t il : range(lines)) {
                    p[il+1] = interpolate(flx, lam, lines[il].lambda[0]*(1.0 + zs[iz]))*
                        sqrt(2.0*dpi)*p[0]*1e4;
                }
            } else {
                for (uint_t il : range(lines)) {
                    if (is_finite(fix_width)) {
                        p[2*il+1] = fix_width;
                    } else {
                        p[2*il+1] = 0.5*(width_max - width_min);
                    }

                    p[2*il+0] = interpolate(flx, lam, lines[il].lambda[0]*(1.0 + zs[iz]))*
                        sqrt(2.0*dpi)*p[2*il+1]*1e4;
                }
            }

            if (is_finite(fix_width)) {
                try_lfit();
            } else {
                try_nlfit();
            }
        }

        if (verbose) progress(pg);
    }

    uint_t ndof = flx.size() - lines.size();
    if (!is_finite(fix_width)) {
        if (same_width) ndof -= 1;
        else            ndof -= lines.size();
    }

    pz = exp(-(pz - chi2)/ndof);

    // Write the result
    std::string filebase = outdir+file::remove_extension(file::get_basename(argv[1]));
    if (verbose) note("write to disk...");
    fits::output_table otbl(filebase+"_slfit.fits");
    otbl.write_columns(ftable(chi2, z, flux, flux_err, width, width_err));
    otbl.write_columns("lines", tlines, "pzx", zs, "pzy", pz);

    if (save_model) {
        fits::output_image ospec(filebase+"_slfit_model.fits");
        ospec.write(vec1d(0)); // empty primary
        ospec.reach_hdu(1);
        ospec.write(best_model);
        ospec.write_header(fimg.read_header());
    }

    return 0;
}

void print_help(const std::map<std::string,line_t>& db) {
    using namespace format;

    print("slinefit v1.0");
    print("usage: slinefit <kmos_spectrum.fits> z0=... dz=... lines=... [options]");
    print("");
    print("Main parameters:");
    paragraph("'kmos_spectrum.fits' must be a valid KMOS spectrum file, i.e., a FITS "
        "file containing 3 extensions: the first must be empty (KMOS convention), the "
        "second contains the flux, while the third contains the uncertainty. Using this "
        "spectrum, the tool searches for lines around a \"first-guess\" redshift of 'z0', "
        "within 'z0-dz' and 'z0+dz' (by default, dz=0.3). It will try to identify spectral "
        "features with the emission lines you provide in the 'lines' list (see below for a "
        "list of available lines and their code names). You can specify as many lines as "
        "you wish, and the order is irrelevant. For example: 'lines=\"[hbeta,o3]\"'. The "
        "program will adjust the redshift, as well as both the line fluxes and widths to "
        "best match the observed spectrum. It will output a column-oriented FITS table "
        "containing the best-fit redshift (and its probability distribution), the chi2, "
        "and the fluxes, widths and uncertainties for each lines. This catalog is saved "
        "into the '*_slfit.fits' file. You can open it within IDL with 'mrdfits(\"...\", 1)'.");
    print("Available lines:");
    print_available_lines(db);
    print("\nAvailable options (in order of importance):");
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
    bullet("same_width", "Set this flag if you want to force all lines to have the same "
        "width, rather than fit them independently. This can help solving "
        "fit instability issues if some lines are very low S/N and the fitted line widths "
        "diverge to unreasonable values. The default is to let each width vary freely and "
        "independently in the fit.");
    bullet("fix_width=...", "Must be a number. If provided, the program will force the "
        "line widths to be equal to this value for all the lines. This can help solving "
        "fit instability issues if some lines are very low S/N and the fitted line widths "
        "diverge to unreasonable values. The default is to let each width vary freely in "
        "the fit.");
    bullet("brute_force_width", "Set this flag to use a brute force approach to compute "
        "the line's widths. This can be used as a fallback solution if the default "
        "algorithm (a non-linear fit using the Levenberg-Marquardt technique) fails or "
        "gives garbage results. Using the brute force approach will provide more stable "
        "results and a guaranteed convergence, but may require more computation time. In "
        "practice the difference in performance is not very large because we make the "
        "assumptions that there is no covariance between the widths of individual lines, "
        "so they are varied one after the other rather than all at once. Also we can use "
        "a simple linear fit to adjust the line fluxes, which is itself much faster.");
    bullet("width_min=...", "Must be a number. Defines the minimum allowed width for a "
        "line in microns. Default is 0.0003. Note that this value is only used if "
        "'brute_force_width' is set. If not, then 'width_min' and 'width_max' are simply "
        "used to estimate the initial value of the line width for the non-linear fit "
        "(which is taken as the average between the two).");
    bullet("width_max=...", "See above.");
    bullet("save_model", "Set this flag to also output the best-fit model spectrum. The "
        "spectrum will be saved into the '*_slfit_model.fits' file as a regular KMOS "
        "spectrum: the first extension is empty, the second contains the flux.");
    bullet("outdir", "Name of the directory into which the output files should be created. "
        "Default is the current directory..");
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