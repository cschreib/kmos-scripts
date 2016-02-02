#include <phypp.hpp>
#include <phypp/mpfit.hpp>

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

// Local functions, defined at the end of the file
void print_help(const std::map<std::string,line_t>& db);
void print_available_lines(const std::map<std::string,line_t>& db);

int main(int argc, char* argv[]) {
    // Build the line data base (you can add your own there!)
    std::map<std::string,line_t> linedb = {
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
        {"palpha",  line_t("palpha",  {1.875},          {1.0})}
    };

    if (argc < 2) {
        print_help(linedb);
        return 0;
    }

    double z0 = dnan;
    double dz = 0.3;
    double width_min = 50.0;
    double width_max = 500.0;
    double delta_width = 0.2;
    double delta_z = 0.2;
    double fix_width = dnan;
    uint_t lambda_pad = 5;
    bool same_width = false;
    bool use_mpfit = false;
    bool verbose = false;
    bool save_model = false;
    std::string outdir;
    bool ascii = false;
    vec1s tlines;

    // Read command line arguments
    read_args(argc-1, argv+1, arg_list(z0, dz, name(tlines, "lines"), width_min, width_max,
        verbose, same_width, save_model, fix_width, use_mpfit, ascii, outdir, delta_width,
        delta_z, lambda_pad
    ));

    if (!outdir.empty()) {
        outdir = file::directorize(outdir);
        file::mkdir(outdir);
    }

    if (is_finite(fix_width) && use_mpfit) {
        warning("'use_mpfit' has no effect if 'fix_width' is set");
        use_mpfit = false;
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
            if (l.find(':') != l.npos) continue;
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
        if (l.find(':') != l.npos) {
            vec1s spl = split(l, ":");
            if (spl.size() < 2 || (spl.size() > 2 && spl.size()%2 != 1)) {
                error("ill-formed line declaration '", l, "'");
                error("custom line declaration must be of the form 'name:lambda' or "
                    "'name:lambda1:lambda2:...:ratio1:ratio2,...'");
                return 1;
            }

            line_t nl;
            nl.name = spl[0];

            vec1d nums;
            if (count(!from_string(spl[1-_], nums)) != 0) {
                error("could not convert line wavelengths and ratios in '", l, "' into a "
                    "list of numbers");
                return 1;
            }

            if (nums.size() == 1) {
                nl.lambda = nums;
            } else {
                nl.lambda = nums[uindgen(nums.size()/2)];
                nl.ratio  = nums[uindgen(nums.size()/2) + nums.size()/2];
            }

            lines.push_back(nl);
        } else {
            lines.push_back(linedb.find(l)->second);
        }

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
    uint_t orig_nlam = lam.size();

    // Identify good regions of the spectrum
    vec1b goodspec = is_finite(flx) && is_finite(err) && err > 0;
    vec1u idl = where(goodspec);
    if (idl.size() <= lambda_pad*2) {
        error("the spectrum does not contain any valid point");
        return 1;
    }

    // Flag out the pixels at the border of the spectrum
    goodspec[idl.front()-_-(idl.front()+lambda_pad)] = false;
    goodspec[(idl.back()-lambda_pad)-_-idl.back()] = false;

    // Select a wavelength domain centered on the line(s)
    idl = where(lam > lambda_min*(1.0+z0-2*dz)
        && lam < lambda_max*(1.0+z0+2*dz)
        && goodspec);

    if (idl.empty()) {
        error("the chosen lines are not covered by the provided cube at z=", z0, " +/- ", dz);
        idl = where(goodspec);
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

    // Define redshift grid so as to have the requested number of samples per wavelength element
    double tdz = delta_z*cdelt/(lines[0].lambda[0]*(1.0+z0));
    uint_t nz = ceil(2*dz/tdz);
    vec1d zs = rgen(z0-dz, z0+dz, nz);

    // Define width grid so as to have the requested number of samples per wavelength element
    double dwidth = 3e5*delta_width*cdelt/(lines[0].lambda[0]*(1.0+z0));
    uint_t nwidth = ceil((width_max - width_min)/dwidth);
    vec1d ws = rgen(width_min, width_max, nwidth);

    // Perform a redshift search
    if (verbose) {
        if (!use_mpfit && !is_finite(fix_width)) {
            note("redshift search (", nz, " redshifts, step = ", zs[1] - zs[0], ", ",
                nwidth, " line widths, step = ", ws[1] - ws[0], ")...");
        } else {
            note("redshift search (", nz, " redshifts, step = ", zs[1] - zs[0], ")...");
        }
    }

    double chi2 = dinf;
    double z = dnan;
    vec1d pz(zs.size());
    vec1d flux(lines.size());
    vec1d flux_err(lines.size());
    vec1d width(lines.size());
    vec1d width_err(lines.size());

    vec1d best_model;

    // Renormalize flux and errors in units of 1e-17 erg/s/cm2 to avoid
    // numerical imprecision (it is always better to deal with numbers close to unity).
    // This helps prevent mpfit getting stuck and not varying the line widths.
    flx *= 1e17; err *= 1e17;

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
                        double lw = (tp[0]/3e5)*lines[il].lambda[isl]*(1.0+zs[iz]);
                        m += (tp[il+1]*lines[il].ratio[isl]/(sqrt(2.0*dpi)*lw*1e4))*exp(
                            -sqr(lines[il].lambda[isl]*(1.0+zs[iz]) - l)/(2.0*sqr(lw))
                        );
                    }
                } else {
                    for (uint_t il : range(lines))
                    for (uint_t isl : range(lines[il].lambda)) {
                        double lw = (tp[2*il+1]/3e5)*lines[il].lambda[isl]*(1.0+zs[iz]);
                        m += (tp[2*il]*lines[il].ratio[isl]/(sqrt(2.0*dpi)*lw*1e4))*exp(
                            -sqr(lines[il].lambda[isl]*(1.0+zs[iz]) - l)/(2.0*sqr(lw))
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
                    double lw = (p[0]/3e5)*lines[il].lambda[isl]*(1.0+zs[iz]);
                    m(il,_) += (lines[il].ratio[isl]/(sqrt(2.0*dpi)*lw*1e4))*exp(
                        -sqr(lines[il].lambda[isl]*(1.0+zs[iz]) - lam)/(2.0*sqr(lw))
                    );
                }
            } else {
                for (uint_t il : range(lines))
                for (uint_t isl : range(lines[il].lambda)) {
                    double lw = (p[2*il+1]/3e5)*lines[il].lambda[isl]*(1.0+zs[iz]);
                    m(il,_) += (lines[il].ratio[isl]/(sqrt(2.0*dpi)*lw*1e4))*exp(
                        -sqr(lines[il].lambda[isl]*(1.0+zs[iz]) - lam)/(2.0*sqr(lw))
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
        if (!use_mpfit && !is_finite(fix_width)) {
            if (same_width) {
                for (uint_t iw : range(ws)) {
                    p[0] = ws[iw];
                    try_lfit();
                }
            } else {
                for (uint_t il : range(lines)) {
                    p[2*il+1] = 0.5*(width_max - width_min);
                }

                for (uint_t il : range(lines)) {
                    for (uint_t iw : range(ws)) {
                        p[2*il+1] = ws[iw];
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

                    for (uint_t il : range(lines)) {
                        double lw = (p[0]/3e5)*lines[il].lambda[0]*(1.0+zs[iz]);
                        p[il+1] = interpolate(flx, lam, lines[il].lambda[0]*(1.0 + zs[iz]))*
                            sqrt(2.0*dpi)*lw*1e4;
                    }
                }
            } else {
                for (uint_t il : range(lines)) {
                    if (is_finite(fix_width)) {
                        p[2*il+1] = fix_width;
                    } else {
                        p[2*il+1] = 0.5*(width_max - width_min);

                        double lw = (p[2*il+1]/3e5)*lines[il].lambda[0]*(1.0+zs[iz]);
                        p[2*il+0] = interpolate(flx, lam, lines[il].lambda[0]*(1.0 + zs[iz]))*
                            sqrt(2.0*dpi)*lw*1e4;
                    }
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

    // Compute number of degrees of freedom for reduced chi2
    uint_t ndof = flx.size() - lines.size();
    if (!is_finite(fix_width)) {
        if (same_width) ndof -= 1;
        else            ndof -= lines.size();
    }

    // Compute redshift 'probability'
    pz = exp(-(pz - chi2)/ndof);
    pz -= min(pz);
    pz /= integrate(zs, pz);

    vec1d cpz = cumul(zs, pz); cpz /= max(cpz);
    double zlow = z - interpolate(zs, cpz, 0.16);
    double zup = interpolate(zs, cpz, 0.84) - z;

    if (verbose) {
        print("best redshift: ", z, " + ", zup, " - ", zlow,
            " (chi2: ", chi2, ", reduced: ", chi2/ndof, ")");
    }

    // Rescale fluxes and uncertainties
    flux *= 1e-17; flux_err *= 1e-17;

    // Build llambda
    vec1d llambda;
    for (uint_t il : range(lines)) {
        llambda.push_back(lines[il].lambda[0]*(1.0 + z));
    }

    // Ungroup line groups
    for (uint_t il : range(lines)) {
        auto& l = lines[il];

        if (l.lambda.size() == 1) continue;
        for (uint_t i : range(1, l.lambda.size())) {
            flux.push_back(flux[il]*l.ratio[i]);
            flux_err.push_back(flux_err[il]*l.ratio[i]);
            width.push_back(width[il]);
            width_err.push_back(width_err[il]);
            tlines.push_back(tlines[il]+"-"+strn(i+1));
            llambda.push_back(l.lambda[i]*(1.0 + z));
        }

        tlines[il] = tlines[il]+"-1";
    }

    // Sort by wavelength
    vec1u ids = sort(llambda);
    flux = flux[ids]; flux_err = flux_err[ids]; width = width[ids]; width_err = width_err[ids];
    llambda = llambda[ids]; tlines = tlines[ids];

    // Write the result
    std::string filebase = outdir+file::remove_extension(file::get_basename(argv[1]));
    if (verbose) note("write to disk...");

    if (ascii) {
        vec1s hdr = {"line", "flux [erg/s/cm2]", "error", "width [um]", "error"};
        file::write_table_hdr(filebase+"_slfit_lines.cat", 18, hdr,
            tlines, llambda, strna_sci(flux), strna_sci(flux_err), strna_sci(width),
            strna_sci(width_err)
        );

        hdr = {"redshift", "P(z)"};
        file::write_table_hdr(filebase+"_slfit_pz.cat", 18, hdr, zs, strna_sci(pz));
    } else {
        fits::output_table otbl(filebase+"_slfit.fits");
        otbl.write_columns(ftable(chi2, z, zup, zlow, flux, flux_err, width, width_err));
        otbl.write_columns("lambda", llambda, "lines", tlines, "pzx", zs, "pzy", pz);
    }

    if (save_model) {
        // Rescale model
        best_model *= 1e-17;

        // First bring back the model into the original wavelength grid
        vec1d nmodel = replicate(dnan, orig_nlam);
        nmodel[idl] = best_model;

        // Then save it
        fits::output_image ospec(filebase+"_slfit_model.fits");
        ospec.write(vec1d(0)); // empty primary
        ospec.reach_hdu(1);
        ospec.write(nmodel);
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
    print("\nNote: you can add your own lines either by modifying the source code of the "
        "program, or directly into the command line arguments. In the 'lines=[...]' "
        "parameter, you can indeed create a new line (or set of lines) with the synthax "
        "'name:lambda' (for a single line) or 'name:lambda1:lambda2:...:ratio1:ratio2:...' "
        "(for a group of lines). In this case, 'name' can be whatever you want (should "
        "not contain spaces), 'lambda' must be the rest-frame wavelength of the line in "
        "microns, and 'ratioX' must be the fixed flux ratio of the the line 'X' and the "
        "first line of the group (i.e., it should be '1' for the first line). For example, "
        "to fit the [SII] doublet: lines=[s2:0.67183:0.67327:1:0.75]. In the fit, the "
        "flux of the [SII]6733 line will be forced to be a factor 0.75 lower than that of "
        "[SII]6718.");
    print("\nAvailable options (in order of importance):");
    bullet("delta_z=...", "Must be a number. Defines the size of a step in the grid of "
        "redshifts, as the fraction of the size of a wavelength element of the spectrum. "
        "In other words, given the spectral resolution R of your spectrum, the redshift "
        "step will be equal to delta_z/R. Default is 0.2, which corresponds to 0.00005 "
        "at R=3800 (H+K) and 0.00003 at R=7100 (K). There is no much need to user smaller "
        "steps since this is already hitting the limits of the spectral resolution, "
        "however you may wish to increase the size of the step if you need more "
        "performance.");
    bullet("width_min=...", "Must be a number. Defines the minimum allowed width for a "
        "line in km/s. Default is 50. Note that this value is not used if 'use_mpfit' is "
        "set (in this case, 'width_min' and 'width_max' are simply used to estimate the "
        "initial value of the line width for the non-linear fit, taken as the average "
        "between the two).");
    bullet("width_max=...", "See above. Default is 500 km/s.");
    bullet("delta_width=...", "Must be a number. Defines the size of a step in the grid of "
        "line widths, as the fraction of the size of a wavelength element of the spectrum. "
        "In other words, given the spectral resolution R of your spectrum, the width step "
        "will be equal to c*delta_width/R. Default is 0.2, which at corresponds to 15 km/s "
        "at R=3800 (H+K) and 8 km/s at R=7100 (K).");
    bullet("same_width", "Set this flag if you want to force all lines to have the same "
        "width, rather than fit them independently. This can help solving "
        "fit instability issues if some lines are very low S/N and the fitted line widths "
        "diverge to unreasonable values. The default is to let each width vary freely and "
        "independently in the fit.");
    bullet("fix_width=...", "Must be a line width (in km/s). If provided, the program will "
        "force the line widths to be equal to this value for all the lines. This can help "
        "solving fit instability issues if some lines are very low S/N and the fitted line "
        "widths diverge to unreasonable values. The default is to let each width vary "
        "freely in the fit.");
    bullet("use_mpfit", "Set this flag to use a non-linear fitting approach to fit the line "
        "profiles. This method uses the Levenberg-Marquardt technique to fit non linear "
        "models, which is more flexible and correct since it allows simultaneous fit of "
        "the fluxes and line widths of all the lines. However these algorithms are more "
        "unstable and can often not converge. The default method, if this flag is not set, "
        "is therefore to use a brute force approach, which is certain to converge, but "
        "may require more computation time. In practice the difference in performance is "
        "not so bad (it may even be faster), because we make the assumption that there is "
        "no covariance between the widths of individual lines, so they are varied one "
        "after the other rather than all at once. Also we can use a simple linear fit to "
        "adjust the line fluxes, which is itself much faster. So, use this flag as an "
        "experiment, but double check that the fit results make sense. Note that if "
        "'fix_width' is used, the fit will always be done with the default approach, since "
        "there is no need for a non-linear fit in this case.");
    bullet("lambda_pad", "Must be an integer. It defines the number of wavelength element "
        "that are ignored both at the beginning and end of the spectrum. Default is 5 "
        "elements. This is used to flag out invalid and poorly covered spectral regions "
        "which could drive the fit toward unrealistic values.");
    bullet("save_model", "Set this flag to also output the best-fit model spectrum. The "
        "spectrum will be saved into the '*_slfit_model.fits' file as a regular KMOS "
        "spectrum: the first extension is empty, the second contains the flux.");
    bullet("outdir", "Name of the directory into which the output files should be created. "
        "Default is the current directory..");
    bullet("ascii", "Set this flag if you want the output catalog to be saved in ASCII "
        "format rather than FITS. In this case, the lines and their fluxes will be saved "
        "in the '*_slfit_lines.cat' file, while the redshift probability distribution "
        "will be saved in '*_slfit_pz.cat'.");
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
