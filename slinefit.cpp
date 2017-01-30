#include <phypp.hpp>
#include <phypp/math/mpfit.hpp>

// Structure to define a line or line group to be fitted simultaneously
struct line_t {
    line_t() = default;
    line_t(std::string n, std::string pn, vec1d lam, vec1d ra) : name(n), pretty_name(pn), lambda(lam), ratio(ra) {
        ratio /= ratio[0];
    }

    void define_continuum(double width) {
        continuum_lmin = mean(lambda)*(1.0 - 0.5*width/3e5);
        continuum_lmax = mean(lambda)*(1.0 + 0.5*width/3e5);
    }

    std::string name; // code identifier of the line
    std::string pretty_name; // verbose description of the line
    vec1d lambda;     // wavelengths of the lines
    vec1d ratio;      // flux ratios of the lines relative to the first
    double continuum_lmin, continuum_lmax;
};

// Structure holding best fit parameters
struct fit_result_t {
    double chi2 = dinf;
    double z = dnan;
    double z_err = dnan;
    vec1d chi2_grid;

    vec1d lambda;
    vec1d lambda_err;
    vec1d flux;
    vec1d flux_err;
    vec1b free_width;
    vec1d width;
    vec1d width_err;
    vec1b free_offset;
    vec1d offset;
    vec1d offset_err;
    vec1d cont;
    vec1d cont_err;
    vec1d ew;
    vec1d ew_err;

    vec1d model;
    vec1d model_continuum;

    fit_result_t(uint_t tnz, uint_t nline) {
        chi2_grid = replicate(dinf, tnz);

        lambda.resize(nline);
        lambda_err.resize(nline);
        flux.resize(nline);
        flux_err.resize(nline);
        free_width.resize(nline);
        width.resize(nline);
        width_err.resize(nline);
        free_offset.resize(nline);
        offset.resize(nline);
        offset_err.resize(nline);
        cont.resize(nline);
        cont_err.resize(nline);
        ew.resize(nline);
        ew_err.resize(nline);
    }
};

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

// Grow a vector back into a larger vector from given size 'n' and IDs 'ids'
// New elements will take the value 'def'.
template <typename T, typename U>
vec<1,meta::rtype_t<T>> reshape(const vec<1,T>& v, const vec1u& ids, uint_t n, const U& def = 0.0) {
    vec<1,meta::rtype_t<T>> nv = replicate(meta::rtype_t<T>{def}, n);
    nv[ids] = v;
    return nv;
}

// Compute MC uncertainties from a set of repeated measurements given by function 'getp(i)'
template <typename F>
double get_mc_error(uint_t nmc, F&& getp) {
    vec1d mc(nmc);
    for (uint_t i : range(nmc)) {
        mc[i] = getp(i);
    }

    vec1u idf = where(is_finite(mc));
    if (!idf.empty()) {
        return stddev(mc[idf]);
    } else {
        return dnan;
    }
}

// Filter a boolean mask to only keep 'false' values when they form a contiguous group of at least
// 'mincount' elements, then expand these groups by 'padding' elements on either side.
// Edges fo the mask (first and last element) are assumed to be gap boundaries in all cases.
vec1b keep_gaps_and_expand(vec1b flags, uint_t mincount, uint_t padding) {
    if (flags.size() < 2*padding) {
        flags[_] = false;
        return flags;
    }

    // Build gaps
    bool is_gap = false;
    uint_t last_gap = npos;
    vec1u gi0, gi1;
    for (uint_t i : range(flags)) {
        if (is_gap) {
            if (flags[i]) {
                if (i-last_gap < mincount && last_gap != 0) {
                    flags[last_gap-_-(i-1)] = true;
                } else {
                    gi0.push_back(last_gap);
                    gi1.push_back(i-1);
                }

                is_gap = false;
                last_gap = npos;
            }
        } else {
            if (!flags[i]) {
                is_gap = true;
                last_gap = i;
            }
        }
    }

    if (is_gap) {
        gi0.push_back(last_gap);
        gi1.push_back(flags.size()-1);
    }

    // Expand them
    if (padding != 0) {
        for (uint_t g : range(gi0)) {
            uint_t i0 = (gi0[g] > padding                ? gi0[g] - padding : 0);
            uint_t i1 = (gi1[g] < flags.size()-1-padding ? gi1[g] + padding : flags.size()-1);
            flags[i0-_-i1] = false;
        }

        // Treat boundaries
        flags[_-(padding-1)] = false;
        flags[(flags.size()-padding)-_] = false;
    }

    return flags;
}

// Local functions, defined at the end of the file
void print_help(const std::map<std::string,line_t>& db);
void print_available_lines(const std::map<std::string,line_t>& db);

int phypp_main(int argc, char* argv[]) {
    // Prepare the line database
    std::map<std::string,line_t> linedb;
    auto add_line = [&](line_t l) {
        auto iter = linedb.find(l.name);
        if (iter != linedb.end()) {
            warning("a line with name '", l.name, "' already exists in the database");
        } else {
            linedb.insert(std::make_pair(l.name, l));
        }
    };

    // Add default lines to the database (you can add your own!)
    //
    // Some references:
    // http://astronomy.nmsu.edu/drewski/tableofemissionlines.html
    // http://physics.nist.gov/PhysRefData/ASD/lines_form.html
    // http://adsabs.harvard.edu/abs/2011AJ....141...37L
    //
    // NB: for doublets of unknown flux ratios, 1 was assumed
    //
    // ---------------------------------------------------------
    // Standard galaxy/AGN UV-NIR emission lines
    add_line({"em_lyalpha",   "Lyman alpha",                      {0.12157},          {1.0}      });
    add_line({"em_n5_1240",   "Nitrogen V -- NV (doublet)",       {0.12388, 0.12428}, {1.0, 1.0} });
    add_line({"em_si4_1400",  "Silicium IV -- SiIV (doublet)",    {0.13938, 0.14028}, {1.0, 1.0} });
    add_line({"em_c4_1550",   "Carbon IV -- CIV (doublet)",       {0.15482, 0.15508}, {1.0, 1.0} });
    add_line({"em_he2_1640",  "Helium II -- HeII",                {0.16404},          {1.0}      });
    add_line({"em_c3_1909",   "Carbon III -- CIII]",              {0.19087},          {1.0}      });
    add_line({"em_c2_2326",   "Carbon II -- CII] (doublet)",      {0.23235, 0.23247}, {1.0, 1.0} });
    add_line({"em_ne4_2422",  "Neon IV -- NeIV",                  {0.24218},          {1.0}      });
    add_line({"em_mg2_2799",  "Magnesium II  -- MgII (doublet)",  {0.27964, 0.28035}, {1.0, 1.0} });
    add_line({"em_ne5_3426",  "Neon V -- NeV",                    {0.34259},          {1.0}      });
    add_line({"em_o2_3727",   "Oxygen II -- [OII] (doublet)",     {0.37260, 0.37288}, {1.0, 1.0} });
    add_line({"em_ne3_3869",  "Neon III -- NeIII",                {0.38688},          {1.0}      });
    add_line({"em_hdelta",    "Hydrogen delta",                   {0.41017},          {1.0}      });
    add_line({"em_hgamma",    "Hydrogen gamma",                   {0.43405},          {1.0}      });
    add_line({"em_hbeta",     "Hydrogen beta",                    {0.48613},          {1.0}      });
    add_line({"em_o3_5007",   "Oxygen III -- [OIII] (doublet)",   {0.50068, 0.49589}, {1.0, 0.3} });
    add_line({"em_halpha",    "Hydrogen alpha",                   {0.65628},          {1.0}      });
    add_line({"em_n2_6583",   "Nitrogen II -- [NII] (doublet)",   {0.65835, 0.65480}, {1.0, 0.3} });
    add_line({"em_s2_6717",   "Sulfur II -- [SII] (doublet)",     {0.67164, 0.67308}, {1.0, 0.75}});
    add_line({"em_palpha",    "Pashen alpha",                     {1.87513},          {1.0}      });

    // Standard galaxy/AGN FIR emission lines
    add_line({"em_o1_145",    "Oxygen I --- [OI]",                {145.525},          {1.0}      });
    add_line({"em_c2_157",    "Carbon II --- [CII]",              {157.741},          {1.0}      });
    add_line({"em_c1_370",    "Carbon I -- CI",                   {370.42},           {1.0}      });
    add_line({"em_c1_609",    "Carbon I -- CI",                   {609.14},           {1.0}      });
    add_line({"em_co98",      "Carbon Monoxyde -- CO(9-8)",       {298.12},           {1.0}      });
    add_line({"em_co87",      "Carbon Monoxyde -- CO(8-7)",       {325.23},           {1.0}      });
    add_line({"em_co76",      "Carbon Monoxyde -- CO(7-6)",       {371.65},           {1.0}      });
    add_line({"em_co65",      "Carbon Monoxyde -- CO(6-5)",       {433.57},           {1.0}      });
    add_line({"em_co54",      "Carbon Monoxyde -- CO(5-4)",       {520.23},           {1.0}      });
    add_line({"em_co43",      "Carbon Monoxyde -- CO(4-3)",       {650.25},           {1.0}      });
    add_line({"em_co32",      "Carbon Monoxyde -- CO(3-2)",       {866.96},           {1.0}      });
    add_line({"em_co21",      "Carbon Monoxyde -- CO(2-1)",       {1300.40},          {1.0}      });
    add_line({"em_co10",      "Carbon Monoxyde -- CO(1-0)",       {2600.75},          {1.0}      });

    // Standard UV absorption lines (from list in Leitherer+11)
    add_line({"abs_c3_1176",  "Carbon III -- CIII",               {0.11755},          {1.0}      });
    add_line({"abs_si2_1260", "Silicium II -- SiII",              {0.12604},          {1.0}      });
    add_line({"abs_o1_1302",  "Oxygen I -- OI",                   {0.13022},          {1.0}      });
    add_line({"abs_si2_1304", "Silicium II -- SiII",              {0.13043},          {1.0}      });
    add_line({"abs_c2_1335",  "Carbon II -- CII",                 {0.13345},          {1.0}      });
    add_line({"abs_o4_1342",  "Oxygen IV -- OIV",                 {0.13416},          {1.0}      });
    add_line({"abs_s5_1502",  "Sulfur V -- SV",                   {0.15018},          {1.0}      });
    add_line({"abs_si2_1526", "Silicium II -- SiII",              {0.15267},          {1.0}      });
    add_line({"abs_fe2_1608", "Iron II -- FeII (doublet)",        {0.16085, 0.16112}, {1.0, 1.0} });
    add_line({"abs_al2_1671", "Aluminum II  -- AlII",             {0.16708},          {1.0}      });
    add_line({"abs_al3_1855", "Aluminum III -- AlIII (doublet)",  {0.18547, 0.18628}, {1.0, 1.0} });
    add_line({"abs_fe2_2344", "Iron II -- FeII",                  {0.23442},          {1.0}      });
    add_line({"abs_fe2_2380", "Iron II -- FeII (doublet)",        {0.23745, 0.23828}, {1.0, 1.0} });
    add_line({"abs_fe2_2600", "Iron II -- FeII (doublet)",        {0.25867, 0.26002}, {1.0, 1.0} });
    // ---------------------------------------------------------

    // Read command line arguments
    if (argc < 2) {
        print_help(linedb);
        return 0;
    }

    double z0 = dnan;
    double dz = 0.3;
    bool forbid_absorption = true;
    bool allow_offsets = false;
    double offset_snr_min = 5.0;
    double width_min = 50.0;
    double width_max = 500.0;
    double offset_max = 1000.0;
    double delta_width = 0.2;
    double delta_z = 0.2;
    double delta_offset = 0.2;
    bool use_global_chi2 = false;
    bool full_range = false;
    bool local_continuum = false;
    double local_continuum_width = 8000.0;
    bool fit_continuum_template = false;
    bool residual_rescale = false;
    std::string template_dir = "templates/";
    double fix_width = dnan;
    uint_t lambda_pad = 5;
    bool same_width = false;
    bool use_mpfit = false;
    bool mc_errors = false;
    uint_t num_mc = 200;
    uint_t nthread = 1;
    bool verbose = false;
    bool save_model = false;
    std::string outdir;
    bool ascii = false;
    std::string tcosmo = "std";
    uint_t tseed = 42;
    uint_t flux_hdu = 1;
    uint_t error_hdu = 2;
    vec1s tlines;

    read_args(argc-1, argv+1, arg_list(z0, dz, name(tlines, "lines"), width_min, width_max,
        verbose, same_width, save_model, fix_width, use_mpfit, ascii, outdir, delta_width,
        delta_z, lambda_pad, local_continuum, local_continuum_width, flux_hdu, error_hdu,
        fit_continuum_template, name(tcosmo, "cosmo"), template_dir, use_global_chi2,
        allow_offsets, offset_max, offset_snr_min, delta_offset, residual_rescale,
        mc_errors, num_mc, name(tseed, "seed"), name(nthread, "threads"), full_range,
        forbid_absorption
    ));

    // Setup cosmological parameters
    auto cosmo = get_cosmo(tcosmo);

    // Check validity of input and create output directories if needed
    if (!outdir.empty()) {
        outdir = file::directorize(outdir);
        file::mkdir(outdir);
    }

    if (!template_dir.empty()) {
        template_dir = file::directorize(template_dir);
    }

    if (is_finite(fix_width) && use_mpfit) {
        warning("'use_mpfit' is disabled when 'fix_width' is set");
        note("will use the linear solver instead (which is faster and more robust)");
        use_mpfit = false;
    }

    if (local_continuum && fit_continuum_template) {
        warning("'local_continuum' has no effect if 'fit_continuum_template' is set");
        local_continuum = false;
    }

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
            nl.pretty_name = "added through command line";

            vec1d nums;
            if (count(!from_string(spl[1-_], nums)) != 0) {
                error("could not convert line wavelengths and ratios in '", l, "' into a "
                    "list of numbers");
                return 1;
            }

            if (nums.size() == 1) {
                nl.lambda = nums;
                nl.ratio = replicate(1.0, nums.size());
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

    // Update continuum ranges
    for (auto& l : lines) {
        l.define_continuum(local_continuum_width);
    }

    // Read spectrum
    if (verbose) note("read input spectrum...");
    std::string infile = argv[1];
    vec1d flx, err;
    fits::input_image fimg(infile);
    fimg.reach_hdu(flux_hdu);
    fimg.read(flx);
    fimg.reach_hdu(error_hdu);
    fimg.read(err);

    // Come back to flux HDU to make sure we read the WCS from there
    fimg.reach_hdu(flux_hdu);

    // Build wavelength axis
    uint_t nlam = flx.size();
    double cdelt = 1, crpix = 1, crval = 1;
    if (!fimg.read_keyword("CDELT1", cdelt) ||
        !fimg.read_keyword("CRPIX1", crpix) ||
        !fimg.read_keyword("CRVAL1", crval)) {
        error("could not read WCS information for wavelength axis");
        return 1;
    }

    // Handle wavelength axis units
    bool frequency = false;
    std::string cunit; {
        if (!fimg.read_keyword("CUNIT1", cunit)) {
            warning("could not find unit of wavelength axis");
            note("assuming wavelengths are given in microns");
        } else {
            cunit = tolower(cunit);
            double conv = 1.0;
            if (cunit == "angstrom") {
                conv = 1e-4;
            } else if (cunit == "nm") {
                conv = 1e-3;
            } else if (cunit == "um" || cunit == "micron") {
                conv = 1.0;
            } else if (cunit == "mm") {
                conv = 1e3;
            } else if (cunit == "cm") {
                conv = 1e4;
            } else if (cunit == "m") {
                conv = 1e6;
            } else if (cunit == "hz") {
                frequency = true;
                conv = 1.0;
            } else if (cunit == "khz") {
                frequency = true;
                conv = 1e3;
            } else if (cunit == "mhz") {
                frequency = true;
                conv = 1e6;
            } else if (cunit == "ghz") {
                frequency = true;
                conv = 1e9;
            } else {
                error("unrecognized wavelength/frequency unit '", cunit, "'");
                return 1;
            }

            crval *= conv;
            cdelt *= conv;
        }
    }

    vec1d lam, laml, lamu;
    if (!frequency) {
        lam = crval + cdelt*(findgen(nlam) + (1 - crpix));
        laml = lam - 0.5*cdelt;
        lamu = lam + 0.5*cdelt;
    } else {
        // x-axis is frequency, convert that to a wavelength
        // and do not forget to reverse the data so that wavelenths are
        // strictly increasing
        vec1d freq = crval + cdelt*(findgen(nlam) + (1 - crpix));
        vec1d freql = freq - 0.5*cdelt;
        vec1d frequ = freq + 0.5*cdelt;
        lam = reverse(1e6*2.99792e8/freq);
        laml = reverse(1e6*2.99792e8/frequ);
        lamu = reverse(1e6*2.99792e8/freql);
        flx = reverse(flx);
        err = reverse(err);
        cdelt = median(lamu-laml);
        if (verbose) {
            note("input spectrum in frequency unit, converting to wavelength");
            note("new coverage: ", lam.front(), " to ", lam.back(), " microns (average cdelt = ", cdelt, ")");
        }
    }

    uint_t orig_nlam = lam.size();

    // Identify good regions of the spectrum
    vec1b goodspec = is_finite(flx) && is_finite(err) && err > 0;
    if (count(goodspec) <= lambda_pad*2) {
        error("the spectrum does not contain any valid point");
        return 1;
    }

    // Flag out the pixels at the border of the spectrum
    vec1b goodspec_flag = goodspec;
    goodspec = keep_gaps_and_expand(goodspec, 10, lambda_pad);

    // Identify regions of the spectrum that are NaN and give them zero weight
    {
        vec1u id_flagged = where(goodspec && !goodspec_flag);
        flx[id_flagged] = 0;
        err[id_flagged] = 1e20*median(err[where(goodspec_flag)]);
    }

    // Select a wavelength domain centered on the line(s)
    const vec1u idl = (full_range ? where(goodspec) : where(goodspec &&
        lam > lambda_min*(1.0+z0-2*dz) && lam < lambda_max*(1.0+z0+2*dz)));

    if (idl.empty()) {
        error("none of the chosen lines are covered by the provided spectrum at z=", z0, " +/- ", dz);
        vec1u tidl = where(goodspec);
        note("your redshift search for these lines requires a range within ",
            lambda_min*(1.0+z0-2*dz), " to ", lambda_max*(1.0+z0+2*dz));
        note("but the spectrum only covers ", min(lam[tidl]), " to ", max(lam[tidl]));
        return 1;
    }

    flx = flx[idl];
    err = err[idl];
    lam = lam[idl];
    laml = laml[idl];
    lamu = lamu[idl];
    goodspec_flag = goodspec_flag[idl];

    // Define redshift grid so as to have the requested number of samples per wavelength element
    double tdz = delta_z*cdelt/mean(lam);
    uint_t nz = 2*ceil(dz/tdz)+1;
    vec1d z_grid = rgen(z0-dz, z0+dz, nz);

    // Define width grid so as to have the requested number of samples per wavelength element
    double dwidth = 3e5*delta_width*cdelt/mean(lam);
    uint_t nwidth = ceil((width_max - width_min)/dwidth);
    vec1d ws = rgen(width_min, width_max, nwidth);

    // Define offset grid so as to have the requested number of samples per wavelength element
    double doffset = 3e5*delta_offset*cdelt/mean(lam);
    uint_t noffset = 2*ceil(offset_max/doffset)+1;
    vec1d os = rgen(-offset_max, offset_max, noffset);

    // Find which lines are covered by this spectrum
    vec1b lcov(lines.size());
    vec1b scov(lam.size());
    for (uint_t il : range(lines)) {
        vec1b tcov = replicate(false, lam.size());
        for (uint_t isl : range(lines[il].lambda)) {
            tcov = tcov || abs(lam/(lines[il].lambda[isl]*(1.0 + z0)) - 1.0) < 2*width_max/3e5 + dz;
        }

        scov = scov || tcov;
        lcov[il] = count(tcov) > 0;
    }

    auto all_lines = lines;
    vec1u id_linecov = where(lcov);
    lines = lines[id_linecov];

    vec1u id_chi2 = where(scov && goodspec_flag);

    if (verbose) {
        note("fitting ", flx.dims[0], " spectral elements between ",
           lam.front(), " and ", lam.back(), " um...");
        note("with ", lines.size(), "/", all_lines.size(), " lines from the database");
    }

    // Load template library
    struct template_t {
        vec1d lam, sed; // intrinsic spectrum of the template
    };

    vec<1,template_t> templates;
    if (fit_continuum_template) {
        vec1s tpls = file::list_files(template_dir+"*.dat");
        for (uint_t it : range(tpls)) {
            std::string fname = template_dir+tpls[it];
            template_t tp;
            ascii::read_table(fname, ascii::find_skip(fname), tp.lam, tp.sed);
            tp.lam *= 1e-4;
            templates.push_back(tp);
        }

        if (verbose) {
            note("and ", templates.size(), " galaxy templates for the continuum level");
        }
    }

    // Perform a redshift search
    if (verbose) {
        note("average resolution of spectrum: ", mean(lam)/cdelt);
        std::string grid_msg = "redshift search ("+strn(nz)+
            " redshifts, step = "+strn(float(z_grid[1] - z_grid[0]));

        if (!use_mpfit && !is_finite(fix_width)) {
            grid_msg += ", "+strn(nwidth)+" line widths, step = "+strn(float(ws[1] - ws[0]));
        }
        if (allow_offsets) {
            grid_msg += ", "+strn(noffset)+" line offsets, step = "+strn(float(os[1] - os[0]));
        }

        grid_msg += ")...";
        note(grid_msg);
    }

    // Renormalize flux and errors in units of 1e-17 erg/s/cm2/A to avoid
    // numerical imprecision (it is always better to deal with numbers close to unity).
    // This helps prevent mpfit getting stuck and not varying the line widths.
    if (!frequency) {
        flx *= 1e17; err *= 1e17;
    }

    // Function to perform the fit!
    // -----------------------------------------------------------------
    auto dofit = [=](fit_result_t& fres, const vec1d& tflx, const vec1d& terr, bool silent = false) {
        // Initialize starting conditions and specify fitting constraints
        // --------------------------------------------------------------
        vec1d p;       // parameter array
        vec1b p_fixed; // do we vary this parameter?
        vec1d p_min;   // do we have a lower limit on this parameter?
        vec1d p_max;   // do we have an upper limit on this parameter?

        vec1u id_amp;    // location in 'p' of amplitude variables
        vec1u id_width;  // location in 'p' of width variables
        vec1u id_offset; // location in 'p' of line offset variables
        vec1u id_cont;   // location in 'p' of continuum variables
        uint_t nparam;   // number of fit parameters in 'p'

        vec1u id_mline;  // location in 'm' of line models (only for linear fit)
        vec1u id_mcont;  // location in 'm' of continuum models (only for linear fit)
        uint_t nmodel;   // number of independent models in 'm' (only for linear fit)

        // Line amplitudes
        id_mline = uindgen(lines.size());
        id_amp = uindgen(lines.size());
        nparam = lines.size();
        nmodel = lines.size();
        append(p,       replicate(0.0,   lines.size()));
        append(p_fixed, replicate(false, lines.size()));
        append(p_max,   replicate(dnan,  lines.size()));
        if (forbid_absorption) {
            append(p_min, replicate(0, lines.size()));
        } else {
            append(p_min, replicate(dnan, lines.size()));
        }

        // Line widths
        if (same_width) {
            id_width = replicate(nparam, lines.size());
            nparam += 1;
            append(p,       replicate(0.0,   1u));
            append(p_fixed, replicate(false, 1u));
            append(p_min,   replicate(dnan,  1u));
            append(p_max,   replicate(dnan,  1u));
        } else {
            id_width = nparam + uindgen(lines.size());
            nparam += lines.size();
            append(p,       replicate(0.0,   lines.size()));
            append(p_fixed, replicate(false, lines.size()));
            append(p_min,   replicate(dnan,  lines.size()));
            append(p_max,   replicate(dnan,  lines.size()));
        }

        // Line offsets
        id_offset = nparam + uindgen(lines.size());
        nparam += lines.size();
        append(p,       replicate(0.0,         lines.size()));
        append(p_fixed, replicate(true,        lines.size())); // fixed at first
        append(p_min,   replicate(-offset_max, lines.size()));
        append(p_max,   replicate(offset_max,  lines.size()));

        // Continuum level
        if (local_continuum) {
            id_cont = nparam + uindgen(lines.size());
            id_mcont = nmodel + uindgen(lines.size());
            nparam += lines.size();
            nmodel += lines.size();
            append(p,       replicate(0.0,   lines.size()));
            append(p_fixed, replicate(false, lines.size()));
            append(p_min,   replicate(dnan,  lines.size()));
            append(p_max,   replicate(dnan,  lines.size()));
        } else if (fit_continuum_template) {
            id_cont = nparam + uindgen(templates.size());
            id_mcont = nmodel + uindgen(templates.size());
            nparam += templates.size();
            nmodel += templates.size();
            append(p,       replicate(0.0,   templates.size()));
            append(p_fixed, replicate(false, templates.size()));
            append(p_min,   replicate(dnan,  templates.size()));
            append(p_max,   replicate(dnan,  templates.size()));
        }

        // Cached computations
        // -------------------
        // Regions of the spectrum defined as continuum around lines
        vec<1,vec1u> line_idc(lines.size());
        vec<1,vec1d> tpl_flux(templates.size());

        // Function to perform a non-linear fit (if varying the line widths)
        // -----------------------------------------------------------------
        auto try_nlfit = [&](double tz, double& gchi2) {
            struct lam_t {
                lam_t(vec1d tl, vec1d tu) : ll(std::move(tl)), lu(std::move(tu)) {}
                vec1d ll, lu;
            };

            auto model = [&](const lam_t& l, const vec1d& tp) {
                vec1d m(l.ll.dims);
                for (uint_t il : range(lines)) {
                    if (local_continuum) {
                        m[line_idc[il]] += tp[id_cont[il]];
                    }

                    for (uint_t isl : range(lines[il].lambda)) {
                        double lw = (tp[id_width[il]]/3e5)*lines[il].lambda[isl]*(1.0+tz);
                        double l0 = lines[il].lambda[isl]*(1.0+tz)*(1.0+tp[id_offset[il]]/3e5);
                        double amp = 1e-4*lines[il].ratio[isl]*tp[id_amp[il]];
                        m += gauss_integral(l.ll, l.lu, l0, lw, amp);
                    }
                }

                if (fit_continuum_template) {
                    for (uint_t it : range(templates)) {
                        m += tp[id_cont[it]]*tpl_flux[it];
                    }
                }

                return m;
            };

            mpfit_options opts(p.size());
            opts.upper_limit = p_max;
            opts.lower_limit = p_min;
            opts.frozen      = p_fixed;

            lam_t lt{laml, lamu};
            mpfit_result res = mpfitfun(tflx, terr, lt, model, p, opts);

            if (!use_global_chi2) {
                // Compute local chi2 (only counting pixels around the lines, not the continuum)
                vec1d tmodel = model(lam_t{laml[id_chi2], lamu[id_chi2]}, res.params);
                res.chi2 = total(sqr((tflx[id_chi2] - tmodel)/terr[id_chi2]));
            }

            bool better = res.success && res.chi2 < fres.chi2;
            if (better) {
                fres.chi2 = res.chi2;
                fres.z = tz;
                fres.flux        = res.params[id_amp];
                fres.flux_err    = res.errors[id_amp];
                fres.width       = res.params[id_width];
                fres.width_err   = res.errors[id_width];
                fres.offset      = res.params[id_offset];
                fres.offset_err  = res.errors[id_offset];

                fres.model = model(lt, res.params);
                vec1d tp = res.params; tp[id_amp] = 0;
                fres.model_continuum = model(lt, tp);
            }

            if (res.success && res.chi2 < gchi2) {
                gchi2 = res.chi2;
            }

            return better;
        };
        // -----------------------------------------------------------------

        // Function to perform a linear fit (fixed line widths)
        // -----------------------------------------------------------------
        auto try_lfit = [&](double tz, double& gchi2) {
            vec2d m(nmodel, lam.size());
            for (uint_t il : range(lines)) {
                if (local_continuum) {
                    m(id_mcont[il],line_idc[il]) = 1.0;
                }

                for (uint_t isl : range(lines[il].lambda)) {
                    double lw = (p[id_width[il]]/3e5)*lines[il].lambda[isl]*(1.0+tz);
                    double l0 = lines[il].lambda[isl]*(1.0+tz)*(1.0+p[id_offset[il]]/3e5);
                    double amp = 1e-4*lines[il].ratio[isl];
                    m(id_mline[il],_) += gauss_integral(laml, lamu, l0, lw, amp);
                }
            }

            if (fit_continuum_template) {
                for (uint_t it : range(templates)) {
                    m(id_mcont[it],_) = tpl_flux[it];
                }
            }

            // Exclude models that are zero from the fit
            vec1u idne = where(partial_total(1, sqr(m)) > 0.0);
            m = m(idne,_);

            linfit_result res = linfit_pack(tflx, terr, m);

            // If asked, ignore lines that would be fitted in absorption
            if (forbid_absorption) {
                vec1u id1, id2;
                match(idne, id_amp, id1, id2);
                res.params[id1[where(res.params[id1] < 0)]] = 0;
                res.errors[id1[where(res.params[id1] < 0)]] = 0;

                // Recompute the chi2 if using the global chi2
                // Else it will be recomputed just below
                if (use_global_chi2) {
                    vec1d tmodel(tflx.size());
                    for (uint_t il : range(idne)) {
                        tmodel += res.params[il]*m(il,_);
                    }

                    res.chi2 = total(sqr((tflx - tmodel)/terr));
                }
            }

            if (!use_global_chi2) {
                // Compute local chi2 (only counting pixels around the lines, not the continuum)
                vec1d tmodel(id_chi2.size());
                for (uint_t il : range(idne)) {
                    tmodel += res.params[il]*m(il,id_chi2);
                }

                res.chi2 = total(sqr((tflx[id_chi2] - tmodel)/terr[id_chi2]));
            }

            bool better = res.success && res.chi2 < fres.chi2;
            if (better) {
                fres.chi2 = res.chi2;
                fres.z = tz;

                // Reform initial model grid
                vec1d rflux = replicate(dnan, nmodel);
                vec1d rerr = replicate(dnan, nmodel);
                rflux[idne] = res.params;
                rerr[idne] = res.errors;

                // Get fluxes from the lines
                fres.flux       = rflux[id_mline];
                fres.flux_err   = rerr[id_mline];

                // Width and offset are given by current grid
                fres.width      = p[id_width];
                fres.width_err  = replicate(0.0, id_amp.size());
                fres.offset     = p[id_offset];
                fres.offset_err = replicate(0.0, id_amp.size());

                fres.model = vec1d(lam.dims);
                fres.model_continuum = vec1d(lam.dims);
                vec1b isline = replicate(false, nmodel);
                isline[id_mline] = true;
                for (uint_t il : range(idne)) {
                    fres.model += res.params[il]*m(il,_);
                    if (!isline[idne[il]]) {
                        fres.model_continuum += res.params[il]*m(il,_);
                    }
                }
            }

            if (res.success && res.chi2 < gchi2) {
                gchi2 = res.chi2;
            }

            return better;
        };
        // -----------------------------------------------------------------

        // Function to build cache variables used in fits for a given redshift
        // -----------------------------------------------------------------
        auto make_cache = [&](double tz) {
            // Cache continuum wavelength regions for each line
            for (uint_t il : range(lines)) {
                line_idc[il] = where(
                    lam >= (1.0 + tz)*lines[il].continuum_lmin &&
                    lam <= (1.0 + tz)*lines[il].continuum_lmax
                );
            }

            // Redshift and rebin galaxy templates to the resolution of the spectrum
            for (uint_t it : range(templates)) {
                auto& tpl = templates[it];
                auto& ttflx = tpl_flux[it];

                vec1d tlam = tpl.lam*(1.0 + tz);
                ttflx = interpolate(tpl.sed, tlam, lam);
                ttflx[where(!is_finite(ttflx) || lam > max(tlam) || lam < min(tlam))] = 0.0;
            }
        };
        // -----------------------------------------------------------------

        auto pg = progress_start(z_grid.size());
        for (uint_t iz : range(z_grid)) {
            // Build caches
            make_cache(z_grid[iz]);

            // Estimate starting parameters and/or fix the parameter values
            // and do the fit
            if (!use_mpfit) {
                // Use best offset to start with
                p[id_offset] = fres.offset;

                if (is_finite(fix_width)) {
                    // Fixed width for all lines
                    fres.free_width[_] = false;
                    p[id_width] = fix_width;
                    try_lfit(z_grid[iz], fres.chi2_grid[iz]);
                } else if (same_width) {
                    // All lines have same width, let it vary
                    fres.free_width[_] = true;
                    for (uint_t iw : range(ws)) {
                        p[id_width] = ws[iw];
                        try_lfit(z_grid[iz], fres.chi2_grid[iz]);
                    }
                } else {
                    // Each line can have a different width, let them vary one by one
                    fres.free_width[_] = true;
                    p[id_width] = 0.5*(width_max - width_min);

                    for (uint_t il : range(lines)) {
                        for (uint_t iw : range(ws)) {
                            p[id_width[il]] = ws[iw];
                            try_lfit(z_grid[iz], fres.chi2_grid[iz]);
                        }

                        p[id_width[il]] = fres.width[il];
                    }
                }
            } else {
                // Each line has its own width, guess start value from allowed range
                fres.free_width[_] = true;
                p[id_width] = 0.5*(width_max - width_min);
                // Set offset to current best fit
                p[id_offset] = fres.offset;

                // Guess amplitude by integrating the spectrum over the guessed width
                for (uint_t il : range(lines)) {
                    double lw = (p[id_width[il]]/3e5)*lines[il].lambda[0]*(1.0+z_grid[iz]);
                    p[id_amp[il]] = interpolate(tflx, lam, lines[il].lambda[0]*(1.0 + z_grid[iz]))*
                        sqrt(2.0*dpi)*lw*1e4;
                }

                // Background level is set to zero at first
                if (local_continuum || fit_continuum_template) {
                    p[id_cont] = 0.0;
                }

                try_nlfit(z_grid[iz], fres.chi2_grid[iz]);
            }

            if (verbose && !silent) progress(pg);
        }

        // Adjust line offsets if asked
        if (allow_offsets) {
            // Rebuild cache for best fit redshift
            make_cache(fres.z);

            if (!use_mpfit) {
                // Start from current best fit offsets
                p[id_offset] = fres.offset;

                // Setup grid from best fit
                p[id_width] = fres.width;

                // Vary the offset of each line
                for (uint_t il : range(lines)) {
                    // Only fit offsets for high SNR lines of for Lyalpha
                    if (abs((fres.flux/fres.flux_err)[il]) < offset_snr_min && lines[il].name != "em_lyalpha") continue;
                    fres.free_offset[il] = true;

                    for (uint_t io : range(os)) {
                        p[id_offset[il]] = os[io];
                        try_lfit(fres.z, fres.chi2);
                    }

                    p[id_offset[il]] = fres.offset[il];
                }
            } else {
                // Set starting parameters from best fit
                p[id_amp] = fres.flux;
                p[id_width] = fres.width;
                p[id_offset] = fres.offset;
                if (local_continuum || fit_continuum_template) {
                    p[id_cont] = 0; // Set this back to zero
                }
                // Fix widths
                p_fixed[id_width] = true;
                // Free the offsets of high SNR lines
                p_fixed[id_offset] = abs(fres.flux/fres.flux_err) > offset_snr_min;
                fres.free_offset = !p_fixed[id_offset];
                // Make sure Lyalpha is allowed all the time
                for (uint_t il : range(lines)) {
                    if (lines[il].name == "em_lyalpha") p_fixed[id_offset[il]] = false;
                }

                // Refit
                try_nlfit(fres.z, fres.chi2);
            }
        } else {
            fres.free_offset[_] = false;
        }

        // Compute centroid wavelength of each line
        for (uint_t il : range(lines)) {
            double ddz = (is_finite(fres.offset[il]) ? fres.offset[il]/3e5 : 0.0);
            fres.lambda[il] = (1.0 + fres.z)*(1.0 + ddz)*lines[il].lambda[0];
        }

        // Compute equivalent widths and continuum fluxes
        if (fit_continuum_template || local_continuum) {
            // Make a new continuum model which is better behaved (strictly positive, no invalid data)
            vec1d ew_cont = fres.model_continuum;
            ew_cont[where(ew_cont < 0)] = 0;
            vec1u idf = where(is_finite(ew_cont));
            vec1u idnf = where(!is_finite(ew_cont));
            ew_cont[idnf] = interpolate(ew_cont[idf], lam[idf], lam[idnf]);

            // Compute EW for each line separately
            for (uint_t il : range(fres.ew)) {
                if (!is_finite(fres.flux[il])) continue;

                double ddz = (is_finite(fres.offset[il]) ? fres.offset[il]/3e5 : 0.0);

                double l0;
                double l1 = lines[il].continuum_lmax*(fres.lambda[il]/lines[il].lambda[0]);
                if (lines[il].name == "em_lyalpha") {
                    // Treat Lyalpha separately: the EW is determined redward of the line only
                    l0 = fres.lambda[il];
                } else {
                    // For all other lines, the EW is determined symmetrically around the line
                    l0 = lines[il].continuum_lmin*(fres.lambda[il]/lines[il].lambda[0]);
                }

                vec1u idc = where(lam >= l0 && lam <= l1);
                fres.cont[il] = mean(ew_cont[idc]);
                fres.cont_err[il] = sqrt(total(sqr(terr[idc])))/idc.size();
                fres.ew[il] = fres.flux[il]/fres.cont[il]/(1.0 + fres.z);
                fres.ew_err[il] = sqrt(sqr(fres.flux_err[il]/fres.cont[il]) +
                    sqr(fres.cont_err[il]*fres.flux[il]/sqr(fres.cont[il])))/(1.0 + fres.z);
            }
        }
    };
    // -----------------------------------------------------------------

    // Do the fit
    fit_result_t best_fit(z_grid.size(), lines.size());
    dofit(best_fit, flx, err);

    // Define whether or not we should do a second pass
    auto should_second_pass = [&](const fit_result_t& fres) {
        return allow_offsets && count(fres.free_offset) > 0;
    };

    bool force_second_pass = false;

    // Rescale uncertainties based on the residuals, if asked
    if (residual_rescale) {
        // Make residual spectrum
        vec1d resi = (flx - best_fit.model)/err;

        vec1d rescale = replicate(dnan, lines.size());
        vec1d rlam(lines.size());
        for (uint_t il : range(lines)) {
            // Find spectral elements of the continuum around that line
            vec1u idg = where(
                lam >= (best_fit.lambda[il]/lines[il].lambda[0])*lines[il].continuum_lmin &&
                lam <= (best_fit.lambda[il]/lines[il].lambda[0])*lines[il].continuum_lmax &&
                is_finite(resi)
            );

            if (idg.size() < 10) continue; // not enough valid data points

            rescale[il] = max(1.0, 1.48*mad(resi[idg]));
            best_fit.flux_err[il]   *= rescale[il];
            best_fit.width_err[il]  *= rescale[il];
            best_fit.offset_err[il] *= rescale[il];
            best_fit.lambda_err[il] *= rescale[il];
            best_fit.ew_err[il]     *= rescale[il];
            best_fit.cont_err[il]   *= rescale[il];
            rlam[il] = mean(lines[il].lambda);
        }

        // For lines that had too few valid data points to estimate RMS, assume the worst
        vec1u idff = where(is_finite(rescale));
        if (!idff.empty()) {
            rescale[where(!is_finite(rescale))] = max(rescale[idff]);

            // Rescale whole error spectrum
            vec1u ids = sort(rlam); rlam = rlam[ids]; rescale = rescale[ids];
            vec1d rsc = interpolate(rescale, rlam, lam);
            rsc[where(lam < min(rlam))] = rescale.front();
            rsc[where(lam > max(rlam))] = rescale.back();
            err *= rsc;

            if (verbose) {
                note("rescaling uncertainties from residual (min=", min(rescale), ", max=", max(rescale), ")");
                note("fitting again...");
            }

            // Remove offsets for lines that went below the SNR threshold
            best_fit.offset[where(abs(best_fit.flux/best_fit.flux_err) < offset_snr_min)] = 0;
            best_fit.offset_err[where(abs(best_fit.flux/best_fit.flux_err) < offset_snr_min)] = 0;

            // Enable second pass fitting
            force_second_pass = true;
        }
    }

    if (should_second_pass(best_fit) || force_second_pass) {
        // Reset chi2
        best_fit.chi2 = dinf;
        best_fit.chi2_grid[_] = dinf;

        // Fit again...
        dofit(best_fit, flx, err);
    }

    // If asked, randomly perturb the spectrum within the uncertainty to derive reliable
    // error bars (this can take some time!)
    if (mc_errors) {
        std::vector<fit_result_t> mc_fits;
        for (uint_t m : range(num_mc)) {
            mc_fits.emplace_back(z_grid.size(), lines.size());
        }

        if (verbose) {
            note("fitting ", num_mc, " MonteCarlo realizations of noise to get uncertainties...");
        }

        auto do_mc = [=,&mc_fits](uint_t m) {
            auto& mc_bfit = mc_fits[m];
            auto seed = make_seed(tseed + m);

            // Perturb the flux spectrum
            vec1d tflx = flx + err*randomn(seed, flx.size());

            // Do the fit
            dofit(mc_bfit, tflx, err, true);

            if (should_second_pass(mc_bfit)) {
                // Reset chi2
                mc_bfit.chi2 = dinf;
                mc_bfit.chi2_grid[_] = dinf;

                // Fit again...
                dofit(mc_bfit, tflx, err, true);
            }
        };

        if (nthread <= 1) {
            auto pg = progress_start(num_mc);
            for (uint_t m : range(num_mc)) {
                do_mc(m);
                if (verbose) progress(pg);
            }
        } else {
            std::atomic<uint_t> iter(0);
            auto pool = thread::pool(nthread);
            uint_t m0 = 0;
            uint_t dm = num_mc/nthread;

            // Start all threads
            for (uint_t it : range(pool)) {
                uint_t m1 = (it == nthread-1 ? num_mc : m0 + dm);

                pool[it].start([&,m0,m1]() {
                    for (uint_t m : range(m0, m1)) {
                        do_mc(m);
                        ++iter;
                    }
                });

                m0 = m1;
            }

            // Wait for them to finish
            auto pg = progress_start(num_mc);
            while (iter < num_mc) {
                thread::sleep_for(0.2);
                if (verbose) print_progress(pg, iter);
            }

            // By now, all threads should have ended their tasks.
            // We must ask them to terminate nicely.
            for (auto& t : pool) {
                t.join();
            }
        }

        // Collect all fits and derive uncertainties from stddev
        for (uint_t il : range(lines)) {
            best_fit.lambda_err[il] = get_mc_error(num_mc, [&](uint_t m) { return mc_fits[m].lambda[il]; });
            best_fit.flux_err[il]   = get_mc_error(num_mc, [&](uint_t m) { return mc_fits[m].flux[il];   });
            best_fit.width_err[il]  = get_mc_error(num_mc, [&](uint_t m) { return mc_fits[m].width[il];  });
            best_fit.offset_err[il] = get_mc_error(num_mc, [&](uint_t m) { return mc_fits[m].offset[il]; });
            best_fit.cont_err[il]   = get_mc_error(num_mc, [&](uint_t m) { return mc_fits[m].cont[il];   });
            best_fit.ew_err[il]     = get_mc_error(num_mc, [&](uint_t m) { return mc_fits[m].ew[il];     });

            // If uncertainty is zero it means our grid is too coarse, just use that as a floor on uncertainty
            if (best_fit.width_err[il] == 0.0  && best_fit.free_width[il]) {
                best_fit.width_err[il]  = 0.5*(ws[1]-ws[0]);
            }
            if (best_fit.offset_err[il] == 0.0 && best_fit.free_offset[il]) {
                best_fit.offset_err[il] = 0.5*(os[1]-os[0]);
            }
        }

        best_fit.z_err = get_mc_error(num_mc, [&](uint_t m) { return mc_fits[m].z; });
    }

    // Compute number of degrees of freedom for reduced chi2
    uint_t ndof = (use_global_chi2 ? flx.size() : id_chi2.size());
    ndof -= lines.size();
    if (!is_finite(fix_width)) {
        if (same_width) ndof -= 1;
        else            ndof -= lines.size();
    }
    if (fit_continuum_template) {
        ndof -= templates.size();
    } else if (local_continuum) {
        ndof -= lines.size();
    }

    // Compute reduced chi2 grid
    best_fit.chi2_grid /= ndof;
    best_fit.chi2      /= ndof;
    vec1d pz = exp(-(best_fit.chi2_grid - best_fit.chi2));
    pz -= min(pz);
    pz /= integrate(z_grid, pz);

    // Estimate redshift uncertainties
    double zlow, zup;
    if (mc_errors) {
        zlow = best_fit.z_err;
        zup  = best_fit.z_err;
    } else {
        vec1d cpz = cumul(z_grid, pz); cpz /= max(cpz);
        zlow = best_fit.z - interpolate(z_grid, cpz, 0.16);
        zup  = interpolate(z_grid, cpz, 0.84) - best_fit.z;
    }

    if (verbose) {
        if (is_finite(best_fit.z)) {
            print("best redshift: ", best_fit.z, " + ", zup, " - ", zlow,
                " (chi2: ", best_fit.chi2*ndof, ", reduced: ", best_fit.chi2, ")");
        } else {
            print("could not fit any redshift...");
        }
    }

    // Reshape the output variables to include the non-covered lines
    lines = all_lines;
    best_fit.lambda      = reshape(best_fit.lambda,      id_linecov, all_lines.size(), dnan);
    best_fit.lambda_err  = reshape(best_fit.lambda_err,  id_linecov, all_lines.size(), dnan);
    best_fit.flux        = reshape(best_fit.flux,        id_linecov, all_lines.size(), dnan);
    best_fit.flux_err    = reshape(best_fit.flux_err,    id_linecov, all_lines.size(), dnan);
    best_fit.free_width  = reshape(best_fit.free_width,  id_linecov, all_lines.size(), false);
    best_fit.width       = reshape(best_fit.width,       id_linecov, all_lines.size(), dnan);
    best_fit.width_err   = reshape(best_fit.width_err,   id_linecov, all_lines.size(), dnan);
    best_fit.free_offset = reshape(best_fit.free_offset, id_linecov, all_lines.size(), false);
    best_fit.offset      = reshape(best_fit.offset,      id_linecov, all_lines.size(), dnan);
    best_fit.offset_err  = reshape(best_fit.offset_err,  id_linecov, all_lines.size(), dnan);
    best_fit.cont        = reshape(best_fit.cont,        id_linecov, all_lines.size(), dnan);
    best_fit.cont_err    = reshape(best_fit.cont_err,    id_linecov, all_lines.size(), dnan);
    best_fit.ew          = reshape(best_fit.ew,          id_linecov, all_lines.size(), dnan);
    best_fit.ew_err      = reshape(best_fit.ew_err,      id_linecov, all_lines.size(), dnan);

    // Build fgroup, line_id
    vec1u fgroup = replicate(0u, lines.size());
    vec1u line_id = uindgen(lines.size());

    // Ungroup line groups
    for (uint_t il : range(lines)) {
        auto& l = lines[il];

        // Set default wavelength for lines that were not covered
        if (!is_finite(best_fit.lambda[il])) {
            best_fit.lambda[il] = (1.0 + best_fit.z)*lines[il].lambda[0];
        }

        if (l.lambda.size() == 1) continue;

        uint_t idg = max(fgroup)+1;
        fgroup[il] = idg;

        for (uint_t i : range(1, l.lambda.size())) {
            best_fit.lambda.push_back(     best_fit.lambda[il]*(l.lambda[i]/l.lambda[0]));
            best_fit.lambda_err.push_back( best_fit.lambda_err[il]*(l.lambda[i]/l.lambda[0]));
            best_fit.flux.push_back(       best_fit.flux[il]*l.ratio[i]);
            best_fit.flux_err.push_back(   best_fit.flux_err[il]*l.ratio[i]);
            best_fit.free_width.push_back( best_fit.free_width[il]);
            best_fit.width.push_back(      best_fit.width[il]);
            best_fit.width_err.push_back(  best_fit.width_err[il]);
            best_fit.free_offset.push_back(best_fit.free_offset[il]);
            best_fit.offset.push_back(     best_fit.offset[il]);
            best_fit.offset_err.push_back( best_fit.offset_err[il]);
            best_fit.cont.push_back(       best_fit.cont[il]);
            best_fit.cont_err.push_back(   best_fit.cont_err[il]);
            best_fit.ew.push_back(         best_fit.ew[il]*l.ratio[i]);
            best_fit.ew_err.push_back(     best_fit.ew_err[il]*l.ratio[i]);

            tlines.push_back(tlines[il]+"-"+strn(i+1));

            fgroup.push_back(idg);
            line_id.push_back(il);
        }

        tlines[il] = tlines[il]+"-1";
    }

    // Sort by wavelength
    vec1u ids = sort(best_fit.lambda);
    best_fit.lambda = best_fit.lambda[ids]; best_fit.lambda_err = best_fit.lambda_err[ids];
    best_fit.flux = best_fit.flux[ids];     best_fit.flux_err = best_fit.flux_err[ids];
    best_fit.width = best_fit.width[ids];   best_fit.width_err = best_fit.width_err[ids];
    best_fit.offset = best_fit.offset[ids]; best_fit.offset_err = best_fit.offset_err[ids];
    best_fit.cont = best_fit.cont[ids];     best_fit.cont_err = best_fit.cont_err[ids];
    best_fit.ew = best_fit.ew[ids];         best_fit.ew_err = best_fit.ew_err[ids];
    best_fit.free_width = best_fit.free_width[ids]; best_fit.free_offset = best_fit.free_offset[ids];
    tlines = tlines[ids]; fgroup = fgroup[ids]; line_id = line_id[ids];

    // Rescale fluxes and uncertainties to original units
    std::string lunit, lfunit, funit, eunit, lamname, lamname2;
    if (!frequency) {
        lamname = "lambda";
        lamname2 = "wavelength";
        lunit = "[um]";
        lfunit = "[erg/s/cm2]";
        funit = "[erg/s/cm2/A]";
        eunit = "[A]";
        best_fit.flux  *= 1e-17; best_fit.flux_err        *= 1e-17;
        best_fit.cont  *= 1e-17; best_fit.cont_err        *= 1e-17;
        best_fit.model *= 1e-17; best_fit.model_continuum *= 1e-17;
    } else {
        lamname = "nu";
        lamname2 = "frequency";
        lunit = "[GHz]";
        lfunit = "[Jy km/s]";
        funit = "[mJy]";
        eunit = "[km/s]";
        best_fit.flux     *= 2.99792e5/(1e4*best_fit.lambda);
        best_fit.flux_err *= 2.99792e5/(1e4*best_fit.lambda);
        best_fit.cont     *= 1e3;
        best_fit.cont_err *= 1e3;
        best_fit.ew       *= 2.99792e5/(1e4*best_fit.lambda);
        best_fit.ew_err   *= 2.99792e5/(1e4*best_fit.lambda);
        best_fit.lambda_err /= best_fit.lambda;
        best_fit.lambda = 2.99792e5/best_fit.lambda;
        best_fit.lambda_err *= best_fit.lambda;
    }

    // Write the result
    std::string filebase = outdir+file::remove_extension(file::get_basename(argv[1]));
    if (verbose) note("write to disk...");

    if (ascii) {
        vec1s hdr = {"line", lamname2+" "+lunit, "error", "fit group",
            "flux "+lfunit, "error", "free width?", "width [km/s]",
            "error", "free offset?", "offset [km/s]", "error",
            "cont "+funit, "error", "EW "+eunit, "error"};

        ascii::write_table_hdr(filebase+"_slfit_lines.cat", 20, hdr,
            tlines, best_fit.lambda, best_fit.lambda_err, fgroup,
            strna_sci(best_fit.flux),   strna_sci(best_fit.flux_err),
            best_fit.free_width,  strna_sci(best_fit.width),  strna_sci(best_fit.width_err),
            best_fit.free_offset, strna_sci(best_fit.offset), strna_sci(best_fit.offset_err),
            strna_sci(best_fit.cont),   strna_sci(best_fit.cont_err),
            strna_sci(best_fit.ew),     strna_sci(best_fit.ew_err)
        );

        hdr = {"redshift", "P(z)", "red.chi2"};
        ascii::write_table_hdr(filebase+"_slfit_pz.cat", 18, hdr,
            z_grid, strna_sci(pz), strna_sci(best_fit.chi2_grid)
        );
    }

    fits::output_table otbl(filebase+"_slfit.fits");
    otbl.write_columns(ftable(
        best_fit.chi2, best_fit.z, zup, zlow, fgroup,
        best_fit.flux,   best_fit.flux_err,
        best_fit.free_width,  best_fit.width,  best_fit.width_err,
        best_fit.free_offset, best_fit.offset, best_fit.offset_err,
        best_fit.cont,   best_fit.cont_err,   best_fit.ew,     best_fit.ew_err
    ));
    otbl.write_columns(
        "lines", tlines, "grid_z", z_grid, "grid_prob", pz,
        lamname, best_fit.lambda, lamname+"_err", best_fit.lambda_err,
        "grid_chi2", best_fit.chi2_grid
    );

    if (save_model && !best_fit.model.empty()) {
        // First bring back the model into the original wavelength grid
        vec1d bmodel  = replicate(dnan, orig_nlam);
        vec1d bmodelc = replicate(dnan, orig_nlam);
        bmodel[idl]  = best_fit.model;
        bmodelc[idl] = best_fit.model_continuum;

        if (frequency) {
            // Reverse back the array to frequency ordering
            bmodel = reverse(bmodel);
            bmodelc = reverse(bmodelc);
        }

        // Then save it
        fits::output_image ospec(filebase+"_slfit_model.fits");
        ospec.write(vec1d(0)); // empty primary

        ospec.reach_hdu(1);
        ospec.write(bmodel);
        ospec.write_header(fimg.read_header());
        ospec.write_keyword("MODEL", "full");
        ospec.write_keyword("BESTZ", best_fit.z);

        ospec.reach_hdu(2);
        ospec.write(bmodelc);
        ospec.write_header(fimg.read_header());
        ospec.write_keyword("MODEL", "continuum");
        ospec.write_keyword("BESTZ", best_fit.z);
    }

    return 0;
}

void print_help(const std::map<std::string,line_t>& db) {
    using namespace format;

    print("slinefit v2.0");
    print("usage: slinefit <spectrum.fits> z0=... dz=... lines=... [options]");
    print("");
    print("Main parameters:");
    paragraph("'spectrum.fits' must be a valid 1D spectrum file (ESO convention), i.e., a FITS "
        "file containing 3 extensions: the first must be empty, the "
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
    paragraph("\nNote: you can add your own lines either by modifying the source code of the "
        "program, or directly into the command line arguments. In the 'lines=[...]' "
        "parameter, you can indeed create a new line (or set of lines) with the syntax "
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
    bullet("fix_width=...", "Must be a line width (in km/s), defined as the 'sigma' in the "
        "Gaussian profile (FWHM/2.335). If provided, the program will force the line "
        "widths to be equal to this value for all the lines. This can help solving fit "
        "instability issues if some lines are very low S/N and the fitted line widths "
        "diverge to unreasonable values. The default is to let each width vary freely in "
        "the fit.");
    bullet("allow_offsets", "Set this flag if you want to allow each line to have its own "
        "velocity offset with respect to the overall redshift. The offset can be positive "
        "or negative, and is limited by 'offset_max'. Only lines with a SNR larger than "
        "'offset_snr_min' will be allowed to move, in addition to Lyman alpha which is "
        "is always allowed to be shifted (since this is a bright resonant line). Enabling "
        "this option will trigger a second fitting pass in the fitting process if "
        "'use_mpfit' is not set.");
    bullet("offset_max=...", "Must be a velocity in km/s. It defines the maximum allowed "
        "velocity offset (when 'allow_offsets' is set) for lines, such that the lines can "
        "move by -offset_max to +offset_max. Default is 1000 km/s.");
    bullet("delta_offset=...", "Must be a number. Defines the size of a step in the grid of "
        "line offsets, as the fraction of the size of a wavelength element of the spectrum. "
        "In other words, given the spectral resolution R of your spectrum, the offset step "
        "will be equal to c*delta_offset/R. Default is 0.2, which at corresponds to 15 km/s "
        "at R=3800 (H+K) and 8 km/s at R=7100 (K).");
    bullet("offset_snr_min=...", "Must be a positive number. It defines the minimum SNR "
        "a line should have to be allowed to shift its centroid owing to velocity offsets "
        "(if 'allow_offsets' is set). Default is SNR=5.");
    bullet("forbid_absorption", "Set this flag to forbid the fit from using absorption lines. "
        "When this option is enabled, all lines are assumed to be emission lines. In cases where "
        "the linear fit would favor an absorption feature, the line is set to zero flux. The "
        "default behavior is to allow each line to be seen either in absorption or emission.");
    bullet("fit_continuum_template", "Set this flag to fit for the continuum emission using "
        "a set of templates provided in the 'template_dir' directory. All these templates "
        "will be linearly combined to best fit the continuum, and the fit is done jointly "
        "with the lines for best accuracy. This flag cannot be combined with "
        "'local_continuum'.");
    bullet("template_dir=...", "Must be the path to a directory containing an unspecified "
        "number of galaxy templates in ASCII format. Each file (with extension '.dat') must "
        "contain one galaxy template defined with two columns: the rest wavelength in "
        "Angstrom, and the flux in F_lambda (normalization is unimportant). Default is "
        "'./templates/'.");
    bullet("local_continuum", "Set this flag to allow a free level of continuum emission "
        "under each line (a constant level). This will decrease the S/N on the line "
        "properties, so be sure to only use it when it is possible to have continuum "
        "emission.");
    bullet("local_continuum_width=...", "Must be a velocity in km/s. It defines the wavelength "
        "region around each line over which the continuum level will be determined to compute "
        "the equivalent widths (EW). In addition, if 'local_continuum' is set, this value is "
        "used to determine the wavelength range around the line over which the continuum level "
        "will be fitted. This range should be large enough that the fit can separate continuum "
        "and line emission. The default value is 3000 km/s, which is substantially larger than "
        "typical line widths, but you can lower it if you need a more precise continuum "
        "determination and you know your lines are relatively narrow.");
    bullet("residual_rescale", "Set this flag to renormalize the uncertainties based on "
        "the residual of the spectrum after subtracting the best-fit lines and continuum "
        "models. This allows to cope for underestimated error spectrum. For each line, "
        "the local RMS of the residual is computed and compared to the expected RMS from "
        "the error spectrum. The uncertainties are rescaled by the ratio of the two (only "
        "if larger than 1, so that errors cannot be reduced). Then, the whole error spectrum "
        "is renormalized accordingly by interpolating between the rescaling factor of the "
        "lines, and the whole fit is performed a second time in a second pass.");
    bullet("mc_errors", "Set this flag to estimate uncertainties on each fitted quantity "
        "by repeatedly and randomly perturbing the spectrum according to the error spectrum "
        "and re-doing the fit each time (\"Monte Carlo\" uncertainties). The uncertainties "
        "are then computed from the standard deviation of their corresponding value among all "
        "these random realizations. This happens after 'residual_rescale', so if this option is "
        "enabled it will influence the amplitude of these random perturbations. The execution "
        "time of the program will be increased by a large factor (see 'num_mc'), so you may "
        "consider enabling multi-threading to compensate.");
    bullet("num_mc=...", "Must be a positive integer. It defines how many random realizations "
        "of the spectrum will be produced to compute the Monte Carlo uncertainties. Default is "
        "200, which should be just enough to get stable 1-sigma confidence intervals. Using "
        "lower values is not advisable unless you only care about a rough estimate of the "
        "uncertainty (in which case a few tens of realizations will suffice). Note that this "
        "value scales linearly with the execution time of the program: the default value of "
        "200 will make the program run 200 times slower! You can enable multi-threading to "
        "compensate, see 'threads'.");
    bullet("threads=...", "Must be a positive integer. It defines how many concurrent threads "
        "the program is allowed to use to speed up the computations. Default is 1, meaning "
        "that multi-threading is not used. The adequate value of this parameter may vary "
        "depending on how much resources are available at run time. A good rule of thumbs is "
        "to allow as many threads as available CPUs. Currently, only the Monte Carlo step "
        "(see 'mc_errors') takes advantage of multi-threading.");
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
    bullet("use_global_chi2", "Set this flag to define the best redshift from the chi2 of "
        "the entire spectrum, rather than the default behavior which is to only use the "
        "spectral elements close to the lines (which is ). This will only make sense if "
        "there is continuum emission and you are using the 'fit_continuum_template' option, "
        "in which case it will increase the weight of the continuum fit in the chi2.");
    bullet("full_range", "Set this flag to fit the entire spectrum. By default the program "
        "will ignore the ranges of the spectrum that are not covered by any line given the "
        "searched redshift range. This speeds up computation, but also ignores some of the "
        "observations, so you can disable this optimization with this flag.");
    bullet("lambda_pad", "Must be an integer. It defines the number of wavelength element "
        "that are ignored both at the beginning and end of the spectrum. Default is 5 "
        "elements. This is used to flag out invalid and poorly covered spectral regions "
        "which could drive the fit toward unrealistic values.");
    bullet("save_model", "Set this flag to also output the best-fit model spectrum. The "
        "spectrum will be saved into the '*_slfit_model.fits' file as a regular 1D "
        "spectrum: the first extension is empty, the second contains the flux.");
    bullet("outdir", "Name of the directory into which the output files should be created. "
        "Default is the current directory.");
    bullet("ascii", "Set this flag if you want the output catalog to be saved in ASCII "
        "format in addition to the default FITS tables. The lines and their fluxes will be "
        "saved in the '*_slfit_lines.cat' file, while the redshift probability distribution "
        "will be saved in '*_slfit_pz.cat'.");
    bullet("flux_hdu=...", "HDU index of the FITS extension containing the flux.");
    bullet("error_hdu=...", "HDU index of the FITS extension containing the uncertainty.");
    bullet("verbose", "Set this flag to print the progress of the detection process in "
        "the terminal. Can be useful if something goes wrong, or just to understand what "
        "is going on.");
}

void print_available_lines(const std::map<std::string,line_t>& db) {
    for (auto& l : db) {
        auto& line = l.second;
        if (line.lambda.size() == 1) {
            print("  - ", line.name, ", lambda=", line.lambda[0], " (", line.pretty_name, ")");
        } else {
            print("  - ", line.name, ", lambda=", line.lambda, ", ratios=", line.ratio, " (", line.pretty_name, ")");
        }
    }
}
