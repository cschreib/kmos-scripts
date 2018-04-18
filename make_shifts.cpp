#include <phypp.hpp>

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print("usage: make_shifts <ob_filter> [options]");
        return 0;
    }

    vec1u exclude;
    std::string helper;

    read_args(argc-1, argv+1, arg_list(exclude, helper));

    if (helper.empty()) {
        error("please provide the name of (one of) the helper target you used to "
            "calibrate the shifts (helper=...)");
        return 1;
    }

    helper = to_lower(helper);

    std::string scis = argv[1];

    vec1d cent_ra, cent_dec;
    ascii::read_table("centroid_helper.txt", 0, cent_ra, cent_dec);

    vec1u ids = uindgen(cent_ra.size())+1;
    vec1u idex = where(is_any_of(ids, exclude));
    inplace_remove(ids, idex);
    inplace_remove(cent_ra, idex);
    inplace_remove(cent_dec, idex);

    std::ofstream cmb("combine.sof");

    vec1d shx, shy;

    vec1d x0, y0;
    bool first_line = true;
    uint_t nexp = 0;

    for (uint_t i : range(ids)) {
        std::string dir = scis+align_right(strn(ids[i]), 2, '0')+"/";
        print(dir);

        vec1s files = dir+file::list_files(dir+"sci_reconstructed*-sci.fits");
        if (files.empty()) {
            files = dir+file::list_files(dir+"SCI_RECONSTRUCTED*-sci.fits");
        }
        inplace_sort(files);

        // Find out which exposures contain the helper target from which the
        // shifts were calibrated
        vec1u ignore;
        for (uint_t k : range(files)) {
            std::string f = files[k];
            fits::generic_file fcubes(f);

            vec1s arms(24);
            bool badfile = false;
            for (uint_t u : range(24)) {
                if (!fcubes.read_keyword("ESO OCS ARM"+strn(u+1)+" NAME", arms[u])) {
                    note("ignoring invalid file '", f, "'");
                    note("missing keyword 'ESO OCS ARM"+strn(u+1)+" NAME'");
                    ignore.push_back(k);
                    badfile = true;
                    break;
                }
            }

            if (badfile) continue;

            arms = to_lower(trim(arms));

            vec1u ida = where(arms == helper);
            if (ida.empty()) {
                ignore.push_back(k);
            } else if (x0.empty()) {
                // Get astrometry of IFUs from first exposure
                // NB: assumes the rotation is the same for all exposures,
                // which is anyway what kmos_combine does later on.
                fcubes.reach_hdu(ida[0]+1);
                astro::ad2xy(astro::wcs(fcubes.read_header()), cent_ra, cent_dec, x0, y0);
            }
        }

        inplace_remove(files, ignore);

        if (files.empty()) {
            warning("folder ", dir, " does not contain any usable file");
            continue;
        }

        for (std::string f : files) {
            cmb << f << " COMMAND_LINE\n";
            ++nexp;
        }

        double dox = x0[0] - x0[i], doy = y0[0] - y0[i];

        if (first_line) {
            // Ommit first line which has, by definition, no offset
            first_line = false;
        } else {
            shx.push_back(dox);
            shy.push_back(doy);
        }

        if (file::exists(dir+"helpers/shifts.txt")) {
            vec1d tx, ty;
            ascii::read_table(dir+"helpers/shifts.txt", 0, tx, ty);
            for (uint_t j : range(tx)) {
                shx.push_back(dox+tx[j]);
                shy.push_back(doy+ty[j]);
            }
        } else {
            warning("could not find dither in ", dir, "helpers/shifts.txt");
        }
    }

    note("found ", nexp, " exposures");

    cmb.close();

    auto truncate_decimals = vectorize_lambda([](double v, uint_t nd) {
        return long(v*e10(nd))/e10(nd);
    });

    shx = truncate_decimals(shx, 2);
    shy = truncate_decimals(shy, 2);

    ascii::write_table("shifts.txt", 10, shx, shy);

    return 0;
}
