#include <phypp.hpp>

int phypp_main(int argc, char* argv[]) {
    vec1s  files;                 // 1D spectra to analyze
    bool   verbose = false;       // print progress in the standard output
    vec1s  filters;               // filters to use to build broadband fluxes

    bool do_sigma_clip = true;          // enable/disable sigma clipping of outliers
    double sigma_clip_threshold = 5.0;  // significance threshold for rejecting outliers
    uint_t sigma_clip_width = 1;        // width (in pixels) of the wavelength bin in which to define outliers

    // Flux catalog used to rescale exposures to the right total flux
    std::string catalog_file;
    std::string catalog_id;

    read_args(argc, argv, arg_list(
        files, filters, name(catalog_file, "catalog"), catalog_id,
        do_sigma_clip, sigma_clip_threshold, sigma_clip_width, verbose
    ));

    if (files.empty()) {
        error("no spectrum to stack (files=...)");
        return 1;
    }

    vec<1,filter_t> fbb;
    if (!filters.empty()) {
        auto db = read_filter_db(data_dir+"fits/filter-db/db.dat");

        if (!get_filters(db, filters, fbb)) {
            return 1;
        }
    }

    vec1d blam(fbb.size());
    for (uint_t f : range(fbb)) {
        blam[f] = fbb[f].rlam;
    }

    vec2d cflx, cerr;

    // Read spectra
    vec1d lam;
    for (uint_t i : range(files)) {
        vec1d tflx, terr;

        fits::input_image iimg(files[i]);
        iimg.reach_hdu(1);
        iimg.read(tflx);

        if (cflx.empty()) {
            double crpix, cdelt, crval;
            if (iimg.read_keyword("CRPIX1", crpix) && iimg.read_keyword("CRVAL1", crval) &&
                iimg.read_keyword("CDELT1", cdelt)) {
                lam = cdelt*(findgen(tflx.size())+1 - crpix) + crval;
            } else {
                std::string ctype;
                if (!iimg.read_keyword("CTYPE1", ctype)) {
                    error("could not read WCS of wavelength axis");
                    return 1;
                }

                if (begins_with(ctype, "TAB")) {
                    // Tabulated axis, read from other extensions
                    uint_t lowext = npos, upext = npos;
                    std::string axis = erase_begin(ctype, "TAB");
                    if (iimg.read_keyword(axis+"LOWEXT", lowext) &&
                        iimg.read_keyword(axis+"UPEXT", upext)) {
                        vec1d xl, xu;
                        iimg.reach_hdu(lowext);
                        iimg.read(xl);
                        iimg.reach_hdu(upext);
                        iimg.read(xu);
                        lam = 0.5*(xl + xu);
                        iimg.reach_hdu(1);
                    }
                } else {
                    error("could not read WCS of wavelength axis");
                    return 1;
                }
            }

            std::string cunit;
            if (iimg.read_keyword("CUNIT1", cunit)) {
                cunit = to_lower(cunit);
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
                } else {
                    error("unrecognized wavelength/frequency unit '", cunit, "'");
                    return false;
                }

                lam *= conv;
            }
        }

        iimg.reach_hdu(2);
        iimg.read(terr);

        if (cflx.empty()) {
            cflx.resize(tflx.dims, files.size());
            cerr.resize(cflx.dims);
        }

        cflx(_,i) = tflx;
        cerr(_,i) = terr;
    }

    // Define weights
    vec2d cwei = 1/sqr(cerr);
    // Adjust global normalization to avoid numerical errors
    cwei /= max(cwei[where(is_finite(cwei))]);

    // Apply sigma clipping (if asked)
    vec2b crej(cflx.dims);
    if (do_sigma_clip) {
        uint_t d1 = (sigma_clip_width-1)/2;
        uint_t d2 = sigma_clip_width-1 - d1;

        for (uint_t l : range(cflx.dims[0])) {
            // Define wavelength region
            uint_t l0 = (l > d1 ?              l - d1 : 0);
            uint_t l1 = (l < cflx.dims[0]-d2 ? l + d2 : cflx.dims[0]-1);

            // First compute the weighted median (which we assume is unbiased)
            vec2d med(l1-l0+1, cflx.dims[1]);
            for (uint_t ll : range(med.dims[0])) {
                med(ll,_) = weighted_median(cflx(l0+ll,_), cwei(l0+ll,_));
            }

            // Compute the absolute residuals
            vec2d res = abs(cflx(l0-_-l1,_) - med);

            // Compute the RMS of these using the MAD (which we assume is unbiased)
            double rms = 1.48*median(res);

            // Select significant outliers in this wavelength element
            crej(l,_) = res(d1,_) > sigma_clip_threshold*rms;
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

    // Resample the filters to the grid of the spectra
    for (auto& f : fbb) {
        f.res = interpolate(f.res, f.lam, lam);
        uint_t i0 = where_first(lam > f.lam.front());
        uint_t i1 = where_last(lam < f.lam.back());
        if (i0 == npos || i1 == npos) {
            error("filter is not covered by spectrum (", mean(f.lam), " vs. ", min(lam), " - ", max(lam), ")");
            return 1;
        }

        f.res[_-i0] = 0;
        f.res[i1-_] = 0;
        f.lam = lam;
    }

    // Read catalog fluxes
    vec2f cflux, cflux_err;
    vec1s csid, cbands;
    uint_t cid;

    if (!catalog_file.empty()) {
        fits::input_table itbl(catalog_file);
        itbl.read_columns("flux", cflux, "flux_err", cflux_err, "bands", cbands);
        if (!itbl.read_column("id", csid)) {
            vec1u cuid;
            if (itbl.read_column("id", cuid)) {
                csid = to_string_vector(cuid);
            } else {
                error("could not find column 'ID' in ", catalog_file);
                return 1;
            }
        }

        cid = where_first(csid == catalog_id);
        if (cid == npos) {
            error("could not find source with ID=", catalog_id, " in ", catalog_file);
            return 1;
        }
    }

    // Compute broad band fluxes for each exposure
    for (uint_t i : range(cflx.dims[1])) {
        vec1f bflx(filters.size());
        vec1f berr(filters.size());
        vec1f cat_flx = replicate(fnan, filters.size());
        vec1f cat_err = replicate(fnan, filters.size());

        for (uint_t b : range(filters)) {
            vec1f tw = cwei(_,i)*fbb[b].res;
            double w = total(tw);

            bflx[b] = total(tw*cflx(_,i))/w;
            berr[b] = sqrt(total(sqr(tw*cerr(_,i))))/w;

            if (!catalog_file.empty()) {
                uint_t ib = where_first(cbands == filters[b]);
                if (ib != npos) {
                    cat_flx[b] = uJy2cgs(fbb[b].rlam, cflux(cid,ib));
                    cat_err[b] = uJy2cgs(fbb[b].rlam, cflux_err(cid,ib));
                }
            }
        }

        if (!catalog_file.empty()) {
            fits::write_table(file::remove_extension(files[i])+"_broadband.fits",
                "flux", bflx, "flux_err", berr, "bands", filters, "lambda", blam,
                "catalog_flux", cat_flx, "catalog_flux_err", cat_err
            );
        } else {
            fits::write_table(file::remove_extension(files[i])+"_broadband.fits",
                "flux", bflx, "flux_err", berr, "bands", filters, "lambda", blam
            );
        }
    }

    return 0;
}
