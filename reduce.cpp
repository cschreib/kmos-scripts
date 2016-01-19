#include <phypp.hpp>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        print("reduce (calib|stdstar|sci) <raw_dir> options=...");
        print("reduce combine <sci_base_name> options=...");
    }

    std::string kmos_calib_dir = file::directorize(system_var("KMOS_CALIB_DIR", ""));
    if (kmos_calib_dir.empty()) {
        error("KMOS_CALIB_DIR variable is not defined!");
        return 1;
    }

    std::string task = argv[1];
    std::string raw_dir = argv[2];
    if (task != "combine") {
        raw_dir = file::directorize(raw_dir);
    }

    vec1s calib;
    std::string stdstar;
    std::string grating;
    vec1s helpers;
    read_args(argc-2, argv+2, arg_list(calib, stdstar, grating, helpers));

    std::string band = tolower(grating.substr(0, grating.size()/3));

    std::ofstream main_file("reduce.sh");
    std::ofstream sof;

    calib = file::directorize(calib);
    if (!stdstar.empty()) stdstar = file::directorize(stdstar);

    auto add_files = [&](std::ofstream& file, std::string suffix, std::string type) {
        vec1s files = file::list_files(raw_dir+"*-"+suffix+".fits");
        inplace_sort(files);

        if (files.empty()) {
            note("no '", suffix, "' frame in ", raw_dir);
            return false;
        }

        for (auto& f : files) {
            file << raw_dir << f << " " << type << "\n";
        }

        return true;
    };

    auto add_stop_fail = [](std::ofstream& file) {
        file << "if [ $? -ne 0 ]; then\n";
        file << "    exit\n";
        file << "fi\n\n";
    };

    if (task == "calib") {
        print("prepare reduction of calibration data in ", raw_dir);

        // DARK
        sof.open("dark.sof");
        if (!add_files(sof, "dark", "DARK")) {
            error("cannot proceed");
            return 1;
        }
        sof.close();

        main_file << "# DARK\n";
        main_file << "esorex --suppress-prefix=TRUE kmos_dark dark.sof\n";
        add_stop_fail(main_file);

        // FLAT
        sof.open("flat.sof");
        if (!add_files(sof, "flat-off", "FLAT_OFF")) {
            error("cannot proceed");
            return 1;
        }
        if (!add_files(sof, "flat-lamp", "FLAT_ON")) {
            error("cannot proceed");
            return 1;
        }
        sof << "badpixel_dark.fits BADPIXEL_DARK\n";
        sof.close();

        main_file << "# FLAT\n";
        main_file << "esorex --suppress-prefix=TRUE kmos_flat flat.sof\n";
        add_stop_fail(main_file);

        // WAVE_CAL
        sof.open("wave_cal.sof");
        if (!add_files(sof, "wave-off", "ARC_OFF")) {
            error("cannot proceed");
            return 1;
        }
        if (!add_files(sof, "wave-lamp", "ARC_ON")) {
            error("cannot proceed");
            return 1;
        }
        sof << kmos_calib_dir+"kmos_wave_ref_table.fits REF_LINES\n";
        sof << kmos_calib_dir+"kmos_wave_band.fits      WAVE_BAND\n";
        sof << kmos_calib_dir+"kmos_ar_ne_list_"+band+".fits  ARC_LIST\n";
        sof << "flat_edge_"+grating+".fits FLAT_EDGE\n";
        sof << "xcal_"+grating+".fits      XCAL\n";
        sof << "ycal_"+grating+".fits      YCAL\n";
        sof.close();

        main_file << "# WAVE_CAL\n";
        main_file << "esorex --suppress-prefix=TRUE kmos_wave_cal wave_cal.sof\n";
        add_stop_fail(main_file);

        // ILLUM
        sof.open("illum.sof");
        if (!add_files(sof, "flat-sky", "FLAT_SKY")) {
            warning("no illumination correction in this set");
            sof.close();
            file::remove("illum.sof");
        } else {
            sof << kmos_calib_dir+"kmos_wave_band.fits WAVE_BAND\n";
            sof << "master_dark.fits        MASTER_DARK\n";
            sof << "master_flat_"+grating+".fits MASTER_FLAT\n";
            sof << "xcal_"+grating+".fits        XCAL\n";
            sof << "ycal_"+grating+".fits        YCAL\n";
            sof << "lcal_"+grating+".fits        LCAL\n";
            sof << "flat_edge_"+grating+".fits   FLAT_EDGE\n";
            sof.close();

            main_file << "# ILLUM\n";
            main_file << "esorex kmos_illumination illum.sof\n";
            add_stop_fail(main_file);
        }
    } else if (task == "stdstar" || task == "sci") {
        if (calib.empty()) {
            error("please provide the reduced calibration directory(ies): calib=...");
            return 1;
        }

        auto add_calib = [&](std::ofstream& file) {
            vec1s calfiles = {"xcal_"+grating+".fits", "ycal_"+grating+".fits", "lcal_"+grating+".fits", "master_flat_"+grating+".fits", "illum_corr_"+grating+".fits"};
            vec1s calcat = {"XCAL", "YCAL", "LCAL", "MASTER_FLAT", "ILLUM_CORR"};
            vec1b found = replicate(false, calfiles.size());
            for (auto& c : calib) {
                vec1u idb = where(!found);
                bool anyfound = false;
                for (uint_t i : idb) {
                    if (file::exists(c+calfiles[i])) {
                        found[i] = true;
                        anyfound = true;
                        file << c+calfiles[i] << " " << calcat[i] << "\n";
                    }
                }

                if (!anyfound) {
                    warning("no calibration data taken from calibration set ", c);
                }
            }

            if (count(!found) != 0) {
                error("missing calibration files:");
                vec1u idb = where(!found);
                for (auto& f : calfiles[idb]) {
                    error(" - ", f);
                }

                return false;
            }

            return true;
        };

        if (task == "stdstar") {
            print("prepare reduction of standard star in ", raw_dir);

            sof.open("stdstar.sof");
            if (!add_files(sof, "object-sky-std-flux", "STD")) {
                error("cannot proceed");
                return 1;
            }

            std::string solar_wave;
            if (band == "h") {
                solar_wave = "2400";
            } else if (band == "hk") {
                solar_wave = "1100";
            } else if (band == "k") {
                solar_wave = "1700";
            }

            sof << kmos_calib_dir+"kmos_wave_band.fits     WAVE_BAND\n";
            sof << kmos_calib_dir+"kmos_spec_type.fits     SPEC_TYPE_LOOKUP\n";
            sof << kmos_calib_dir+"kmos_atmos_"+band+".fits      ATMOS_MODEL\n";
            sof << kmos_calib_dir+"kmos_solar_"+band+"_"+solar_wave+".fits SOLAR_SPEC\n";
            if (!add_calib(sof)) {
                error("cannot proceed");
                return 1;
            }
            sof.close();

            main_file << "# STD_STAR\n";
            main_file << "esorex kmos_std_star -save_cubes stdstar.sof\n";
            add_stop_fail(main_file);
        } else if (task == "sci") {
            print("prepare reduction of science frames in ", raw_dir);

            sof.open("sci.sof");
            if (!add_files(sof, "sci", "SCIENCE")) {
                error("cannot proceed");
                return 1;
            }
            sof << kmos_calib_dir+"kmos_wave_band.fits  WAVE_BAND\n";
            sof << kmos_calib_dir+"kmos_oh_spec_"+band+".fits OH_SPEC\n";
            if (!add_calib(sof)) {
                error("cannot proceed");
                return 1;
            }
            if (!stdstar.empty()) {
                if (!file::exists(stdstar+"telluric_"+grating+".fits")) {
                    error("missing standard star calibration ("+stdstar+"telluric_"+grating+".fits)");
                    return 1;
                }
                sof << stdstar+"telluric_"+grating+".fits TELLURIC\n";
            }
            sof.close();

            main_file << "# SCI\n";
            main_file << "esorex kmos_sci_red -no_combine -background sci.sof\n";
            add_stop_fail(main_file);
        }
    } else if (task == "helpers") {
        print("prepare reduction of helper targets in ", raw_dir);

        vec1s dithers = raw_dir+file::list_files(raw_dir+"*.fits");
        inplace_sort(dithers);

        for (uint_t i : range(dithers)) {
            sof.open("cont"+strn(i+1)+".sof");
            sof << dithers[i] << "\n";
            sof << kmos_calib_dir+"kmos_oh_spec_"+band+".fits\n";
            sof.close();

            std::string out_file = file::remove_extension(file::get_basename(dithers[i]))+
                "_img_cont.fits";

            main_file << "# Dither " << i+1 << "\n";
            main_file << "esorex kmo_make_image cont" << i+1 << ".sof\n";
            main_file << "mv make_image.fits " << out_file << "\n";
            main_file << "../../extract_ifu " << out_file << " names=[" <<
                collapse(helpers, ",") << "]\n";
            main_file << "rm " << out_file << "\n\n";
        }
    } else if (task == "collapse") {
        vec1s cubes = raw_dir+file::list_files(raw_dir+"combine_sci_reconstructed_*.fits");
        inplace_sort(cubes);

        for (uint_t i : range(cubes)) {
            std::string out_file = file::remove_extension(file::get_basename(cubes[i]))+
                "_img_cont.fits";

            std::string sid = erase_begin(
                erase_end(out_file, "_img_cont.fits"),
                "combine_sci_reconstructed_");

            sof.open("cont_"+sid+".sof");
            sof << cubes[i] << "\n";
            sof << kmos_calib_dir+"kmos_oh_spec_"+band+".fits\n";
            sof.close();

            main_file << "# " << sid << "\n";
            main_file << "esorex kmo_make_image cont_" << sid << ".sof\n";
            main_file << "mv make_image.fits " << out_file << "\n\n";
        }
    } else if (task == "combine") {
        std::string base_dir = file::get_directory(raw_dir);
        vec1s dirs = file::list_directories(raw_dir+"*");
        vec1s files;

        for (auto& d : dirs) {
            vec1s tf = file::list_files(base_dir+d+"/sci_reconstructed_*-sci.fits");
            append(files, base_dir+d+"/"+tf);
        }

        sof.open("combine.sof");
        for (auto& f : files) {
            sof << f << "\n";
        }
        sof.close();

        main_file << "# COMBINE\n";
        main_file << "esorex kmos_combine -edge_nan -method='header' combine.sof\n";
        add_stop_fail(main_file);
    } else {
        error("unknown task '", task, "'");
        return 1;
    }

    main_file.close();
    spawn("chmod +x reduce.sh");

    return 0;
}

