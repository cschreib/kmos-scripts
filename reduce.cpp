#include <phypp.hpp>

int main(int argc, char* argv[]) {
    // If too few arguments provided, print usage and exit
    if (argc < 3) {
        print("reduce (calib|stdstar|sci) <raw_dir> options=...");
        print("reduce combine <sci_base_name> options=...");
    }

    // Get static calibration directory from environment variable
    std::string kmos_calib_dir = file::directorize(system_var("KMOS_CALIB_DIR", ""));
    if (kmos_calib_dir.empty()) {
        error("KMOS_CALIB_DIR variable is not defined!");
        return 1;
    }

    // Read main parameters from the command line: task and directory
    std::string task = argv[1];
    std::string raw_dir = argv[2];
    if (task != "combine") {
        raw_dir = file::directorize(raw_dir);
    }

    // Read secondary parameters from the command line
    vec1s calib;
    std::string stdstar;
    std::string grating;
    vec1s helpers;
    read_args(argc-2, argv+2, arg_list(calib, stdstar, grating, helpers));

    // Extract 'simple' band ID from the grating: HKHKHK -> HK
    std::string band = tolower(grating.substr(0, grating.size()/3));


    // Function to add raw files of a given type to a SOF file
    // file: output SOF file
    // suffix: FITS file suffix to filter (e.g., "dark" or "flat-off")
    // type: pipeline type corresponding to this FITS file (e.g., "DARK" or "FLAT_OFF")
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

    // Function to add a "exit if previous command failed" check in bash scripts
    auto add_stop_fail = [](std::ofstream& file) {
        file << "if [ $? -ne 0 ]; then\n";
        file << "    exit\n";
        file << "fi\n\n";
    };

    // Common variables for all tasks
    std::ofstream main_file("reduce.sh");
    std::ofstream sof;

    calib = file::directorize(calib);
    if (!stdstar.empty()) stdstar = file::directorize(stdstar);

    // Task specific code
    if (task == "calib") {
        // --------------------------------------------
        // Reduce calibration
        // --------------------------------------------

        // A calibration set must contain the following frames:
        //  - dark
        //  - flat-off and flat-lamp
        //  - wave-off and wave-lamp
        //  - flat-sky (optional)
        //
        // Each of the bullets above corresponds to a different step in the calibration
        // reduction process:
        //  - "dark" frames are used to obtain a first badpixel and "zero flux" image
        //    that will be used for the calibration itself, but not for the science
        //    reduction ("kmos_dark" pipeline).
        //  - "flat-off" and "flat-lamp" are used to determine the borders of the IFUs on
        //    the detector, and obtain the final badpixel, "zero flux" and "flat" images
        //    ("kmos_flat" pipeline).
        //  - "wave-off" and "wave-lamp" are used for the wavelength calibration
        //    ("kmos_wave_cal" pipeline).
        //  - "flat-sky" is used to obtain illumination correction
        //    ("kmos_illumination" pipeline).
        //
        // Each of these steps requires a different "SOF" file, containing the list of
        // raw data and existing calibration to use. This is documented in the KMOS
        // pipeline manual.

        if (grating.empty()) {
            error("please indicate the observing band: grating=... (example: grating=HHH)");
            return 1;
        }

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
        // --------------------------------------------
        // Reduce standard stars and science frames
        // --------------------------------------------

        // The reduction in both cases is very similar. The only difference is that, for
        // standard stars, the pipeline will derive the telluric correction (absolute
        // flux calibration) from the reduced cubes ("kmos_std_star" pipeline). For
        // science targets, only IFU cubes will be produced ("kmos_sci_red" pipeline).
        // Here we enforce that the individual exposures within a given OB are not
        // collapsed into a single cube, because we will combine multiple OBs afterwards
        // and it is more efficient to leave them uncombined for now. It also allows
        // checking the astronmetry and quality of each exposure.

        if (calib.empty()) {
            error("please provide the reduced calibration directory(ies): calib=...");
            return 1;
        }

        if (grating.empty()) {
            error("please indicate the observing band: grating=... (example: grating=HHH)");
            return 1;
        }

        auto add_calib = [&](std::ofstream& file) {
            vec1s calfiles = {"xcal_"+grating+".fits", "ycal_"+grating+".fits",
                "lcal_"+grating+".fits", "master_flat_"+grating+".fits",
                "illum_corr_"+grating+".fits"};
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
                    error("missing standard star calibration ("+stdstar+
                        "telluric_"+grating+".fits)");
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
        // --------------------------------------------
        // Obtain continuum images of helper targets
        // --------------------------------------------

        // This task works on the output of the task "sci", i.e., uncollapsed IFU cubes
        // that are stored into a single FITS file.

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
    } else if (task == "combine") {
        // --------------------------------------------
        // Combine multiple OBs into a single master data set
        // --------------------------------------------

        // The individual exposures reduced by the task "sci" are combined into a master
        // cube for all IFUs. Identification of targets is made through the FITS headers
        // and should be reliable. If a dithering pattern is used, the pipeline recognizes
        // it and applies the corresponding expected position shifts. The real shift may
        // be different than the expected one though, and this is not taken care of here.
        // You have to correct for it yourself, manually.

        std::string base_dir = file::get_directory(raw_dir);
        vec1s dirs = file::list_directories(raw_dir+"*");
        inplace_sort(dirs);
        vec1s files;

        for (auto& d : dirs) {
            vec1s tf = file::list_files(base_dir+d+"/sci_reconstructed_*-sci.fits");
            inplace_sort(tf);
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
    } else if (task == "collapse") {
        // --------------------------------------------
        // Obtain continuum images of all science targets
        // --------------------------------------------

        // This task works on the output of "combine".

        if (grating.empty()) {
            error("please indicate the observing band: grating=... (example: grating=HHH)");
            return 1;
        }

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
    } else {
        error("unknown task '", task, "'");
        return 1;
    }

    main_file.close();
    spawn("chmod +x reduce.sh");

    return 0;
}

