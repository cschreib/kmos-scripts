#include <phypp.hpp>

void print_help();

int phypp_main(int argc, char* argv[]) {
    // If too few arguments provided, print usage and exit
    if (argc < 3) {
        print_help();
        return 0;
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
    vec1s options;
    vec1s helpers;
    std::string pipeline_version_str;
    read_args(argc-2, argv+2, arg_list(
        calib, stdstar, grating, helpers, options, name(pipeline_version_str, "pipeline_version")
    ));

    // Read pipeline version (default 1.3.19)
    uint_t pipeline_version = 1;
    float  pipeline_version_rev = 3.19;
    if (!pipeline_version_str.empty()) {
        vec1s spl = split(pipeline_version_str, ".");
        if (spl.size() != 3) {
            error("pipeline version must be of the form 'x.y.z' (got '", pipeline_version_str, "')");
            return 1;
        }

        uint_t pipeline_version_major, pipeline_version_minor;
        if (!from_string(spl[0], pipeline_version) ||
            !from_string(spl[1], pipeline_version_major) ||
            !from_string(spl[1], pipeline_version_minor)) {
            error("failed to read pipeline version from '", pipeline_version_str, "'");
            return 1;
        }

        pipeline_version_rev = pipeline_version_major;
        pipeline_version_rev += pipeline_version_minor*e10(-log10(pipeline_version_minor));
    }

    if (pipeline_version != 1 && pipeline_version_rev < 3.0) {
        error("this script was built for pipeline versions 1.x.y, with x > 3");
        return 1;
    }

    // Extract 'simple' band ID from the grating: HKHKHK -> HK
    std::string band = to_lower(grating.substr(0, grating.size()/3));


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
        main_file << "esorex --log-file=esorex_dark.log --suppress-prefix=TRUE kmos_dark dark.sof\n";
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
        if (pipeline_version_rev >= 3.19) {
            sof << "BADPIXEL_DARK.fits BADPIXEL_DARK\n";
        } else {
            sof << "badpixel_dark.fits BADPIXEL_DARK\n";
        }
        sof.close();

        main_file << "# FLAT\n";
        main_file << "esorex --log-file=esorex_flat.log --suppress-prefix=TRUE kmos_flat flat.sof\n";
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
        if (pipeline_version_rev >= 3.19) {
            sof << "FLAT_EDGE_"+grating+".fits FLAT_EDGE\n";
            sof << "XCAL_"+grating+".fits      XCAL\n";
            sof << "YCAL_"+grating+".fits      YCAL\n";
        } else {
            sof << "flat_edge_"+grating+".fits FLAT_EDGE\n";
            sof << "xcal_"+grating+".fits      XCAL\n";
            sof << "ycal_"+grating+".fits      YCAL\n";
        }
        sof.close();

        main_file << "# WAVE_CAL\n";
        main_file << "esorex --log-file=esorex_wave.log --suppress-prefix=TRUE kmos_wave_cal wave_cal.sof\n";
        add_stop_fail(main_file);

        // ILLUM
        sof.open("illum.sof");
        if (!add_files(sof, "flat-sky", "FLAT_SKY")) {
            warning("no illumination correction in this set");
            sof.close();
            file::remove("illum.sof");
        } else {
            sof << kmos_calib_dir+"kmos_wave_band.fits WAVE_BAND\n";
            if (pipeline_version_rev >= 3.19) {
                sof << "MASTER_DARK.fits        MASTER_DARK\n";
                sof << "MASTER_FLAT_"+grating+".fits MASTER_FLAT\n";
                sof << "XCAL_"+grating+".fits        XCAL\n";
                sof << "YCAL_"+grating+".fits        YCAL\n";
                sof << "LCAL_"+grating+".fits        LCAL\n";
                sof << "FLAT_EDGE_"+grating+".fits   FLAT_EDGE\n";
            } else {
                sof << "master_dark.fits        MASTER_DARK\n";
                sof << "master_flat_"+grating+".fits MASTER_FLAT\n";
                sof << "xcal_"+grating+".fits        XCAL\n";
                sof << "ycal_"+grating+".fits        YCAL\n";
                sof << "lcal_"+grating+".fits        LCAL\n";
                sof << "flat_edge_"+grating+".fits   FLAT_EDGE\n";
            }
            sof.close();

            main_file << "# ILLUM\n";
            main_file << "esorex  --log-file=esorex_illum.log kmos_illumination illum.sof\n";
            add_stop_fail(main_file);
        }
    } else if (task == "stdstar" || task == "sci" || task == "acq") {
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
        // checking the astrometry and quality of each exposure.

        if (calib.empty()) {
            error("please provide the reduced calibration directory(ies): calib=...");
            return 1;
        }

        if (grating.empty()) {
            error("please indicate the observing band: grating=... (example: grating=HHH)");
            return 1;
        }

        auto add_calib = [&](std::ofstream& file) {
            vec1s calfiles = {"xcal_"+grating, "ycal_"+grating,
                "lcal_"+grating, "master_flat_"+grating,
                "illum_corr_"+grating};

            if (pipeline_version_rev >= 3.19) {
                calfiles = to_upper(calfiles);
            }

            calfiles += ".fits";

            vec1s calcat = {"XCAL", "YCAL", "LCAL", "MASTER_FLAT", "ILLUM_CORR"};
            vec1b required = {true, true,   true,   true,          false};
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

            if (count(!found && !required) != 0) {
                warning("missing optional calibration files:");
                vec1u idb = where(!found);
                for (auto& f : calfiles[idb]) {
                    warning(" - ", f);
                }
            }

            if (count(!found && required) != 0) {
                error("missing mandatory calibration files:");
                vec1u idb = where(!found);
                for (auto& f : calfiles[idb]) {
                    error(" - ", f);
                }

                return false;
            }

            return true;
        };

        std::string baked_options = collapse(options, " ");

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
            if (!solar_wave.empty()) {
                sof << kmos_calib_dir+"kmos_solar_"+band+"_"+solar_wave+".fits SOLAR_SPEC\n";
            }

            if (!add_calib(sof)) {
                error("cannot proceed");
                return 1;
            }
            sof.close();

            main_file << "# STD_STAR\n";
            main_file << "esorex  --log-file=esorex_stdstar.log kmos_std_star -save_cubes "
                << baked_options << " stdstar.sof\n";
            add_stop_fail(main_file);

            // We add an extra step to strip the empty IFUs from the image file
            if (pipeline_version_rev >= 3.19) {
                // main_file << "esorex  --log-file=esorex_strip.log kmo_fits_strip -empty STD_IMAGE_"
                //     << grating << ".fits\n";
                // main_file << "mv strip.fits STD_IMAGE_" << grating << ".fits\n";
            } else {
                main_file << "esorex  --log-file=esorex_strip.log kmo_fits_strip -empty std_image_"
                    << grating << ".fits\n";
                main_file << "mv strip.fits std_image_" << grating << ".fits\n";
            }
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
                std::string telluric;
                if (pipeline_version_rev >= 3.19) {
                    telluric = stdstar+"TELLURIC_"+grating+".fits";
                } else {
                    telluric = stdstar+"telluric_"+grating+".fits";
                }
                if (!file::exists(telluric)) {
                    error("missing standard star calibration ("+telluric+")");
                    return 1;
                }
                sof << telluric+" TELLURIC\n";
            }
            sof.close();

            main_file << "# SCI\n";
            main_file << "# Cubic spline interp., to get best flux\n";
            main_file << "mkdir -p CS\n";
            main_file << "esorex  --log-file=esorex_sci.log kmos_sci_red -no_combine "
                "-background " << baked_options << " sci.sof\n";
            main_file << "mv *.fits CS/\n\n";
            main_file << "# Nearest neighbor interp., to get uncertainty\n";
            main_file << "mkdir -p NN\n";
            main_file << "esorex  --log-file=esorex_sci.log kmos_sci_red -no_combine "
                "-background " << baked_options << " sci.sof\n";
            main_file << "mv *.fits NN/\n\n";
            main_file << "# Combine the two interpolations\n";
            main_file << "${KMOS_SCRIPTS_DIR}/merge_cs_nn CS NN\n\n";
            add_stop_fail(main_file);
        } else if (task == "acq") {
            print("prepare reduction of acquisition frames in ", raw_dir);

            sof.open("acq.sof");
            if (!add_files(sof, "acq", "SCIENCE")) {
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
                std::string telluric;
                if (pipeline_version_rev >= 3.19) {
                    telluric = stdstar+"TELLURIC_"+grating+".fits";
                } else {
                    telluric = stdstar+"telluric_"+grating+".fits";
                }
                if (!file::exists(telluric)) {
                    error("missing standard star calibration ("+telluric+")");
                    return 1;
                }
                sof << telluric+" TELLURIC\n";
            }
            sof.close();

            main_file << "# ACQ\n";
            main_file << "esorex  --log-file=esorex_acq.log kmos_sci_red -no_combine "
                "-background " << baked_options << " acq.sof\n";
            add_stop_fail(main_file);
        }
    } else if (task == "helpers") {
        // --------------------------------------------
        // Obtain continuum images of helper targets
        // --------------------------------------------

        // This task does two things works on the output of the task "sci", i.e.,
        // uncollapsed IFU cubes that are stored into a single FITS file. It does two
        // things:
        // 1) Extract the continuum image of helper targets for each exposure. This allows
        //    you to check for systematic position offsets between exposures.
        // 2) Extract the continuum image of helper targets for the combined OB, merging
        //    all the exposures. This allows you to identify systematic position offsets
        //    between OBs.

        if (grating.empty()) {
            error("please indicate the observing band: grating=... (example: grating=HHH)");
            return 1;
        }

        print("prepare reduction of helper targets in ", raw_dir);

        vec1s dithers = raw_dir+file::list_files(raw_dir+"*.fits");
        inplace_sort(dithers);

        sof.open("combine.sof");

        for (uint_t i : range(dithers)) {
            std::ofstream sof2("cont"+strn(i+1)+".sof");
            sof2 << dithers[i] << " COMMAND_LINE\n";
            sof2 << kmos_calib_dir+"kmos_oh_spec_"+band+".fits COMMAND_LINE\n";
            sof2.close();

            sof << dithers[i] << " COMMAND_LINE\n";

            std::string out_file = file::remove_extension(file::get_basename(dithers[i]))+
                "_img_cont.fits";

            main_file << "# Dither " << i+1 << "\n";
            main_file << "esorex kmo_make_image cont" << i+1 << ".sof\n";
            if (pipeline_version_rev >= 3.19) {
                main_file << "mv MAKE_IMAGE.fits " << out_file << "\n";
            } else {
                main_file << "mv make_image.fits " << out_file << "\n";
            }
            main_file << "${KMOS_SCRIPTS_DIR}/extract_ifu " << out_file << " names=[" <<
                collapse(helpers, ",") << "]\n";
            main_file << "rm " << out_file << "\n\n";
        }

        sof.close();

        for (std::string helper : helpers) {
            sof.open("image"+helper+".sof");
            if (pipeline_version_rev >= 3.19) {
                sof << "COMBINE_SCI_RECONSTRUCTED_" << helper << ".fits COMMAND_LINE\n";
            } else {
                sof << "combine_sci_reconstructed_" << helper << ".fits COMMAND_LINE\n";
            }
            sof << kmos_calib_dir+"kmos_oh_spec_"+band+".fits COMMAND_LINE\n";
            sof.close();

            main_file << "# Full " << helper << "\n";
            main_file << "esorex --log-file=esorex_combine_" << helper <<
                ".log kmos_combine -method='header' -cmethod='median' "
                "-name='" << to_upper(helper) << "' combine.sof\n";
            main_file << "esorex kmo_make_image image"+helper+".sof\n";
            if (pipeline_version_rev >= 3.19) {
                main_file << "rm COMBINE_SCI_RECONSTRUCTED_" << helper << ".fits\n";
                main_file << "rm EXP_MASK_SCI_RECONSTRUCTED_" << helper << ".fits\n";
                main_file << "mv MAKE_IMAGE.fits COMBINE_SCI_RECONSTRUCTED_" << helper
                    << "_img_cont.fits\n\n";
            } else {
                main_file << "rm combine_sci_reconstructed_" << helper << ".fits\n";
                main_file << "rm exp_mask_sci_reconstructed_" << helper << ".fits\n";
                main_file << "mv make_image.fits combine_sci_reconstructed_" << helper
                    << "_img_cont.fits\n\n";
            }
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
            vec1s tf;
            if (pipeline_version_rev >= 3.19) {
                tf = file::list_files(base_dir+d+"/SCI_RECONSTRUCTED_*-sci.fits");
            } else {
                tf = file::list_files(base_dir+d+"/sci_reconstructed_*-sci.fits");
            }
            inplace_sort(tf);
            append(files, base_dir+d+"/"+tf);
        }

        sof.open("combine.sof");
        for (auto& f : files) {
            sof << f << " COMMAND_LINE\n";
        }
        sof.close();

        main_file << "# COMBINE\n";
        main_file << "esorex  --log-file=esorex_combine.log kmos_combine -edge_nan -method='header' combine.sof\n";
        add_stop_fail(main_file);
        if (pipeline_version_rev >= 3.19) {
            main_file << "for f in COMBINE*.fits; do echo $f; ${KMOS_SCRIPTS_DIR}/fill_nan $f; done\n";
        } else {
            main_file << "for f in combine*.fits; do echo $f; ${KMOS_SCRIPTS_DIR}/fill_nan $f; done\n";
        }
    } else if (task == "collapse") {
        // --------------------------------------------
        // Obtain continuum images of all science targets
        // --------------------------------------------

        // This task works on the output of "combine".

        if (grating.empty()) {
            error("please indicate the observing band: grating=... (example: grating=HHH)");
            return 1;
        }

        std::string prefix;
        if (pipeline_version_rev >= 3.19) {
            prefix = "COMBINE_SCI_RECONSTRUCTED_";
        } else {
            prefix = "combine_sci_reconstructed_";
        }

        vec1s cubes = raw_dir+file::list_files(raw_dir+prefix+"*.fits");
        inplace_sort(cubes);

        for (uint_t i : range(cubes)) {
            std::string out_file = file::remove_extension(file::get_basename(cubes[i]))+"_img_cont.fits";
            std::string sid = erase_begin(erase_end(out_file, "_img_cont.fits"), prefix);

            sof.open("cont_"+sid+".sof");
            sof << cubes[i] << " COMMAND_LINE\n";
            sof << kmos_calib_dir+"kmos_oh_spec_"+band+".fits COMMAND_LINE\n";
            sof.close();

            main_file << "# " << sid << "\n";
            main_file << "esorex kmo_make_image -cmethod='median' cont_" << sid << ".sof\n";
            if (pipeline_version_rev >= 3.19) {
                main_file << "mv MAKE_IMAGE.fits " << out_file << "\n\n";
            } else {
                main_file << "mv make_image.fits " << out_file << "\n\n";
            }
        }
    } else {
        error("unknown task '", task, "'");
        return 1;
    }

    main_file.close();
    spawn("chmod +x reduce.sh");

    return 0;
}

void print_help() {
    using namespace format;

    print("reduce v1.0");
    print("usage: reduce <task> <target> grating=... [options]");
    print("");
    paragraph("This program will help you create the files needed by the KMOS pipeline "
        "(EsoRex) to run a full reduction of raw data into spectral cubes. It decomposes "
        "the process into various \"tasks\" which are listed below. In all cases, the "
        "output files and scripts are created in the current directory, including the SOF "
        "files needed by the pipeline and a 'reduce.sh' script that will run the pipeline "
        "reduction steps. You can start the reduction by calling this newly created "
        "script. Also, for most tasks you have to specify the chosen grating for your "
        "observations in 'grating=...'. For example, for observations in the K band, "
        "'grating=KKK', while H+K band would be 'grating=HKHKHK'.");
    print("'task' must be one of the following:");
    bullet("calib", "Prepare the reduction of a KMOS calibration set, including darks, "
        "flats, illumination, spatial and spectral calibration. The parameter 'target' "
        "must be the directory containing the raw frames of this calibration set (no "
        "other calibration sets shall be present in this same directory).");
    bullet("stdstar", "Prepare the reduction of a KMOS standard star calibration set, to "
        "obtain absolute flux calibration. The parameter 'target' must be the directory "
        "containing the raw frames of this standard star calibration set (no other "
        "calibration shall be present in this same directory). You also have to provide "
        "the 'calib=...' parameter (see below) to indicate which calibration files to use "
        "for the reduction of these stars.");
    bullet("sci", "Prepare the reduction of a KMOS science observing block (OB), to "
        "obtain spectral cubes for all targets in each exposure of this OB. The "
        "parameter 'target' must be the directory containing the raw frames of this OB "
        "(no other OB shall be present in this same directory). You also have to provide "
        "the 'calib=...' parameter (see below) to indicate which calibration files to use "
        "for the reduction of these targets, as well as the 'stdstar=...' parameter (see "
        "below) to specify which standard stars should be used for the flux calibration.");
    bullet("helpers", "Prepare the necessary files to create continuum images of 'helper' "
        "targets; bright continuum sources in your OBs. These targets are useful to check "
        "the quality of the reduction, the astrometry, and the absolute flux calibration. "
        "This task will create images of these targets for a given OB, for each exposure, "
        "and also combining all the exposures of the OB. The parameter 'target' must be "
        "the directory containing the reduced cubes, created by the task 'sci'. You also "
        "have to provide the 'helpers=...' parameter (see below) to list the names of "
        "your helper targets (name as given in the KARMA catalog).");
    bullet("combine", "Prepare the reduction of multiple KMOS OBs into a single data set, "
        "combining all the exposures previously reduced using the 'sci' task. The paramter "
        "'target' must give the pattern to use to look for directories containing the "
        "reduced OBs. For example, if you have two OBs in the 'sci-01' and 'sci-02' "
        "directories, then 'target' must be 'sci-*'. The pipeline will create one FITS "
        "file per target observed in these OBs, containing the combined cube.");
    bullet("collapse", "Prepare the necessary files to create continuum images of all "
        "the observed targets in a combined data set. The parameter 'target' must be "
        "the directory containing the combined cubes, created by the task 'combine'. "
        "The pipeline will create one FITS image for each target.");
    print("\nAvailable options:");
    bullet("calib=[...]", "Needed for tasks 'stdstar' and 'sci'. This parameter lists "
        "the directories containing the reduced calibration files to be used for the "
        "reduction of spectral cubes. Multiple directories can be listed, in case "
        "not all required calibration files are found into a single calibration set "
        "(i.e., in particular illumination correction is often missing). In that case, "
        "list the different directories in order of decreasing priority: all the files of "
        "the first directory will be picked, and only the files that could not be found "
        "there will be searched for in the next directories.");
    bullet("stdstar=...", "Needed for task 'sci'. This parameter gives the directory "
        "that contains the reduced flux calibration files to be used for the reduction of "
        "fully calibrated spectral cubes. Only a single directory should be listed.");
    bullet("helpers=[...]", "Needed for task 'helpers'. Lists the names of the helper "
        "targets observed in this OB. The names should be the ones given in the KARMA "
        "catalog when you prepared the Phase 2 before the observations started. These "
        "names are case-insensitive. If a helper target is not found, it will simply be "
        "ignored.");
    bullet("options=[...]", "Optional, available for all tasks. This parameter lists all "
        "additional options you wish to pass to the pipeline for the reduction. These "
        "options may not contain commas (,) or quotes (\"). If they contain spaces, be "
        "sure to write the parameter as: options=\"[-option1=foo, -option2='foo bar'].\"");
    bullet("pipeline_version=...", "Optional, available for all tasks. This parameter "
        "allows you to specify which version of the KMOS pipeline you have installed "
        "on your computer. The reduction scripts sometimes need to be adjusted depending "
        "on which version you are running, e.g., because the names of the intermediate "
        "FITS files have changed (which happened between v1.3.17 and v1.3.19).");
    print("");
}
