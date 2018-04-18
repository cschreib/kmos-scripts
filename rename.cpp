#include <phypp.hpp>

void print_help();

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    std::string dir = file::directorize(argv[1]);

    vec1s files = file::list_files(dir+"KMOS*.fits");
    print(files.size(), " FITS files to analyze");

    for (auto& file : files) {
        fits::header hdr = fits::read_header_hdu(dir+file, 0);

        std::string obj;
        if (!fits::getkey(hdr, "OBJECT", obj)) {
            error("missing OBJECT keyword in '", file, "'");
            continue;
        }

        obj = trim(to_lower(replace(obj, ",", "-")));

        if (obj == "object") obj = "acq";
        else if (obj == "dark" || obj == "flat-off" || obj == "flat-lamp" ||
            obj == "wave-off" || obj == "wave-lamp" || obj == "flat-sky" ||
            obj == "sky" || obj == "object-sky-std-flux") {
            // nothing to do
        } else {
            // This must be a science frame
            // Make sure
            std::string tpl;
            if (!fits::getkey(hdr, "HIERARCH ESO TPL ID", tpl)) {
                error("missing HIERARCH ESO TPL ID keyword in '", file, "'");
                continue;
            }

            tpl = trim(to_lower(tpl));

            if (!begins_with(tpl, "kmos_spec_obs")) {
                warning("unknown frame type '", obj, "' with template '", tpl, "'");
                continue;
            }

            obj = "sci";
        }

        if (!ends_with(file, obj+".fits")) {
            spawn("mv "+dir+file+" "+dir+file::remove_extension(file)+"-"+obj+".fits");
        }
    }

    return 0;
}

void print_help() {
    using namespace format;

    print("rename v1.0");
    print("usage: rename <directory>");
    print("");
    print("Main parameters:");
    paragraph("'directory' must contain an abitrary number of KMOS raw data frames in "
        "FITS format. This program will rename these files one by one, appending a suffix "
        "to the file's name to allow you to identify the file's purpose.");
}
