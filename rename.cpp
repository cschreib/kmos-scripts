#include <phypp.hpp>

int main(int argc, char* argv[]) {
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

        obj = trim(tolower(replace(obj, ",", "-")));

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

            tpl = trim(tolower(tpl));

            if (!start_with(tpl, "kmos_spec_obs")) {
                warning("unknown frame type '", obj, "' with template '", tpl, "'");
                continue;
            }

            obj = "sci";
        }

        spawn("mv "+dir+file+" "+dir+file::remove_extension(file)+"-"+obj+".fits");
    }

    return 0;
}

