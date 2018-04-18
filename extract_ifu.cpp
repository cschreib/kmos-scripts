#include <phypp.hpp>

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print("usage: extract_ifu <file> names=[...]");
        return 0;
    }

    fits::image fimg(argv[1]);
    fimg.reach_hdu(0);

    fits::header phdr = fimg.read_header();

    vec1s names;
    read_args(argc-1, argv+1, arg_list(names));

    names = to_lower(names);

    vec1b found(names.size());
    for (uint_t i : range(24)) {
        fimg.reach_hdu(0);

        std::string name;
        fimg.read_keyword("HIERARCH ESO OCS ARM"+strn(i+1)+" NAME", name);
        name = to_lower(trim(name, " '"));

        uint_t k = where_first(name == names);
        if (k == npos) continue;

        fimg.reach_hdu(i+1);
        if (fimg.axis_count() == 0) continue;

        fits::image out(file::remove_extension(argv[1])+"-"+name+".fits");
        fits::header hdr = fimg.read_header();

        if (fimg.axis_count() == 3) {
            vec3d cube;
            fimg.read(cube);
            out.write(cube);
            out.write_header(hdr);
        } else if (fimg.axis_count() == 2) {
            vec2d img;
            fimg.read(img);
            out.write(img);
            out.write_header(hdr);
        }

        found[k] = true;
    }

    for (uint_t i : range(found)) {
        if (!found[i]) warning("could not find '", names[i], "' in '", argv[1], "'");
    }

    return 0;
}
