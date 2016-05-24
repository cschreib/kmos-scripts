#include <phypp.hpp>

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print("usage: fill_nan <cube>");
        return 0;
    }

    vec3d cube;
    fits::image img(argv[1]);
    for (uint_t i : range(img.hdu_count())) {
        img.reach_hdu(i);
        if (img.axis_count() == 3 && total(img.image_dims()) != 0) {
            img.read(cube);
            break;
        }
    }

    for (uint_t l : range(cube.dims[0])) {
        auto mima = minmax(cube(l,_,_));
        if (mima.first == mima.second) {
            cube(l,_,_) = dnan;
        }
    }

    img.update(cube);

    return 0;
}
