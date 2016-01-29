#include <phypp.hpp>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print("usage: fill_nan <cube>");
        return 0;
    }

    vec3f cube;
    fits::image img(argv[1]);
    img.read(cube);

    for (uint_t l : range(cube.dims[0])) {
        auto mima = minmax(cube(l,_,_));
        if (mima.first == mima.second) {
            cube(l,_,_) = fnan;
        }
    }

    img.update(cube);

    return 0;
}
