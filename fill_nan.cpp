#include <phypp.hpp>

int main(int argc, char* argv[]) {
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
