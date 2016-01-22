#include <phypp.hpp>

int main(int argc, char* argv[]) {
    std::string dir = file::directorize(argv[1]);
    vec1s files = dir+file::list_files(dir+"sci_reconstructed*-sci.fits");

    for (auto oname : files) {
        std::string newfile = file::get_basename(oname);
        file::copy(oname, newfile);
        fits::image fimg(newfile);
        for (uint_t i : range(1, fimg.hdu_count())) {
            fimg.reach_hdu(i);
            if (fimg.axis_count() != 3) continue;

            vec3d cube;
            fimg.read(cube);
            for (uint_t l : range(cube.dims[0])) {
                cube(l,_,_) -= median(cube(l,_,_));
            }

            fimg.update(cube);
        }
    }

    return 0;
}
