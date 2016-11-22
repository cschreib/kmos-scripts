#include <phypp.hpp>

int phypp_main(int argc, char* argv[]) {
    std::string cs_dir = file::directorize(argv[1]);
    std::string nn_dir = file::directorize(argv[2]);

    vec1s files = file::list_files(cs_dir+"*.fits");

    for (auto& f : files) {
        print(f);
        file::copy(nn_dir+f, f);
        fits::image img(f);
        fits::input_image iimg(cs_dir+f);

        for (uint_t i : range(24)) {
            iimg.reach_hdu(i+1);

            if (iimg.image_dims().size() != 3) continue;

            vec3d cube;
            iimg.read(cube);

            img.reach_hdu(2*i+1);
            img.update(cube);
        }
    }

    return 0;
}
