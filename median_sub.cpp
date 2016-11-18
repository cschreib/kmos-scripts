#include <phypp.hpp>

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print("usage: median_sub <directory> <mask_catalog>");
        return 0;
    }

    std::string dir = file::directorize(argv[1]);
    vec1s files = dir+file::list_files(dir+"sci_reconstructed*-sci.fits");
    if (files.empty()) {
        files = dir+file::list_files(dir+"SCI_RECONSTRUCTED*-sci.fits");
    }

    vec1d ra, dec, size;
    if (argc > 2) {
        vec1s id;
        ascii::read_table(argv[2], ascii::find_skip(argv[2]), id, ra, dec, size);
    }

    auto ksigma = [](vec1d data) {
        double v = median(data);
        double m = 1.48*median(abs(data - v));
        vec1u idg = where(abs(data - v) < 7*m);
        return mean(data[idg]);
    };

    for (auto oname : files) {
        std::string newfile = file::get_basename(oname);
        file::copy(oname, newfile);
        fits::image fimg(newfile);
        for (uint_t i : range(1, fimg.hdu_count())) {
            fimg.reach_hdu(i);
            if (fimg.axis_count() != 3) continue;

            vec3d cube;
            fimg.read(cube);

            uint_t ny = cube.dims[1];
            uint_t nx = cube.dims[2];

            vec2b mask = replicate(true, ny, nx);

            // Mask nearby sources from the provided catalog (if any)
            if (!ra.empty()) {
                astro::wcs w(fimg.read_header());
                vec1d x, y;
                astro::ad2xy(w, ra, dec, x, y);
                x -= 1.0; y -= 1.0;
                double aspix;
                astro::get_pixel_size(w, aspix);
                vec1d r = size/aspix;

                vec2d ix = generate_img(mask.dims, [](int_t,    int_t tx) { return tx; });
                vec2d iy = generate_img(mask.dims, [](int_t ty, int_t)    { return ty; });

                for (uint_t s : range(ra)) {
                    mask = mask && sqr(x[s] - ix) + sqr(y[s] - iy) > sqr(r[s]);
                }
            }

            // Mask borders
            mask(0,_) = mask(ny-1,_) = mask(_,0) = mask(_,nx-1) = false;

            // Subtract whole IFU
            vec1u idg = where(mask);
            for (uint_t l : range(cube.dims[0])) {
                cube(l,_,_) -= ksigma(cube(l,_,_)[idg]);
            }

            fimg.update(cube);
        }
    }

    return 0;
}
