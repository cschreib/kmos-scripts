#include <phypp.hpp>

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print("usage: median_sub <directory> [options]");
        return 0;
    }

    std::string dir = file::directorize(argv[1]);
    vec1s files = dir+file::list_files(dir+"sci_reconstructed*-sci.fits");
    if (files.empty()) {
        files = dir+file::list_files(dir+"SCI_RECONSTRUCTED*-sci.fits");
    }

    bool mask_center = false;
    double mask_radius = 0.5; // arcsec
    std::string masks;
    double cx = dnan, cy = dnan;
    read_args(argc-1, argv+1, arg_list(mask_center, masks, mask_radius, cx, cy));

    vec1d ra, dec, size;
    if (!masks.empty()) {
        ascii::read_table(masks, ascii::find_skip(masks), _, ra, dec, size);
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
            std::string extname;
            if (fimg.read_keyword("EXTNAME", extname) && ends_with(extname, ".NOISE")) continue;

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
                double aspix = 1.0;
                if (!astro::get_pixel_size(w, aspix)) {
                    error("could not read pixel size from cube");
                    return 1;
                }

                vec1d r = size/aspix;

                vec2d ix = generate_img(mask.dims, [](int_t,    int_t tx) { return tx; });
                vec2d iy = generate_img(mask.dims, [](int_t ty, int_t)    { return ty; });

                for (uint_t s : range(ra)) {
                    mask = mask && sqr(x[s] - ix) + sqr(y[s] - iy) > sqr(r[s]);
                }
            } else if (mask_center) {
                astro::wcs w(fimg.read_header());
                double aspix = 1.0;
                if (!astro::get_pixel_size(w, aspix)) {
                    error("could not read pixel size from cube");
                    return 1;
                }

                vec2d ix = generate_img(mask.dims, [](int_t,    int_t tx) { return tx; });
                vec2d iy = generate_img(mask.dims, [](int_t ty, int_t)    { return ty; });

                double tcx = cx, tcy = cy;
                if (!is_finite(tcx)) tcx = nx/2;
                if (!is_finite(tcy)) tcy = ny/2;

                mask = mask && sqr(tcx - ix) + sqr(tcy - iy) > sqr(mask_radius/aspix);
            }

            // Mask borders
            mask(0,_) = mask(ny-1,_) = mask(_,0) = mask(_,nx-1) = false;

            // Subtract median per wavelength slice
            vec1u idg = where(mask);
            for (uint_t l : range(cube.dims[0])) {
                cube(l,_,_) -= ksigma(cube(l,_,_)[idg]);
            }

            // Subtract median of total cube
            cube -= median(reform(cube, cube.dims[0], cube.dims[1]*cube.dims[2])(_,idg));

            // Save file
            fimg.update(cube);
        }
    }

    return 0;
}
