#include <phypp.hpp>

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print("usage: fill_nan <cube>");
        return 0;
    }

    vec1d lines = {
        // H band
        1.4564, 1.4606, 1.4666, 1.4699, 1.4784, 1.4832, 1.4864, 1.4888,
        1.4909, 1.4932, 1.4803, 1.5054, 1.5057, 1.5069, 1.5088, 1.5188,
        1.5241, 1.5287, 1.5332, 1.5395, 1.5433, 1.5505, 1.5515, 1.5541,
        1.5544, 1.5570, 1.5598, 1.5632, 1.5656, 1.5702, 1.5832, 1.5848,
        1.5869, 1.5973, 1.6031, 1.6080, 1.6129, 1.6195, 1.6236, 1.6317,
        1.6352, 1.6388, 1.6415, 1.6444, 1.6477, 1.6503, 1.6555, 1.6611,
        1.6690, 1.6693, 1.6706, 1.6709, 1.6733, 1.6841, 1.6904, 1.6955,
        1.7009, 1.7079, 1.7124, 1.7211, 1.7249, 1.7283, 1.7330, 1.7350,
        1.7358, 1.7386, 1.7428, 1.7450, 1.7505, 1.7652, 1.7810, 1.7879,
        1.7991, 1.7995, 1.8067, 1.8120, 1.8209, 1.8216,

        // H+K gap
        1.8645, 1.8809, 1.8843,

        // K band
        1.9200, 1.9207, 1.9246, 1.9250, 1.9277, 1.9283, 1.9309, 1.9350,
        1.9399, 1.9560, 1.9642, 1.9701, 1.9753, 1.9774, 1.9841, 2.0010,
        2.0034, 2.0278, 2.0342, 2.0414, 2.0500, 2.0566, 2.0731, 2.0862,
        2.0909, 2.1177, 2.1233, 2.1252, 2.1509, 2.1541, 2.1710, 2.1804,
        2.1874, 2.1957, 2.2125
    };

    vec1d bands_low = {1.8256, 1.8678, 1.8898, 1.8981, 1.9109};
    vec1d bands_up  = {1.8609, 1.8779, 1.8939, 1.9072, 1.9169};

    vec3d cube;
    fits::image img(argv[1]);
    for (uint_t i : range(img.hdu_count())) {
        img.reach_hdu(i);
        if (img.axis_count() == 3 && total(img.image_dims()) != 0) {
            img.read(cube);
            break;
        }
    }

    astro::wcs wcs(img.read_header());
    vec1d lam = astro::build_axis(wcs, 0, astro::axis_unit::wave_um);
    double cdelt = abs(lam[1]-lam[0]);

    vec1b flagged(lam.dims);
    for (uint_t l : range(lines)) {
        flagged = flagged || (lam >= lines[l]-0.0006 && lam <= lines[l]+0.0006);
    }
    for (uint_t l : range(bands_low)) {
        flagged = flagged || (lam >= bands_low[l]-cdelt && lam <= bands_up[l]+cdelt);
    }

    for (uint_t l : range(cube.dims[0])) {
        if (flagged[l]) {
            cube(l,_,_) = dnan;
        } else {
            auto mima = minmax(cube(l,_,_));
            if (mima.first == mima.second) {
                cube(l,_,_) = dnan;
            }
        }
    }

    img.update(cube);

    return 0;
}
