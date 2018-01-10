#include <phypp.hpp>

void print_help();

int phypp_main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    uint_t continuum_width = 350; // in spectral element

    read_args(argc-1, argv+1, arg_list(continuum_width));

    // Copy and open cube
    std::string infile = argv[1];
    std::string outfile = file::remove_extension(file::get_basename(infile))+"_contsub.fits";
    file::copy(infile, outfile);

    vec3d cflx, cerr;
    fits::image fimg(outfile);
    fimg.reach_hdu(1);
    fimg.read(cflx);
    fimg.reach_hdu(1);
    fimg.read(cerr);

    // Estimate and subtract the continuum emission
    for (uint_t y : range(cflx.dims[1]))
    for (uint_t x : range(cflx.dims[2])) {
        vec1d tflx = cflx.safe(_,y,x);
        // vec1d twei = 1.0/cerr(_,y,x);
        vec1d twei = 0.0*cerr.safe(_,y,x) + 1.0;

        for (uint_t l : range(cflx.dims[0])) {
            uint_t l0 = max(0, int_t(l)-int_t(continuum_width/2));
            uint_t l1 = min(cflx.dims[0]-1, l+continuum_width/2);
            cflx.safe(l,y,x) -= weighted_median(tflx.safe[l0-_-l1], twei.safe[l0-_-l1]);
        }
    }

    fimg.reach_hdu(1);
    fimg.update(cflx);

    return 0;
}

void print_help() {
    using namespace format;

    print("contsub v1.0");
    print("usage: contsub <kmos_cube.fits> continuum_width=...");
    print("");
    print("Main parameter:");
    paragraph("'kmos_cube.fits' should be a cube created by the KMOS pipeline, with 3 "
        "extensions: the first is empty (KMOS convention), the second contains the flux "
        "and the third contains the uncertainty. This program will then compute the "
        "weighted median flux of each pixel along a wide wavelength range to estimate "
        "the flux of the continuum, and subtract it from the cube.");
    print("Available option:");
    bullet("continuum_width=...", "Must be an integer. It defines the width of the "
        "wavelength range within which the continuum level is estimated, for each pixel "
        "and each wavelength element. It must be given in number of wavelength elements. "
        "Default is 350. Be cautious not to pick too small a value, else you may start to "
        "subtract part of the flux of your emission lines. You will, however, obtain a "
        "higher resolution spectrum of the continuum.");
}
