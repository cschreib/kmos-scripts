#include <phypp.hpp>
#include <phypp/math/mpfit.hpp>

int phypp_main(int argc, char* argv[]) {
    vec1d lam, sed;
    ascii::read_table(argv[1], 0, lam, sed);

    lam *= 1e-4;

    vec1f lines = {0.41017, 0.43405, 0.48613};
    double dl = 0.02;
    double width1 = 200.0;

    read_args(argc-1, argv+1, arg_list(lines, dl, width1));

    for (uint_t l : range(lines)) {
        vec1u idl = where(abs(lam - lines[l])/lines[l] < dl);

        auto model = [](const vec1d& x, const vec1d& p) {
            return p[0] + p[1]*(x - p[2]) + p[3]*exp(-sqr(x - p[2])/(2.0*sqr(p[4]))) +  p[5]*exp(-sqr(x - p[2])/(2.0*sqr(p[6])));
        };

        vec1d start_p = {mean(sed[idl]), 0.0, lines[l], -mean(sed[idl]), lines[l]*400.0/3e5, -mean(sed[idl]), lines[l]*100.0/3e5};

        auto p = mpfitfun(sed[idl], 1.0, lam[idl], model, start_p);

        p.params[0] = 0; p.params[1] = 0;
        sed -= model(lam, p.params);

        print(lines[l], ": " , p.params);
    }

    lam /= 1e-4;

    file::mkdir("../templates_noabs");
    ascii::write_table("../templates_noabs/"+file::remove_extension(argv[1])+".dat", 18, lam, sed);

    return 0;
}
