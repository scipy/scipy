#pragma once

#include <complex>
#include "npy_2_npymathcompat.h"

#include "special/binom.h"
#include "special/lambertw.h"

inline double binom(double n, double k) {
    return special::binom(n, k);
}

inline npy_cdouble lambertw_scalar(npy_cdouble zp, long k, double tol) {
    std::complex<double> z(npymath_creal(zp), npymath_cimag(zp));
    std::complex<double> w = special::lambertw(z, k, tol);
    return npymath_cpack(real(w), imag(w));
}
