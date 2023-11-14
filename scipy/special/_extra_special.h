#pragma once

#include <complex>
#include <numpy/npy_math.h>

#include "_binom.h"
#include "extra_special/lambertw.h"

inline double binom(double n, double k) {
    return extra_special::binom(n, k);
}

inline npy_cdouble lambertw_scalar(npy_cdouble zp, long k, double tol) {
    std::complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = extra_special::lambertw(z, k, tol);
    return npy_cpack(real(w), imag(w));
}
