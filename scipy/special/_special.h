#pragma once

#include <complex>
#include <numpy/npy_math.h>

#include "special/binom.h"
#include "special/lambertw.h"
#include "_trig.h"


inline double binom(double n, double k) {
    return special::binom(n, k);
}

inline npy_cdouble lambertw_scalar(npy_cdouble zp, long k, double tol) {
    std::complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = special::lambertw(z, k, tol);
    return npy_cpack(real(w), imag(w));
}

inline npy_cdouble csinpi(npy_cdouble zp) {
    std::complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = scipy::special::csinpi(z);
    return npy_cpack(real(w), imag(w));
}


inline npy_cdouble ccospi(npy_cdouble zp) {
    std::complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = scipy::special::ccospi(z);
    return npy_cpack(real(w), imag(w));
}
						    
