#pragma once

#include <numpy/npy_math.h>
#include <complex>


#include "special/binom.h"
#include "special/hyp2f1.h"
#include "special/lambertw.h"
#include "special/loggamma.h"
#include "special/trig.h"
#include "special/digamma.h"
#include "special/wright_bessel.h"


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
    std::complex<double> w = special::sinpi(z);
    return npy_cpack(real(w), imag(w));
}


inline npy_cdouble ccospi(npy_cdouble zp) {
    std::complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = special::cospi(z);
    return npy_cpack(real(w), imag(w));
}


inline npy_cdouble cgamma(npy_cdouble zp) {
    std::complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = special::gamma(z);
    return npy_cpack(real(w), imag(w));
}

inline npy_cdouble crgamma(npy_cdouble zp) {
    std::complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = special::rgamma(z);
    return npy_cpack(real(w), imag(w));
}

inline npy_cdouble loggamma(npy_cdouble zp) {
    std::complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = special::loggamma(z);
    return npy_cpack(real(w), imag(w));
}

inline double loggamma_real(double x) {
    return special::loggamma(x);
}

						    
inline double digamma(double z) {
    return special::digamma(z);
}

inline npy_cdouble cdigamma(npy_cdouble zp) {
    std::complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = special::digamma(z);
    return npy_cpack(real(w), imag(w));
}

inline npy_cdouble hyp2f1_complex(double a, double b, double c, npy_cdouble zp) {
    std::complex<double> z(npy_creal(zp), npy_cimag(zp));
    std::complex<double> w = special::hyp2f1(a, b, c, z);
    return npy_cpack(real(w), imag(w));
}

inline double wright_bessel_scalar(double a, double b, double x) {
    return special::wright_bessel(a, b, x);
}

inline double rgamma(double x) {
    return special::cephes::rgamma(x);
}
