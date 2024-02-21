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

#include "special/cephes/besselpoly.h"
#include "special/cephes/beta.h"
#include "special/cephes/gamma.h"
#include "special/cephes/hyp2f1.h"
#include "special/cephes/i0.h"
#include "special/cephes/i1.h"
#include "special/cephes/lanczos.h"
#include "special/cephes/poch.h"
#include "special/cephes/rgamma.h"
#include "special/cephes/trig.h"
#include "special/cephes/zeta.h"


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


// Special functions from cephes

inline double besselpoly(double a, double lambda, double nu) {
    return special::cephes::besselpoly(a, lambda, nu);
}

inline double beta(double a, double b) {
    return special::cephes::beta(a, b);
}

inline double lbeta(double a, double b) {
    return special::cephes::lbeta(a, b);
}

inline double sinpi(double x) {
    return special::cephes::sinpi(x);
}

inline double cospi(double x) {
    return special::cephes::cospi(x);
}

inline double Gamma(double x) {
    return special::cephes::Gamma(x);
}

inline double gammasgn(double x) {
    return special::cephes::gammasgn(x);
}

inline double lgam(double x) {
    return special::cephes::lgam(x);
}

inline double hyp2f1_real(double a, double b, double c, double x) {
    return special::cephes::hyp2f1(a, b, c, x);
}

inline double i0(double x) {
    return special::cephes::i0(x);
}

inline double i0e(double x) {
    return special::cephes::i0e(x);
}

inline double i1(double x) {
    return special::cephes::i1(x);
}

inline double i1e(double x) {
    return special::cephes::i1e(x);
}

inline double lanczos_sum_expg_scaled(double x) {
    return special::cephes::lanczos_sum_expg_scaled(x);
}

inline double poch(double x, double m) {
    return special::cephes::poch(x, m);
}

inline double rgamma(double x) {
    return special::cephes::rgamma(x);
}

inline double zeta(double x, double q) {
    return special::cephes::zeta(x, q);
}
