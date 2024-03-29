#pragma once

#include <numpy/npy_math.h>
#include <complex>


#include "special/binom.h"
#include "special/ellipk.h"
#include "special/hyp2f1.h"
#include "special/lambertw.h"
#include "special/loggamma.h"
#include "special/trig.h"
#include "special/digamma.h"
#include "special/wright_bessel.h"

#include "special/cephes/besselpoly.h"
#include "special/cephes/beta.h"
#include "special/cephes/cbrt.h"
#include "special/cephes/chdtr.h"
#include "special/cephes/ellie.h"
#include "special/cephes/ellik.h"
#include "special/cephes/ellpe.h"
#include "special/cephes/ellpk.h"
#include "special/cephes/ellpj.h"
#include "special/cephes/expn.h"
#include "special/cephes/gamma.h"
#include "special/cephes/hyp2f1.h"
#include "special/cephes/hyperg.h"
#include "special/cephes/i0.h"
#include "special/cephes/i1.h"
#include "special/cephes/scipy_iv.h"
#include "special/cephes/j0.h"
#include "special/cephes/j1.h"
#include "special/cephes/jv.h"
#include "special/cephes/k0.h"
#include "special/cephes/k1.h"
#include "special/cephes/igam.h"
#include "special/cephes/igami.h"
#include "special/cephes/lanczos.h"
#include "special/cephes/ndtr.h"
#include "special/cephes/poch.h"
#include "special/cephes/rgamma.h"
#include "special/cephes/trig.h"
#include "special/cephes/unity.h"
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

inline double ellipk(double m) {
    return special::ellipk(m);
}


// Special functions from cephes

inline double cephes_besselpoly(double a, double lambda, double nu) {
    return special::cephes::besselpoly(a, lambda, nu);
}

inline double cephes_beta(double a, double b) {
    return special::cephes::beta(a, b);
}

inline double cephes_chdtr(double df, double x) {
    return special::cephes::chdtr(df, x);
}

inline double cephes_chdtrc(double df, double x) {
    return special::cephes::chdtrc(df, x);
}

inline double cephes_chdtri(double df, double y) {
    return special::cephes::chdtri(df, y);
}

inline double cephes_lbeta(double a, double b) {
    return special::cephes::lbeta(a, b);
}

inline double cephes_sinpi(double x) {
    return special::cephes::sinpi(x);
}

inline double cephes_cospi(double x) {
    return special::cephes::cospi(x);
}

inline double cephes_cbrt(double x) {
    return special::cephes::detail::cbrt(x);
}

inline double cephes_Gamma(double x) {
    return special::cephes::Gamma(x);
}

inline double cephes_gammasgn(double x) {
    return special::cephes::gammasgn(x);
}

inline double cephes_lgam(double x) {
    return special::cephes::lgam(x);
}

inline double cephes_hyp2f1(double a, double b, double c, double x) {
    return special::cephes::hyp2f1(a, b, c, x);
}

inline double cephes_i0(double x) {
    return special::cephes::i0(x);
}

inline double cephes_i0e(double x) {
    return special::cephes::i0e(x);
}

inline double cephes_i1(double x) {
    return special::cephes::i1(x);
}

inline double cephes_i1e(double x) {
    return special::cephes::i1e(x);
}

inline double cephes_iv(double v, double x) {
    return special::cephes::iv(v, x);
}

inline double cephes_j0(double x) {
    return special::cephes::j0(x);
}

inline double cephes_j1(double x) {
    return special::cephes::j1(x);
}

inline double cephes_k0(double x) {
    return special::cephes::k0(x);
}

inline double cephes_k0e(double x) {
    return special::cephes::k0e(x);
}

inline double cephes_k1(double x) {
    return special::cephes::k1(x);
}

inline double cephes_k1e(double x) {
    return special::cephes::k1e(x);
}

inline double cephes_igam(double a, double x) {
    return special::cephes::igam(a, x);
}

inline double cephes_igamc(double a, double x) {
    return special::cephes::igamc(a, x);
}

inline double cephes_igami(double a, double p) {
    return special::cephes::igami(a, p);
}

inline double cephes_igamci(double a, double p) {
    return special::cephes::igamci(a, p);
}

inline double cephes_igam_fac(double a, double x) {
    return special::cephes::detail::igam_fac(a, x);
}

inline double cephes_lanczos_sum_expg_scaled(double x) {
    return special::cephes::lanczos_sum_expg_scaled(x);
}

inline double cephes_ndtr(double x) {
    return special::cephes::ndtr(x);
}

inline double cephes_erf(double x) {
    return special::cephes::erf(x);
}

inline double cephes_erfc(double x) {
    return special::cephes::erfc(x);
}

inline double cephes_poch(double x, double m) {
    return special::cephes::poch(x, m);
}

inline double cephes_rgamma(double x) {
    return special::cephes::rgamma(x);
}

inline double cephes_zeta(double x, double q) {
    return special::cephes::zeta(x, q);
}

inline double cephes_log1p(double x) {
    return special::cephes::log1p(x);
}

inline double cephes_log1pmx(double x) {
    return special::cephes::log1pmx(x);
}

inline double cephes_lgam1p(double x) {
    return special::cephes::lgam1p(x);
}

inline double cephes_expm1(double x) {
    return special::cephes::expm1(x);
}

inline double cephes_cosm1(double x) {
    return special::cephes::cosm1(x);
}

inline double cephes_expn(int n, double x) {
    return special::cephes::expn(n, x);
}

inline double cephes_ellpe(double x) {
    return special::cephes::ellpe(x);
}

inline double cephes_ellpk(double x) {
    return special::cephes::ellpk(x);
}

inline double cephes_ellie(double phi, double m) {
    return special::cephes::ellie(phi, m);
}

inline double cephes_ellik(double phi, double m) {
    return special::cephes::ellik(phi, m);
}
