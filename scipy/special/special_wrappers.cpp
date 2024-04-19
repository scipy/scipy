#include "special_wrappers.h"
#include "special.h"
#include "special/airy.h"
#include "special/amos.h"
#include "special/bessel.h"
#include "special/expint.h"
#include "special/fresnel.h"
#include "special/gamma.h"
#include "special/hyp2f1.h"
#include "special/kelvin.h"
#include "special/lambertw.h"
#include "special/log_exp.h"
#include "special/loggamma.h"
#include "special/mathieu.h"
#include "special/par_cyl.h"
#include "special/specfun.h"
#include "special/sph_bessel.h"
#include "special/sph_harm.h"
#include "special/sphd_wave.h"
#include "special/struve.h"
#include "special/trig.h"
#include "special/wright_bessel.h"

using namespace std;

namespace {

complex<double> to_complex(npy_cdouble z) {
    union {
        npy_cdouble cvalue;
        complex<double> value;
    } z_union{z};
    return z_union.value;
}

npy_cdouble to_ccomplex(complex<double> z) {
    union {
        complex<double> value;
        npy_cdouble cvalue;
    } z_union{z};
    return z_union.cvalue;
}

} // namespace

npy_cdouble chyp2f1_wrap(double a, double b, double c, npy_cdouble z) {
    return to_ccomplex(special::chyp2f1(a, b, c, to_complex(z)));
}

npy_cdouble chyp1f1_wrap(double a, double b, npy_cdouble z) {
    return to_ccomplex(special::chyp1f1(a, b, to_complex(z)));
}

double hypU_wrap(double a, double b, double x) { return special::hypu(a, b, x); }

double hyp1f1_wrap(double a, double b, double x) { return special::hyp1f1(a, b, x); }

void special_itairy(double x, double *apt, double *bpt, double *ant, double *bnt) {
    special::itairy(x, *apt, *bpt, *ant, *bnt);
}

double special_exp1(double x) { return special::exp1(x); }

npy_cdouble special_cexp1(npy_cdouble z) { return to_ccomplex(special::exp1(to_complex(z))); }

double special_expi(double x) { return special::expi(x); }

npy_cdouble special_cexpi(npy_cdouble z) { return to_ccomplex(special::expi(to_complex(z))); }

npy_double special_exprel(npy_double x) { return special::exprel(x); }

npy_cdouble cerf_wrap(npy_cdouble z) { return to_ccomplex(special::cerf(to_complex(z))); }

double special_itstruve0(double x) { return special::itstruve0(x); }

double special_it2struve0(double x) { return special::it2struve0(x); }

double special_itmodstruve0(double x) { return special::itmodstruve0(x); }

double special_ber(double x) { return special::ber(x); }

double special_bei(double x) { return special::bei(x); }

double special_ker(double x) { return special::ker(x); }

double special_kei(double x) { return special::kei(x); }

double special_berp(double x) { return special::berp(x); }

double special_beip(double x) { return special::beip(x); }

double special_kerp(double x) { return special::kerp(x); }

double special_keip(double x) { return special::keip(x); }

void special_ckelvin(double x, npy_cdouble *Be, npy_cdouble *Ke, npy_cdouble *Bep, npy_cdouble *Kep) {
    special::kelvin(
        x, *reinterpret_cast<complex<double> *>(Be), *reinterpret_cast<complex<double> *>(Ke),
        *reinterpret_cast<complex<double> *>(Bep), *reinterpret_cast<complex<double> *>(Kep)
    );
}

npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble z) {
    return to_ccomplex(special::hyp2f1(a, b, c, to_complex(z)));
}

void it1j0y0_wrap(double x, double *j0int, double *y0int) { special::it1j0y0(x, *j0int, *y0int); }

void it2j0y0_wrap(double x, double *j0int, double *y0int) { special::it2j0y0(x, *j0int, *y0int); }

void it1i0k0_wrap(double x, double *i0int, double *k0int) { special::it1i0k0(x, *i0int, *k0int); }

void it2i0k0_wrap(double x, double *i0int, double *k0int) { special::it2i0k0(x, *i0int, *k0int); }

int cfresnl_wrap(npy_cdouble z, npy_cdouble *zfs, npy_cdouble *zfc) {
    special::cfresnl(to_complex(z), reinterpret_cast<complex<double> *>(zfs), reinterpret_cast<complex<double> *>(zfc));
    return 0;
}

double cem_cva_wrap(double m, double q) { return special::cem_cva(m, q); }

double sem_cva_wrap(double m, double q) { return special::sem_cva(m, q); }

void cem_wrap(double m, double q, double x, double *csf, double *csd) { special::cem(m, q, x, *csf, *csd); }

void sem_wrap(double m, double q, double x, double *csf, double *csd) { special::sem(m, q, x, *csf, *csd); }

void mcm1_wrap(double m, double q, double x, double *f1r, double *d1r) { special::mcm1(m, q, x, *f1r, *d1r); }

void msm1_wrap(double m, double q, double x, double *f1r, double *d1r) { special::msm1(m, q, x, *f1r, *d1r); }

void mcm2_wrap(double m, double q, double x, double *f2r, double *d2r) { special::mcm2(m, q, x, *f2r, *d2r); }

void msm2_wrap(double m, double q, double x, double *f2r, double *d2r) { special::msm2(m, q, x, *f2r, *d2r); }

double pmv_wrap(double m, double v, double x) { return special::pmv(m, v, x); }

void pbwa_wrap(double a, double x, double *wf, double *wd) { special::pbwa(a, x, *wf, *wd); }

void pbdv_wrap(double v, double x, double *pdf, double *pdd) { special::pbdv(v, x, *pdf, *pdd); }

void pbvv_wrap(double v, double x, double *pvf, double *pvd) { special::pbvv(v, x, *pvf, *pvd); }

double prolate_segv_wrap(double m, double n, double c) { return special::prolate_segv(m, n, c); }

double oblate_segv_wrap(double m, double n, double c) { return special::oblate_segv(m, n, c); }

double prolate_aswfa_nocv_wrap(double m, double n, double c, double x, double *s1d) {
    double s1f;
    special::prolate_aswfa_nocv(m, n, c, x, s1f, *s1d);

    return s1f;
}

double oblate_aswfa_nocv_wrap(double m, double n, double c, double x, double *s1d) {
    double s1f;
    special::oblate_aswfa_nocv(m, n, c, x, s1f, *s1d);

    return s1f;
}

void prolate_aswfa_wrap(double m, double n, double c, double cv, double x, double *s1f, double *s1d) {
    special::prolate_aswfa(m, n, c, cv, x, *s1f, *s1d);
}

void oblate_aswfa_wrap(double m, double n, double c, double cv, double x, double *s1f, double *s1d) {
    special::oblate_aswfa(m, n, c, cv, x, *s1f, *s1d);
}

double prolate_radial1_nocv_wrap(double m, double n, double c, double x, double *r1d) {
    double r1f;
    special::prolate_radial1_nocv(m, n, c, x, r1f, *r1d);

    return r1f;
}

double prolate_radial2_nocv_wrap(double m, double n, double c, double x, double *r2d) {
    double r2f;
    special::prolate_radial2_nocv(m, n, c, x, r2f, *r2d);

    return r2f;
}

void prolate_radial1_wrap(double m, double n, double c, double cv, double x, double *r1f, double *r1d) {
    special::prolate_radial1(m, n, c, cv, x, *r1f, *r1d);
}

void prolate_radial2_wrap(double m, double n, double c, double cv, double x, double *r2f, double *r2d) {
    special::prolate_radial2(m, n, c, cv, x, *r2f, *r2d);
}

double oblate_radial1_nocv_wrap(double m, double n, double c, double x, double *r1d) {
    double r1f;
    special::oblate_radial1_nocv(m, n, c, x, r1f, *r1d);

    return r1f;
}

double oblate_radial2_nocv_wrap(double m, double n, double c, double x, double *r2d) {
    double r2f;
    special::oblate_radial2_nocv(m, n, c, x, r2f, *r2d);

    return r2f;
}

void oblate_radial1_wrap(double m, double n, double c, double cv, double x, double *r1f, double *r1d) {
    special::oblate_radial1(m, n, c, cv, x, *r1f, *r1d);
}

void oblate_radial2_wrap(double m, double n, double c, double cv, double x, double *r2f, double *r2d) {
    special::oblate_radial2(m, n, c, cv, x, *r2f, *r2d);
}

void modified_fresnel_plus_wrap(double x, npy_cdouble *Fplus, npy_cdouble *Kplus) {
    special::modified_fresnel_plus(
        x, *reinterpret_cast<complex<double> *>(Fplus), *reinterpret_cast<complex<double> *>(Kplus)
    );
}

void modified_fresnel_minus_wrap(double x, npy_cdouble *Fminus, npy_cdouble *Kminus) {
    special::modified_fresnel_minus(
        x, *reinterpret_cast<complex<double> *>(Fminus), *reinterpret_cast<complex<double> *>(Kminus)
    );
}

double special_sinpi(double x) { return special::sinpi(x); }

npy_cdouble special_csinpi(npy_cdouble z) { return to_ccomplex(special::sinpi(to_complex(z))); }

double special_cospi(double x) { return special::cospi(x); }

void special_airy(double x, double *ai, double *aip, double *bi, double *bip) {
    special::airy(x, *ai, *aip, *bi, *bip);
}

void special_cairy(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip) {
    special::airy(
        to_complex(z), *reinterpret_cast<complex<double> *>(ai), *reinterpret_cast<complex<double> *>(aip),
        *reinterpret_cast<complex<double> *>(bi), *reinterpret_cast<complex<double> *>(bip)
    );
}

void special_airye(double z, double *ai, double *aip, double *bi, double *bip) {
    special::airye(z, *ai, *aip, *bi, *bip);
}

void special_cairye(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip) {
    special::airye(
        to_complex(z), *reinterpret_cast<complex<double> *>(ai), *reinterpret_cast<complex<double> *>(aip),
        *reinterpret_cast<complex<double> *>(bi), *reinterpret_cast<complex<double> *>(bip)
    );
}

double special_cyl_bessel_j(double v, double x) { return special::cyl_bessel_j(v, x); }

npy_cdouble special_ccyl_bessel_j(double v, npy_cdouble z) {
    return to_ccomplex(special::cyl_bessel_j(v, to_complex(z)));
}

double special_cyl_bessel_je(double v, double z) { return special::cyl_bessel_je(v, z); }

npy_cdouble special_ccyl_bessel_je(double v, npy_cdouble z) {
    return to_ccomplex(special::cyl_bessel_je(v, to_complex(z)));
}

double special_cyl_bessel_y(double v, double x) { return special::cyl_bessel_y(v, x); }

npy_cdouble special_ccyl_bessel_y(double v, npy_cdouble z) {
    return to_ccomplex(special::cyl_bessel_y(v, to_complex(z)));
}

double special_cyl_bessel_ye(double v, double z) { return special::cyl_bessel_ye(v, z); }

npy_cdouble special_ccyl_bessel_ye(double v, npy_cdouble z) {
    return to_ccomplex(special::cyl_bessel_ye(v, to_complex(z)));
}

double special_cyl_bessel_i(double v, double z) { return special::cyl_bessel_i(v, z); }

npy_cdouble special_ccyl_bessel_i(double v, npy_cdouble z) {
    return to_ccomplex(special::cyl_bessel_i(v, to_complex(z)));
}

double special_cyl_bessel_ie(double v, double z) { return special::cyl_bessel_ie(v, z); }

npy_cdouble special_ccyl_bessel_ie(double v, npy_cdouble z) {
    return to_ccomplex(special::cyl_bessel_ie(v, to_complex(z)));
}

double special_cyl_bessel_k_int(int n, double z) { return special::cyl_bessel_k(static_cast<double>(n), z); }

double special_cyl_bessel_k(double v, double z) { return special::cyl_bessel_k(v, z); }

npy_cdouble special_ccyl_bessel_k(double v, npy_cdouble z) {
    return to_ccomplex(special::cyl_bessel_k(v, to_complex(z)));
}

double special_cyl_bessel_ke(double v, double z) { return special::cyl_bessel_ke(v, z); }

npy_cdouble special_ccyl_bessel_ke(double v, npy_cdouble z) {
    return to_ccomplex(special::cyl_bessel_ke(v, to_complex(z)));
}

npy_cdouble special_ccyl_hankel_1(double v, npy_cdouble z) {
    return to_ccomplex(special::cyl_hankel_1(v, to_complex(z)));
}

npy_cdouble special_ccyl_hankel_1e(double v, npy_cdouble z) {
    return to_ccomplex(special::cyl_hankel_1e(v, to_complex(z)));
}

npy_cdouble special_ccyl_hankel_2(double v, npy_cdouble z) {
    return to_ccomplex(special::cyl_hankel_2(v, to_complex(z)));
}

npy_cdouble special_ccyl_hankel_2e(double v, npy_cdouble z) {
    return to_ccomplex(special::cyl_hankel_2e(v, to_complex(z)));
}

double special_binom(double n, double k) { return special::binom(n, k); }

double special_digamma(double z) { return special::digamma(z); }

npy_cdouble special_cdigamma(npy_cdouble z) { return to_ccomplex(special::digamma(to_complex(z))); }

double special_gamma(double x) { return special::gamma(x); }

npy_cdouble special_cgamma(npy_cdouble z) { return to_ccomplex(special::gamma(to_complex(z))); }

double special_rgamma(double x) { return special::rgamma(x); }

npy_cdouble special_crgamma(npy_cdouble z) { return to_ccomplex(special::rgamma(to_complex(z))); }

float special_expitf(float x) { return special::expit(x); };

double special_expit(double x) { return special::expit(x); };

npy_longdouble special_expitl(npy_longdouble x) { return special::expit(x); };

float special_log_expitf(float x) { return special::log_expit(x); };

double special_log_expit(double x) { return special::log_expit(x); };

npy_longdouble special_log_expitl(npy_longdouble x) { return special::log_expit(x); };

float special_logitf(float x) { return special::logit(x); };

double special_logit(double x) { return special::logit(x); };

npy_longdouble special_logitl(npy_longdouble x) { return special::logit(x); };

double special_loggamma(double x) { return special::loggamma(x); }

npy_cdouble special_cloggamma(npy_cdouble z) { return to_ccomplex(special::loggamma(to_complex(z))); }

double special_hyp2f1(double a, double b, double c, double z) { return special::hyp2f1(a, b, c, z); }

npy_cdouble special_chyp2f1(double a, double b, double c, npy_cdouble z) {
    return to_ccomplex(special::hyp2f1(a, b, c, to_complex(z)));
}

npy_cdouble special_lambertw(npy_cdouble z, long k, double tol) {
    return to_ccomplex(special::lambertw(to_complex(z), k, tol));
}

npy_cdouble special_sph_harm(long m, long n, double theta, double phi) {
    return to_ccomplex(::sph_harm(m, n, theta, phi));
}

npy_cdouble special_sph_harm_unsafe(double m, double n, double theta, double phi) {
    return to_ccomplex(::sph_harm(static_cast<long>(m), static_cast<long>(n), theta, phi));
}

double special_wright_bessel(double a, double b, double x) { return special::wright_bessel(a, b, x); }

double special_scaled_exp1(double x) { return special::scaled_exp1(x); }

double special_sph_bessel_j(long n, double x) { return special::sph_bessel_j(n, x); }

npy_cdouble special_csph_bessel_j(long n, npy_cdouble z) {
    return to_ccomplex(special::sph_bessel_j(n, to_complex(z)));
}

double special_sph_bessel_j_jac(long n, double x) { return special::sph_bessel_j_jac(n, x); }

npy_cdouble special_csph_bessel_j_jac(long n, npy_cdouble z) {
    return to_ccomplex(special::sph_bessel_j_jac(n, to_complex(z)));
}

double special_sph_bessel_y(long n, double x) { return special::sph_bessel_y(n, x); }

npy_cdouble special_csph_bessel_y(long n, npy_cdouble z) {
    return to_ccomplex(special::sph_bessel_y(n, to_complex(z)));
}

double special_sph_bessel_y_jac(long n, double x) { return special::sph_bessel_y_jac(n, x); }

npy_cdouble special_csph_bessel_y_jac(long n, npy_cdouble z) {
    return to_ccomplex(special::sph_bessel_y_jac(n, to_complex(z)));
}

double special_sph_bessel_i(long n, double x) { return special::sph_bessel_i(n, x); }

npy_cdouble special_csph_bessel_i(long n, npy_cdouble z) {
    return to_ccomplex(special::sph_bessel_i(n, to_complex(z)));
}

double special_sph_bessel_i_jac(long n, double x) { return special::sph_bessel_i_jac(n, x); }

npy_cdouble special_csph_bessel_i_jac(long n, npy_cdouble z) {
    return to_ccomplex(special::sph_bessel_i_jac(n, to_complex(z)));
}

double special_sph_bessel_k(long n, double x) { return special::sph_bessel_k(n, x); }

npy_cdouble special_csph_bessel_k(long n, npy_cdouble z) {
    return to_ccomplex(special::sph_bessel_k(n, to_complex(z)));
}

double special_sph_bessel_k_jac(long n, double x) { return special::sph_bessel_k_jac(n, x); }

npy_cdouble special_csph_bessel_k_jac(long n, npy_cdouble z) {
    return to_ccomplex(special::sph_bessel_k_jac(n, to_complex(z)));
}
