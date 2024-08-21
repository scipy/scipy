#include "xsf_wrappers.h"
#include "xsf_special.h"
#include "xsf/airy.h"
#include "xsf/amos.h"
#include "xsf/bessel.h"
#include "xsf/binom.h"
#include "xsf/expint.h"
#include "xsf/fresnel.h"
#include "xsf/gamma.h"
#include "xsf/hyp2f1.h"
#include "xsf/kelvin.h"
#include "xsf/lambertw.h"
#include "xsf/log_exp.h"
#include "xsf/loggamma.h"
#include "xsf/mathieu.h"
#include "xsf/par_cyl.h"
#include "xsf/specfun.h"
#include "xsf/sph_bessel.h"
#include "xsf/sph_harm.h"
#include "xsf/sphd_wave.h"
#include "xsf/struve.h"
#include "xsf/trig.h"
#include "xsf/wright_bessel.h"

#include "xsf/binom.h"
#include "xsf/digamma.h"
#include "xsf/ellipk.h"
#include "xsf/gamma.h"
#include "xsf/hyp2f1.h"
#include "xsf/lambertw.h"
#include "xsf/loggamma.h"
#include "xsf/trig.h"
#include "xsf/wright_bessel.h"

#include "xsf/cephes/bdtr.h"
#include "xsf/cephes/besselpoly.h"
#include "xsf/cephes/beta.h"
#include "xsf/cephes/cbrt.h"
#include "xsf/cephes/chdtr.h"
#include "xsf/cephes/ellie.h"
#include "xsf/cephes/ellik.h"
#include "xsf/cephes/ellpe.h"
#include "xsf/cephes/ellpj.h"
#include "xsf/cephes/ellpk.h"
#include "xsf/cephes/erfinv.h"
#include "xsf/cephes/exp10.h"
#include "xsf/cephes/exp2.h"
#include "xsf/cephes/expn.h"
#include "xsf/cephes/fdtr.h"
#include "xsf/cephes/gamma.h"
#include "xsf/cephes/gdtr.h"
#include "xsf/cephes/hyp2f1.h"
#include "xsf/cephes/hyperg.h"
#include "xsf/cephes/i0.h"
#include "xsf/cephes/i1.h"
#include "xsf/cephes/igam.h"
#include "xsf/cephes/igami.h"
#include "xsf/cephes/incbet.h"
#include "xsf/cephes/incbi.h"
#include "xsf/cephes/j0.h"
#include "xsf/cephes/j1.h"
#include "xsf/cephes/jv.h"
#include "xsf/cephes/k0.h"
#include "xsf/cephes/k1.h"
#include "xsf/cephes/kolmogorov.h"
#include "xsf/cephes/lanczos.h"
#include "xsf/cephes/nbdtr.h"
#include "xsf/cephes/ndtr.h"
#include "xsf/cephes/ndtri.h"
#include "xsf/cephes/owens_t.h"
#include "xsf/cephes/pdtr.h"
#include "xsf/cephes/poch.h"
#include "xsf/cephes/rgamma.h"
#include "xsf/cephes/round.h"
#include "xsf/cephes/scipy_iv.h"
#include "xsf/cephes/sindg.h"
#include "xsf/cephes/spence.h"
#include "xsf/cephes/struve.h"
#include "xsf/cephes/tandg.h"
#include "xsf/cephes/trig.h"
#include "xsf/cephes/tukey.h"
#include "xsf/cephes/unity.h"
#include "xsf/cephes/yn.h"
#include "xsf/cephes/zeta.h"
#include "xsf/cephes/zetac.h"

#include "xsf/cephes/airy.h"
#include "xsf/cephes/bdtr.h"
#include "xsf/cephes/beta.h"
#include "xsf/cephes/ellpj.h"
#include "xsf/cephes/ellpk.h"
#include "xsf/cephes/expn.h"
#include "xsf/cephes/fresnl.h"
#include "xsf/cephes/gamma.h"
#include "xsf/cephes/hyp2f1.h"
#include "xsf/cephes/jv.h"
#include "xsf/cephes/kolmogorov.h"
#include "xsf/cephes/nbdtr.h"
#include "xsf/cephes/ndtr.h"
#include "xsf/cephes/ndtri.h"
#include "xsf/cephes/pdtr.h"
#include "xsf/cephes/shichi.h"
#include "xsf/cephes/sici.h"

using namespace std;

namespace {

complex<double> to_complex(npy_cdouble z) { return {npy_creal(z), npy_cimag(z)}; }

npy_cdouble to_ccomplex(complex<double> z) { return {z.real(), z.imag()}; }

} // namespace

npy_cdouble chyp2f1_wrap(double a, double b, double c, npy_cdouble z) {
    return to_ccomplex(xsf::chyp2f1(a, b, c, to_complex(z)));
}

npy_cdouble chyp1f1_wrap(double a, double b, npy_cdouble z) { return to_ccomplex(xsf::chyp1f1(a, b, to_complex(z))); }

double hypU_wrap(double a, double b, double x) { return xsf::hypu(a, b, x); }

double hyp1f1_wrap(double a, double b, double x) { return xsf::hyp1f1(a, b, x); }

void special_itairy(double x, double *apt, double *bpt, double *ant, double *bnt) {
    xsf::itairy(x, *apt, *bpt, *ant, *bnt);
}

double special_exp1(double x) { return xsf::exp1(x); }

npy_cdouble special_cexp1(npy_cdouble z) { return to_ccomplex(xsf::exp1(to_complex(z))); }

double special_expi(double x) { return xsf::expi(x); }

npy_cdouble special_cexpi(npy_cdouble z) { return to_ccomplex(xsf::expi(to_complex(z))); }

npy_double special_exprel(npy_double x) { return xsf::exprel(x); }

npy_cdouble cerf_wrap(npy_cdouble z) { return to_ccomplex(xsf::cerf(to_complex(z))); }

double special_itstruve0(double x) { return xsf::itstruve0(x); }

double special_it2struve0(double x) { return xsf::it2struve0(x); }

double special_itmodstruve0(double x) { return xsf::itmodstruve0(x); }

double special_ber(double x) { return xsf::ber(x); }

double special_bei(double x) { return xsf::bei(x); }

double special_ker(double x) { return xsf::ker(x); }

double special_kei(double x) { return xsf::kei(x); }

double special_berp(double x) { return xsf::berp(x); }

double special_beip(double x) { return xsf::beip(x); }

double special_kerp(double x) { return xsf::kerp(x); }

double special_keip(double x) { return xsf::keip(x); }

void special_ckelvin(double x, npy_cdouble *Be, npy_cdouble *Ke, npy_cdouble *Bep, npy_cdouble *Kep) {
    xsf::kelvin(x, *reinterpret_cast<complex<double> *>(Be), *reinterpret_cast<complex<double> *>(Ke),
                *reinterpret_cast<complex<double> *>(Bep), *reinterpret_cast<complex<double> *>(Kep));
}

npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble z) {
    return to_ccomplex(xsf::hyp2f1(a, b, c, to_complex(z)));
}

void it1j0y0_wrap(double x, double *j0int, double *y0int) { xsf::it1j0y0(x, *j0int, *y0int); }

void it2j0y0_wrap(double x, double *j0int, double *y0int) { xsf::it2j0y0(x, *j0int, *y0int); }

void it1i0k0_wrap(double x, double *i0int, double *k0int) { xsf::it1i0k0(x, *i0int, *k0int); }

void it2i0k0_wrap(double x, double *i0int, double *k0int) { xsf::it2i0k0(x, *i0int, *k0int); }

int cfresnl_wrap(npy_cdouble z, npy_cdouble *zfs, npy_cdouble *zfc) {
    xsf::cfresnl(to_complex(z), reinterpret_cast<complex<double> *>(zfs), reinterpret_cast<complex<double> *>(zfc));
    return 0;
}

double cem_cva_wrap(double m, double q) { return xsf::cem_cva(m, q); }

double sem_cva_wrap(double m, double q) { return xsf::sem_cva(m, q); }

void cem_wrap(double m, double q, double x, double *csf, double *csd) { xsf::cem(m, q, x, *csf, *csd); }

void sem_wrap(double m, double q, double x, double *csf, double *csd) { xsf::sem(m, q, x, *csf, *csd); }

void mcm1_wrap(double m, double q, double x, double *f1r, double *d1r) { xsf::mcm1(m, q, x, *f1r, *d1r); }

void msm1_wrap(double m, double q, double x, double *f1r, double *d1r) { xsf::msm1(m, q, x, *f1r, *d1r); }

void mcm2_wrap(double m, double q, double x, double *f2r, double *d2r) { xsf::mcm2(m, q, x, *f2r, *d2r); }

void msm2_wrap(double m, double q, double x, double *f2r, double *d2r) { xsf::msm2(m, q, x, *f2r, *d2r); }

double pmv_wrap(double m, double v, double x) { return xsf::pmv(m, v, x); }

void pbwa_wrap(double a, double x, double *wf, double *wd) { xsf::pbwa(a, x, *wf, *wd); }

void pbdv_wrap(double v, double x, double *pdf, double *pdd) { xsf::pbdv(v, x, *pdf, *pdd); }

void pbvv_wrap(double v, double x, double *pvf, double *pvd) { xsf::pbvv(v, x, *pvf, *pvd); }

double prolate_segv_wrap(double m, double n, double c) { return xsf::prolate_segv(m, n, c); }

double oblate_segv_wrap(double m, double n, double c) { return xsf::oblate_segv(m, n, c); }

double prolate_aswfa_nocv_wrap(double m, double n, double c, double x, double *s1d) {
    double s1f;
    xsf::prolate_aswfa_nocv(m, n, c, x, s1f, *s1d);

    return s1f;
}

double oblate_aswfa_nocv_wrap(double m, double n, double c, double x, double *s1d) {
    double s1f;
    xsf::oblate_aswfa_nocv(m, n, c, x, s1f, *s1d);

    return s1f;
}

void prolate_aswfa_wrap(double m, double n, double c, double cv, double x, double *s1f, double *s1d) {
    xsf::prolate_aswfa(m, n, c, cv, x, *s1f, *s1d);
}

void oblate_aswfa_wrap(double m, double n, double c, double cv, double x, double *s1f, double *s1d) {
    xsf::oblate_aswfa(m, n, c, cv, x, *s1f, *s1d);
}

double prolate_radial1_nocv_wrap(double m, double n, double c, double x, double *r1d) {
    double r1f;
    xsf::prolate_radial1_nocv(m, n, c, x, r1f, *r1d);

    return r1f;
}

double prolate_radial2_nocv_wrap(double m, double n, double c, double x, double *r2d) {
    double r2f;
    xsf::prolate_radial2_nocv(m, n, c, x, r2f, *r2d);

    return r2f;
}

void prolate_radial1_wrap(double m, double n, double c, double cv, double x, double *r1f, double *r1d) {
    xsf::prolate_radial1(m, n, c, cv, x, *r1f, *r1d);
}

void prolate_radial2_wrap(double m, double n, double c, double cv, double x, double *r2f, double *r2d) {
    xsf::prolate_radial2(m, n, c, cv, x, *r2f, *r2d);
}

double oblate_radial1_nocv_wrap(double m, double n, double c, double x, double *r1d) {
    double r1f;
    xsf::oblate_radial1_nocv(m, n, c, x, r1f, *r1d);

    return r1f;
}

double oblate_radial2_nocv_wrap(double m, double n, double c, double x, double *r2d) {
    double r2f;
    xsf::oblate_radial2_nocv(m, n, c, x, r2f, *r2d);

    return r2f;
}

void oblate_radial1_wrap(double m, double n, double c, double cv, double x, double *r1f, double *r1d) {
    xsf::oblate_radial1(m, n, c, cv, x, *r1f, *r1d);
}

void oblate_radial2_wrap(double m, double n, double c, double cv, double x, double *r2f, double *r2d) {
    xsf::oblate_radial2(m, n, c, cv, x, *r2f, *r2d);
}

void modified_fresnel_plus_wrap(double x, npy_cdouble *Fplus, npy_cdouble *Kplus) {
    xsf::modified_fresnel_plus(x, *reinterpret_cast<complex<double> *>(Fplus),
                               *reinterpret_cast<complex<double> *>(Kplus));
}

void modified_fresnel_minus_wrap(double x, npy_cdouble *Fminus, npy_cdouble *Kminus) {
    xsf::modified_fresnel_minus(x, *reinterpret_cast<complex<double> *>(Fminus),
                                *reinterpret_cast<complex<double> *>(Kminus));
}

double special_sinpi(double x) { return xsf::sinpi(x); }

npy_cdouble special_csinpi(npy_cdouble z) { return to_ccomplex(xsf::sinpi(to_complex(z))); }

double special_cospi(double x) { return xsf::cospi(x); }

void special_airy(double x, double *ai, double *aip, double *bi, double *bip) { xsf::airy(x, *ai, *aip, *bi, *bip); }

void special_cairy(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip) {
    xsf::airy(to_complex(z), *reinterpret_cast<complex<double> *>(ai), *reinterpret_cast<complex<double> *>(aip),
              *reinterpret_cast<complex<double> *>(bi), *reinterpret_cast<complex<double> *>(bip));
}

void special_airye(double z, double *ai, double *aip, double *bi, double *bip) { xsf::airye(z, *ai, *aip, *bi, *bip); }

void special_cairye(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip) {
    xsf::airye(to_complex(z), *reinterpret_cast<complex<double> *>(ai), *reinterpret_cast<complex<double> *>(aip),
               *reinterpret_cast<complex<double> *>(bi), *reinterpret_cast<complex<double> *>(bip));
}

double special_cyl_bessel_j(double v, double x) { return xsf::cyl_bessel_j(v, x); }

npy_cdouble special_ccyl_bessel_j(double v, npy_cdouble z) { return to_ccomplex(xsf::cyl_bessel_j(v, to_complex(z))); }

double special_cyl_bessel_je(double v, double z) { return xsf::cyl_bessel_je(v, z); }

npy_cdouble special_ccyl_bessel_je(double v, npy_cdouble z) {
    return to_ccomplex(xsf::cyl_bessel_je(v, to_complex(z)));
}

double special_cyl_bessel_y(double v, double x) { return xsf::cyl_bessel_y(v, x); }

npy_cdouble special_ccyl_bessel_y(double v, npy_cdouble z) { return to_ccomplex(xsf::cyl_bessel_y(v, to_complex(z))); }

double special_cyl_bessel_ye(double v, double z) { return xsf::cyl_bessel_ye(v, z); }

npy_cdouble special_ccyl_bessel_ye(double v, npy_cdouble z) {
    return to_ccomplex(xsf::cyl_bessel_ye(v, to_complex(z)));
}

double special_cyl_bessel_i(double v, double z) { return xsf::cyl_bessel_i(v, z); }

npy_cdouble special_ccyl_bessel_i(double v, npy_cdouble z) { return to_ccomplex(xsf::cyl_bessel_i(v, to_complex(z))); }

double special_cyl_bessel_ie(double v, double z) { return xsf::cyl_bessel_ie(v, z); }

npy_cdouble special_ccyl_bessel_ie(double v, npy_cdouble z) {
    return to_ccomplex(xsf::cyl_bessel_ie(v, to_complex(z)));
}

double special_cyl_bessel_k_int(Py_ssize_t n, double z) {
    return xsf::cyl_bessel_k(static_cast<double>(n), z);
}

double special_cyl_bessel_k(double v, double z) { return xsf::cyl_bessel_k(v, z); }

npy_cdouble special_ccyl_bessel_k(double v, npy_cdouble z) { return to_ccomplex(xsf::cyl_bessel_k(v, to_complex(z))); }

double special_cyl_bessel_ke(double v, double z) { return xsf::cyl_bessel_ke(v, z); }

npy_cdouble special_ccyl_bessel_ke(double v, npy_cdouble z) {
    return to_ccomplex(xsf::cyl_bessel_ke(v, to_complex(z)));
}

npy_cdouble special_ccyl_hankel_1(double v, npy_cdouble z) { return to_ccomplex(xsf::cyl_hankel_1(v, to_complex(z))); }

npy_cdouble special_ccyl_hankel_1e(double v, npy_cdouble z) {
    return to_ccomplex(xsf::cyl_hankel_1e(v, to_complex(z)));
}

npy_cdouble special_ccyl_hankel_2(double v, npy_cdouble z) { return to_ccomplex(xsf::cyl_hankel_2(v, to_complex(z))); }

npy_cdouble special_ccyl_hankel_2e(double v, npy_cdouble z) {
    return to_ccomplex(xsf::cyl_hankel_2e(v, to_complex(z)));
}

double binom_wrap(double n, double k) { return xsf::binom(n, k); }

double special_binom(double n, double k) { return xsf::binom(n, k); }

double special_digamma(double z) { return xsf::digamma(z); }

npy_cdouble special_cdigamma(npy_cdouble z) { return to_ccomplex(xsf::digamma(to_complex(z))); }

double special_gamma(double x) { return xsf::gamma(x); }

npy_cdouble special_cgamma(npy_cdouble z) { return to_ccomplex(xsf::gamma(to_complex(z))); }

double special_rgamma(double x) { return xsf::rgamma(x); }

npy_cdouble special_crgamma(npy_cdouble z) { return to_ccomplex(xsf::rgamma(to_complex(z))); }

float special_expitf(float x) { return xsf::expit(x); };

double special_expit(double x) { return xsf::expit(x); };

npy_longdouble special_expitl(npy_longdouble x) { return xsf::expit(x); };

float special_log_expitf(float x) { return xsf::log_expit(x); };

double special_log_expit(double x) { return xsf::log_expit(x); };

npy_longdouble special_log_expitl(npy_longdouble x) { return xsf::log_expit(x); };

float special_logitf(float x) { return xsf::logit(x); };

double special_logit(double x) { return xsf::logit(x); };

npy_longdouble special_logitl(npy_longdouble x) { return xsf::logit(x); };

double special_loggamma(double x) { return xsf::loggamma(x); }

npy_cdouble special_cloggamma(npy_cdouble z) { return to_ccomplex(xsf::loggamma(to_complex(z))); }

double special_hyp2f1(double a, double b, double c, double z) { return xsf::hyp2f1(a, b, c, z); }

npy_cdouble special_chyp2f1(double a, double b, double c, npy_cdouble z) {
    return to_ccomplex(xsf::hyp2f1(a, b, c, to_complex(z)));
}

npy_cdouble special_lambertw(npy_cdouble z, long k, double tol) {
    return to_ccomplex(xsf::lambertw(to_complex(z), k, tol));
}

npy_cdouble special_sph_harm(long m, long n, double theta, double phi) {
    return to_ccomplex(::sph_harm(m, n, theta, phi));
}

npy_cdouble special_sph_harm_unsafe(double m, double n, double theta, double phi) {
    return to_ccomplex(::sph_harm(static_cast<long>(m), static_cast<long>(n), theta, phi));
}

double cephes_hyp2f1_wrap(double a, double b, double c, double x) { return xsf::cephes::hyp2f1(a, b, c, x); }

double cephes_airy_wrap(double x, double *ai, double *aip, double *bi, double *bip) {
    return xsf::cephes::airy(x, ai, aip, bi, bip);
}

double cephes_beta_wrap(double a, double b) { return xsf::cephes::beta(a, b); }

double cephes_lbeta_wrap(double a, double b) { return xsf::cephes::lbeta(a, b); }

double cephes_bdtr_wrap(double k, Py_ssize_t n, double p) {
    return xsf::cephes::bdtr(k, static_cast<int>(n), p);
}

double cephes_bdtri_wrap(double k, Py_ssize_t n, double y) {
    return xsf::cephes::bdtri(k, static_cast<int>(n), y);
}

double cephes_bdtrc_wrap(double k, Py_ssize_t n, double p) {
    return xsf::cephes::bdtrc(k, static_cast<int>(n), p);
}

double cephes_cosm1_wrap(double x) { return xsf::cephes::cosm1(x); }

double cephes_expm1_wrap(double x) { return xsf::cephes::expm1(x); }

double cephes_expn_wrap(Py_ssize_t n, double x) {
    return xsf::cephes::expn(static_cast<int>(n), x);
}

double cephes_log1p_wrap(double x) { return xsf::cephes::log1p(x); }

double cephes_gamma_wrap(double x) { return xsf::cephes::Gamma(x); }

double cephes_gammasgn_wrap(double x) { return xsf::cephes::gammasgn(x); }

double cephes_lgam_wrap(double x) { return xsf::cephes::lgam(x); }

double cephes_iv_wrap(double v, double x) { return xsf::cephes::iv(v, x); }

double cephes_jv_wrap(double v, double x) { return xsf::cephes::jv(v, x); }

int cephes_ellpj_wrap(double u, double m, double *sn, double *cn, double *dn, double *ph) {
    return xsf::cephes::ellpj(u, m, sn, cn, dn, ph);
}

double cephes_ellpk_wrap(double x) { return xsf::cephes::ellpk(x); }

int cephes_fresnl_wrap(double xxa, double *ssa, double *cca) { return xsf::cephes::fresnl(xxa, ssa, cca); }

double cephes_nbdtr_wrap(Py_ssize_t k, Py_ssize_t n, double p) {
    return xsf::cephes::nbdtr(static_cast<int>(k), static_cast<int>(n), p);
}

double cephes_nbdtrc_wrap(Py_ssize_t k, Py_ssize_t n, double p) {
    return xsf::cephes::nbdtrc(static_cast<int>(k), static_cast<int>(n), p);
}

double cephes_nbdtri_wrap(Py_ssize_t k, Py_ssize_t n, double p) {
    return xsf::cephes::nbdtri(static_cast<int>(k), static_cast<int>(n), p);
}

double cephes_ndtr_wrap(double x) { return xsf::cephes::ndtr(x); }

double cephes_ndtri_wrap(double x) { return xsf::cephes::ndtri(x); }

double cephes_pdtri_wrap(Py_ssize_t k, double y) {
    return xsf::cephes::pdtri(static_cast<int>(k), y);
}

double cephes_poch_wrap(double x, double m) { return xsf::cephes::poch(x, m); }

int cephes_sici_wrap(double x, double *si, double *ci) { return xsf::cephes::sici(x, si, ci); }

int cephes_shichi_wrap(double x, double *si, double *ci) { return xsf::cephes::shichi(x, si, ci); }

double cephes_smirnov_wrap(Py_ssize_t n, double x) {
    return xsf::cephes::smirnov(static_cast<int>(n), x);
}

double cephes_smirnovc_wrap(Py_ssize_t n, double x) {
    return xsf::cephes::smirnovc(static_cast<int>(n), x);
}

double cephes_smirnovi_wrap(Py_ssize_t n, double x) {
    return xsf::cephes::smirnovi(static_cast<int>(n), x);
}

double cephes_smirnovci_wrap(Py_ssize_t n, double x) {
    return xsf::cephes::smirnovci(static_cast<int>(n), x);
}

double cephes_smirnovp_wrap(Py_ssize_t n, double x) {
    return xsf::cephes::smirnovp(static_cast<int>(n), x);
}

double cephes__struve_asymp_large_z(double v, double z, Py_ssize_t is_h, double *err) {
    return xsf::cephes::detail::struve_asymp_large_z(v, z, static_cast<int>(is_h), err);
}

double cephes__struve_bessel_series(double v, double z, Py_ssize_t is_h, double *err) {
    return xsf::cephes::detail::struve_bessel_series(v, z, static_cast<int>(is_h), err);
}

double cephes__struve_power_series(double v, double z, Py_ssize_t is_h, double *err) {
    return xsf::cephes::detail::struve_power_series(v, z, static_cast<int>(is_h), err);
}

double cephes_yn_wrap(Py_ssize_t n, double x) {
    return xsf::cephes::yn(static_cast<int>(n), x);
}

double cephes_polevl_wrap(double x, const double coef[], int N) { return xsf::cephes::polevl(x, coef, N); }

double cephes_p1evl_wrap(double x, const double coef[], int N) { return xsf::cephes::p1evl(x, coef, N); }

double gammaln_wrap(double x) { return xsf::gammaln(x); }
double special_wright_bessel(double a, double b, double x) { return xsf::wright_bessel(a, b, x); }
double special_log_wright_bessel(double a, double b, double x) { return xsf::log_wright_bessel(a, b, x); }

double special_scaled_exp1(double x) { return xsf::scaled_exp1(x); }

double special_sph_bessel_j(long n, double x) { return xsf::sph_bessel_j(n, x); }

npy_cdouble special_csph_bessel_j(long n, npy_cdouble z) { return to_ccomplex(xsf::sph_bessel_j(n, to_complex(z))); }

double special_sph_bessel_j_jac(long n, double x) { return xsf::sph_bessel_j_jac(n, x); }

npy_cdouble special_csph_bessel_j_jac(long n, npy_cdouble z) {
    return to_ccomplex(xsf::sph_bessel_j_jac(n, to_complex(z)));
}

double special_sph_bessel_y(long n, double x) { return xsf::sph_bessel_y(n, x); }

npy_cdouble special_csph_bessel_y(long n, npy_cdouble z) { return to_ccomplex(xsf::sph_bessel_y(n, to_complex(z))); }

double special_sph_bessel_y_jac(long n, double x) { return xsf::sph_bessel_y_jac(n, x); }

npy_cdouble special_csph_bessel_y_jac(long n, npy_cdouble z) {
    return to_ccomplex(xsf::sph_bessel_y_jac(n, to_complex(z)));
}

double special_sph_bessel_i(long n, double x) { return xsf::sph_bessel_i(n, x); }

npy_cdouble special_csph_bessel_i(long n, npy_cdouble z) { return to_ccomplex(xsf::sph_bessel_i(n, to_complex(z))); }

double special_sph_bessel_i_jac(long n, double x) { return xsf::sph_bessel_i_jac(n, x); }

npy_cdouble special_csph_bessel_i_jac(long n, npy_cdouble z) {
    return to_ccomplex(xsf::sph_bessel_i_jac(n, to_complex(z)));
}

double special_sph_bessel_k(long n, double x) { return xsf::sph_bessel_k(n, x); }

npy_cdouble special_csph_bessel_k(long n, npy_cdouble z) { return to_ccomplex(xsf::sph_bessel_k(n, to_complex(z))); }

double special_sph_bessel_k_jac(long n, double x) { return xsf::sph_bessel_k_jac(n, x); }

npy_cdouble special_csph_bessel_k_jac(long n, npy_cdouble z) {
    return to_ccomplex(xsf::sph_bessel_k_jac(n, to_complex(z)));
}

double special_ellipk(double m) { return xsf::ellipk(m); }

double cephes_besselpoly(double a, double lambda, double nu) { return xsf::cephes::besselpoly(a, lambda, nu); }

double cephes_beta(double a, double b) { return xsf::cephes::beta(a, b); }

double cephes_chdtr(double df, double x) { return xsf::cephes::chdtr(df, x); }

double cephes_chdtrc(double df, double x) { return xsf::cephes::chdtrc(df, x); }

double cephes_chdtri(double df, double y) { return xsf::cephes::chdtri(df, y); }

double cephes_lbeta(double a, double b) { return xsf::cephes::lbeta(a, b); }

double cephes_sinpi(double x) { return xsf::cephes::sinpi(x); }

double cephes_cospi(double x) { return xsf::cephes::cospi(x); }

double cephes_cbrt(double x) { return xsf::cephes::detail::cbrt(x); }

double cephes_Gamma(double x) { return xsf::cephes::Gamma(x); }

double cephes_gammasgn(double x) { return xsf::cephes::gammasgn(x); }

double cephes_hyp2f1(double a, double b, double c, double x) { return xsf::cephes::hyp2f1(a, b, c, x); }

double cephes_i0(double x) { return xsf::cephes::i0(x); }

double cephes_i0e(double x) { return xsf::cephes::i0e(x); }

double cephes_i1(double x) { return xsf::cephes::i1(x); }

double cephes_i1e(double x) { return xsf::cephes::i1e(x); }

double cephes_iv(double v, double x) { return xsf::cephes::iv(v, x); }

double cephes_j0(double x) { return xsf::cephes::j0(x); }

double cephes_j1(double x) { return xsf::cephes::j1(x); }

double cephes_k0(double x) { return xsf::cephes::k0(x); }

double cephes_k0e(double x) { return xsf::cephes::k0e(x); }

double cephes_k1(double x) { return xsf::cephes::k1(x); }

double cephes_k1e(double x) { return xsf::cephes::k1e(x); }

double cephes_y0(double x) { return xsf::cephes::y0(x); }

double cephes_y1(double x) { return xsf::cephes::y1(x); }

double cephes_yn(int n, double x) { return xsf::cephes::yn(n, x); }

double cephes_igam(double a, double x) { return xsf::cephes::igam(a, x); }

double cephes_igamc(double a, double x) { return xsf::cephes::igamc(a, x); }

double cephes_igami(double a, double p) { return xsf::cephes::igami(a, p); }

double cephes_igamci(double a, double p) { return xsf::cephes::igamci(a, p); }

double cephes_igam_fac(double a, double x) { return xsf::cephes::detail::igam_fac(a, x); }

double cephes_lanczos_sum_expg_scaled(double x) { return xsf::cephes::lanczos_sum_expg_scaled(x); }

double cephes_kolmogorov(double x) { return xsf::cephes::kolmogorov(x); }

double cephes_kolmogc(double x) { return xsf::cephes::kolmogc(x); }

double cephes_kolmogi(double x) { return xsf::cephes::kolmogi(x); }

double cephes_kolmogci(double x) { return xsf::cephes::kolmogci(x); }

double cephes_kolmogp(double x) { return xsf::cephes::kolmogp(x); }

double cephes_smirnov(int n, double x) { return xsf::cephes::smirnov(n, x); }

double cephes_smirnovc(int n, double x) { return xsf::cephes::smirnovc(n, x); }

double cephes_smirnovi(int n, double x) { return xsf::cephes::smirnovi(n, x); }

double cephes_smirnovci(int n, double x) { return xsf::cephes::smirnovci(n, x); }

double cephes_smirnovp(int n, double x) { return xsf::cephes::smirnovp(n, x); }

double cephes_ndtr(double x) { return xsf::cephes::ndtr(x); }

double cephes_erf(double x) { return xsf::cephes::erf(x); }

double cephes_erfc(double x) { return xsf::cephes::erfc(x); }

double cephes_poch(double x, double m) { return xsf::cephes::poch(x, m); }

double cephes_rgamma(double x) { return xsf::cephes::rgamma(x); }

double cephes_zeta(double x, double q) { return xsf::cephes::zeta(x, q); }

double cephes_zetac(double x) { return xsf::cephes::zetac(x); }

double cephes_riemann_zeta(double x) { return xsf::cephes::riemann_zeta(x); }

double cephes_log1p(double x) { return xsf::cephes::log1p(x); }

double cephes_log1pmx(double x) { return xsf::cephes::log1pmx(x); }

double cephes_lgam1p(double x) { return xsf::cephes::lgam1p(x); }

double cephes_expm1(double x) { return xsf::cephes::expm1(x); }

double cephes_cosm1(double x) { return xsf::cephes::cosm1(x); }

double cephes_expn(int n, double x) { return xsf::cephes::expn(n, x); }

double cephes_ellpe(double x) { return xsf::cephes::ellpe(x); }

double cephes_ellpk(double x) { return xsf::cephes::ellpk(x); }

double cephes_ellie(double phi, double m) { return xsf::cephes::ellie(phi, m); }

double cephes_ellik(double phi, double m) { return xsf::cephes::ellik(phi, m); }

double cephes_sindg(double x) { return xsf::cephes::sindg(x); }

double cephes_cosdg(double x) { return xsf::cephes::cosdg(x); }

double cephes_tandg(double x) { return xsf::cephes::tandg(x); }

double cephes_cotdg(double x) { return xsf::cephes::cotdg(x); }

double cephes_radian(double d, double m, double s) { return xsf::cephes::radian(d, m, s); }

double cephes_ndtri(double x) { return xsf::cephes::ndtri(x); }

double cephes_bdtr(double k, int n, double p) { return xsf::cephes::bdtr(k, n, p); }

double cephes_bdtri(double k, int n, double y) { return xsf::cephes::bdtri(k, n, y); }

double cephes_bdtrc(double k, int n, double p) { return xsf::cephes::bdtrc(k, n, p); }

double cephes_btdtri(double aa, double bb, double yy0) { return xsf::cephes::incbi(aa, bb, yy0); }

double cephes_btdtr(double a, double b, double x) { return xsf::cephes::incbet(a, b, x); }

double cephes_erfcinv(double y) { return xsf::cephes::erfcinv(y); }

double cephes_exp10(double x) { return xsf::cephes::exp10(x); }

double cephes_exp2(double x) { return xsf::cephes::exp2(x); }

double cephes_fdtr(double a, double b, double x) { return xsf::cephes::fdtr(a, b, x); }

double cephes_fdtrc(double a, double b, double x) { return xsf::cephes::fdtrc(a, b, x); }

double cephes_fdtri(double a, double b, double y) { return xsf::cephes::fdtri(a, b, y); }

double cephes_gdtr(double a, double b, double x) { return xsf::cephes::gdtr(a, b, x); }

double cephes_gdtrc(double a, double b, double x) { return xsf::cephes::gdtrc(a, b, x); }

double cephes_owens_t(double h, double a) { return xsf::cephes::owens_t(h, a); }

double cephes_nbdtr(int k, int n, double p) { return xsf::cephes::nbdtr(k, n, p); }

double cephes_nbdtrc(int k, int n, double p) { return xsf::cephes::nbdtrc(k, n, p); }

double cephes_nbdtri(int k, int n, double p) { return xsf::cephes::nbdtri(k, n, p); }

double cephes_pdtr(double k, double m) { return xsf::cephes::pdtr(k, m); }

double cephes_pdtrc(double k, double m) { return xsf::cephes::pdtrc(k, m); }

double cephes_pdtri(int k, double y) { return xsf::cephes::pdtri(k, y); }

double cephes_round(double x) { return xsf::cephes::round(x); }

double cephes_spence(double x) { return xsf::cephes::spence(x); }

double cephes_tukeylambdacdf(double x, double lmbda) { return xsf::cephes::tukeylambdacdf(x, lmbda); }

double cephes_struve_h(double v, double z) { return xsf::cephes::struve_h(v, z); }

double cephes_struve_l(double v, double z) { return xsf::cephes::struve_l(v, z); }
