#include "special_wrappers.h"
#include "special.h"
#include "special/airy.h"
#include "special/amos.h"
#include "special/bessel.h"
#include "special/binom.h"
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

#include "special/binom.h"
#include "special/digamma.h"
#include "special/ellipk.h"
#include "special/gamma.h"
#include "special/hyp2f1.h"
#include "special/lambertw.h"
#include "special/loggamma.h"
#include "special/trig.h"
#include "special/wright_bessel.h"

#include "special/cephes/bdtr.h"
#include "special/cephes/besselpoly.h"
#include "special/cephes/beta.h"
#include "special/cephes/cbrt.h"
#include "special/cephes/chdtr.h"
#include "special/cephes/ellie.h"
#include "special/cephes/ellik.h"
#include "special/cephes/ellpe.h"
#include "special/cephes/ellpj.h"
#include "special/cephes/ellpk.h"
#include "special/cephes/erfinv.h"
#include "special/cephes/exp10.h"
#include "special/cephes/exp2.h"
#include "special/cephes/expn.h"
#include "special/cephes/fdtr.h"
#include "special/cephes/gamma.h"
#include "special/cephes/gdtr.h"
#include "special/cephes/hyp2f1.h"
#include "special/cephes/hyperg.h"
#include "special/cephes/i0.h"
#include "special/cephes/i1.h"
#include "special/cephes/igam.h"
#include "special/cephes/igami.h"
#include "special/cephes/incbet.h"
#include "special/cephes/incbi.h"
#include "special/cephes/j0.h"
#include "special/cephes/j1.h"
#include "special/cephes/jv.h"
#include "special/cephes/k0.h"
#include "special/cephes/k1.h"
#include "special/cephes/kolmogorov.h"
#include "special/cephes/lanczos.h"
#include "special/cephes/nbdtr.h"
#include "special/cephes/ndtr.h"
#include "special/cephes/ndtri.h"
#include "special/cephes/owens_t.h"
#include "special/cephes/pdtr.h"
#include "special/cephes/poch.h"
#include "special/cephes/rgamma.h"
#include "special/cephes/round.h"
#include "special/cephes/scipy_iv.h"
#include "special/cephes/sindg.h"
#include "special/cephes/spence.h"
#include "special/cephes/struve.h"
#include "special/cephes/tandg.h"
#include "special/cephes/trig.h"
#include "special/cephes/tukey.h"
#include "special/cephes/unity.h"
#include "special/cephes/yn.h"
#include "special/cephes/zeta.h"
#include "special/cephes/zetac.h"

#include "special/cephes/airy.h"
#include "special/cephes/bdtr.h"
#include "special/cephes/beta.h"
#include "special/cephes/ellpj.h"
#include "special/cephes/ellpk.h"
#include "special/cephes/expn.h"
#include "special/cephes/fresnl.h"
#include "special/cephes/gamma.h"
#include "special/cephes/hyp2f1.h"
#include "special/cephes/jv.h"
#include "special/cephes/kolmogorov.h"
#include "special/cephes/nbdtr.h"
#include "special/cephes/ndtr.h"
#include "special/cephes/ndtri.h"
#include "special/cephes/pdtr.h"
#include "special/cephes/shichi.h"
#include "special/cephes/sici.h"

using namespace std;

namespace {

complex<double> to_complex(npy_cdouble z) {
    return {npy_creal(z), npy_cimag(z)};
}

npy_cdouble to_ccomplex(complex<double> z) {
    return {z.real(), z.imag()};
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

double binom_wrap(double n, double k) { return special::binom(n, k); }

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

double cephes_hyp2f1_wrap(double a, double b, double c, double x) { return special::cephes::hyp2f1(a, b, c, x); }

double cephes_airy_wrap(double x, double *ai, double *aip, double *bi, double *bip) {
    return special::cephes::airy(x, ai, aip, bi, bip);
}

double cephes_beta_wrap(double a, double b) { return special::cephes::beta(a, b); }

double cephes_lbeta_wrap(double a, double b) { return special::cephes::lbeta(a, b); }

double cephes_bdtr_wrap(double k, int n, double p) { return special::cephes::bdtr(k, n, p); }

double cephes_bdtri_wrap(double k, int n, double y) { return special::cephes::bdtri(k, n, y); }

double cephes_bdtrc_wrap(double k, int n, double p) { return special::cephes::bdtrc(k, n, p); }

double cephes_cosm1_wrap(double x) { return special::cephes::cosm1(x); }

double cephes_expm1_wrap(double x) { return special::cephes::expm1(x); }

double cephes_expn_wrap(int n, double x) { return special::cephes::expn(n, x); }

double cephes_log1p_wrap(double x) { return special::cephes::log1p(x); }

double cephes_gamma_wrap(double x) { return special::cephes::Gamma(x); }

double cephes_gammasgn_wrap(double x) { return special::cephes::gammasgn(x); }

double cephes_lgam_wrap(double x) { return special::cephes::lgam(x); }

double cephes_iv_wrap(double v, double x) { return special::cephes::iv(v, x); }

double cephes_jv_wrap(double v, double x) { return special::cephes::jv(v, x); }

int cephes_ellpj_wrap(double u, double m, double *sn, double *cn, double *dn, double *ph) {
    return special::cephes::ellpj(u, m, sn, cn, dn, ph);
}

double cephes_ellpk_wrap(double x) { return special::cephes::ellpk(x); }

int cephes_fresnl_wrap(double xxa, double *ssa, double *cca) { return special::cephes::fresnl(xxa, ssa, cca); }

double cephes_nbdtr_wrap(int k, int n, double p) { return special::cephes::nbdtr(k, n, p); }

double cephes_nbdtrc_wrap(int k, int n, double p) { return special::cephes::nbdtrc(k, n, p); }

double cephes_nbdtri_wrap(int k, int n, double p) { return special::cephes::nbdtri(k, n, p); }

double cephes_ndtr_wrap(double x) { return special::cephes::ndtr(x); }

double cephes_ndtri_wrap(double x) { return special::cephes::ndtri(x); }

double cephes_pdtri_wrap(int k, double y) { return special::cephes::pdtri(k, y); }

double cephes_poch_wrap(double x, double m) { return special::cephes::poch(x, m); }

int cephes_sici_wrap(double x, double *si, double *ci) { return special::cephes::sici(x, si, ci); }

int cephes_shichi_wrap(double x, double *si, double *ci) { return special::cephes::shichi(x, si, ci); }

double cephes_smirnov_wrap(int n, double x) { return special::cephes::smirnov(n, x); }

double cephes_smirnovc_wrap(int n, double x) { return special::cephes::smirnovc(n, x); }

double cephes_smirnovi_wrap(int n, double x) { return special::cephes::smirnovi(n, x); }

double cephes_smirnovci_wrap(int n, double x) { return special::cephes::smirnovci(n, x); }

double cephes_smirnovp_wrap(int n, double x) { return special::cephes::smirnovp(n, x); }

double cephes__struve_asymp_large_z(double v, double z, int is_h, double *err) {
    return special::cephes::detail::struve_asymp_large_z(v, z, is_h, err);
}

double cephes__struve_bessel_series(double v, double z, int is_h, double *err) {
    return special::cephes::detail::struve_bessel_series(v, z, is_h, err);
}

double cephes__struve_power_series(double v, double z, int is_h, double *err) {
    return special::cephes::detail::struve_power_series(v, z, is_h, err);
}

double cephes_yn_wrap(int n, double x) { return special::cephes::yn(n, x); }

double cephes_polevl_wrap(double x, const double coef[], int N) { return special::cephes::polevl(x, coef, N); }

double cephes_p1evl_wrap(double x, const double coef[], int N) { return special::cephes::p1evl(x, coef, N); }

double gammaln_wrap(double x) { return special::gammaln(x); }
double special_wright_bessel(double a, double b, double x) { return special::wright_bessel(a, b, x); }
double special_log_wright_bessel(double a, double b, double x) { return special::log_wright_bessel(a, b, x); }

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

double special_ellipk(double m) { return special::ellipk(m); }

double cephes_besselpoly(double a, double lambda, double nu) { return special::cephes::besselpoly(a, lambda, nu); }

double cephes_beta(double a, double b) { return special::cephes::beta(a, b); }

double cephes_chdtr(double df, double x) { return special::cephes::chdtr(df, x); }

double cephes_chdtrc(double df, double x) { return special::cephes::chdtrc(df, x); }

double cephes_chdtri(double df, double y) { return special::cephes::chdtri(df, y); }

double cephes_lbeta(double a, double b) { return special::cephes::lbeta(a, b); }

double cephes_sinpi(double x) { return special::cephes::sinpi(x); }

double cephes_cospi(double x) { return special::cephes::cospi(x); }

double cephes_cbrt(double x) { return special::cephes::detail::cbrt(x); }

double cephes_Gamma(double x) { return special::cephes::Gamma(x); }

double cephes_gammasgn(double x) { return special::cephes::gammasgn(x); }

double cephes_hyp2f1(double a, double b, double c, double x) { return special::cephes::hyp2f1(a, b, c, x); }

double cephes_i0(double x) { return special::cephes::i0(x); }

double cephes_i0e(double x) { return special::cephes::i0e(x); }

double cephes_i1(double x) { return special::cephes::i1(x); }

double cephes_i1e(double x) { return special::cephes::i1e(x); }

double cephes_iv(double v, double x) { return special::cephes::iv(v, x); }

double cephes_j0(double x) { return special::cephes::j0(x); }

double cephes_j1(double x) { return special::cephes::j1(x); }

double cephes_k0(double x) { return special::cephes::k0(x); }

double cephes_k0e(double x) { return special::cephes::k0e(x); }

double cephes_k1(double x) { return special::cephes::k1(x); }

double cephes_k1e(double x) { return special::cephes::k1e(x); }

double cephes_y0(double x) { return special::cephes::y0(x); }

double cephes_y1(double x) { return special::cephes::y1(x); }

double cephes_yn(int n, double x) { return special::cephes::yn(n, x); }

double cephes_igam(double a, double x) { return special::cephes::igam(a, x); }

double cephes_igamc(double a, double x) { return special::cephes::igamc(a, x); }

double cephes_igami(double a, double p) { return special::cephes::igami(a, p); }

double cephes_igamci(double a, double p) { return special::cephes::igamci(a, p); }

double cephes_igam_fac(double a, double x) { return special::cephes::detail::igam_fac(a, x); }

double cephes_lanczos_sum_expg_scaled(double x) { return special::cephes::lanczos_sum_expg_scaled(x); }

double cephes_kolmogorov(double x) { return special::cephes::kolmogorov(x); }

double cephes_kolmogc(double x) { return special::cephes::kolmogc(x); }

double cephes_kolmogi(double x) { return special::cephes::kolmogi(x); }

double cephes_kolmogci(double x) { return special::cephes::kolmogci(x); }

double cephes_kolmogp(double x) { return special::cephes::kolmogp(x); }

double cephes_smirnov(int n, double x) { return special::cephes::smirnov(n, x); }

double cephes_smirnovc(int n, double x) { return special::cephes::smirnovc(n, x); }

double cephes_smirnovi(int n, double x) { return special::cephes::smirnovi(n, x); }

double cephes_smirnovci(int n, double x) { return special::cephes::smirnovci(n, x); }

double cephes_smirnovp(int n, double x) { return special::cephes::smirnovp(n, x); }

double cephes_ndtr(double x) { return special::cephes::ndtr(x); }

double cephes_erf(double x) { return special::cephes::erf(x); }

double cephes_erfc(double x) { return special::cephes::erfc(x); }

double cephes_poch(double x, double m) { return special::cephes::poch(x, m); }

double cephes_rgamma(double x) { return special::cephes::rgamma(x); }

double cephes_zeta(double x, double q) { return special::cephes::zeta(x, q); }

double cephes_zetac(double x) { return special::cephes::zetac(x); }

double cephes_riemann_zeta(double x) { return special::cephes::riemann_zeta(x); }

double cephes_log1p(double x) { return special::cephes::log1p(x); }

double cephes_log1pmx(double x) { return special::cephes::log1pmx(x); }

double cephes_lgam1p(double x) { return special::cephes::lgam1p(x); }

double cephes_expm1(double x) { return special::cephes::expm1(x); }

double cephes_cosm1(double x) { return special::cephes::cosm1(x); }

double cephes_expn(int n, double x) { return special::cephes::expn(n, x); }

double cephes_ellpe(double x) { return special::cephes::ellpe(x); }

double cephes_ellpk(double x) { return special::cephes::ellpk(x); }

double cephes_ellie(double phi, double m) { return special::cephes::ellie(phi, m); }

double cephes_ellik(double phi, double m) { return special::cephes::ellik(phi, m); }

double cephes_sindg(double x) { return special::cephes::sindg(x); }

double cephes_cosdg(double x) { return special::cephes::cosdg(x); }

double cephes_tandg(double x) { return special::cephes::tandg(x); }

double cephes_cotdg(double x) { return special::cephes::cotdg(x); }

double cephes_radian(double d, double m, double s) { return special::cephes::radian(d, m, s); }

double cephes_ndtri(double x) { return special::cephes::ndtri(x); }

double cephes_bdtr(double k, int n, double p) { return special::cephes::bdtr(k, n, p); }

double cephes_bdtri(double k, int n, double y) { return special::cephes::bdtri(k, n, y); }

double cephes_bdtrc(double k, int n, double p) { return special::cephes::bdtrc(k, n, p); }

double cephes_btdtri(double aa, double bb, double yy0) { return special::cephes::incbi(aa, bb, yy0); }

double cephes_btdtr(double a, double b, double x) { return special::cephes::incbet(a, b, x); }

double cephes_erfcinv(double y) { return special::cephes::erfcinv(y); }

double cephes_exp10(double x) { return special::cephes::exp10(x); }

double cephes_exp2(double x) { return special::cephes::exp2(x); }

double cephes_fdtr(double a, double b, double x) { return special::cephes::fdtr(a, b, x); }

double cephes_fdtrc(double a, double b, double x) { return special::cephes::fdtrc(a, b, x); }

double cephes_fdtri(double a, double b, double y) { return special::cephes::fdtri(a, b, y); }

double cephes_gdtr(double a, double b, double x) { return special::cephes::gdtr(a, b, x); }

double cephes_gdtrc(double a, double b, double x) { return special::cephes::gdtrc(a, b, x); }

double cephes_owens_t(double h, double a) { return special::cephes::owens_t(h, a); }

double cephes_nbdtr(int k, int n, double p) { return special::cephes::nbdtr(k, n, p); }

double cephes_nbdtrc(int k, int n, double p) { return special::cephes::nbdtrc(k, n, p); }

double cephes_nbdtri(int k, int n, double p) { return special::cephes::nbdtri(k, n, p); }

double cephes_pdtr(double k, double m) { return special::cephes::pdtr(k, m); }

double cephes_pdtrc(double k, double m) { return special::cephes::pdtrc(k, m); }

double cephes_pdtri(int k, double y) { return special::cephes::pdtri(k, y); }

double cephes_round(double x) { return special::cephes::round(x); }

double cephes_spence(double x) { return special::cephes::spence(x); }

double cephes_tukeylambdacdf(double x, double lmbda) { return special::cephes::tukeylambdacdf(x, lmbda); }

double cephes_struve_h(double v, double z) { return special::cephes::struve_h(v, z); }

double cephes_struve_l(double v, double z) { return special::cephes::struve_l(v, z); }
