#include "specfun_wrappers.h"
#include "special/airy.h"
#include "special/bessel.h"
#include "special/fresnel.h"
#include "special/kelvin.h"
#include "special/mathieu.h"
#include "special/par_cyl.h"
#include "special/specfun.h"
#include "special/sphd_wave.h"
#include "special/struve.h"

extern "C" {

npy_cdouble chyp2f1_wrap(double a, double b, double c, npy_cdouble z) {
    std::complex<double> res = special::chyp2f1(a, b, c, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

npy_cdouble chyp1f1_wrap(double a, double b, npy_cdouble z) {
    std::complex<double> res = special::chyp1f1(a, b, {npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

double hypU_wrap(double a, double b, double x) { return special::hypu(a, b, x); }

double hyp1f1_wrap(double a, double b, double x) { return special::hyp1f1(a, b, x); }

int itairy_wrap(double x, double *apt, double *bpt, double *ant, double *bnt) {
    special::itairy(x, apt, bpt, ant, bnt);
    return 0;
}

double exp1_wrap(double x) { return special::exp1(x); }

npy_cdouble cexp1_wrap(npy_cdouble z) {
    std::complex<double> res = special::exp1(std::complex<double>{npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

double expi_wrap(double x) { return special::expi(x); }

npy_cdouble cexpi_wrap(npy_cdouble z) {
    std::complex<double> res = special::expi(std::complex<double>{npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

npy_cdouble cerf_wrap(npy_cdouble z) {
    std::complex<double> res = special::cerf({npy_creal(z), npy_cimag(z)});
    return npy_cpack(res.real(), res.imag());
}

double itstruve0_wrap(double x) { return special::itstruve0(x); }

double it2struve0_wrap(double x) { return special::it2struve0(x); }

double itmodstruve0_wrap(double x) { return special::itmodstruve0(x); }

double ber_wrap(double x) { return special::ber(x); }

double bei_wrap(double x) { return special::bei(x); }

double ker_wrap(double x) { return special::ker(x); }

double kei_wrap(double x) { return special::kei(x); }

double berp_wrap(double x) { return special::berp(x); }

double beip_wrap(double x) { return special::beip(x); }

double kerp_wrap(double x) { return special::kerp(x); }

double keip_wrap(double x) { return special::keip(x); }

int kelvin_wrap(double x, npy_cdouble *Be, npy_cdouble *Ke, npy_cdouble *Bep, npy_cdouble *Kep) {
    special::kelvin(x, reinterpret_cast<std::complex<double> *>(Be), reinterpret_cast<std::complex<double> *>(Ke),
                    reinterpret_cast<std::complex<double> *>(Bep), reinterpret_cast<std::complex<double> *>(Kep));
    return 0;
}

int it1j0y0_wrap(double x, double *j0int, double *y0int) {
    special::it1j0y0(x, j0int, y0int);
    return 0;
}

int it2j0y0_wrap(double x, double *j0int, double *y0int) {
    special::it2j0y0(x, j0int, y0int);
    return 0;
}

int it1i0k0_wrap(double x, double *i0int, double *k0int) {
    special::it1i0k0(x, i0int, k0int);
    return 0;
}

int it2i0k0_wrap(double x, double *i0int, double *k0int) {
    special::it2i0k0(x, i0int, k0int);
    return 0;
}

int cfresnl_wrap(npy_cdouble z, npy_cdouble *zfs, npy_cdouble *zfc) {
    special::cfresnl({npy_creal(z), npy_cimag(z)}, reinterpret_cast<std::complex<double> *>(zfs),
                     reinterpret_cast<std::complex<double> *>(zfc));
    return 0;
}

double cem_cva_wrap(double m, double q) { return special::cem_cva(m, q); }

double sem_cva_wrap(double m, double q) { return special::sem_cva(m, q); }

int cem_wrap(double m, double q, double x, double *csf, double *csd) {
    special::cem(m, q, x, csf, csd);
    return 0;
}

int sem_wrap(double m, double q, double x, double *csf, double *csd) {
    special::sem(m, q, x, csf, csd);
    return 0;
}

int mcm1_wrap(double m, double q, double x, double *f1r, double *d1r) {
    special::mcm1(m, q, x, f1r, d1r);
    return 0;
}

int msm1_wrap(double m, double q, double x, double *f1r, double *d1r) {
    special::msm1(m, q, x, f1r, d1r);
    return 0;
}

int mcm2_wrap(double m, double q, double x, double *f2r, double *d2r) {
    special::mcm2(m, q, x, f2r, d2r);
    return 0;
}

int msm2_wrap(double m, double q, double x, double *f2r, double *d2r) {
    special::msm2(m, q, x, f2r, d2r);
    return 0;
}

double pmv_wrap(double m, double v, double x) { return special::pmv(m, v, x); }

int pbwa_wrap(double a, double x, double *wf, double *wd) {
    special::pbwa(a, x, wf, wd);
    return 0;
}

int pbdv_wrap(double v, double x, double *pdf, double *pdd) {
    special::pbdv(v, x, pdf, pdd);
    return 0;
}

int pbvv_wrap(double v, double x, double *pvf, double *pvd) {
    special::pbvv(v, x, pvf, pvd);
    return 0;
}

double prolate_segv_wrap(double m, double n, double c) { return special::prolate_segv(m, n, c); }

double oblate_segv_wrap(double m, double n, double c) { return special::oblate_segv(m, n, c); }

double prolate_aswfa_nocv_wrap(double m, double n, double c, double x, double *s1d) {
    double s1f;
    special::prolate_aswfa_nocv(m, n, c, x, &s1f, s1d);

    return s1f;
}

double oblate_aswfa_nocv_wrap(double m, double n, double c, double x, double *s1d) {
    double s1f;
    special::oblate_aswfa_nocv(m, n, c, x, &s1f, s1d);

    return s1f;
}

int prolate_aswfa_wrap(double m, double n, double c, double cv, double x, double *s1f, double *s1d) {
    special::prolate_aswfa(m, n, c, cv, x, s1f, s1d);
    return 0;
}

int oblate_aswfa_wrap(double m, double n, double c, double cv, double x, double *s1f, double *s1d) {
    special::oblate_aswfa(m, n, c, cv, x, s1f, s1d);
    return 0;
}

double prolate_radial1_nocv_wrap(double m, double n, double c, double x, double *r1d) {
    double r1f;
    special::prolate_radial1_nocv(m, n, c, x, &r1f, r1d);

    return r1f;
}

double prolate_radial2_nocv_wrap(double m, double n, double c, double x, double *r2d) {
    double r2f;
    special::prolate_radial2_nocv(m, n, c, x, &r2f, r2d);

    return r2f;
}

int prolate_radial1_wrap(double m, double n, double c, double cv, double x, double *r1f, double *r1d) {
    special::prolate_radial1(m, n, c, cv, x, r1f, r1d);
    return 0;
}

int prolate_radial2_wrap(double m, double n, double c, double cv, double x, double *r2f, double *r2d) {
    special::prolate_radial2(m, n, c, cv, x, r2f, r2d);
    return 0;
}

double oblate_radial1_nocv_wrap(double m, double n, double c, double x, double *r1d) {
    double r1f;
    special::oblate_radial1_nocv(m, n, c, x, &r1f, r1d);

    return r1f;
}

double oblate_radial2_nocv_wrap(double m, double n, double c, double x, double *r2d) {
    double r2f;
    special::oblate_radial2_nocv(m, n, c, x, &r2f, r2d);

    return r2f;
}

int oblate_radial1_wrap(double m, double n, double c, double cv, double x, double *r1f, double *r1d) {
    special::oblate_radial1(m, n, c, cv, x, r1f, r1d);
    return 0;
}

int oblate_radial2_wrap(double m, double n, double c, double cv, double x, double *r2f, double *r2d) {
    special::oblate_radial2(m, n, c, cv, x, r2f, r2d);
    return 0;
}

int modified_fresnel_plus_wrap(double x, npy_cdouble *Fplus, npy_cdouble *Kplus) {
    special::modified_fresnel_plus(x, reinterpret_cast<std::complex<double> *>(Fplus),
                                   reinterpret_cast<std::complex<double> *>(Kplus));
    return 0;
}

int modified_fresnel_minus_wrap(double x, npy_cdouble *Fminus, npy_cdouble *Kminus) {
    special::modified_fresnel_minus(x, reinterpret_cast<std::complex<double> *>(Fminus),
                                    reinterpret_cast<std::complex<double> *>(Kminus));
    return 0;
}
}
