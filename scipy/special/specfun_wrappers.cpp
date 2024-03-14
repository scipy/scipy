#include "specfun_wrappers.h"
#include "special/specfun.h"

extern "C" {

npy_cdouble chyp2f1_wrap(double a, double b, double c, npy_cdouble z) {
  std::complex<double> res = special::chyp2f1_wrap(
      a, b, c, *reinterpret_cast<std::complex<double> *>(&z));
  return *reinterpret_cast<npy_cdouble *>(&res);
}

npy_cdouble chyp1f1_wrap(double a, double b, npy_cdouble z) {
  std::complex<double> res = special::chyp1f1_wrap(
      a, b, *reinterpret_cast<std::complex<double> *>(&z));
  return *reinterpret_cast<npy_cdouble *>(&res);
}

double hypU_wrap(double a, double b, double x) {
  return special::hypU_wrap(a, b, x);
}

double hyp1f1_wrap(double a, double b, double x) {
  return special::hyp1f1_wrap(a, b, x);
}

int itairy_wrap(double x, double *apt, double *bpt, double *ant, double *bnt) {
  return special::itairy_wrap(x, apt, bpt, ant, bnt);
}

double exp1_wrap(double x) {
  return special::exp1_wrap(x);
}

npy_cdouble cexp1_wrap(npy_cdouble z) {
  std::complex<double> res = special::cexp1_wrap(*reinterpret_cast<std::complex<double> *>(&z));
  return *reinterpret_cast<npy_cdouble *>(&res);
}

double expi_wrap(double x) {
  return special::expi_wrap(x);
}

npy_cdouble cexpi_wrap(npy_cdouble z) {
  std::complex<double> res = special::cexpi_wrap(*reinterpret_cast<std::complex<double> *>(&z));
  return *reinterpret_cast<npy_cdouble *>(&res);
}

npy_cdouble cerf_wrap(npy_cdouble z) {
  std::complex<double> res = special::cerf_wrap(*reinterpret_cast<std::complex<double> *>(&z));
  return *reinterpret_cast<npy_cdouble *>(&res);
}

double itstruve0_wrap(double x) {
  return special::itstruve0_wrap(x);
}

double it2struve0_wrap(double x) {
  return special::it2struve0_wrap(x);
}

double itmodstruve0_wrap(double x) {
  return special::itmodstruve0_wrap(x);
}

double ber_wrap(double x) {
  return special::ber_wrap(x);
}

double bei_wrap(double x) {
  return special::bei_wrap(x);
}

double ker_wrap(double x) {
  return special::ker_wrap(x);
}

double kei_wrap(double x) {
  return special::kei_wrap(x);
}

double berp_wrap(double x) {
  return special::berp_wrap(x);
}

double beip_wrap(double x) {
  return special::beip_wrap(x);
}

double kerp_wrap(double x) {
  return special::kerp_wrap(x);
}

double keip_wrap(double x) {
  return special::keip_wrap(x);
}

int kelvin_wrap(double x, npy_cdouble *Be, npy_cdouble *Ke, npy_cdouble *Bep,
                npy_cdouble *Kep) {
  return special::kelvin_wrap(x, reinterpret_cast<std::complex<double> *>(Be), reinterpret_cast<std::complex<double> *>(Ke),
    reinterpret_cast<std::complex<double> *>(Bep), reinterpret_cast<std::complex<double> *>(Kep));
}

int it1j0y0_wrap(double x, double *j0int, double *y0int) {
  return special::it1j0y0_wrap(x, j0int, y0int);
}

int it2j0y0_wrap(double x, double *j0int, double *y0int) {
  return special::it2j0y0_wrap(x, j0int, y0int);
}

int it1i0k0_wrap(double x, double *i0int, double *k0int) {
  return special::it1i0k0_wrap(x, i0int, k0int);
}

int it2i0k0_wrap(double x, double *i0int, double *k0int) {
  return special::it2i0k0_wrap(x, i0int, k0int);
}

int cfresnl_wrap(npy_cdouble z, npy_cdouble *zfs, npy_cdouble *zfc) {
  return special::cfresnl_wrap(*reinterpret_cast<std::complex<double> *>(&z),
    reinterpret_cast<std::complex<double> *>(zfs), reinterpret_cast<std::complex<double> *>(zfc));
}

double cem_cva_wrap(double m, double q) {
  return special::cem_cva_wrap(m, q);
}

double sem_cva_wrap(double m, double q) {
  return special::sem_cva_wrap(m, q);
}

int cem_wrap(double m, double q, double x, double *csf, double *csd) {
  return special::cem_wrap(m, q, x, csf, csd);
}

int sem_wrap(double m, double q, double x, double *csf, double *csd) {
  return special::sem_wrap(m, q, x, csf, csd);
}

int mcm1_wrap(double m, double q, double x, double *f1r, double *d1r) {
  return special::mcm1_wrap(m, q, x, f1r, d1r);
}

int msm1_wrap(double m, double q, double x, double *f1r, double *d1r) {
  return special::msm1_wrap(m, q, x, f1r, d1r);
}

int mcm2_wrap(double m, double q, double x, double *f2r, double *d2r) {
  return special::mcm2_wrap(m, q, x, f2r, d2r);
}

int msm2_wrap(double m, double q, double x, double *f2r, double *d2r) {
  return special::msm2_wrap(m, q, x, f2r, d2r);
}

double pmv_wrap(double m, double v, double x) {
  return special::pmv_wrap(m, v, x);
}

int pbwa_wrap(double a, double x, double *wf, double *wd) {
  return special::pbwa_wrap(a, x, wf, wd);
}

int pbdv_wrap(double v, double x, double *pdf, double *pdd) {
  return special::pbdv_wrap(v, x, pdf, pdd);
}

int pbvv_wrap(double v, double x, double *pvf, double *pvd) {
  return special::pbvv_wrap(v, x, pvf, pvd);
}

double prolate_segv_wrap(double m, double n, double c) {
  return special::prolate_segv_wrap(m, n, c);
}

double oblate_segv_wrap(double m, double n, double c) {
  return special::oblate_segv_wrap(m, n, c);
}

double prolate_aswfa_nocv_wrap(double m, double n, double c, double x,
                               double *s1d) {
  return special::prolate_aswfa_nocv_wrap(m, n, c, x, s1d);
}

double oblate_aswfa_nocv_wrap(double m, double n, double c, double x,
                              double *s1d) {
  return special::oblate_aswfa_nocv_wrap(m, n, c, x, s1d);
}

int prolate_aswfa_wrap(double m, double n, double c, double cv, double x,
                       double *s1f, double *s1d) {
  return special::prolate_aswfa_wrap(m, n, c, cv, x, s1f, s1d);
}

int oblate_aswfa_wrap(double m, double n, double c, double cv, double x,
                      double *s1f, double *s1d) {
  return special::oblate_aswfa_wrap(m, n, c, cv, x, s1f, s1d);
}

double prolate_radial1_nocv_wrap(double m, double n, double c, double x,
                                 double *r1d) {
  return special::prolate_radial1_nocv_wrap(m, n, c, x, r1d);
}

double prolate_radial2_nocv_wrap(double m, double n, double c, double x,
                                 double *r2d) {
  return special::prolate_radial2_nocv_wrap(m, n, c, x, r2d);
}

int prolate_radial1_wrap(double m, double n, double c, double cv, double x,
                         double *r1f, double *r1d) {
  return special::prolate_radial1_wrap(m, n, c, cv, x, r1f, r1d);
}

int prolate_radial2_wrap(double m, double n, double c, double cv, double x,
                         double *r2f, double *r2d) {
  return special::prolate_radial2_wrap(m, n, c, cv, x, r2f, r2d);
}

double oblate_radial1_nocv_wrap(double m, double n, double c, double x,
                                double *r1d) {
  return special::oblate_radial1_nocv_wrap(m, n, c, x, r1d);
}

double oblate_radial2_nocv_wrap(double m, double n, double c, double x,
                                double *r2d) {
  return special::oblate_radial2_nocv_wrap(m, n, c, x, r2d);
}

int oblate_radial1_wrap(double m, double n, double c, double cv, double x,
                        double *r1f, double *r1d) {
  return special::oblate_radial1_wrap(m, n, c, cv, x, r1f, r1d);
}

int oblate_radial2_wrap(double m, double n, double c, double cv, double x,
                        double *r2f, double *r2d) {
  return special::oblate_radial2_wrap(m, n, c, cv, x, r2f, r2d);
}

int modified_fresnel_plus_wrap(double x, npy_cdouble *Fplus,
                               npy_cdouble *Kplus) {
  return special::modified_fresnel_minus_wrap(
      x, reinterpret_cast<std::complex<double> *>(Fplus),
      reinterpret_cast<std::complex<double> *>(Kplus));
}

int modified_fresnel_minus_wrap(double x, npy_cdouble *Fminus,
                                npy_cdouble *Kminus) {
  return special::modified_fresnel_minus_wrap(
      x, reinterpret_cast<std::complex<double> *>(Fminus),
      reinterpret_cast<std::complex<double> *>(Kminus));
}
}