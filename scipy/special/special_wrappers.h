/* This file is a collection of wrappers around the
 *  Special Function  Fortran library of functions
 *  to be compiled with the other special functions in cephes
 *
 * Functions written by Shanjie Zhang and Jianming Jin.
 * Interface by
 *  Travis E. Oliphant
 */

#pragma once

#include "Python.h"
#include "npy_2_complexcompat.h"
#include "sf_error.h"
#include <math.h>
#include <numpy/npy_math.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

npy_cdouble clngamma_wrap(npy_cdouble z);
npy_cdouble chyp2f1_wrap(double a, double b, double c, npy_cdouble z);
npy_cdouble chyp1f1_wrap(double a, double b, npy_cdouble z);
double hyp1f1_wrap(double a, double b, double x);
double hypU_wrap(double a, double b, double x);
npy_cdouble cerf_wrap(npy_cdouble z);

double special_exp1(double x);
npy_cdouble special_cexp1(npy_cdouble z);

double special_expi(double x);
npy_cdouble special_cexpi(npy_cdouble z);

double struve_wrap(double v, double x);
double special_itstruve0(double x);
double special_it2struve0(double x);

double modstruve_wrap(double v, double x);
double special_itmodstruve0(double x);

double special_ber(double x);
double special_bei(double x);
double special_ker(double x);
double special_kei(double x);
double special_berp(double x);
double special_beip(double x);
double special_kerp(double x);
double special_keip(double x);

void special_ckelvin(double x, npy_cdouble *Be, npy_cdouble *Ke, npy_cdouble *Bep, npy_cdouble *Kep);

void it1j0y0_wrap(double x, double *, double *);
void it2j0y0_wrap(double x, double *, double *);
void it1i0k0_wrap(double x, double *, double *);
void it2i0k0_wrap(double x, double *, double *);

int cfresnl_wrap(npy_cdouble x, npy_cdouble *sf, npy_cdouble *cf);
double cem_cva_wrap(double m, double q);
double sem_cva_wrap(double m, double q);
void cem_wrap(double m, double q, double x, double *csf, double *csd);
void sem_wrap(double m, double q, double x, double *csf, double *csd);
void mcm1_wrap(double m, double q, double x, double *f1r, double *d1r);
void msm1_wrap(double m, double q, double x, double *f1r, double *d1r);
void mcm2_wrap(double m, double q, double x, double *f2r, double *d2r);
void msm2_wrap(double m, double q, double x, double *f2r, double *d2r);
double pmv_wrap(double, double, double);
void pbwa_wrap(double, double, double *, double *);
void pbdv_wrap(double, double, double *, double *);
void pbvv_wrap(double, double, double *, double *);

void prolate_aswfa_wrap(double, double, double, double, double, double *, double *);
void prolate_radial1_wrap(double, double, double, double, double, double *, double *);
void prolate_radial2_wrap(double, double, double, double, double, double *, double *);
void oblate_aswfa_wrap(double, double, double, double, double, double *, double *);
void oblate_radial1_wrap(double, double, double, double, double, double *, double *);
void oblate_radial2_wrap(double, double, double, double, double, double *, double *);
double prolate_aswfa_nocv_wrap(double, double, double, double, double *);
double prolate_radial1_nocv_wrap(double, double, double, double, double *);
double prolate_radial2_nocv_wrap(double, double, double, double, double *);
double oblate_aswfa_nocv_wrap(double, double, double, double, double *);
double oblate_radial1_nocv_wrap(double, double, double, double, double *);
double oblate_radial2_nocv_wrap(double, double, double, double, double *);
double prolate_segv_wrap(double, double, double);
double oblate_segv_wrap(double, double, double);

void modified_fresnel_plus_wrap(double x, npy_cdouble *F, npy_cdouble *K);
void modified_fresnel_minus_wrap(double x, npy_cdouble *F, npy_cdouble *K);

void special_airy(double x, double *ai, double *aip, double *bi, double *bip);
void special_cairy(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip);

void special_airye(double z, double *ai, double *aip, double *bi, double *bip);
void special_cairye(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip);

void special_itairy(double x, double *apt, double *bpt, double *ant, double *bnt);

double special_cyl_bessel_j(double v, double z);
npy_cdouble special_ccyl_bessel_j(double v, npy_cdouble z);

double special_cyl_bessel_je(double v, double z);
npy_cdouble special_ccyl_bessel_je(double v, npy_cdouble z);

double special_cyl_bessel_y(double v, double x);
npy_cdouble special_ccyl_bessel_y(double v, npy_cdouble z);

double special_cyl_bessel_ye(double v, double z);
npy_cdouble special_ccyl_bessel_ye(double v, npy_cdouble z);

double special_cyl_bessel_i(double v, double z);
npy_cdouble special_ccyl_bessel_i(double v, npy_cdouble z);

double special_cyl_bessel_ie(double v, double z);
npy_cdouble special_ccyl_bessel_ie(double v, npy_cdouble z);

double special_cyl_bessel_k_int(int n, double z);
double special_cyl_bessel_k(double v, double z);
npy_cdouble special_ccyl_bessel_k(double v, npy_cdouble z);

double special_cyl_bessel_ke(double v, double z);
npy_cdouble special_ccyl_bessel_ke(double v, npy_cdouble z);

npy_cdouble special_ccyl_hankel_1(double v, npy_cdouble z);
npy_cdouble special_ccyl_hankel_2(double v, npy_cdouble z);
npy_cdouble special_ccyl_hankel_1e(double v, npy_cdouble z);
npy_cdouble special_ccyl_hankel_2e(double v, npy_cdouble z);

npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble zp);
double sin_pi(double x);
double gammaln_wrap(double x);

double binom_wrap(double n, double k);

double special_binom(double n, double k);

double special_cospi(double x);

double special_sinpi(double x);
npy_cdouble special_csinpi(npy_cdouble z);

double special_digamma(double z);
npy_cdouble special_cdigamma(npy_cdouble z);

float special_expitf(float x);
double special_expit(double x);
npy_longdouble special_expitl(npy_longdouble x);

npy_double special_exprel(npy_double x);

float special_log_expitf(float x);
double special_log_expit(double x);
npy_longdouble special_log_expitl(npy_longdouble x);

float special_logitf(float x);
double special_logit(double x);
npy_longdouble special_logitl(npy_longdouble x);

double special_loggamma(double x);
npy_cdouble special_cloggamma(npy_cdouble z);

double special_gamma(double x);
npy_cdouble special_cgamma(npy_cdouble z);

double special_hyp2f1(double a, double b, double c, double z);
npy_cdouble special_chyp2f1(double a, double b, double c, npy_cdouble z);

npy_cdouble special_lambertw(npy_cdouble z, long k, double tol);

double special_rgamma(double x);
npy_cdouble special_crgamma(npy_cdouble z);

npy_cdouble special_sph_harm(long m, long n, double theta, double phi);
npy_cdouble special_sph_harm_unsafe(double m, double n, double theta, double phi);

double binom_wrap(double n, double k);
npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble zp);
double cephes_hyp2f1_wrap(double a, double b, double c, double x);
double cephes_airy_wrap(double x, double *ai, double *aip, double *bi, double *bip);
double cephes_beta_wrap(double a, double b);
double cephes_lbeta_wrap(double a, double b);
double cephes_bdtr_wrap(double k, int n, double p);
double cephes_bdtri_wrap(double k, int n, double y);
double cephes_bdtrc_wrap(double k, int n, double p);
double cephes_cosm1_wrap(double x);
double cephes_expm1_wrap(double x);
double cephes_expn_wrap(int n, double x);
double cephes_log1p_wrap(double x);
double cephes_gamma_wrap(double x);
double cephes_gammasgn_wrap(double x);
double cephes_lgam_wrap(double x);
double cephes_iv_wrap(double v, double x);
double cephes_jv_wrap(double v, double x);
double cephes_ellpk_wrap(double x);
int cephes_ellpj_wrap(double u, double m, double *sn, double *cn, double *dn, double *ph);
int cephes_fresnl_wrap(double xxa, double *ssa, double *cca);
double cephes_nbdtr_wrap(int k, int n, double p);
double cephes_nbdtrc_wrap(int k, int n, double p);
double cephes_nbdtri_wrap(int k, int n, double p);
double cephes_ndtr_wrap(double x);
double cephes_ndtri_wrap(double x);
double cephes_pdtri_wrap(int k, double y);
double cephes_poch_wrap(double x, double m);
int cephes_sici_wrap(double x, double *si, double *ci);
int cephes_shichi_wrap(double x, double *si, double *ci);
double cephes__struve_asymp_large_z(double v, double z, int is_h, double *err);
double cephes__struve_bessel_series(double v, double z, int is_h, double *err);
double cephes__struve_power_series(double v, double z, int is_h, double *err);
double cephes_smirnov_wrap(int n, double x);
double cephes_smirnovc_wrap(int n, double x);
double cephes_smirnovi_wrap(int n, double x);
double cephes_smirnovci_wrap(int n, double x);
double cephes_smirnovp_wrap(int n, double x);
double cephes_yn_wrap(int n, double x);
double cephes_polevl_wrap(double x, const double coef[], int N);
double cephes_p1evl_wrap(double x, const double coef[], int N);
double special_wright_bessel(double a, double b, double x);

double special_scaled_exp1(double x);

double special_sph_bessel_j(long n, double x);
npy_cdouble special_csph_bessel_j(long n, npy_cdouble z);

double special_sph_bessel_j_jac(long n, double x);
npy_cdouble special_csph_bessel_j_jac(long n, npy_cdouble z);

double special_sph_bessel_y(long n, double x);
npy_cdouble special_csph_bessel_y(long n, npy_cdouble z);

double special_sph_bessel_y_jac(long n, double x);
npy_cdouble special_csph_bessel_y_jac(long n, npy_cdouble z);

double special_sph_bessel_i(long n, double x);
npy_cdouble special_csph_bessel_i(long n, npy_cdouble z);

double special_sph_bessel_i_jac(long n, double x);
npy_cdouble special_csph_bessel_i_jac(long n, npy_cdouble z);

double special_sph_bessel_k(long n, double x);
npy_cdouble special_csph_bessel_k(long n, npy_cdouble z);

double special_sph_bessel_k_jac(long n, double x);
npy_cdouble special_csph_bessel_k_jac(long n, npy_cdouble z);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */
