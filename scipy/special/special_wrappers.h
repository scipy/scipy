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
double exp1_wrap(double x);
double expi_wrap(double x);
npy_cdouble cexp1_wrap(npy_cdouble z);
npy_cdouble cexpi_wrap(npy_cdouble z);
npy_cdouble cerf_wrap(npy_cdouble z);
void itairy_wrap(double x, double *apt, double *bpt, double *ant, double *bnt);

double struve_wrap(double v, double x);
double itstruve0_wrap(double x);
double it2struve0_wrap(double x);

double modstruve_wrap(double v, double x);
double itmodstruve0_wrap(double x);

double ber_wrap(double x);
double bei_wrap(double x);
double ker_wrap(double x);
double kei_wrap(double x);
double berp_wrap(double x);
double beip_wrap(double x);
double kerp_wrap(double x);
double keip_wrap(double x);

void kelvin_wrap(double x, npy_cdouble *Be, npy_cdouble *Ke, npy_cdouble *Bep, npy_cdouble *Kep);

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

void airy_wrap(double x, double *ai, double *aip, double *bi, double *bip);
void cairy_wrap(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip);
void cairy_wrap_e(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip);
void cairy_wrap_e_real(double z, double *ai, double *aip, double *bi, double *bip);
npy_cdouble cbesi_wrap(double v, npy_cdouble z);
npy_cdouble cbesi_wrap_e(double v, npy_cdouble z);
double cbesi_wrap_e_real(double v, double z);
npy_cdouble cbesj_wrap(double v, npy_cdouble z);
npy_cdouble cbesj_wrap_e(double v, npy_cdouble z);
double cbesj_wrap_real(double v, double z);
double cbesj_wrap_e_real(double v, double z);
npy_cdouble cbesy_wrap(double v, npy_cdouble z);
double cbesy_wrap_real(double v, double x);
npy_cdouble cbesy_wrap_e(double v, npy_cdouble z);
double cbesy_wrap_e_real(double v, double z);
npy_cdouble cbesk_wrap(double v, npy_cdouble z);
npy_cdouble cbesk_wrap_e(double v, npy_cdouble z);
double cbesk_wrap_real(double v, double z);
double cbesk_wrap_e_real(double v, double z);
double cbesk_wrap_real_int(int n, double z);
npy_cdouble cbesh_wrap1(double v, npy_cdouble z);
npy_cdouble cbesh_wrap1_e(double v, npy_cdouble z);
npy_cdouble cbesh_wrap2(double v, npy_cdouble z);
npy_cdouble cbesh_wrap2_e(double v, npy_cdouble z);
npy_cdouble hyp2f1_complex_wrap(double a, double b, double c, npy_cdouble zp);

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

double special_wright_bessel(double a, double b, double x);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */
