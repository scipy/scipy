#include "specfun_wrappers.h"
#include "sf_error.h"

#define CADDR(z) ((double *) (&(z))), (&(((double *) (&(z)))[1]))

#ifndef CMPLX
#define CMPLX(x, y) ((double complex)((double)(x) + I * (double)(y)))
#endif /* CMPLX */


npy_cdouble chyp2f1_wrap( double a, double b, double c, npy_cdouble z) {
  npy_cdouble outz;
  double complex z99;
  double complex outz99;
  int l1, l0, isfer = 0;

  l0 = ((c == floor(c)) && (c < 0));
  l1 = ((fabs(1-npy_creal(z)) < 1e-15) && (npy_cimag(z) == 0) && (c-a-b <= 0));
  if (l0 || l1) {
    sf_error("chyp2f1", SF_ERROR_OVERFLOW, NULL);
    NPY_CSETREAL(&outz, INFINITY);
    NPY_CSETIMAG(&outz, 0.0);
    return outz;
  }
  z99 = CMPLX(npy_creal(z), npy_cimag(z));
  outz99 = specfun_hygfz(a, b, c, z99, &isfer);
  NPY_CSETREAL(&outz, creal(outz99));
  NPY_CSETIMAG(&outz, cimag(outz99));
  if (isfer == 3) {
    sf_error("chyp2f1", SF_ERROR_OVERFLOW, NULL);
    NPY_CSETREAL(&outz, INFINITY);
    NPY_CSETIMAG(&outz, 0.0);
  } else if (isfer == 5) {
    sf_error("chyp2f1", SF_ERROR_LOSS, NULL);
  } else if (isfer != 0) {
    sf_error("chyp2f1", isfer, NULL);
    NPY_CSETREAL(&outz, NAN);
    NPY_CSETIMAG(&outz, NAN);
  }
  return outz;
}

npy_cdouble chyp1f1_wrap(double a, double b, npy_cdouble z) {
  npy_cdouble outz;
  double complex outz99;
  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));

  outz99 = specfun_cchg(a, b, z99);
  NPY_CSETREAL(&outz, creal(outz99));
  NPY_CSETIMAG(&outz, cimag(outz99));

  if (npy_creal(outz) == 1e300) {
    sf_error("chyp1f1", SF_ERROR_OVERFLOW, NULL);
    NPY_CSETREAL(&outz, INFINITY);
  }
  return outz;
}


double hypU_wrap(double a, double b, double x) {
  double out;
  int md; /* method code --- not returned */
  int isfer = 0;

  out = specfun_chgu(x, a, b, &md, &isfer);
  if (out == 1e300) {
      sf_error("hypU", SF_ERROR_OVERFLOW, NULL);
      out = INFINITY;
  }
  if (isfer == 6) {
    sf_error("hypU", SF_ERROR_NO_RESULT, NULL);
    out = NAN;
  } else if (isfer != 0) {
    sf_error("hypU", isfer, NULL);
    out = NAN;
  }
  return out;
}

double hyp1f1_wrap(double a, double b, double x) {
   double outy;

   outy = specfun_chgm(x, a, b);
   return outy;
}

int itairy_wrap(double x, double *apt, double *bpt, double *ant, double *bnt) {
  double tmp;
  int flag = 0;

  if (x < 0) {
    x = -x;
    flag = 1;
  }
  specfun_itairy(x, apt, bpt, ant, bnt);
  if (flag) {  /* negative limit -- switch signs and roles */
    tmp = *apt;
    *apt = -*ant;
    *ant = -tmp;
    tmp = *bpt;
    *bpt = -*bnt;
    *bnt = -tmp;
  }
  return 0;
}


double exp1_wrap(double x) {
  double out;

  out = specfun_e1xb(x);
  CONVINF("exp1", out);
  return out;
}

npy_cdouble cexp1_wrap(npy_cdouble z) {
  npy_cdouble outz;
  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));

  double complex outz99 = specfun_e1z(z99);
  NPY_CSETREAL(&outz, creal(outz99));
  NPY_CSETIMAG(&outz, cimag(outz99));
  ZCONVINF("cexp1", outz);
  return outz;
}

double expi_wrap(double x) {
  double out;

  out = specfun_eix(x);
  CONVINF("expi", out);
  return out;
}

npy_cdouble cexpi_wrap(npy_cdouble z) {
  npy_cdouble outz;
  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));

  double complex outz99 = specfun_eixz(z99);
  NPY_CSETREAL(&outz, creal(outz99));
  NPY_CSETIMAG(&outz, cimag(outz99));
  ZCONVINF("cexpi", outz);
  return outz;
}

npy_cdouble cerf_wrap(npy_cdouble z) {
  npy_cdouble outz;
  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));

  double complex outz99 = specfun_cerror(z99);
  NPY_CSETREAL(&outz, creal(outz99));
  NPY_CSETIMAG(&outz, cimag(outz99));
  return outz;
}

double itstruve0_wrap(double x) {
  double out;

  if (x<0) x=-x;
  out = specfun_itsh0(x);
  CONVINF("itstruve0", out);
  return out;
}

double it2struve0_wrap(double x) {
  double out;
  int flag=0;

  if (x < 0) {
      x=-x;
      flag=1;
  }
  out = specfun_itth0(x);
  CONVINF("it2struve0", out);
  if (flag) {
    out = M_PI - out;
  }
  return out;
}

double itmodstruve0_wrap(double x) {
  double out;

  if (x<0) {
      x=-x;
  }
  out = specfun_itsl0(x);
  CONVINF("itmodstruve0", out);
  return out;
}


double ber_wrap(double x)
{
  npy_cdouble Be;
  double ber, bei, ger, gei, der, dei, her, hei;

  if (x<0) {
      x=-x;
  }
  specfun_klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
  NPY_CSETREAL(&Be, ber);
  NPY_CSETIMAG(&Be, bei);
  ZCONVINF("ber", Be);
  return npy_creal(Be);
}

double bei_wrap(double x)
{
  npy_cdouble Be;
  double ber, bei, ger, gei, der, dei, her, hei;
  
  if (x<0) {
      x=-x;
  }
  specfun_klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
  NPY_CSETREAL(&Be, ber);
  NPY_CSETIMAG(&Be, bei);
  ZCONVINF("bei", Be);
  return npy_cimag(Be);
}

double ker_wrap(double x)
{
  npy_cdouble Ke;
  double ber, bei, ger, gei, der, dei, her, hei;
  if (x<0) {
      return NAN;
  }
  specfun_klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
  NPY_CSETREAL(&Ke, ger);
  NPY_CSETIMAG(&Ke, gei);
  ZCONVINF("ker", Ke);
  return npy_creal(Ke);
}

double kei_wrap(double x)
{
  npy_cdouble Ke;
  double ber, bei, ger, gei, der, dei, her, hei;
  if (x<0) {
      return NAN;
  }
  specfun_klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
  NPY_CSETREAL(&Ke, ger);
  NPY_CSETIMAG(&Ke, gei);
  ZCONVINF("kei", Ke);
  return npy_cimag(Ke);
}

double berp_wrap(double x)
{
  npy_cdouble Bep;
  double ber, bei, ger, gei, der, dei, her, hei;
  int flag = 0;

  if (x<0) {
      x=-x;
      flag=1;
  }
  specfun_klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
  NPY_CSETREAL(&Bep, der);
  NPY_CSETIMAG(&Bep, dei);
  ZCONVINF("berp", Bep);
  if (flag) {
      return -npy_creal(Bep);
  }
  return npy_creal(Bep);
}

double beip_wrap(double x)
{
  npy_cdouble Bep;
  double ber, bei, ger, gei, der, dei, her, hei;
  int flag = 0;

  if (x<0) {
      x=-x;
      flag=1;
  }
  specfun_klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
  NPY_CSETREAL(&Bep, der);
  NPY_CSETIMAG(&Bep, dei);
  ZCONVINF("beip", Bep);
  if (flag) {
      return -npy_cimag(Bep);
  }
  return npy_cimag(Bep);
}

double kerp_wrap(double x)
{
  npy_cdouble Kep;
  double ber, bei, ger, gei, der, dei, her, hei;

  if (x<0) {
      return NAN;
  }
  specfun_klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
  NPY_CSETREAL(&Kep, her);
  NPY_CSETIMAG(&Kep, hei);
  ZCONVINF("kerp", Kep);
  return npy_creal(Kep);
}

double keip_wrap(double x)
{
  npy_cdouble Kep;
  double ber, bei, ger, gei, der, dei, her, hei;

  if (x<0) {
      return NAN;
  }
  specfun_klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
  NPY_CSETREAL(&Kep, her);
  NPY_CSETIMAG(&Kep, hei);
  ZCONVINF("keip", Kep);
  return npy_cimag(Kep);
}


int kelvin_wrap(double x, npy_cdouble *Be, npy_cdouble *Ke, npy_cdouble *Bep, npy_cdouble *Kep) {
  int flag = 0;
  double ber, bei, ger, gei, der, dei, her, hei;
  if (x<0) {
      x=-x;
      flag=1;
  }

  specfun_klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
  NPY_CSETREAL(Be, ber);
  NPY_CSETIMAG(Be, bei);
  NPY_CSETREAL(Ke, ger);
  NPY_CSETIMAG(Ke, gei);
  NPY_CSETREAL(Bep, der);
  NPY_CSETIMAG(Bep, dei);
  NPY_CSETREAL(Kep, her);
  NPY_CSETIMAG(Kep, hei);

  ZCONVINF("klvna", *Be);
  ZCONVINF("klvna", *Ke);
  ZCONVINF("klvna", *Bep);
  ZCONVINF("klvna", *Kep);
  if (flag) {
    NPY_CSETREAL(Bep, -npy_creal(*Bep));
    NPY_CSETIMAG(Bep, -npy_cimag(*Bep));
    NPY_CSETREAL(Ke, NAN);
    NPY_CSETIMAG(Ke, NAN);
    NPY_CSETREAL(Kep, NAN);
    NPY_CSETIMAG(Kep, NAN);
  }
  return 0;
}

/* Integrals of bessel functions */

/* int(j0(t),t=0..x) */
/* int(y0(t),t=0..x) */

int it1j0y0_wrap(double x, double *j0int, double *y0int)
{
  int flag = 0;

  if (x < 0) {
      x = -x;
      flag=1;
  }
  specfun_itjya(x, j0int, y0int);
  if (flag) {
    *j0int = -(*j0int);
    *y0int = NAN;    /* domain error */
  }
  return 0;
}

/* int((1-j0(t))/t,t=0..x) */
/* int(y0(t)/t,t=x..inf) */

int it2j0y0_wrap(double x, double *j0int, double *y0int)
{
  int flag = 0;

  if (x < 0) {
      x=-x;
      flag=1;
  }
  specfun_ittjya(x, j0int, y0int);
  if (flag) {
    *y0int = NAN;  /* domain error */
  }
  return 0;
}

/* Integrals of modified bessel functions */

int it1i0k0_wrap(double x, double *i0int, double *k0int)
{
  int flag = 0;

  if (x < 0) {
      x = -x;
      flag=1;
  }
  specfun_itika(x, i0int, k0int);
  if (flag) {
    *i0int = -(*i0int);
    *k0int = NAN;    /* domain error */
  }
  return 0;
}

int it2i0k0_wrap(double x, double *i0int, double *k0int)
{
  int flag = 0;

  if (x < 0) {
      x=-x;
      flag=1;
  }
  specfun_ittika(x, i0int, k0int);
  if (flag) {
    *k0int = NAN;  /* domain error */
  }
  return 0;
}


/* Fresnel integrals of complex numbers */

int cfresnl_wrap(npy_cdouble z, npy_cdouble *zfs, npy_cdouble *zfc)
{
  double complex z99 = CMPLX(npy_creal(z), npy_cimag(z));
  double complex zfs99, zfc99, zfd;

  specfun_cfs(z99, &zfs99, &zfd);
  specfun_cfc(z99, &zfc99, &zfd);

  NPY_CSETREAL(zfs, creal(zfs99));
  NPY_CSETIMAG(zfs, cimag(zfs99));
  NPY_CSETREAL(zfc, creal(zfc99));
  NPY_CSETIMAG(zfc, cimag(zfc99));  

  return 0;
}

/* Mathieu functions */
/* Characteristic values */
double cem_cva_wrap(double m, double q) {
  int int_m, kd=1;
  double out;

  if ((m < 0) || (m != floor(m))) {
    sf_error("cem_cva", SF_ERROR_DOMAIN, NULL);
    return NAN;
  }
  int_m = (int)m;
  if (q < 0) {
    /* https://dlmf.nist.gov/28.2#E26 */
    if (int_m % 2 == 0) {
      return cem_cva_wrap(m, -q);
    }
    else {
      return sem_cva_wrap(m, -q);
    }
  }

  if (int_m % 2) {
      kd=2;
  }
  out = specfun_cva2(kd, int_m, q);
  return out;
}

double sem_cva_wrap(double m, double q) {
  int int_m, kd=4;
  double out;

  if ((m <= 0) || (m != floor(m))) {
    sf_error("cem_cva", SF_ERROR_DOMAIN, NULL);
    return NAN;
  }
  int_m = (int)m;
  if (q < 0) {
    /* https://dlmf.nist.gov/28.2#E26 */
    if (int_m % 2 == 0) {
      return sem_cva_wrap(m, -q);
    }
    else {
      return cem_cva_wrap(m, -q);
    }
  }
  if (int_m % 2) {
      kd=3;
  }
  out = specfun_cva2(kd, int_m, q);
  return out;
}

/* Mathieu functions */
int cem_wrap(double m, double q, double x, double *csf, double *csd)
{
  int int_m, kf=1, sgn;
  double f = 0.0, d = 0.0;
  if ((m < 0) || (m != floor(m))) {
    *csf = NAN;
    *csd = NAN;
    sf_error("cem", SF_ERROR_DOMAIN, NULL);
    return -1;
  }
  int_m = (int)m;
  if (q < 0) {
      /* https://dlmf.nist.gov/28.2#E34 */
      if (int_m % 2 == 0) {
        sgn = ((int_m/2) % 2 == 0) ? 1 : -1;
        cem_wrap(m, -q, 90 - x, &f, &d);
        *csf =  sgn * f;
        *csd = -sgn * d;
        return 0;
      }
      else {
        sgn = ((int_m/2) % 2 == 0) ? 1 : -1;
        sem_wrap(m, -q, 90 - x, &f, &d);
        *csf =  sgn * f;
        *csd = -sgn * d;
        return 0;
      }
  }
  specfun_mtu0(kf, int_m, q, x, csf, csd);
  return 0;
}

int sem_wrap(double m, double q, double x, double *csf, double *csd)
{
  int int_m, kf=2, sgn;
  double f = 0.0, d = 0.0;
  if ((m < 0) || (m != floor(m))) {
    *csf = NAN;
    *csd = NAN;
    sf_error("sem", SF_ERROR_DOMAIN, NULL);
    return -1;
  }
  int_m = (int)m;
  if (int_m == 0) {
    *csf = 0;
    *csd = 0;
    return 0;
  }
  if (q < 0) {
      /* https://dlmf.nist.gov/28.2#E34 */
      if (int_m % 2 == 0) {
        sgn = ((int_m/2) % 2 == 0) ? -1 : 1;
        sem_wrap(m, -q, 90 - x, &f, &d);
        *csf =  sgn * f;
        *csd = -sgn * d;
        return 0;
      }
      else {
        sgn = ((int_m/2) % 2 == 0) ? 1 : -1;
        cem_wrap(m, -q, 90 - x, &f, &d);
        *csf =  sgn * f;
        *csd = -sgn * d;
        return 0;
      }
  }
  specfun_mtu0(kf, int_m, q, x, csf, csd);
  return 0;
}


int mcm1_wrap(double m, double q, double x, double *f1r, double *d1r)
{
  int int_m, kf=1, kc=1;
  double f2r = 0.0, d2r = 0.0;

  if ((m < 0) || (m != floor(m)) || (q<0)) {
    *f1r = NAN;
    *d1r = NAN;
    sf_error("mcm1", SF_ERROR_DOMAIN, NULL);
    return -1;
  }
  int_m = (int)m;
  specfun_mtu12(kf, kc, int_m, q, x, f1r, d1r, &f2r, &d2r);
  return 0;
}

int msm1_wrap(double m, double q, double x, double *f1r, double *d1r)
{
  int int_m, kf=2, kc=1;
  double f2r = 0.0, d2r = 0.0;

  if ((m < 1) || (m != floor(m)) || (q<0)) {
    *f1r = NAN;
    *d1r = NAN;
    sf_error("msm1", SF_ERROR_DOMAIN, NULL);
    return -1;
  }
  int_m = (int)m;
  specfun_mtu12(kf, kc, int_m, q, x, f1r, d1r, &f2r, &d2r);
  return 0;
}

int mcm2_wrap(double m, double q, double x, double *f2r, double *d2r)
{
  int int_m, kf=1, kc=2;
  double f1r = 0.0, d1r = 0.0;

  if ((m < 0) || (m != floor(m)) || (q<0)) {
    *f2r = NAN;
    *d2r = NAN;
    sf_error("mcm2", SF_ERROR_DOMAIN, NULL);
    return -1;
  }
  int_m = (int)m;
  specfun_mtu12(kf, kc, int_m, q, x, &f1r, &d1r, f2r, d2r);
  return 0;
}

int msm2_wrap(double m, double q, double x, double *f2r, double *d2r)
{
  int int_m, kf=2, kc=2;
  double f1r = 0.0, d1r = 0.0;

  if ((m < 1) || (m != floor(m)) || (q<0)) {
    *f2r = NAN;
    *d2r = NAN;
    sf_error("msm2", SF_ERROR_DOMAIN, NULL);
    return -1;
  }
  int_m = (int)m;
  specfun_mtu12(kf, kc, int_m, q, x, &f1r, &d1r, f2r, d2r);
  return 0;
}


double pmv_wrap(double m, double v, double x){
  int int_m;
  double out;

  if (m != floor(m)) {
      return NAN;
  }
  int_m = (int)m;
  out = specfun_lpmv(x, int_m, v);
  CONVINF("pmv", out);
  return out;
}


/*
 * If x > 0 return w1f and w1d. Otherwise set x = abs(x) and return
 * w2f and -w2d.
*/
int pbwa_wrap(double a, double x, double *wf, double *wd) {
  int flag = 0;
  double w1f = 0.0, w1d = 0.0, w2f = 0.0, w2d = 0.0;

  if (x < -5 || x > 5 || a < -5 || a > 5) {
    /*
     * The Zhang and Jin implementation only uses Taylor series;
     * return NaN outside of the range which they are accurate.
     */
    *wf = NAN;
    *wd = NAN;
    sf_error("pbwa", SF_ERROR_LOSS, NULL);
    return 0;
  }

  if (x < 0) {
    x = -x;
    flag = 1;
  }
  specfun_pbwa(a, x, &w1f, &w1d, &w2f, &w2d);
  if (flag) {
    *wf = w2f;
    *wd = -w2d;
  }
  else {
    *wf = w1f;
    *wd = w1d;
  }
  return 0;
}

int pbdv_wrap(double v, double x, double *pdf, double *pdd) {

  double *dv;
  double *dp;
  int num;

  if (isnan(v) || isnan(x)) {
    *pdf = NAN;
    *pdd = NAN;
    return 0;
  }
  /* NB. Indexing of DV/DP in specfun.f:PBDV starts from 0, hence +2 */
  num = ABS((int)v) + 2;
  dv = (double *)PyMem_Malloc(sizeof(double)*2*num);
  if (dv==NULL) {
    sf_error("pbdv", SF_ERROR_OTHER, "memory allocation error");
    *pdf = NAN;
    *pdd = NAN;
    return -1;
  }
  dp = dv + num;
  specfun_pbdv(x, v, dv, dp, pdf, pdd);
  PyMem_Free(dv);
  return 0;
}

int pbvv_wrap(double v, double x, double *pvf, double *pvd) {
  double *vv;
  double *vp;
  int num;

  if (isnan(v) || isnan(x)) {
    *pvf = NAN;
    *pvd = NAN;
    return 0;
  }
  /* NB. Indexing of DV/DP in specfun.f:PBVV starts from 0, hence +2 */
  num = ABS((int)v) + 2;
  vv = (double *)PyMem_Malloc(sizeof(double)*2*num);
  if (vv==NULL) {
    sf_error("pbvv", SF_ERROR_OTHER, "memory allocation error");
    *pvf = NAN;
    *pvd = NAN;
    return -1;
  }
  vp = vv + num;
  specfun_pbvv(x, v, vv, vp, pvf, pvd);
  PyMem_Free(vv);
  return 0;
}

double prolate_segv_wrap(double m, double n, double c)
{
  int kd=1;
  int int_m, int_n;
  double cv = 0.0, *eg;

  if ((m<0) || (n<m) || (m!=floor(m)) || (n!=floor(n)) || ((n-m)>198)) {
    return NAN;
  }
  int_m = (int) m;
  int_n = (int) n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    sf_error("prolate_segv", SF_ERROR_OTHER, "memory allocation error");
    return NAN;
  }
  specfun_segv(int_m, int_n, c, kd, &cv, eg);
  PyMem_Free(eg);
  return cv;
}

double oblate_segv_wrap(double m, double n, double c)
{
  int kd=-1;
  int int_m, int_n;
  double cv= 0.0, *eg;

  if ((m<0) || (n<m) || (m!=floor(m)) || (n!=floor(n)) || ((n-m)>198)) {
    return NAN;
  }
  int_m = (int) m;
  int_n = (int) n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    sf_error("oblate_segv", SF_ERROR_OTHER, "memory allocation error");
    return NAN;
  }
  specfun_segv(int_m, int_n, c, kd, &cv,eg);
  PyMem_Free(eg);
  return cv;
}


double prolate_aswfa_nocv_wrap(double m, double n, double c, double x, double *s1d)
{
  int kd = 1;
  int int_m, int_n;
  double cv = 0.0, s1f, *eg;

  if ((x >=1) || (x <=-1) || (m<0) || (n<m) || \
      (m!=floor(m)) || (n!=floor(n)) || ((n-m)>198)) {
    sf_error("prolate_aswfa_nocv", SF_ERROR_DOMAIN, NULL);
    *s1d = NAN;
    return NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    sf_error("prolate_aswfa_nocv", SF_ERROR_OTHER, "memory allocation error");
    *s1d = NAN;
    return NAN;
  }
  specfun_segv(int_m, int_n, c, kd, &cv, eg);
  specfun_aswfa(x, int_m, int_n, c, kd, cv, &s1f, s1d);
  PyMem_Free(eg);
  return s1f;
}


double oblate_aswfa_nocv_wrap(double m, double n, double c, double x, double *s1d)
{
  int kd = -1;
  int int_m, int_n;
  double cv = 0.0, s1f = 0.0, *eg;

  if ((x >=1) || (x <=-1) || (m<0) || (n<m) || \
      (m!=floor(m)) || (n!=floor(n)) || ((n-m)>198)) {
    sf_error("oblate_aswfa_nocv", SF_ERROR_DOMAIN, NULL);
    *s1d = NAN;
    return NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    sf_error("oblate_aswfa_nocv", SF_ERROR_OTHER, "memory allocation error");
    *s1d = NAN;
    return NAN;
  }
  specfun_segv(int_m, int_n, c, kd, &cv, eg);
  specfun_aswfa(x, int_m, int_n, c, kd, cv, &s1f, s1d);
  PyMem_Free(eg);
  return s1f;
}


int prolate_aswfa_wrap(double m, double n, double c, double cv, double x, double *s1f, double *s1d)
{
  int kd = 1;
  int int_m, int_n;

  if ((x >=1) || (x <=-1) || (m<0) || (n<m) || \
      (m!=floor(m)) || (n!=floor(n))) {
    sf_error("prolate_aswfa", SF_ERROR_DOMAIN, NULL);
    *s1f = NAN;
    *s1d = NAN;
    return 0;
  }
  int_m = (int )m;
  int_n = (int )n;
  specfun_aswfa(x, int_m, int_n, c, kd, cv, s1f, s1d);
  return 0;
}


int oblate_aswfa_wrap(double m, double n, double c, double cv, double x, double *s1f, double *s1d)
{
  int kd = -1;
  int int_m, int_n;

  if ((x >=1) || (x <=-1) || (m<0) || (n<m) || \
      (m!=floor(m)) || (n!=floor(n))) {
    sf_error("oblate_aswfa", SF_ERROR_DOMAIN, NULL);
    *s1f = NAN;
    *s1d = NAN;
    return 0;
  }
  int_m = (int )m;
  int_n = (int )n;
  specfun_aswfa(x, int_m, int_n, c, kd, cv, s1f, s1d);
  return 0;
}


double prolate_radial1_nocv_wrap(double m, double n, double c, double x, double *r1d)
{
  int kf=1, kd=1;
  double r2f = 0.0, r2d = 0.0, r1f = 0.0, cv = 0.0, *eg;
  int int_m, int_n;

  if ((x <=1.0) || (m<0) || (n<m) || \
     (m!=floor(m)) || (n!=floor(n)) || ((n-m)>198)) {
    sf_error("prolate_radial1_nocv", SF_ERROR_DOMAIN, NULL);
    *r1d = NAN;
    return NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    sf_error("prolate_radial1_nocv", SF_ERROR_OTHER, "memory allocation error");
    *r1d = NAN;
    return NAN;
  }
  specfun_segv(int_m, int_n, c, kd, &cv, eg);
  specfun_rswfp(int_m, int_n, c, x, cv, kf, &r1f, r1d, &r2f, &r2d);
  PyMem_Free(eg);
  return r1f;
}

double prolate_radial2_nocv_wrap(double m, double n, double c, double x, double *r2d)
{
  int kf=2, kd=1;
  double r1f = 0.0, r1d = 0.0, r2f = 0.0, cv = 0.0, *eg;
  int int_m, int_n;

  if ((x <=1.0) || (m<0) || (n<m) || \
     (m!=floor(m)) || (n!=floor(n)) || ((n-m)>198)) {
    sf_error("prolate_radial2_nocv", SF_ERROR_DOMAIN, NULL);
    *r2d = NAN;
    return NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    sf_error("prolate_radial2_nocv", SF_ERROR_OTHER, "memory allocation error");
    *r2d = NAN;
    return NAN;
  }
  specfun_segv(int_m, int_n, c, kd, &cv, eg);
  specfun_rswfp(int_m, int_n, c, x, cv, kf, &r1f, &r1d, &r2f, r2d);
  PyMem_Free(eg);
  return r2f;
}

int prolate_radial1_wrap(double m, double n, double c, double cv, double x, double *r1f, double *r1d)
{
  int kf=1;
  double r2f = 0.0, r2d = 0.0;
  int int_m, int_n;

  if ((x <= 1.0) || (m<0) || (n<m) || \
     (m!=floor(m)) || (n!=floor(n))) {
    sf_error("prolate_radial1", SF_ERROR_DOMAIN, NULL);
    *r1f = NAN;
    *r1d = NAN;
    return 0;
  }
  int_m = (int )m;
  int_n = (int )n;
  specfun_rswfp(int_m, int_n, c, x, cv, kf, r1f, r1d, &r2f, &r2d);
  return 0;
}

int prolate_radial2_wrap(double m, double n, double c, double cv, double x, double *r2f, double *r2d)
{
  int kf=2;
  double r1f = 0.0, r1d = 0.0;
  int int_m, int_n;

  if ((x <= 1.0) || (m<0) || (n<m) || \
     (m!=floor(m)) || (n!=floor(n))) {
    sf_error("prolate_radial2", SF_ERROR_DOMAIN, NULL);
    *r2f = NAN;
    *r2d = NAN;
    return 0;
  }
  int_m = (int )m;
  int_n = (int )n;
  specfun_rswfp(int_m, int_n, c, x, cv, kf, &r1f, &r1d, r2f, r2d);
  return 0;
}

double oblate_radial1_nocv_wrap(double m, double n, double c, double x, double *r1d)
{
  int kf=1, kd=-1;
  double r2f = 0.0, r2d = 0.0, r1f = 0.0, cv = 0.0, *eg;
  int int_m, int_n;

  if ((x < 0.0) || (m<0) || (n<m) || \
     (m!=floor(m)) || (n!=floor(n)) || ((n-m)>198)) {
    sf_error("oblate_radial1_nocv", SF_ERROR_DOMAIN, NULL);
    *r1d = NAN;
    return NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    sf_error("oblate_radial1_nocv", SF_ERROR_OTHER, "memory allocation error");
    *r1d = NAN;
    return NAN;
  }
  specfun_segv(int_m, int_n, c, kd, &cv, eg);
  specfun_rswfo(int_m, int_n, c, x, cv, kf, &r1f, r1d, &r2f, &r2d);
  PyMem_Free(eg);
  return r1f;
}

double oblate_radial2_nocv_wrap(double m, double n, double c, double x, double *r2d)
{
  int kf=2, kd=-1;
  double r1f = 0.0, r1d = 0.0, r2f = 0.0, cv = 0.0, *eg;
  int int_m, int_n;

  if ((x < 0.0) || (m<0) || (n<m) || \
     (m!=floor(m)) || (n!=floor(n)) || ((n-m)>198)) {
    sf_error("oblate_radial2_nocv", SF_ERROR_DOMAIN, NULL);
    *r2d = NAN;
    return NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    sf_error("oblate_radial2_nocv", SF_ERROR_OTHER, "memory allocation error");
    *r2d = NAN;
    return NAN;
  }
  specfun_segv(int_m, int_n, c, kd, &cv, eg);
  specfun_rswfo(int_m, int_n, c, x, cv, kf, &r1f, &r1d, &r2f, r2d);
  PyMem_Free(eg);
  return r2f;
}

int oblate_radial1_wrap(double m, double n, double c, double cv, double x, double *r1f, double *r1d)
{
  int kf=1;
  double r2f = 0.0, r2d = 0.0;
  int int_m, int_n;

  if ((x <0.0) || (m<0) || (n<m) || \
     (m!=floor(m)) || (n!=floor(n))) {
    sf_error("oblate_radial1", SF_ERROR_DOMAIN, NULL);
    *r1f = NAN;
    *r1d = NAN;
    return 0;
  }
  int_m = (int )m;
  int_n = (int )n;
  specfun_rswfo(int_m, int_n, c, x, cv, kf, r1f, r1d, &r2f, &r2d);
  return 0;
}

int oblate_radial2_wrap(double m, double n, double c, double cv, double x, double *r2f, double *r2d)
{
  int kf=2;
  double r1f = 0.0, r1d = 0.0;
  int int_m, int_n;

  if ((x <0.0) || (m<0) || (n<m) || \
     (m!=floor(m)) || (n!=floor(n))) {
    sf_error("oblate_radial2", SF_ERROR_DOMAIN, NULL);
    *r2f = NAN;
    *r2d = NAN;
    return 0;
  }
  int_m = (int )m;
  int_n = (int )n;
  specfun_rswfo(int_m, int_n, c, x, cv, kf, &r1f, &r1d, r2f, r2d);
  return 0;
}


int modified_fresnel_plus_wrap(double x, npy_cdouble *Fplus, npy_cdouble *Kplus)
{
  int ks=0;
  double fr = 0.0, gr = 0.0, fi = 0.0, gi = 0.0, fa = 0.0, ga = 0.0, fm = 0.0, gm = 0.0;

  specfun_ffk(ks, x, &fr, &fi, &fm, &fa, &gr, &gi, &gm, &ga);
  NPY_CSETREAL(Fplus, fr);
  NPY_CSETIMAG(Fplus, fi);
  NPY_CSETREAL(Kplus, gr);
  NPY_CSETIMAG(Kplus, gi); 
  return 0;
}

int modified_fresnel_minus_wrap(double x, npy_cdouble *Fminus, npy_cdouble *Kminus)
{
  int ks=1;
  double fr = 0.0, gr = 0.0, fi = 0.0, gi = 0.0, fa = 0.0, ga = 0.0, fm = 0.0, gm = 0.0;

  specfun_ffk(ks, x, &fr, &fi, &fm, &fa, &gr, &gi, &gm, &ga);
  NPY_CSETREAL(Fminus, fr);
  NPY_CSETIMAG(Fminus, fi);
  NPY_CSETREAL(Kminus, gr);
  NPY_CSETIMAG(Kminus, gi); 
  return 0;
}
