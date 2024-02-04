/* This file is a collection of wrappers around the
 *  Specfun Fortran library of functions
 */

#include "specfun_wrappers.h"

#include "sf_error.h"

#define CADDR(z) ((double *) (&(z))), (&(((double *) (&(z)))[1]))
#define F2C_CST(z) ((double *) (z)), (&(((double *) (z))[1]))

#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif

extern void F_FUNC(cgama,CGAMA)(double*,double*,int*,double*,double*);
extern void F_FUNC(hygfz,HYGFZ)(double*,double*,double*,npy_cdouble*,npy_cdouble*,int*);
extern void F_FUNC(cchg,CCHG)(double*,double*,npy_cdouble*,npy_cdouble*);
extern void F_FUNC(chgm,CHGM)(double*,double*,double*,double*);
extern void F_FUNC(chgu,CHGU)(double*,double*,double*,double*,int*,int*);
extern void F_FUNC(itairy,ITAIRY)(double*,double*,double*,double*,double*);
extern void F_FUNC(e1xb,E1XB)(double*,double*);
extern void F_FUNC(e1z,E1Z)(npy_cdouble*,npy_cdouble*);
extern void F_FUNC(eix,EIX)(double*,double*);
extern void F_FUNC(eixz,EIXZ)(npy_cdouble*,npy_cdouble*);
extern void F_FUNC(cerror,CERROR)(npy_cdouble*,npy_cdouble*);
extern void F_FUNC(stvh0,STVH0)(double*,double*);
extern void F_FUNC(stvh1,STVH1)(double*,double*);
extern void F_FUNC(stvhv,STVHV)(double*,double*,double*);
extern void F_FUNC(stvl0,STVL0)(double*,double*);
extern void F_FUNC(stvl1,STVL1)(double*,double*);
extern void F_FUNC(stvlv,STVLV)(double*,double*,double*);
extern void F_FUNC(itsh0,ITSH0)(double*,double*);
extern void F_FUNC(itth0,ITTH0)(double*,double*);
extern void F_FUNC(itsl0,ITSL0)(double*,double*);
extern void F_FUNC(klvna,KLVNA)(double*,double*,double*,double*,double*,double*,double*,double*,double*);
extern void F_FUNC(itjya,ITJYA)(double*,double*,double*);
extern void F_FUNC(ittjya,ITTJYA)(double*,double*,double*);
extern void F_FUNC(itika,ITIKA)(double*,double*,double*);
extern void F_FUNC(ittika,ITTIKA)(double*,double*,double*);
extern void F_FUNC(cfc,CFC)(npy_cdouble*,npy_cdouble*,npy_cdouble*);
extern void F_FUNC(cfs,CFS)(npy_cdouble*,npy_cdouble*,npy_cdouble*);
extern void F_FUNC(cva2,CVA2)(int*,int*,double*,double*);
extern void F_FUNC(mtu0,MTU0)(int*,int*,double*,double*,double*,double*);
extern void F_FUNC(mtu12,MTU12)(int*,int*,int*,double*,double*,double*,double*,double*,double*);
extern void F_FUNC(lpmv,LPMV)(double*,int*,double*,double*);
extern void F_FUNC(pbwa,PBWA)(double*,double*,double*,double*,double*,double*);
extern void F_FUNC(pbdv,PBDV)(double*,double*,double*,double*,double*,double*);
extern void F_FUNC(pbvv,PBVV)(double*,double*,double*,double*,double*,double*);
extern void F_FUNC(segv,SEGV)(int*,int*,double*,int*,double*,double*);
extern void F_FUNC(aswfa,ASWFA)(int*,int*,double*,double*,int*,double*,double*,double*);
extern void F_FUNC(rswfp,RSWFP)(int*,int*,double*,double*,double*,int*,double*,double*,double*,double*);
extern void F_FUNC(rswfo,RSWFO)(int*,int*,double*,double*,double*,int*,double*,double*,double*,double*);
extern void F_FUNC(ffk,FFK)(int*,double*,double*,double*,double*,double*,double*,double*,double*,double*);


/* This must be linked with fortran
 */

npy_cdouble chyp2f1_wrap( double a, double b, double c, npy_cdouble z) {
  npy_cdouble outz;
  int l1, l0, isfer = 0;


  l0 = ((c == floor(c)) && (c < 0));
  l1 = ((fabs(1-npy_creal(z)) < 1e-15) && (npy_cimag(z) == 0) && (c-a-b <= 0));
  if (l0 || l1) {
    sf_error("chyp2f1", SF_ERROR_OVERFLOW, NULL);
    NPY_CSETREAL(&outz, INFINITY);
    NPY_CSETIMAG(&outz, 0.0);
    return outz;
  }
  F_FUNC(hygfz, HYGFZ)(&a, &b, &c, &z, &outz, &isfer);
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

  F_FUNC(cchg,CCHG)(&a, &b, &z, &outz);
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

  F_FUNC(chgu,CHGU)(&a, &b, &x, &out, &md, &isfer);
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

   F_FUNC(chgm,CHGM)(&a, &b, &x, &outy);
   return outy;
}

int itairy_wrap(double x, double *apt, double *bpt, double *ant, double *bnt) {
  double tmp;
  int flag = 0;

  if (x < 0) {
    x = -x;
    flag = 1;
  }
  F_FUNC(itairy,ITAIRY)(&x, apt, bpt, ant, bnt);
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

  F_FUNC(e1xb,E1XB)(&x, &out);
  CONVINF("exp1", out);
  return out;
}

npy_cdouble cexp1_wrap(npy_cdouble z) {
  npy_cdouble outz;

  F_FUNC(e1z,E1Z)(&z, &outz);
  ZCONVINF("cexp1", outz);
  return outz;
}

double expi_wrap(double x) {
  double out;

  F_FUNC(eix,EIX)(&x, &out);
  CONVINF("expi", out);
  return out;
}

npy_cdouble cexpi_wrap(npy_cdouble z) {
  npy_cdouble outz;

  F_FUNC(eixz,EIXZ)(&z, &outz);
  ZCONVINF("cexpi", outz);
  return outz;
}

npy_cdouble cerf_wrap(npy_cdouble z) {
  npy_cdouble outz;

  F_FUNC(cerror,CERROR)(&z, &outz);
  return outz;
}

double itstruve0_wrap(double x) {
  double out;

  if (x<0) x=-x;
  F_FUNC(itsh0,ITSH0)(&x,&out);
  CONVINF("itstruve0", out);
  return out;
}

double it2struve0_wrap(double x) {
  double out;
  int flag=0;

  if (x<0) {x=-x; flag=1;}
  F_FUNC(itth0,ITTH0)(&x,&out);
  CONVINF("it2struve0", out);
  if (flag) {
    out = M_PI - out;
  }
  return out;
}

double itmodstruve0_wrap(double x) {
  double out;

  if (x<0) x=-x;
  F_FUNC(itsl0,ITSL0)(&x,&out);
  CONVINF("itmodstruve0", out);
  return out;
}


double ber_wrap(double x)
{
  npy_cdouble Be, Ke, Bep, Kep;

  if (x<0) x=-x;
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF("ber", Be);
  return npy_creal(Be);
}

double bei_wrap(double x)
{
  npy_cdouble Be, Ke, Bep, Kep;

  if (x<0) x=-x;
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF("bei", Be);
  return npy_cimag(Be);
}

double ker_wrap(double x)
{
  npy_cdouble Be, Ke, Bep, Kep;

  if (x<0) return NAN;
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF("ker", Ke);
  return npy_creal(Ke);
}

double kei_wrap(double x)
{
  npy_cdouble Be, Ke, Bep, Kep;

  if (x<0) return NAN;
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF("kei", Ke);
  return npy_cimag(Ke);
}

double berp_wrap(double x)
{
  npy_cdouble Be, Ke, Bep, Kep;
  int flag = 0;

  if (x<0) {x=-x; flag=1;}
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF("berp", Bep);
  if (flag) return -npy_creal(Bep);
  return npy_creal(Bep);
}

double beip_wrap(double x)
{
  npy_cdouble Be, Ke, Bep, Kep;
  int flag = 0;

  if (x<0) {x=-x; flag=1;}
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF("beip", Bep);
  if (flag) return -npy_cimag(Bep);
  return npy_cimag(Bep);
}

double kerp_wrap(double x)
{
  npy_cdouble Be, Ke, Bep, Kep;

  if (x<0) return NAN;
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF("kerp", Kep);
  return npy_creal(Kep);
}

double keip_wrap(double x)
{
  npy_cdouble Be, Ke, Bep, Kep;

  if (x<0) return NAN;
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF("keip", Kep);
  return npy_cimag(Kep);
}


int kelvin_wrap(double x, npy_cdouble *Be, npy_cdouble *Ke, npy_cdouble *Bep, npy_cdouble *Kep) {
  int flag = 0;

  if (x<0) {x=-x; flag=1;}
  F_FUNC(klvna,KLVNA)(&x, F2C_CST(Be), F2C_CST(Ke), F2C_CST(Bep), F2C_CST(Kep));
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

  if (x < 0) {x = -x; flag=1;}
  F_FUNC(itjya, ITJYA)(&x, j0int, y0int);
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

  if (x < 0) {x=-x; flag=1;}
  F_FUNC(ittjya, ITTJYA)(&x, j0int, y0int);
  if (flag) {
    *y0int = NAN;  /* domain error */
  }
  return 0;
}

/* Integrals of modified bessel functions */

int it1i0k0_wrap(double x, double *i0int, double *k0int)
{
  int flag = 0;

  if (x < 0) {x = -x; flag=1;}
  F_FUNC(itika, ITIKA)(&x, i0int, k0int);
  if (flag) {
    *i0int = -(*i0int);
    *k0int = NAN;    /* domain error */
  }
  return 0;
}

int it2i0k0_wrap(double x, double *i0int, double *k0int)
{
  int flag = 0;

  if (x < 0) {x=-x; flag=1;}
  F_FUNC(ittika, ITTIKA)(&x, i0int, k0int);
  if (flag) {
    *k0int = NAN;  /* domain error */
  }
  return 0;
}


/* Fresnel integrals of complex numbers */

int cfresnl_wrap(npy_cdouble z, npy_cdouble *zfs, npy_cdouble *zfc)
{
  npy_cdouble zfd;
  F_FUNC(cfs,CFS)(&z,zfs,&zfd);
  F_FUNC(cfc,CFC)(&z,zfc,&zfd);
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
  int_m = (int )m;
  if (q < 0) {
    /* https://dlmf.nist.gov/28.2#E26 */
    if (int_m % 2 == 0) {
      return cem_cva_wrap(m, -q);
    }
    else {
      return sem_cva_wrap(m, -q);
    }
  }

  if (int_m % 2) kd=2;
  F_FUNC(cva2,CVA2)(&kd, &int_m, &q, &out);
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
  if (int_m % 2) kd=3;
  F_FUNC(cva2,CVA2)(&kd, &int_m, &q, &out);
  return out;
}

/* Mathieu functions */
int cem_wrap(double m, double q, double x, double *csf, double *csd)
{
  int int_m, kf=1, sgn;
  double f, d;
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
  F_FUNC(mtu0,MTU0)(&kf,&int_m, &q, &x, csf, csd);
  return 0;
}

int sem_wrap(double m, double q, double x, double *csf, double *csd)
{
  int int_m, kf=2, sgn;
  double f, d;
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
  F_FUNC(mtu0,MTU0)(&kf,&int_m, &q, &x, csf, csd);
  return 0;
}


int mcm1_wrap(double m, double q, double x, double *f1r, double *d1r)
{
  int int_m, kf=1, kc=1;
  double f2r, d2r;

  if ((m < 0) || (m != floor(m)) || (q<0)) {
    *f1r = NAN;
    *d1r = NAN;
    sf_error("mcm1", SF_ERROR_DOMAIN, NULL);
    return -1;
  }
  int_m = (int )m;
  F_FUNC(mtu12,MTU12)(&kf,&kc,&int_m, &q, &x, f1r, d1r, &f2r, &d2r);
  return 0;
}

int msm1_wrap(double m, double q, double x, double *f1r, double *d1r)
{
  int int_m, kf=2, kc=1;
  double f2r, d2r;

  if ((m < 1) || (m != floor(m)) || (q<0)) {
    *f1r = NAN;
    *d1r = NAN;
    sf_error("msm1", SF_ERROR_DOMAIN, NULL);
    return -1;
  }
  int_m = (int )m;
  F_FUNC(mtu12,MTU12)(&kf,&kc,&int_m, &q, &x, f1r, d1r, &f2r, &d2r);
  return 0;
}

int mcm2_wrap(double m, double q, double x, double *f2r, double *d2r)
{
  int int_m, kf=1, kc=2;
  double f1r, d1r;

  if ((m < 0) || (m != floor(m)) || (q<0)) {
    *f2r = NAN;
    *d2r = NAN;
    sf_error("mcm2", SF_ERROR_DOMAIN, NULL);
    return -1;
  }
  int_m = (int )m;
  F_FUNC(mtu12,MTU12)(&kf,&kc,&int_m, &q, &x, &f1r, &d1r, f2r, d2r);
  return 0;
}

int msm2_wrap(double m, double q, double x, double *f2r, double *d2r)
{
  int int_m, kf=2, kc=2;
  double f1r, d1r;

  if ((m < 1) || (m != floor(m)) || (q<0)) {
    *f2r = NAN;
    *d2r = NAN;
    sf_error("msm2", SF_ERROR_DOMAIN, NULL);
    return -1;
  }
  int_m = (int )m;
  F_FUNC(mtu12,MTU12)(&kf,&kc,&int_m, &q, &x, &f1r, &d1r, f2r, d2r);
  return 0;
}


double pmv_wrap(double m, double v, double x){
  int int_m;
  double out;

  if (m != floor(m)) return NAN;
  int_m = (int ) m;
  F_FUNC(lpmv,LPMV)(&v, &int_m, &x, &out);
  CONVINF("pmv", out);
  return out;
}


/*
 * If x > 0 return w1f and w1d. Otherwise set x = abs(x) and return
 * w2f and -w2d.
*/
int pbwa_wrap(double a, double x, double *wf, double *wd) {
  int flag = 0;
  double w1f, w1d, w2f, w2d;

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
  F_FUNC(pbwa,PBWA)(&a, &x, &w1f, &w1d, &w2f, &w2d);
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
  F_FUNC(pbdv,PBDV)(&v, &x, dv, dp, pdf, pdd);
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
  F_FUNC(pbvv,PBVV)(&v, &x, vv, vp, pvf, pvd);
  PyMem_Free(vv);
  return 0;
}

double prolate_segv_wrap(double m, double n, double c)
{
  int kd=1;
  int int_m, int_n;
  double cv, *eg;

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
  F_FUNC(segv,SEGV)(&int_m,&int_n,&c,&kd,&cv,eg);
  PyMem_Free(eg);
  return cv;
}

double oblate_segv_wrap(double m, double n, double c)
{
  int kd=-1;
  int int_m, int_n;
  double cv, *eg;

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
  F_FUNC(segv,SEGV)(&int_m,&int_n,&c,&kd,&cv,eg);
  PyMem_Free(eg);
  return cv;
}


double prolate_aswfa_nocv_wrap(double m, double n, double c, double x, double *s1d)
{
  int kd = 1;
  int int_m, int_n;
  double cv, s1f, *eg;

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
  F_FUNC(segv,SEGV)(&int_m,&int_n,&c,&kd,&cv,eg);
  F_FUNC(aswfa,ASWFA)(&int_m,&int_n,&c,&x,&kd,&cv,&s1f,s1d);
  PyMem_Free(eg);
  return s1f;
}


double oblate_aswfa_nocv_wrap(double m, double n, double c, double x, double *s1d)
{
  int kd = -1;
  int int_m, int_n;
  double cv, s1f, *eg;

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
  F_FUNC(segv,SEGV)(&int_m,&int_n,&c,&kd,&cv,eg);
  F_FUNC(aswfa,ASWFA)(&int_m,&int_n,&c,&x,&kd,&cv,&s1f,s1d);
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
  F_FUNC(aswfa,ASWFA)(&int_m,&int_n,&c,&x,&kd,&cv,s1f,s1d);
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
  F_FUNC(aswfa,ASWFA)(&int_m,&int_n,&c,&x,&kd,&cv,s1f,s1d);
  return 0;
}


double prolate_radial1_nocv_wrap(double m, double n, double c, double x, double *r1d)
{
  int kf=1, kd=1;
  double r2f, r2d, r1f, cv, *eg;
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
  F_FUNC(segv,SEGV)(&int_m,&int_n,&c,&kd,&cv,eg);
  F_FUNC(rswfp,RSWFP)(&int_m,&int_n,&c,&x,&cv,&kf,&r1f,r1d,&r2f,&r2d);
  PyMem_Free(eg);
  return r1f;
}

double prolate_radial2_nocv_wrap(double m, double n, double c, double x, double *r2d)
{
  int kf=2, kd=1;
  double r1f, r1d, r2f, cv, *eg;
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
  F_FUNC(segv,SEGV)(&int_m,&int_n,&c,&kd,&cv,eg);
  F_FUNC(rswfp,RSWFP)(&int_m,&int_n,&c,&x,&cv,&kf,&r1f,&r1d,&r2f,r2d);
  PyMem_Free(eg);
  return r2f;
}

int prolate_radial1_wrap(double m, double n, double c, double cv, double x, double *r1f, double *r1d)
{
  int kf=1;
  double r2f, r2d;
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
  F_FUNC(rswfp,RSWFP)(&int_m,&int_n,&c,&x,&cv,&kf,r1f,r1d,&r2f,&r2d);
  return 0;
}

int prolate_radial2_wrap(double m, double n, double c, double cv, double x, double *r2f, double *r2d)
{
  int kf=2;
  double r1f, r1d;
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
  F_FUNC(rswfp,RSWFP)(&int_m,&int_n,&c,&x,&cv,&kf,&r1f,&r1d,r2f,r2d);
  return 0;
}

double oblate_radial1_nocv_wrap(double m, double n, double c, double x, double *r1d)
{
  int kf=1, kd=-1;
  double r2f, r2d, r1f, cv, *eg;
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
  F_FUNC(segv,SEGV)(&int_m,&int_n,&c,&kd,&cv,eg);
  F_FUNC(rswfo,RSWFO)(&int_m,&int_n,&c,&x,&cv,&kf,&r1f,r1d,&r2f,&r2d);
  PyMem_Free(eg);
  return r1f;
}

double oblate_radial2_nocv_wrap(double m, double n, double c, double x, double *r2d)
{
  int kf=2, kd=-1;
  double r1f, r1d, r2f, cv, *eg;
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
  F_FUNC(segv,SEGV)(&int_m,&int_n,&c,&kd,&cv,eg);
  F_FUNC(rswfo,RSWFO)(&int_m,&int_n,&c,&x,&cv,&kf,&r1f,&r1d,&r2f,r2d);
  PyMem_Free(eg);
  return r2f;
}

int oblate_radial1_wrap(double m, double n, double c, double cv, double x, double *r1f, double *r1d)
{
  int kf=1;
  double r2f, r2d;
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
  F_FUNC(rswfo,RSWFO)(&int_m,&int_n,&c,&x,&cv,&kf,r1f,r1d,&r2f,&r2d);
  return 0;
}

int oblate_radial2_wrap(double m, double n, double c, double cv, double x, double *r2f, double *r2d)
{
  int kf=2;
  double r1f, r1d;
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
  F_FUNC(rswfo,RSWFO)(&int_m,&int_n,&c,&x,&cv,&kf,&r1f,&r1d,r2f,r2d);
  return 0;
}


int modified_fresnel_plus_wrap(double x, npy_cdouble *Fplus, npy_cdouble *Kplus)
{
  int ks=0;
  double fm, fa, gm, ga;

  F_FUNC(ffk,FFK)(&ks,&x,F2C_CST(Fplus),&fm,&fa,F2C_CST(Kplus),&gm,&ga);
  return 0;
}

int modified_fresnel_minus_wrap(double x, npy_cdouble *Fminus, npy_cdouble *Kminus)
{
  int ks=1;
  double fm, fa, gm, ga;

  F_FUNC(ffk,FFK)(&ks,&x,F2C_CST(Fminus),&fm,&fa,F2C_CST(Kminus),&gm,&ga);
  return 0;
}
