/* This file is a collection of wrappers around the
 *  Amos Fortran library of functions that take complex
 *  variables (see www.netlib.org) so that they can be called from
 *  the cephes library of corresponding name but work with complex
 *  arguments.
 */

#include "specfun_wrappers.h"
#include <stdio.h>

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

extern double psi(double);
extern double struve(double, double);

/* This must be linked with fortran
 */

Py_complex cgamma_wrap( Py_complex z) {
  int kf = 1;
  Py_complex cy;

  F_FUNC(cgama,CGAMA)(CADDR(z), &kf, CADDR(cy));
  return cy;
}

Py_complex clngamma_wrap( Py_complex z) {
  int kf = 0;
  Py_complex cy;

  F_FUNC(cgama,CGAMA)(CADDR(z), &kf, CADDR(cy));
  return cy;
}

Py_complex cpsi_wrap( Py_complex z) {
  Py_complex cy;
  
  if (IMAG(z)==0.0) {
    REAL(cy) = psi(REAL(z));
    IMAG(cy) = 0.0;
  }
  else {
    F_FUNC(cpsi,CPSI)(CADDR(z), CADDR(cy));
  }
  return cy;
}

Py_complex crgamma_wrap( Py_complex z) {
  int kf = 1;
  Py_complex cy;
  Py_complex cy2;
  double magsq;

  F_FUNC(cgama,CGAMA)(CADDR(z), &kf, CADDR(cy));
  magsq = ABSQ(cy);
  REAL(cy2) = REAL(cy) / magsq;
  IMAG(cy2) = -IMAG(cy) / magsq;
  return cy2;
}

Py_complex chyp2f1_wrap( double a, double b, double c, Py_complex z) {
  Py_complex outz;
  int l1, l0;
 
 
  l0 = ((c == floor(c)) && (c < 0));
  l1 = ((fabs(1-REAL(z)) < 1e-15) && (IMAG(z) == 0) && (c-a-b <= 0));
  if (l0 || l1) {
    REAL(outz) = INFINITY;
    IMAG(outz) = 0.0;
    return outz;
  }
  F_FUNC(hygfz, HYGFZ)(&a, &b, &c, &z, &outz);
  return outz;
}

Py_complex chyp1f1_wrap(double a, double b, Py_complex z) {
  Py_complex outz;

  F_FUNC(cchg,CCHG)(&a, &b, &z, &outz);
  if (REAL(outz) == 1e300) {
    REAL(outz) = INFINITY;
  }
  return outz;
}


double hypU_wrap(double a, double b, double x) {
  double out;
  int md; /* method code --- not returned */

  F_FUNC(chgu,CHGU)(&a, &b, &x, &out, &md);
  if (out == 1e300) out = INFINITY;
  return out;
  
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
  CONVINF(out);
  return out;
}

Py_complex cexp1_wrap(Py_complex z) {
  Py_complex outz;
  
  F_FUNC(e1z,E1Z)(&z, &outz);
  ZCONVINF(outz);
  return outz;
}

double expi_wrap(double x) {
  double out;
  
  F_FUNC(eix,EIX)(&x, &out);
  CONVINF(out);
  return out;
}

Py_complex cerf_wrap(Py_complex z) {
  Py_complex outz;
  
  F_FUNC(cerror,CERROR)(&z, &outz);
  return outz;
}

double struve_wrap(double v, double x) {
  double out;
  int flag=0;

  if ((v<-8.0) || (v>12.5)) {
    return struve(v, x);  /* from cephes */
  }
  if (v==0.0) {
    if (x < 0) {x = -x; flag=1;}
    F_FUNC(stvh0,STVH0)(&x,&out);
    CONVINF(out);
    if (flag) out = -out;
    return out;
  }
  if (v==1.0) {
    if (x < 0) x=-x;
    F_FUNC(stvh1,STVH1)(&x,&out);
    CONVINF(out);
    return out;
  }
  F_FUNC(stvhv,STVHV)(&v,&x,&out);
  CONVINF(out);
  return out;  
}

double modstruve_wrap(double v, double x) {
  double out;
  int flag=0;

  if ((x < 0) & (floor(v)!=v)) return NAN;
  if (v==0.0) {
    if (x < 0) {x = -x; flag=1;}
    F_FUNC(stvl0,STVl0)(&x,&out);
    CONVINF(out);
    if (flag) out = -out;
    return out;
  }
  if (v==1.0) {
    if (x < 0) x=-x;
    F_FUNC(stvl1,STVl1)(&x,&out);
    CONVINF(out);
    return out;
  }
  if (x<0) {
    x = -x;
    flag = 1;
  }
  F_FUNC(stvlv,STVlV)(&v,&x,&out);
  CONVINF(out);
  if (flag && (!((int)floor(v) % 2))) out = -out;
  return out;  
}

double itstruve0_wrap(double x) {
  double out;

  if (x<0) x=-x;
  F_FUNC(itsh0,ITSH0)(&x,&out);
  CONVINF(out);
  return out;
}

double it2struve0_wrap(double x) {
  double out;
  int flag=0;
  
  if (x<0) {x=-x; flag=1;}
  F_FUNC(itth0,ITTH0)(&x,&out);
  CONVINF(out);
  if (flag) {
    out = PI - out;
  }
  return out;
}

double itmodstruve0_wrap(double x) {
  double out;

  if (x<0) x=-x;
  F_FUNC(itsl0,ITSL0)(&x,&out);
  CONVINF(out);
  return out;
}

int kelvin_wrap(double x, Py_complex *Be, Py_complex *Ke, Py_complex *Bep, Py_complex *Kep) {
  Py_complex outz;
  int flag = 0;
  
  if (x<0) {x=-x; flag=1;}
  F_FUNC(klvna,KLVNA)(&x, F2C_CST(Be), F2C_CST(Ke), F2C_CST(Bep), F2C_CST(Kep));
  ZCONVINF(*Be);
  ZCONVINF(*Ke);
  ZCONVINF(*Bep);
  ZCONVINF(*Kep);
  if (flag) {
    REAL(*Bep) = -REAL(*Bep);
    IMAG(*Bep) = -IMAG(*Bep);
    REAL(*Ke) = NAN;
    IMAG(*Ke) = NAN;
    REAL(*Kep) = NAN;
    IMAG(*Kep) = NAN;    
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

int cfresnl_wrap(Py_complex z, Py_complex *zfs, Py_complex *zfc)
{
  Py_complex zfd;
  printf("z = %f + j%f\n",z.real,z.imag);
  F_FUNC(cfs,CFS)(&z,zfs,&zfd);
  F_FUNC(cfc,CFC)(&z,zfc,&zfd); 
  return 0;
}

/* Mathieu functions */
/* Characteristic values */
double cem_cva_wrap(double m, double q) {
  int int_m, kd=1;
  double out;

  if ((m < 0) || (m != floor(m))) 
    return NAN;
  int_m = (int )m;
  if (int_m % 2) kd=2;
  F_FUNC(cva2,CVA2)(&kd, &int_m, &q, &out);
  return out;               
}

double sem_cva_wrap(double m, double q) {
  int int_m, kd=4;
  double out;

  if ((m < 1) || (m != floor(m))) 
    return NAN;
  int_m = (int )m;
  if (int_m % 2) kd=3;
  F_FUNC(cva2,CVA2)(&kd, &int_m, &q, &out);
  return out;               
}

/* Mathieu functions */
int cem_wrap(double m, double q, double x, double *csf, double *csd)
{
  int int_m, kf=1;
  if ((m < 1) || (m != floor(m)) || (q<0)) {
    *csf = NAN;
    *csd = NAN;
  }
  int_m = (int )m;
  F_FUNC(mtu0,MTU0)(&kf,&int_m, &q, &x, csf, csd);
  return 0;  
}

int sem_wrap(double m, double q, double x, double *csf, double *csd)
{
  int int_m, kf=2;
  if ((m < 1) || (m != floor(m)) || (q<0)) {
    *csf = NAN;
    *csd = NAN;
  }
  int_m = (int )m;
  F_FUNC(mtu0,MTU0)(&kf,&int_m, &q, &x, csf, csd);
  return 0;  
}


int mcm1_wrap(double m, double q, double x, double *f1r, double *d1r)
{
  int int_m, kf=1, kc=1;
  double f2r, d2r;

  if ((m < 1) || (m != floor(m)) || (q<0)) {
    *f1r = NAN;
    *d1r = NAN;
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
  }
  int_m = (int )m;
  F_FUNC(mtu12,MTU12)(&kf,&kc,&int_m, &q, &x, f1r, d1r, &f2r, &d2r);
  return 0;  
}

int mcm2_wrap(double m, double q, double x, double *f2r, double *d2r)
{
  int int_m, kf=1, kc=2;
  double f1r, d1r;

  if ((m < 1) || (m != floor(m)) || (q<0)) {
    *f2r = NAN;
    *d2r = NAN;
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
  }
  int_m = (int )m;
  F_FUNC(mtu12,MTU12)(&kf,&kc,&int_m, &q, &x, &f1r, &d1r, f2r, d2r);
  return 0;  
}

/* Stopped here --- below these need to be added to cephes */


double pmv_wrap(double m, double v, double x){
  int int_m;
  double out;

  if (m != floor(m)) return NAN;
  int_m = (int ) m;
  F_FUNC(lpmv,LPMV)(&v, &int_m, &x, &out);
  return out;
}

/* if x > 0 return w1f and w1d.
    otherwise return w2f and w2d (after abs(x))
*/
int pbwa_wrap(double a, double x, double *wf, double *wd) {
  int flag = 0;
  double w1f, w1d, w2f, w2d;
   
  if (x < 0) {x=-x; flag=1;}
  F_FUNC(pbwa,PBWA)(&a, &x, &w1f, &w1d, &w2f, &w2d);
  if (flag) {
    *wf = w2f;
    *wd = w2d;
  }
  else {
    *wf = w1f;
    *wd = w1d;
  }
}


int prolate_aswfa_wrap(double m, double n, double c, double x, double cv, double *s1f, double *s1d)
{
  int kd = 1;
  int int_m, int_n;

  if ((x >=1) || (x <=-1) || (m<0) || (n<m) || \
      (m!=floor(m)) || (n!=floor(n))) {
    *s1f = NAN;
    *s1d = NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  F_FUNC(aswfa,ASWFA)(&int_m,&int_n,&c,&x,&kd,&cv,s1f,s1d);
  return 0;
}

int oblate_aswfa_wrap(double m, double n, double c, double x, double cv, double *s1f, double *s1d)
{
  int kd = -1;
  int int_m, int_n;

  if ((x >=1) || (x <=-1) || (m<0) || (n<m) || \
      (m!=floor(m)) || (n!=floor(n))) {
    *s1f = NAN;
    *s1d = NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  F_FUNC(aswfa,ASWFA)(&int_m,&int_n,&c,&x,&kd,&cv,s1f,s1d);
  return 0;
}

int prolate_radial1_wrap(double m, double n, double c, double x, double cv, double *r1f, double *r1d)
{
  int kf=1;
  double r2f, r2d;
  int int_m, int_n;

  if ((x <=1.0) || (m<0) || (n<m) || \
     (m!=floor(m)) || (n!=floor(n))) {
    *r1f = NAN;
    *r1d = NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  F_FUNC(rswfp,RSWFP)(&int_m,&int_n,&c,&x,&cv,&kf,r1f,r1d,&r2f,&r2d);
  return 0;  
}

int prolate_radial2_wrap(double m, double n, double c, double x, double cv, double *r2f, double *r2d)
{
  int kf=2;
  double r1f, r1d;
  int int_m, int_n;

  if ((x <=1.0) || (m<0) || (n<m) || \
     (m!=floor(m)) || (n!=floor(n))) {
    *r2f = NAN;
    *r2d = NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  F_FUNC(rswfp,RSWFP)(&int_m,&int_n,&c,&x,&cv,&kf,&r1f,&r1d,r2f,r2d);
  return 0;  
}

int oblate_radial1_wrap(double m, double n, double c, double x, double cv, double *r1f, double *r1d)
{
  int kf=1;
  double r2f, r2d;
  int int_m, int_n;

  if ((x <0.0) || (m<0) || (n<m) || \
     (m!=floor(m)) || (n!=floor(n))) {
    *r1f = NAN;
    *r1d = NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  F_FUNC(rswfo,RSWFO)(&int_m,&int_n,&c,&x,&cv,&kf,r1f,r1d,&r2f,&r2d);
  return 0;  
}

int oblate_radial2_wrap(double m, double n, double c, double x, double cv, double *r2f, double *r2d)
{
  int kf=2;
  double r1f, r1d;
  int int_m, int_n;

  if ((x <0.0) || (m<0) || (n<m) || \
     (m!=floor(m)) || (n!=floor(n))) {
    *r2f = NAN;
    *r2d = NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  F_FUNC(rswfo,RSWFO)(&int_m,&int_n,&c,&x,&cv,&kf,&r1f,&r1d,r2f,r2d);
  return 0;  
}

int modified_fresnel_plus_wrap(double x, Py_complex *Fplus, Py_complex *Kplus)
{
  int ks=0;
  double fm, fa, gm, ga;
  
  F_FUNC(ffk,FFK)(&ks,&x,F2C_CST(Fplus),&fm,&fa,F2C_CST(Kplus),gm,ga);
  return 0;
}

int modified_fresnel_minus_wrap(double x, Py_complex *Fminus, Py_complex *Kminus)
{
  int ks=1;
  double fm, fa, gm, ga;
  
  F_FUNC(ffk,FFK)(&ks,&x,F2C_CST(Fminus),&fm,&fa,F2C_CST(Kminus),gm,ga);
  return 0;
}

