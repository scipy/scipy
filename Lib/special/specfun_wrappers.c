/* This file is a collection of wrappers around the
 *  Specfun Fortran library of functions 
 */

#include "specfun_wrappers.h"
#include <stdio.h>

#define CADDR(z) (double *)(&((z).real)), (double*)(&((z).imag))
#define F2C_CST(z) (double *)&((z)->real), (double *)&((z)->imag)

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

extern void F_FUNC(cgama,CGAMA)(double*,double*,int*,double*,double*);
extern void F_FUNC(cpsi,CPSI)(double*,double*,double*,double*);
extern void F_FUNC(hygfz,HYGFZ)(double*,double*,double*,Py_complex*,Py_complex*);
extern void F_FUNC(cchg,CCHG)(double*,double*,Py_complex*,Py_complex*);
extern void F_FUNC(chgu,CHGU)(double*,double*,double*,double*,int*);
extern void F_FUNC(itairy,ITAIRY)(double*,double*,double*,double*,double*);
extern void F_FUNC(e1xb,E1XB)(double*,double*);
extern void F_FUNC(e1z,E1Z)(Py_complex*,Py_complex*);
extern void F_FUNC(eix,EIX)(double*,double*);
extern void F_FUNC(cerror,CERROR)(Py_complex*,Py_complex*);
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
extern void F_FUNC(cfc,CFC)(Py_complex*,Py_complex*,Py_complex*);
extern void F_FUNC(cfs,CFS)(Py_complex*,Py_complex*,Py_complex*);
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


double ber_wrap(double x)
{
  Py_complex Be, Ke, Bep, Kep;

  if (x<0) x=-x;
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF(Be);
  return REAL(Be);
}

double bei_wrap(double x)
{
  Py_complex Be, Ke, Bep, Kep;

  if (x<0) x=-x;
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF(Be);
  return IMAG(Be);
}

double ker_wrap(double x)
{
  Py_complex Be, Ke, Bep, Kep;

  if (x<0) return NAN;
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF(Ke);
  return REAL(Ke);  
}

double kei_wrap(double x)
{
  Py_complex Be, Ke, Bep, Kep;

  if (x<0) return NAN;
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF(Ke);
  return IMAG(Ke);  
}

double berp_wrap(double x)
{
  Py_complex Be, Ke, Bep, Kep;
  int flag = 0;

  if (x<0) {x=-x; flag=1;}
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF(Bep);
  if (flag) return -REAL(Bep);
  return REAL(Bep);
}

double beip_wrap(double x)
{
  Py_complex Be, Ke, Bep, Kep;
  int flag = 0;

  if (x<0) {x=-x; flag=1;}
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF(Bep);
  if (flag) return -IMAG(Bep);
  return IMAG(Bep);
}

double kerp_wrap(double x)
{
  Py_complex Be, Ke, Bep, Kep;

  if (x<0) return NAN;
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF(Kep);
  return REAL(Kep);  
}

double keip_wrap(double x)
{
  Py_complex Be, Ke, Bep, Kep;

  if (x<0) return NAN;
  F_FUNC(klvna,KLVNA)(&x, CADDR(Be), CADDR(Ke), CADDR(Bep), CADDR(Kep));
  ZCONVINF(Kep);
  return IMAG(Kep);  
}


int kelvin_wrap(double x, Py_complex *Be, Py_complex *Ke, Py_complex *Bep, Py_complex *Kep) {
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

int cfresnl_wrap(Py_complex z, Py_complex *zfs, Py_complex *zfc)
{
  Py_complex zfd;
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
  return 0;
}

int pbdv_wrap(double v, double x, double *pdf, double *pdd) {

  double *dv;
  double *dp;
  int num;

  num = ABS((int) v)+1;    
  dv = (double *)PyMem_Malloc(sizeof(double)*2*num);
  if (dv==NULL) {
    printf("Warning: Memory allocation error.\n"); 
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

  num = ABS((int) v)+1;    
  vv = (double *)PyMem_Malloc(sizeof(double)*2*num);
  if (vv==NULL) {
    printf("Warning: Memory allocation error.\n"); 
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
    printf("Warning: Memory allocation error.\n"); 
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
    printf("Warning: Memory allocation error.\n"); 
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
    *s1d = NAN;
    return NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    printf("Warning: Memory allocation error.\n"); 
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
    *s1d = NAN;
    return NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    printf("Warning: Memory allocation error.\n"); 
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
    *r1d = NAN;
    return NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    printf("Warning: Memory allocation error.\n"); 
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
    *r2d = NAN;
    return NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    printf("Warning: Memory allocation error.\n"); 
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
    *r1d = NAN;
    return NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    printf("Warning: Memory allocation error.\n"); 
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
    *r2d = NAN;
    return NAN;
  }
  int_m = (int )m;
  int_n = (int )n;
  eg = (double *)PyMem_Malloc(sizeof(double)*(n-m+2));
  if (eg==NULL) {
    printf("Warning: Memory allocation error.\n"); 
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
    *r2f = NAN;
    *r2d = NAN;
    return 0;
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
  
  F_FUNC(ffk,FFK)(&ks,&x,F2C_CST(Fplus),&fm,&fa,F2C_CST(Kplus),&gm,&ga);
  return 0;
}

int modified_fresnel_minus_wrap(double x, Py_complex *Fminus, Py_complex *Kminus)
{
  int ks=1;
  double fm, fa, gm, ga;
  
  F_FUNC(ffk,FFK)(&ks,&x,F2C_CST(Fminus),&fm,&fa,F2C_CST(Kminus),&gm,&ga);
  return 0;
}

