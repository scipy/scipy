/* This file is a collection of wrappers around the
 *  Amos Fortran library of functions that take complex
 *  variables (see www.netlib.org) so that they can be called from
 *  the cephes library of corresponding name but work with complex
 *  arguments.
 */

#include "specfun_wrappers.h"

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
  if (x<0) {
    x = -x;
    flag = 1;
  }
  F_FUNC(stvhv,STVHV)(&v,&x,&out);
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
  F_FUNC(itsh0,ITSH0)(&x,&out);
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
}



