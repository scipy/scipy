/* This file is a collection of wrappers around the
 *  Special Function  Fortran library of functions 
 *  to be compiled with the other special functions in cephes
 *
 * Functions written by Shanjie Zhang and Jianming Jin.
 * Interface by
 *  Travis E. Oliphant
 */

#ifndef _SPEC_WRAPPERS_H
#define _SPEC_WRAPPERS_H
#include "Python.h"

#undef NAN
#undef INFINITY

extern double NAN;
extern double INFINITY;
extern double PI;

#define CADDR(z) (double *)(&((z).real)), (double*)(&((z).imag))
#define F2C_CST(z) (double *)&((z)->real), (double *)&((z)->imag)
#define REAL(z) (z).real
#define IMAG(z) (z).imag
#define ABSQ(z) (z).real*(z).real + (z).imag*(z).imag;
#define ZCONVINF(z) if (REAL((z))==1.0e300) REAL((z))=INFINITY; if (REAL((z))==-1.0e300) REAL((z))=-INFINITY
#define CONVINF(x) if ((x)==1.0e300) (x)=INFINITY; if ((x)==-1.0e300) (x)=-INFINITY

Py_complex cgamma_wrap( Py_complex z);
Py_complex clngamma_wrap( Py_complex z);
Py_complex cpsi_wrap( Py_complex z);
Py_complex crgamma_wrap( Py_complex z);
Py_complex chyp2f1_wrap( double a, double b, double c, Py_complex z);
Py_complex chyp1f1_wrap( double a, double b, Py_complex z);
double hypU_wrap(double a, double b, double x);
double exp1_wrap(double x);
double expi_wrap(double x);
Py_complex cexp1_wrap( Py_complex z);
Py_complex cerf_wrap( Py_complex z);
int itairy_wrap(double x, double *apt, double *bpt, double *ant, double *bnt);

double struve_wrap(double v, double x);
double itstruve0_wrap(double x);
double it2struve0_wrap(double x);

double modstruve_wrap(double v, double x);
double itmodstruve0_wrap(double x);

int kelvin_wrap(double x, Py_complex *Be, Py_complex *Ke, Py_complex *Bep, Py_complex *Kep);
#endif



  







