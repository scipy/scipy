/* This file is a collection of wrappers around the
 *  Amos Fortran library of functions that take complex
 *  variables (see www.netlib.org) so that they can be called from
 *  the cephes library of corresponding name but work with complex
 *  arguments.
 */

#include "amos_wrappers.h"

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

extern int F_FUNC(zairy,ZAIRY)
     (double*, double*, int*, int*, double*, double*, int*, int*);
extern int F_FUNC(zbiry,ZBIRY)
     (double*, double*, int*, int*, double*, double*, int*, int*);
extern int F_FUNC(zbesi,ZBESI)
     (double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern int F_FUNC(zbesj,ZBESJ)
     (double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern int F_FUNC(zbesk,ZBESK)
     (double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern int F_FUNC(zbesy,ZBESY)
     (double*, double*, double*, int*, int*, double*, double*, int*, double*, double*, int*);
extern int F_FUNC(zbesh,ZBESH)
     (double*, double*, double*, int*, int*, int*, double*, double*, int*, int*);

/* This must be linked with fortran
 */

int ierr_to_mtherr( int nz, int ierr) {
     /* Return mtherr equivalents for ierr values */
  
  if (nz != 0) return UNDERFLOW;

  switch (ierr) {
  case 1:
    return DOMAIN;
  case 2:
    return OVERFLOW;
  case 3:
    return PLOSS;
  case 4:
    return TLOSS;
  case 5:   /* Algorithm termination condition not met */
    return TLOSS;    
  }
  return -1;
}

int cairy_wrap(Py_complex z, Py_complex *ai, Py_complex *aip, Py_complex *bi, Py_complex *bip) {
  int id = 0;
  int ierr = 0;
  int kode = 1;
  int nz;

  F_FUNC(zairy,ZAIRY)(CADDR(z), &id, &kode, F2C_CST(ai), &nz, &ierr);
  DO_MTHERR("airy:");
  F_FUNC(zbiry,ZBIRY)(CADDR(z), &id, &kode, F2C_CST(bi), &nz, &ierr);
  DO_MTHERR("airy:");
  
  id = 1;
  F_FUNC(zairy,ZAIRY)(CADDR(z), &id, &kode, F2C_CST(aip), &nz, &ierr);
  DO_MTHERR("airy:");
  F_FUNC(zbiry,ZBIRY)(CADDR(z), &id, &kode, F2C_CST(bip), &nz, &ierr);
  DO_MTHERR("airy:");
  return 0;
}

int cairy_wrap_e(Py_complex z, Py_complex *ai, Py_complex *aip, Py_complex *bi, Py_complex *bip) {
  int id = 0;
  int kode = 2;        /* Exponential scaling */
  int nz, ierr;

  F_FUNC(zairy,ZAIRY)(CADDR(z), &id, &kode, F2C_CST(ai), &nz, &ierr);
  DO_MTHERR("airye:");
  F_FUNC(zbiry,ZBIRY)(CADDR(z), &id, &kode, F2C_CST(bi), &nz, &ierr);
  DO_MTHERR("airye:");
  
  id = 1;
  F_FUNC(zairy,ZAIRY)(CADDR(z), &id, &kode, F2C_CST(aip), &nz, &ierr);
  DO_MTHERR("airye:");
  F_FUNC(zbiry,ZBIRY)(CADDR(z), &id, &kode, F2C_CST(bip), &nz, &ierr);
  DO_MTHERR("airye:");
  return 0;
}

Py_complex cbesi_wrap( double v, Py_complex z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  Py_complex cy;

  F_FUNC(zbesi,ZBESI)(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("iv:");
  return cy;
}

Py_complex cbesi_wrap_e( double v, Py_complex z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  Py_complex cy;

  F_FUNC(zbesi,ZBESI)(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("ive:");
  return cy;
}

  
Py_complex cbesj_wrap( double v, Py_complex z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  Py_complex cy;

  F_FUNC(zbesj,ZBESJ)(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("jv:");
  return cy;
}

Py_complex cbesj_wrap_e( double v, Py_complex z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  Py_complex cy;

  F_FUNC(zbesj,ZBESJ)(CADDR(z), &v, &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("jve:");
  return cy;
}

  
Py_complex cbesy_wrap( double v, Py_complex z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  Py_complex cy, cwork;

  F_FUNC(zbesy,ZBESY)(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, CADDR(cwork), &ierr);

  DO_MTHERR("yv:");
  return cy;
}

Py_complex cbesy_wrap_e( double v, Py_complex z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  Py_complex cy, cwork;

  F_FUNC(zbesy,ZBESY)(CADDR(z), &v, &kode, &n, CADDR(cy), &nz, CADDR(cwork), &ierr);
  DO_MTHERR("yve:");
  return cy;
}

  
Py_complex cbesk_wrap( double v, Py_complex z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  Py_complex cy;

  F_FUNC(zbesk,ZBESK)(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("kv:");
  return cy;
}

Py_complex cbesk_wrap_e( double v, Py_complex z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  Py_complex cy;

  F_FUNC(zbesk,ZBESK)(CADDR(z), &v, &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("kve:");
  return cy;
}
  
Py_complex cbesh_wrap1( double v, Py_complex z) {
  int n = 1;
  int kode = 1;
  int m = 1;
  int nz, ierr;
  Py_complex cy;

  F_FUNC(zbesh,ZBESH)(CADDR(z), &v,  &kode, &m, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("hankel1:");
  return cy;
}

Py_complex cbesh_wrap1_e( double v, Py_complex z) {
  int n = 1;
  int kode = 2;
  int m = 1;
  int nz, ierr;
  Py_complex cy;

  F_FUNC(zbesh,ZBESH)(CADDR(z), &v, &kode, &m, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("hankel1e:");
  return cy;
}
  
Py_complex cbesh_wrap2( double v, Py_complex z) {
  int n = 1;
  int kode = 1;
  int m = 2;
  int nz, ierr;
  Py_complex cy;

  F_FUNC(zbesh,ZBESH)(CADDR(z), &v,  &kode, &m, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("hankel2:");
  return cy;
}

Py_complex cbesh_wrap2_e( double v, Py_complex z) {
  int n = 1;
  int kode = 2;
  int m = 2;
  int nz, ierr;
  Py_complex cy;

  F_FUNC(zbesh,ZBESH)(CADDR(z), &v, &kode, &m, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("hankel2e:");
  return cy;
}
