/* This file is a collection of wrappers around the
 *  Amos Fortran library of functions that take complex
 *  variables (see www.netlib.org) so that they can be called from
 *  the cephes library of corresponding name but work with complex
 *  arguments.
 */

#include "amos_wrappers.h"

/* This must be linked with g77
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

  zairy_(CADDR(z), &id, &kode, F2C_CST(ai), &nz, &ierr);
  DO_MTHERR("airy:");
  zbiry_(CADDR(z), &id, &kode, F2C_CST(bi), &nz, &ierr);
  DO_MTHERR("airy:");
  
  id = 1;
  zairy_(CADDR(z), &id, &kode, F2C_CST(aip), &nz, &ierr);
  DO_MTHERR("airy:");
  zbiry_(CADDR(z), &id, &kode, F2C_CST(bip), &nz, &ierr);
  DO_MTHERR("airy:");
  return 0;
}

int cairy_wrap_e(Py_complex z, Py_complex *ai, Py_complex *aip, Py_complex *bi, Py_complex *bip) {
  int id = 0;
  int kode = 2;        /* Exponential scaling */
  int nz, ierr;

  zairy_(CADDR(z), &id, &kode, F2C_CST(ai), &nz, &ierr);
  DO_MTHERR("airye:");
  zbiry_(CADDR(z), &id, &kode, F2C_CST(bi), &nz, &ierr);
  DO_MTHERR("airye:");
  
  id = 1;
  zairy_(CADDR(z), &id, &kode, F2C_CST(aip), &nz, &ierr);
  DO_MTHERR("airye:");
  zbiry_(CADDR(z), &id, &kode, F2C_CST(bip), &nz, &ierr);
  DO_MTHERR("airye:");
  return 0;
}

Py_complex cbesi_wrap( double v, Py_complex z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  Py_complex cy;

  zbesi_(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("iv:");
  return cy;
}

Py_complex cbesi_wrap_e( double v, Py_complex z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  Py_complex cy;

  zbesi_(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("ive:");
  return cy;
}

  
Py_complex cbesj_wrap( double v, Py_complex z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  Py_complex cy;

  zbesj_(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("jv:");
  return cy;
}

Py_complex cbesj_wrap_e( double v, Py_complex z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  Py_complex cy;

  zbesj_(CADDR(z), &v, &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("jve:");
  return cy;
}

  
Py_complex cbesy_wrap( double v, Py_complex z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  Py_complex cy, cwork;

  zbesy_(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, CADDR(cwork), &ierr);

  DO_MTHERR("yv:");
  return cy;
}

Py_complex cbesy_wrap_e( double v, Py_complex z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  Py_complex cy, cwork;

  zbesy_(CADDR(z), &v, &kode, &n, CADDR(cy), &nz, CADDR(cwork), &ierr);
  DO_MTHERR("yve:");
  return cy;
}

  
Py_complex cbesk_wrap( double v, Py_complex z) {
  int n = 1;
  int kode = 1;
  int nz, ierr;
  Py_complex cy;

  zbesk_(CADDR(z), &v,  &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("kv:");
  return cy;
}

Py_complex cbesk_wrap_e( double v, Py_complex z) {
  int n = 1;
  int kode = 2;
  int nz, ierr;
  Py_complex cy;

  zbesk_(CADDR(z), &v, &kode, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("kve:");
  return cy;
}
  
Py_complex cbesh_wrap1( double v, Py_complex z) {
  int n = 1;
  int kode = 1;
  int m = 1;
  int nz, ierr;
  Py_complex cy;

  zbesh_(CADDR(z), &v,  &kode, &m, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("hankel1:");
  return cy;
}

Py_complex cbesh_wrap1_e( double v, Py_complex z) {
  int n = 1;
  int kode = 2;
  int m = 1;
  int nz, ierr;
  Py_complex cy;

  zbesh_(CADDR(z), &v, &kode, &m, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("hankel1e:");
  return cy;
}
  
Py_complex cbesh_wrap2( double v, Py_complex z) {
  int n = 1;
  int kode = 1;
  int m = 2;
  int nz, ierr;
  Py_complex cy;

  zbesh_(CADDR(z), &v,  &kode, &m, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("hankel2:");
  return cy;
}

Py_complex cbesh_wrap2_e( double v, Py_complex z) {
  int n = 1;
  int kode = 2;
  int m = 2;
  int nz, ierr;
  Py_complex cy;

  zbesh_(CADDR(z), &v, &kode, &m, &n, CADDR(cy), &nz, &ierr);
  DO_MTHERR("hankel2e:");
  return cy;
}
