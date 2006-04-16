/* This file is a collection of wrappers around the
 *  Amos Fortran library of functions that take complex
 *  variables (see www.netlib.org) so that they can be called from
 *  the cephes library of corresponding name but work with complex
 *  arguments.
 */

#ifndef _AMOS_WRAPPERS_H
#define _AMOS_WRAPPERS_H
#include "Python.h"
#include "cephes/mconf.h"

#define DO_MTHERR(name) if (nz !=0 || ierr !=0) mtherr(name, ierr_to_mtherr(nz, ierr))
int ierr_to_mtherr( int nz, int ierr); 
int cairy_wrap(Py_complex z, Py_complex *ai, Py_complex *aip, Py_complex *bi, Py_complex *bip);
int cairy_wrap_e(Py_complex z, Py_complex *ai, Py_complex *aip, Py_complex *bi, Py_complex *bip);
Py_complex cbesi_wrap( double v, Py_complex z);
Py_complex cbesi_wrap_e( double v, Py_complex z);
Py_complex cbesj_wrap( double v, Py_complex z);
Py_complex cbesj_wrap_e( double v, Py_complex z);
Py_complex cbesy_wrap( double v, Py_complex z);
Py_complex cbesy_wrap_e( double v, Py_complex z);
Py_complex cbesk_wrap( double v, Py_complex z);
Py_complex cbesk_wrap_e( double v, Py_complex z);  
Py_complex cbesh_wrap1( double v, Py_complex z);
Py_complex cbesh_wrap1_e( double v, Py_complex z);  
Py_complex cbesh_wrap2( double v, Py_complex z);
Py_complex cbesh_wrap2_e( double v, Py_complex z);
/* 
int cairy_(double *, int *, int *, doublecomplex *, int *, int *);
int cbiry_(doublecomplex *, int *, int *, doublecomplex *, int *, int *);
int cbesi_(doublecomplex *, double *, int *, int *, doublecomplex *, int *, int *);
int cbesj_(doublecomplex *, double *, int *, int *, doublecomplex *, int *, int *);
int cbesk_(doublecomplex *, double *, int *, int *, doublecomplex *, int *, int *);
int cbesy_(doublecomplex *, double *, int *, int *, doublecomplex *, int *, doublecomplex *, int *);
int cbesh_(doublecomplex *, double *, int *, int *, int *, doublecomplex *, int *, int *);
*/

#endif



  







