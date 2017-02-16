/* This file is a collection of wrappers around the
 *  Amos Fortran library of functions that take complex
 *  variables (see www.netlib.org) so that they can be called from
 *  the cephes library of corresponding name but work with complex
 *  arguments.
 */

#ifndef _AMOS_WRAPPERS_H
#define _AMOS_WRAPPERS_H
#include "Python.h"
#include "sf_error.h"

#include <numpy/npy_math.h>

#define DO_SFERR(name, varp)                          \
    do {                                              \
      if (nz !=0 || ierr != 0) {                      \
        sf_error(name, ierr_to_sferr(nz, ierr), NULL);\
        set_nan_if_no_computation_done(varp, ierr);   \
      }                                               \
    } while (0)

int ierr_to_sferr( int nz, int ierr);
void set_nan_if_no_computation_done(npy_cdouble *var, int ierr);
int airy_wrap(double x, double *ai, double *aip, double *bi, double *bip);
int cairy_wrap(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip);
int cairy_wrap_e(npy_cdouble z, npy_cdouble *ai, npy_cdouble *aip, npy_cdouble *bi, npy_cdouble *bip);
int cairy_wrap_e_real(double z, double *ai, double *aip, double *bi, double *bip);
npy_cdouble cbesi_wrap( double v, npy_cdouble z);
npy_cdouble cbesi_wrap_e( double v, npy_cdouble z);
double cbesi_wrap_e_real( double v, double z);
npy_cdouble cbesj_wrap( double v, npy_cdouble z);
npy_cdouble cbesj_wrap_e( double v, npy_cdouble z);
double cbesj_wrap_real(double v, double z);
double cbesj_wrap_e_real( double v, double z);
npy_cdouble cbesy_wrap( double v, npy_cdouble z);
double cbesy_wrap_real(double v, double x);
npy_cdouble cbesy_wrap_e( double v, npy_cdouble z);
double cbesy_wrap_e_real( double v, double z);
npy_cdouble cbesk_wrap( double v, npy_cdouble z);
npy_cdouble cbesk_wrap_e( double v, npy_cdouble z);  
double cbesk_wrap_real( double v, double z);
double cbesk_wrap_e_real( double v, double z);
double cbesk_wrap_real_int(int n, double z);
npy_cdouble cbesh_wrap1( double v, npy_cdouble z);
npy_cdouble cbesh_wrap1_e( double v, npy_cdouble z);  
npy_cdouble cbesh_wrap2( double v, npy_cdouble z);
npy_cdouble cbesh_wrap2_e( double v, npy_cdouble z);
double sin_pi(double x);
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



  







