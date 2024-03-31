/* This file is a collection of wrappers around the
 *  Amos Fortran library of functions that take complex
 *  variables (see www.netlib.org) so that they can be called from
 *  the cephes library of corresponding name but work with complex
 *  arguments.
 */

#pragma once

#include "Python.h"
#include "sf_error.h"
#include "npy_2_complexcompat.h"
#include <numpy/npy_math.h>

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

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

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

  







