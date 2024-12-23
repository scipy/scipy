#ifndef FADDEEVA_H_
#define FADDEEVA_H_

#ifdef __cplusplus
#define EXTERN_C_START extern "C" {
#define EXTERN_C_END }
#else
#define EXTERN_C_START
#define EXTERN_C_END
#endif

#include <Python.h>
#include <complex>

#include "Faddeeva.hh"

EXTERN_C_START

#include <numpy/npy_math.h>

double faddeeva_log_ndtr(double x);
npy_cdouble faddeeva_log_ndtr_complex(npy_cdouble zp);

EXTERN_C_END

#endif
