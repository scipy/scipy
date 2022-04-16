#ifndef FADDEEVA_H_
#define FADDEEVA_H_

#ifdef __cplusplus
#define EXTERN_C_START extern "C" {
#define EXTERN_C_END }
#else
#define EXTERN_C_START
#define EXTERN_C_END
#endif

#include <complex>
#include <Python.h>

#include "Faddeeva.hh"

EXTERN_C_START

#include <numpy/npy_math.h>

npy_cdouble faddeeva_w(npy_cdouble zp);
npy_cdouble faddeeva_erf(npy_cdouble zp);

double faddeeva_erfc(double x);
npy_cdouble faddeeva_erfc_complex(npy_cdouble zp);

double faddeeva_erfcx(double x);
npy_cdouble faddeeva_erfcx_complex(npy_cdouble zp);

double faddeeva_erfi(double zp);
npy_cdouble faddeeva_erfi_complex(npy_cdouble zp);

double faddeeva_dawsn(double zp);
npy_cdouble faddeeva_dawsn_complex(npy_cdouble zp);

npy_cdouble faddeeva_ndtr(npy_cdouble zp);

double faddeeva_log_ndtr(double x);
npy_cdouble faddeeva_log_ndtr_complex(npy_cdouble zp);

double faddeeva_voigt_profile(double x, double sigma, double gamma);

EXTERN_C_END

#endif
