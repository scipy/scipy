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

extern std::complex<double> Faddeeva_w(std::complex<double> z,double relerr=0);
extern double ImFaddeeva_w(double x); // special-case code for Im[w(x)]
extern double erfcx(double x);

EXTERN_C_START

#include <numpy/npy_math.h>

npy_cdouble faddeeva_w(npy_cdouble zp);
npy_cdouble faddeeva_erf(npy_cdouble zp);
npy_cdouble faddeeva_erfc(npy_cdouble zp);
double faddeeva_erfcx(double zp);
npy_cdouble faddeeva_erfcx_complex(npy_cdouble zp);

EXTERN_C_END

#endif
