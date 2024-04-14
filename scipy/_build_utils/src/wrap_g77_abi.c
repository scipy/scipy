/*
Some linear algebra libraries are built with a Fortran ABI that is incompatible
with the compiler used to build scipy (see gh-11812). This results in segfaults
when calling functions with complex-valued return types.

The wrappers in THIS FILE ensure compatibility by calling either:
1. The ABI-independent CBLAS API (cdotc, cdotu, zdotc, zdotu)
2. Fortran functions without complex-valued args (cladiv, zladiv)

When these wrappers are not necessary, wrap_g77_dummy_abi.c provides wrappers
that have the same name ('w' prefix) and calling convention as those here.

The choice of which wrapper file to compile with is handled at build time by
Meson (g77_abi_wrappers in scipy/meson.build).

On x86 machines, segfaults occur when Cython/F2PY-generated C code calls the
'w'-prefixed wrappers because the wrappers return C99 complex types while
Cython/F2PY use struct complex types (`{float r, i;}`).

Cython/F2PY code should instead call the 'wrp'-suffixed wrappers in this file,
passing a pointer to a variable in which to store the computed result. Unlike
return values, struct complex arguments work without segfaulting.
*/

#include "fortran_defs.h"
#include "npy_cblas.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
/* MSVC uses non-standard names for complex types. Intel compilers do not. */
#if defined(_MSC_VER) && !defined(__INTEL_COMPILER) && !defined(__INTEL_LLVM_COMPILER)
typedef _Dcomplex double_complex;
typedef _Fcomplex float_complex;
#else /* !defined(_MSC_VER) || defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)*/
typedef double _Complex double_complex;
typedef float _Complex float_complex;
#endif

float_complex F_FUNC(wcdotc,WCDOTC)(CBLAS_INT *n, float_complex *cx, \
        CBLAS_INT *incx, float_complex *cy, CBLAS_INT *incy){
    float_complex ret;
    /* Prototype in npy_cblas_base.h which is included in npy_cblas.h */
    CBLAS_FUNC(cblas_cdotc_sub)(*n, cx, *incx, cy, *incy,&ret);
    return ret;
}

double_complex F_FUNC(wzdotc,WZDOTC)(CBLAS_INT *n, double_complex *zx, \
        CBLAS_INT *incx, double_complex *zy, CBLAS_INT *incy){
    double_complex ret;
    /* Prototype in npy_cblas_base.h which is included in npy_cblas.h */
    CBLAS_FUNC(cblas_zdotc_sub)(*n, zx, *incx, zy, *incy,&ret);
    return ret;
}

float_complex F_FUNC(wcdotu,WCDOTU)(CBLAS_INT *n, float_complex *cx, \
        CBLAS_INT *incx, float_complex *cy, CBLAS_INT *incy){
    float_complex ret;
    /* Prototype in npy_cblas_base.h which is included in npy_cblas.h */
    CBLAS_FUNC(cblas_cdotu_sub)(*n, cx, *incx, cy, *incy,&ret);
    return ret;
}

double_complex F_FUNC(wzdotu,WZDOTU)(CBLAS_INT *n, double_complex *zx, \
        CBLAS_INT *incx, double_complex *zy, CBLAS_INT *incy){
    double_complex ret;
    /* Prototype in npy_cblas_base.h which is included in npy_cblas.h */
    CBLAS_FUNC(cblas_zdotu_sub)(*n, zx, *incx, zy, *incy,&ret);
    return ret;
}

void BLAS_FUNC(sladiv)(float *xr, float *xi, float *yr, float *yi, \
    float *retr, float *reti);
float_complex F_FUNC(wcladiv,WCLADIV)(float_complex *x, float_complex *y){
    float_complex ret;
    /* float_complex has the same memory layout as float[2], so we can
       cast and use pointer arithmetic to get the real and imaginary parts */
    BLAS_FUNC(sladiv)((float*)(x), (float*)(x)+1, \
        (float*)(y), (float*)(y)+1, \
        (float*)(&ret), (float*)(&ret)+1);
    return ret;
}

void BLAS_FUNC(dladiv)(double *xr, double *xi, double *yr, double *yi, \
    double *retr, double *reti);
double_complex F_FUNC(wzladiv,WZLADIV)(double_complex *x, double_complex *y){
    double_complex ret;
    /* double_complex has the same memory layout as double[2], so we can
       cast and use pointer arithmetic to get the real and imaginary parts */
    BLAS_FUNC(dladiv)((double*)(x), (double*)(x)+1, \
        (double*)(y), (double*)(y)+1, \
        (double*)(&ret), (double*)(&ret)+1);
    return ret;
}

void F_FUNC(cdotcwrp,WCDOTCWRP)(float_complex *ret, CBLAS_INT *n, float_complex *cx, \
        CBLAS_INT *incx, float_complex *cy, CBLAS_INT *incy){
    *ret = F_FUNC(wcdotc,WCDOTC)(n, cx, incx, cy, incy);
}

void F_FUNC(zdotcwrp,WZDOTCWRP)(double_complex *ret, CBLAS_INT *n, double_complex *zx, \
        CBLAS_INT *incx, double_complex *zy, CBLAS_INT *incy){
    *ret = F_FUNC(wzdotc,WZDOTC)(n, zx, incx, zy, incy);
}

void F_FUNC(cdotuwrp,CDOTUWRP)(float_complex *ret, CBLAS_INT *n, float_complex *cx, \
        CBLAS_INT *incx, float_complex *cy, CBLAS_INT *incy){
    *ret = F_FUNC(wcdotu,WCDOTU)(n, cx, incx, cy, incy);
}

void F_FUNC(zdotuwrp,ZDOTUWRP)(double_complex *ret, CBLAS_INT *n, double_complex *zx, \
        CBLAS_INT *incx, double_complex *zy, CBLAS_INT *incy){
    *ret = F_FUNC(wzdotu,WZDOTU)(n, zx, incx, zy, incy);
}

void F_FUNC(cladivwrp,CLADIVWRP)(float_complex *ret, float_complex *x, float_complex *y){
    *ret = F_FUNC(wcladiv,WCLADIV)(x, y);
}

void F_FUNC(zladivwrp,ZLADIVWRP)(double_complex *ret, double_complex *x, double_complex *y){
    *ret = F_FUNC(wzladiv,WZLADIV)(x, y);
}

#ifdef __cplusplus
}
#endif
