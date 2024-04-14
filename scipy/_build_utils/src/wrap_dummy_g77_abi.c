/*
Some linear algebra libraries are built with a Fortran ABI that is incompatible
with the compiler used to build scipy (see gh-11812). This results in segfaults
when calling functions with complex-valued return types.

The wrappers in wrap_g77_abi.c ensure compatibility by calling either:
1. The ABI-independent CBLAS API (cdotc, cdotu, zdotc, zdotu)
2. Fortran functions without complex-valued args/return type (cladiv, zladiv)

When these wrappers are not necessary, THIS FILE provides wrappers that have
the same name ('w' prefix) and calling convention as those in wrap_g77_abi.c.

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

float_complex BLAS_FUNC(cdotc)(CBLAS_INT *n, float_complex *cx, CBLAS_INT *incx, \
    float_complex *cy, CBLAS_INT *incy);
float_complex F_FUNC(wcdotc,WCDOTC)(CBLAS_INT *n, float_complex *cx, \
        CBLAS_INT *incx, float_complex *cy, CBLAS_INT *incy){
    return BLAS_FUNC(cdotc)(n, cx, incx, cy, incy);
}

double_complex BLAS_FUNC(zdotc)(CBLAS_INT *n, double_complex *zx, CBLAS_INT *incx, \
    double_complex *zy, CBLAS_INT *incy);
double_complex F_FUNC(wzdotc,WZDOTC)(CBLAS_INT *n, double_complex *zx, \
        CBLAS_INT *incx, double_complex *zy, CBLAS_INT *incy){
    return BLAS_FUNC(zdotc)(n, zx, incx, zy, incy);
}

float_complex BLAS_FUNC(cdotu)(CBLAS_INT *n, float_complex *cx, CBLAS_INT *incx, \
    float_complex *cy, CBLAS_INT *incy);
float_complex F_FUNC(wcdotu,WCDOTU)(CBLAS_INT *n, float_complex *cx, \
        CBLAS_INT *incx, float_complex *cy, CBLAS_INT *incy){
    return BLAS_FUNC(cdotu)(n, cx, incx, cy, incy);
}

double_complex BLAS_FUNC(zdotu)(CBLAS_INT *n, double_complex *zx, CBLAS_INT *incx, \
    double_complex *zy, CBLAS_INT *incy);
double_complex F_FUNC(wzdotu,WZDOTU)(CBLAS_INT *n, double_complex *zx, \
        CBLAS_INT *incx, double_complex *zy, CBLAS_INT *incy){
    return BLAS_FUNC(zdotu)(n, zx, incx, zy, incy);
}

float_complex BLAS_FUNC(cladiv)(float_complex *x, float_complex *y);
float_complex F_FUNC(wcladiv,WCLADIV)(float_complex *x, float_complex *y){
    return BLAS_FUNC(cladiv)(x, y);
}

double_complex BLAS_FUNC(zladiv)(double_complex *x, double_complex *y);
double_complex F_FUNC(wzladiv,WZLADIV)(double_complex *x, double_complex *y){
    return BLAS_FUNC(zladiv)(x, y);
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
