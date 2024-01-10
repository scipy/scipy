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

float_complex BLAS_FUNC(cdotc)(int *n, float_complex *cx, int *incx, \
    float_complex *cy, int *incy);
float_complex F_FUNC(wcdotc, WCDOTC)(CBLAS_INT *n, float_complex *cx, \
        CBLAS_INT *incx, float_complex *cy, CBLAS_INT *incy){
    return BLAS_FUNC(cdotc)(n, cx, incx, cy, incy);
}

double_complex BLAS_FUNC(zdotc)(int *n, double_complex *zx, int *incx, \
    double_complex *zy, int *incy);
double_complex F_FUNC(wzdotc, WZDOTC)(CBLAS_INT *n, double_complex *zx, \
        CBLAS_INT *incx, double_complex *zy, CBLAS_INT *incy){
    return BLAS_FUNC(zdotc)(n, zx, incx, zy, incy);
}

float_complex BLAS_FUNC(cdotu)(int *n, float_complex *cx, int *incx, \
    float_complex *cy, int *incy);
float_complex F_FUNC(wcdotu, WCDOTU)(CBLAS_INT *n, float_complex *cx, \
        CBLAS_INT *incx, float_complex *cy, CBLAS_INT *incy){
    return BLAS_FUNC(cdotu)(n, cx, incx, cy, incy);
}

double_complex BLAS_FUNC(zdotu)(int *n, double_complex *zx, int *incx, \
    double_complex *zy, int *incy);
double_complex F_FUNC(wzdotu, WZDOTU)(CBLAS_INT *n, double_complex *zx, \
        CBLAS_INT *incx, double_complex *zy, CBLAS_INT *incy){
    return BLAS_FUNC(zdotu)(n, zx, incx, zy, incy);
}

float_complex BLAS_FUNC(cladiv)(float_complex *x, float_complex *y);
float_complex F_FUNC(wcladiv, WCLADIV)(float_complex *x, float_complex *y){
    return BLAS_FUNC(cladiv)(x, y);
}

double_complex BLAS_FUNC(zladiv)(double_complex *x, double_complex *y);
double_complex F_FUNC(wzladiv, WZLADIV)(double_complex *x, double_complex *y){
    return BLAS_FUNC(zladiv)(x, y);
}

#ifdef __cplusplus
}
#endif
