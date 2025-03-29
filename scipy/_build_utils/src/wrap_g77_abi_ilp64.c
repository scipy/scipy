/*
 * This file is an analog of wrap_g77_abi.c, only for the ILP64 BLAS build.
 * See `wrap_g77_abi.c` for the purpose and scope.
 *
 * The only difference between this and wrap_g77_abi.c is that here we
 * export the wrappers mangled with BLAS_FUNC (i.e. with BLAS_PREFIX) instead
 * of mangling with F_FUNC from fortran_defs.h.
 * This mirrors the usage in scipy/linalg/fblas_64.pyf.src
 *
 */
#include "npy_cblas.h"
#define F_FUNC(f, F) BLAS_FUNC(f)

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
