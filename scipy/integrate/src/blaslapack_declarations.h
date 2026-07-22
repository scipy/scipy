#ifndef BLASLAPACK_DECLARATIONS_H
#define BLASLAPACK_DECLARATIONS_H

#include "scipy_complex_support.h"
#include "scipy_blas_defines.h"

#if defined(_MSC_VER)
    // MSVC definition
    typedef _Dcomplex ZVODE_CPLX_TYPE;
    #define ZVODE_cplx(real, imag) ((_Dcomplex){real, imag})

#else
    // C99 compliant compilers
    typedef double complex ZVODE_CPLX_TYPE;
    #define ZVODE_cplx(real, imag) CMPLX(real, imag)
#endif


void BLAS_FUNC(daxpy)(CBLAS_INT* n, double* a, double* x, CBLAS_INT* incx, double* y, CBLAS_INT* incy);
void BLAS_FUNC(dcopy)(CBLAS_INT* n, double* x, CBLAS_INT* incx, double* y, CBLAS_INT* incy);
void BLAS_FUNC(dscal)(CBLAS_INT* n, double* a, double* x, CBLAS_INT* incx);
void BLAS_FUNC(dgbtrf)(CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* kl, CBLAS_INT* ku, double* ab, CBLAS_INT* ldab, CBLAS_INT* ipiv, CBLAS_INT* info);
void BLAS_FUNC(dgbtrs)(char* trans, CBLAS_INT* n, CBLAS_INT* kl, CBLAS_INT* ku, CBLAS_INT* nrhs, double* ab, CBLAS_INT* ldab, CBLAS_INT* ipiv, double* b, CBLAS_INT* ldb, CBLAS_INT* info);
void BLAS_FUNC(dgetrf)(CBLAS_INT* m, CBLAS_INT* n, double* a, CBLAS_INT* lda, CBLAS_INT* ipiv, CBLAS_INT* info);
void BLAS_FUNC(dgetrs)(char* trans, CBLAS_INT* n, CBLAS_INT* nrhs, double* a, CBLAS_INT* lda, CBLAS_INT* ipiv, double* b, CBLAS_INT* ldb, CBLAS_INT* info);

// For ZVODE complex routines
void BLAS_FUNC(zaxpy)(CBLAS_INT* n, ZVODE_CPLX_TYPE* a, ZVODE_CPLX_TYPE* x, CBLAS_INT* incx, ZVODE_CPLX_TYPE* y, CBLAS_INT* incy);
void BLAS_FUNC(zcopy)(CBLAS_INT* n, ZVODE_CPLX_TYPE* x, CBLAS_INT* incx, ZVODE_CPLX_TYPE* y, CBLAS_INT* incy);
void BLAS_FUNC(zscal)(CBLAS_INT* n, ZVODE_CPLX_TYPE* a, ZVODE_CPLX_TYPE* x, CBLAS_INT* incx);
void BLAS_FUNC(zgbtrf)(CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* kl, CBLAS_INT* ku, ZVODE_CPLX_TYPE* ab, CBLAS_INT* ldab, CBLAS_INT* ipiv, CBLAS_INT* info);
void BLAS_FUNC(zgbtrs)(char* trans, CBLAS_INT* n, CBLAS_INT* kl, CBLAS_INT* ku, CBLAS_INT* nrhs, ZVODE_CPLX_TYPE* ab, CBLAS_INT* ldab, CBLAS_INT* ipiv, ZVODE_CPLX_TYPE* b, CBLAS_INT* ldb, CBLAS_INT* info);
void BLAS_FUNC(zgetrf)(CBLAS_INT* m, CBLAS_INT* n, ZVODE_CPLX_TYPE* a, CBLAS_INT* lda, CBLAS_INT* ipiv, CBLAS_INT* info);
void BLAS_FUNC(zgetrs)(char* trans, CBLAS_INT* n, CBLAS_INT* nrhs, ZVODE_CPLX_TYPE* a, CBLAS_INT* lda, CBLAS_INT* ipiv, ZVODE_CPLX_TYPE* b, CBLAS_INT* ldb, CBLAS_INT* info);

#endif
