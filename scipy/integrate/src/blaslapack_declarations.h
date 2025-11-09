#ifndef BLASLAPACK_DECLARATIONS_H
#define BLASLAPACK_DECLARATIONS_H

#include <complex.h>
#if defined(_MSC_VER)
    // MSVC definition
    typedef _Dcomplex ZVODE_CPLX_TYPE;
    #define ZVODE_cplx(real, imag) ((_Dcomplex){real, imag})

#else
    // C99 compliant compilers
    typedef double complex ZVODE_CPLX_TYPE;
    #define ZVODE_cplx(real, imag) ((real) + (imag)*I)
#endif


void daxpy_(int* n, double* a, double* x, int* incx, double* y, int* incy);
void dcopy_(int* n, double* x, int* incx, double* y, int* incy);
void dscal_(int* n, double* a, double* x, int* incx);
void dgbtrf_(int* m, int* n, int* kl, int* ku, double* ab, int* ldab, int* ipiv, int* info);
void dgbtrs_(char* trans, int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info);
void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
void dgetrs_(char* trans, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);

// For ZVODE complex routines
void zaxpy_(int* n, ZVODE_CPLX_TYPE* a, ZVODE_CPLX_TYPE* x, int* incx, ZVODE_CPLX_TYPE* y, int* incy);
void zcopy_(int* n, ZVODE_CPLX_TYPE* x, int* incx, ZVODE_CPLX_TYPE* y, int* incy);
void zscal_(int* n, ZVODE_CPLX_TYPE* a, ZVODE_CPLX_TYPE* x, int* incx);
void zgbtrf_(int* m, int* n, int* kl, int* ku, ZVODE_CPLX_TYPE* ab, int* ldab, int* ipiv, int* info);
void zgbtrs_(char* trans, int* n, int* kl, int* ku, int* nrhs, ZVODE_CPLX_TYPE* ab, int* ldab, int* ipiv, ZVODE_CPLX_TYPE* b, int* ldb, int* info);
void zgetrf_(int* m, int* n, ZVODE_CPLX_TYPE* a, int* lda, int* ipiv, int* info);
void zgetrs_(char* trans, int* n, int* nrhs, ZVODE_CPLX_TYPE* a, int* lda, int* ipiv, ZVODE_CPLX_TYPE* b, int* ldb, int* info);

#endif
