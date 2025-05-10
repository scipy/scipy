#ifndef _ARPACK_N_DOUBLE_COMPLEX_H
#define _ARPACK_N_DOUBLE_COMPLEX_H

#include "_arpack.h"

// BLAS Routines used
void zaxpy_(int* n, ARPACK_CPLX_TYPE* alpha, ARPACK_CPLX_TYPE* x, int* incx, ARPACK_CPLX_TYPE* y, int* incy);
void zcopy_(int* n, ARPACK_CPLX_TYPE* x, int* incx, ARPACK_CPLX_TYPE* y, int* incy);
ARPACK_CPLX_TYPE zdotc_(int* n, ARPACK_CPLX_TYPE* x, int* incx, ARPACK_CPLX_TYPE* y, int* incy);
void zgeru_(int* m, int* n, ARPACK_CPLX_TYPE* alpha, ARPACK_CPLX_TYPE* x, int* incx, ARPACK_CPLX_TYPE* y, int* incy, ARPACK_CPLX_TYPE* a, int* lda);
double dznrm2_(int* n, ARPACK_CPLX_TYPE* x, int* incx);
void zscal_(int* n, ARPACK_CPLX_TYPE* alpha, ARPACK_CPLX_TYPE* x, int* incx);
void zdscal_(int* n, double* da, ARPACK_CPLX_TYPE* zx, int* incx);
void zgemv_(char* trans, int* m, int* n, ARPACK_CPLX_TYPE* alpha, ARPACK_CPLX_TYPE* a, int* lda, ARPACK_CPLX_TYPE* x, int* incx, ARPACK_CPLX_TYPE* beta, ARPACK_CPLX_TYPE* y, int* incy);
void ztrmm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, ARPACK_CPLX_TYPE* alpha, ARPACK_CPLX_TYPE* a, int* lda, ARPACK_CPLX_TYPE* b, int* ldb);

// LAPACK Routines used
void zgeqr2_(int* m, int* n, ARPACK_CPLX_TYPE* a, int* lda, ARPACK_CPLX_TYPE* tau, ARPACK_CPLX_TYPE* work, int* info);
void zlacpy_(char* uplo, int* m, int* n, ARPACK_CPLX_TYPE* a, int* lda, ARPACK_CPLX_TYPE* b, int* ldb);
void zlahqr_(int* wantt, int* wantz, int* n, int* ilo, int* ihi, ARPACK_CPLX_TYPE* h, int* ldh, ARPACK_CPLX_TYPE* w, int* iloz, int* ihiz, ARPACK_CPLX_TYPE* z, int* ldz, int* info );
double zlanhs_(char* norm, int* n, ARPACK_CPLX_TYPE* a, int* lda, double* work);
void zlarf_(char* side, int* m, int* n, ARPACK_CPLX_TYPE* v, int* incv, ARPACK_CPLX_TYPE* tau, ARPACK_CPLX_TYPE* c, int* ldc, ARPACK_CPLX_TYPE* work);
void zlarfg_(int* n, ARPACK_CPLX_TYPE* alpha, ARPACK_CPLX_TYPE* x, int* incx, ARPACK_CPLX_TYPE* tau);
void zlartg_(ARPACK_CPLX_TYPE* f, ARPACK_CPLX_TYPE* g, double* c, ARPACK_CPLX_TYPE* s, ARPACK_CPLX_TYPE* r);
void zlascl_(char* mtype, int* kl, int* ku, double* cfrom, double* cto, int* m, int* n, ARPACK_CPLX_TYPE* a, int* lda, int* info);
void zlaset_(char* uplo, int* m, int* n, ARPACK_CPLX_TYPE* alpha, ARPACK_CPLX_TYPE* beta, ARPACK_CPLX_TYPE* a, int* lda);
void ztrevc_(char* side, char* howmny, int* select, int* n, ARPACK_CPLX_TYPE* t, int* ldt, ARPACK_CPLX_TYPE* vl, int* ldvl, ARPACK_CPLX_TYPE* vr, int* ldvr, int* mm, int* m, ARPACK_CPLX_TYPE* work, double* rwork, int* info);
void ztrsen_(char* job, char* compq, int* select, int* n, ARPACK_CPLX_TYPE* t, int* ldt, ARPACK_CPLX_TYPE* q, int* ldq, ARPACK_CPLX_TYPE* w, int* m, double* s, double* sep, ARPACK_CPLX_TYPE* work, int* lwork, int* info);
void zunm2r_(char* side, char* trans, int* m, int* n, int* k, ARPACK_CPLX_TYPE* a, int* lda, ARPACK_CPLX_TYPE* tau, ARPACK_CPLX_TYPE* c, int* ldc, ARPACK_CPLX_TYPE* work, int* info);

#if defined(_MSC_VER)
    // MSVC definitions
    #include <complex.h>  // MSVC C++ header
    typedef _Dcomplex ARPACK_CPLX_TYPE;
    #define ARPACK_cplx(real, imag) ((_Dcomplex){real, imag})

#else
    // C99 compliant compilers
    #include <complex.h>
    typedef double complex ARPACK_CPLX_TYPE;
    #define ARPACK_cplx(real, imag) ((real) + (imag)*I)

#endif

#endif
