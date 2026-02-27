#ifndef BLASLAPACK_DECLARATIONS_H
#define BLASLAPACK_DECLARATIONS_H
#include "npy_cblas.h"

void   BLAS_FUNC(daxpy)(CBLAS_INT* n, double* alpha, double* x, CBLAS_INT* incx, double* y, CBLAS_INT* incy);
void   BLAS_FUNC(dscal)(CBLAS_INT* n, double* alpha, double* x, CBLAS_INT* incx);
void   BLAS_FUNC(dcopy)(CBLAS_INT* n, double* x, CBLAS_INT* incx, double* y, CBLAS_INT* incy);
double BLAS_FUNC(dnrm2)(CBLAS_INT* n, double* x, CBLAS_INT* incx);
double BLAS_FUNC(ddot)(CBLAS_INT* n, double* x, CBLAS_INT* incx, double* y, CBLAS_INT* incy);
void   BLAS_FUNC(dgelsy)(CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* nrhs, double* a, CBLAS_INT* lda, double* b, CBLAS_INT* ldb, CBLAS_INT* jpvt, double* rcond, CBLAS_INT* rank, double* work, CBLAS_INT* lwork, CBLAS_INT* info);
void   BLAS_FUNC(dgemv)(char* trans, CBLAS_INT* m, CBLAS_INT* n, double* alpha, double* a, CBLAS_INT* lda, double* x, CBLAS_INT* incx, double* beta, double* y, CBLAS_INT* incy);
void   BLAS_FUNC(dgeqr2)(CBLAS_INT* m, CBLAS_INT* n, double* a, CBLAS_INT* lda, double* tau, double* work, CBLAS_INT* info);
void   BLAS_FUNC(dgeqrf)(CBLAS_INT* m, CBLAS_INT* n, double* a, CBLAS_INT* lda, double* tau, double* work, double* lwork, CBLAS_INT* info);
void   BLAS_FUNC(dgerq2)(CBLAS_INT* m, CBLAS_INT* n, double* a, CBLAS_INT* lda, double* tau, double* work, CBLAS_INT* info);
void   BLAS_FUNC(dlarf)(char* side, CBLAS_INT* m, CBLAS_INT* n, double* v, CBLAS_INT* incv, double* tau, double* c, CBLAS_INT* ldc, double* work);
void   BLAS_FUNC(dlarfgp)(CBLAS_INT* n, double* alpha, double* x, CBLAS_INT* incx, double* tau);
void   BLAS_FUNC(dlartgp)(double* f, double* g, double* cs, double* sn, double* r);
void   BLAS_FUNC(dorm2r)(char* side, char* trans, CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* k, double* a, CBLAS_INT* lda, double* tau, double* c, CBLAS_INT* ldc, double* work, CBLAS_INT* info);
void   BLAS_FUNC(dormr2)(char* side, char* trans, CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* k, double* a, CBLAS_INT* lda, double* tau, double* c, CBLAS_INT* ldc, double* work, CBLAS_INT* info);
void   BLAS_FUNC(dtpmv)(char* uplo, char* trans, char* diag, CBLAS_INT* n, double* ap, double* x, CBLAS_INT* incx);
void   BLAS_FUNC(dtpsv)(char* uplo, char* trans, char* diag, CBLAS_INT* n, double* ap, double* x, CBLAS_INT* incx);
void   BLAS_FUNC(dtrsm)(char* side, char* uplo, char* transa, char* diag, CBLAS_INT* m, CBLAS_INT* n, double* alpha, double* a, CBLAS_INT* lda, double* b, CBLAS_INT* ldb);
void   BLAS_FUNC(dtrsv)(char* uplo, char* trans, char* diag, CBLAS_INT* n, double* a, CBLAS_INT* lda, double* x, CBLAS_INT* incx);
void   BLAS_FUNC(dpotrf)(char* uplo, CBLAS_INT* n, double* a, CBLAS_INT* lda, CBLAS_INT* info);
void   BLAS_FUNC(dtrtrs)(char* uplo, char* trans, char* diag, CBLAS_INT* n, CBLAS_INT* nrhs, double* a, CBLAS_INT* lda, double* b, CBLAS_INT* ldb, CBLAS_INT* info);

#endif
