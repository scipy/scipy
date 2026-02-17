#ifndef BLASLAPACK_DECLARATIONS_H
#define BLASLAPACK_DECLARATIONS_H

void   daxpy_(int* n, double* alpha, double* x, int* incx, double* y, int* incy);
void   dscal_(int* n, double* alpha, double* x, int* incx);
void   dcopy_(int* n, double* x, int* incx, double* y, int* incy);
double dnrm2_(int* n, double* x, int* incx);
double ddot_(int* n, double* x, int* incx, double* y, int* incy);
void   dgelsy_(int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, int* jpvt, double* rcond, int* rank, double* work, int* lwork, int* info);
void   dgemv_(char* trans, int* m, int* n, double* alpha, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
void   dgeqr2_(int* m, int* n, double* a, int* lda, double* tau, double* work, int* info);
void   dgeqrf_(int* m, int* n, double* a, int* lda, double* tau, double* work, double* lwork, int* info);
void   dgerq2_(int* m, int* n, double* a, int* lda, double* tau, double* work, int* info);
void   dlarf_(char* side, int* m, int* n, double* v, int* incv, double* tau, double* c, int* ldc, double* work);
void   dlarfgp_(int* n, double* alpha, double* x, int* incx, double* tau);
void   dlartgp_(double* f, double* g, double* cs, double* sn, double* r);
void   dorm2r_(char* side, char* trans, int* m, int* n, int* k, double* a, int* lda, double* tau, double* c, int* ldc, double* work, int* info);
void   dormr2_(char* side, char* trans, int* m, int* n, int* k, double* a, int* lda, double* tau, double* c, int* ldc, double* work, int* info);
void   dtpmv_(char* uplo, char* trans, char* diag, int* n, double* ap, double* x, int* incx);
void   dtpsv_(char* uplo, char* trans, char* diag, int* n, double* ap, double* x, int* incx);
void   dtrsm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha, double* a, int* lda, double* b, int* ldb);
void   dtrsv_(char* uplo, char* trans, char* diag, int* n, double* a, int* lda, double* x, int* incx);
void   dpotrf_(char* uplo, int* n, double* a, int* lda, int* info);
void   dtrtrs_(char* uplo, char* trans, char* diag, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, int* info);

#endif
