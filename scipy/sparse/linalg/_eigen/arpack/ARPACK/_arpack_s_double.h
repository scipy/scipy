#include "_arpack.h"

// BLAS Routines used
void daxpy_(int* n, double* alpha, double* x, int* incx, double* y, int* incy);
void dcopy_(int* n, double* x, int* incx, double* y, int* incy);
double ddot_(int* n, double* x, int* incx, double* y, int* incy);
void dgemv_(char* trans, int* m, int* n, double* alpha, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
void dger_(int* m, int* n, double* alpha, double* x, int* incx, double* y, int* incy, double* a, int* lda);
double dnrm2_(int* n, double* x, int* incx);
void dscal_(int* n, double* alpha, double* x, int* incx);
void dswap_(int* n, double* x, int* incx, double* y, int* incy);
void dtrmm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha, double* a, int* lda, double* b, int* ldb);


// LAPACK Routines used
void dgeqr2_(int* m, int* n, double* a, int* lda, double* tau, double* work, int* info);
void dlacpy_(char* uplo, int* m, int* n, double* a, int* lda, double* b, int* ldb);
void dlaev2_(double* a, double* b, double* c, double* rt1, double* rt2, double* cs1, double* sn1);
double dlanst_(double* norm, int* n, double* d, double* e);
void dlartg_(double* f, double* g, double* c, double* s, double* r);
void dlascl_(char* mtype, int* kl, int* ku, double* cfrom, double* cto, int* m, int* n, double* a, int* lda, int* info);
void dlasr_(char* side, char* pivot, char* direct, int* m, int* n, double* c, double* s, double* a, int* lda);
void dorm2r_(char* side, char* trans, int* m, int* n, int* k, double* a, int* lda, double* tau, double* c, int* ldc, double* work, int* info);
void dsteqr_(char* compz, int* n, double* d, double* e, double* z, int* ldz, double* work, int* info);