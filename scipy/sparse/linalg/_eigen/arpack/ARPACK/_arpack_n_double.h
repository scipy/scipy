#ifndef _ARPACK_N_DOUBLE_H
#define _ARPACK_N_DOUBLE_H

#include "_arpack.h"

// BLAS Routines used
void daxpy_(int* n, double* alpha, double* x, int* incx, double* y, int* incy);
void dcopy_(int* n, double* x, int* incx, double* y, int* incy);
double ddot_(int* n, double* x, int* incx, double* y, int* incy);
void dger_(int* m, int* n, double* alpha, double* x, int* incx, double* y, int* incy, double* a, int* lda);
double dnrm2_(int* n, double* x, int* incx);
void dscal_(int* n, double* alpha, double* x, int* incx);
void dgemv_(char* trans, int* m, int* n, double* alpha, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
void dtrmm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha, double* a, int* lda, double* b, int* ldb);

// LAPACK Routines used
void dgeqr2_(int* m, int* n, double* a, int* lda, double* tau, double* work, int* info);
void dlacpy_(char* uplo, int* m, int* n, double* a, int* lda, double* b, int* ldb);
void dlahqr_(int* wantt, int* wantz, int* n, int* ilo, int* ihi, double* h, int* ldh, double* wr, double* wi, int* iloz, int* ihiz, double* z, int* ldz, int* info );
double dlanhs_(char* norm, int* n, double* a, int* lda, double* work);
void dlaset_(char* uplo, int* m, int* n, double* alpha, double* beta, double* a, int* lda);
void dlarf_(char* side, int* m, int* n, double* v, int* incv, double* tau, double* c, int* ldc, double* work);
void dlarfg_(int* n, double* alpha, double* x, int* incx, double* tau);
void dlartg_(double* f, double* g, double* c, double* s, double* r);
void dlascl_(char* mtype, int* kl, int* ku, double* cfrom, double* cto, int* m, int* n, double* a, int* lda, int* info);
void dorm2r_(char* side, char* trans, int* m, int* n, int* k, double* a, int* lda, double* tau, double* c, int* ldc, double* work, int* info);
void dtrevc_(char* side, char* howmny, int* select, int* n, double* t, int* ldt, double* vl, int* ldvl, double* vr, int* ldvr, int* mm, int* m, double* work, int* info);
void dtrsen_(char* job, char* compq, int* select, int* n, double* t, int* ldt, double* q, int* ldq, double* wr, double* wi, int* m, double* s, double* sep, double* work, int* lwork, int* iwork, int* liwork, int* info);

#endif
