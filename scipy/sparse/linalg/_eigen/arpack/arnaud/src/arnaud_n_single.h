#ifndef _ARPACK_N_SINGLE_H
#define _ARPACK_N_SINGLE_H

#include "_arpack.h"

// BLAS Routines used
void saxpy_(int* n, float* alpha, float* x, int* incx, float* y, int* incy);
void scopy_(int* n, float* x, int* incx, float* y, int* incy);
float sdot_(int* n, float* x, int* incx, float* y, int* incy);
void sger_(int* m, int* n, float* alpha, float* x, int* incx, float* y, int* incy, float* a, int* lda);
float snrm2_(int* n, float* x, int* incx);
void sscal_(int* n, float* alpha, float* x, int* incx);
void sgemv_(char* trans, int* m, int* n, float* alpha, float* a, int* lda, float* x, int* incx, float* beta, float* y, int* incy);
void srot_(int* n, float* x, int* incx, float* y, int* incy, float* c, float* s);
void strmm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, float* alpha, float* a, int* lda, float* b, int* ldb);

// LAPACK Routines used
void sgeqr2_(int* m, int* n, float* a, int* lda, float* tau, float* work, int* info);
void slacpy_(char* uplo, int* m, int* n, float* a, int* lda, float* b, int* ldb);
void slahqr_(int* wantt, int* wantz, int* n, int* ilo, int* ihi, float* h, int* ldh, float* wr, float* wi, int* iloz, int* ihiz, float* z, int* ldz, int* info );
float slanhs_(char* norm, int* n, float* a, int* lda, float* work);
void slaset_(char* uplo, int* m, int* n, float* alpha, float* beta, float* a, int* lda);
void slarf_(char* side, int* m, int* n, float* v, int* incv, float* tau, float* c, int* ldc, float* work);
void slarfg_(int* n, float* alpha, float* x, int* incx, float* tau);
void slartg_(float* f, float* g, float* c, float* s, float* r);
void slartgp_(float* f, float* g, float* c, float* s, float* r);
void slascl_(char* mtype, int* kl, int* ku, float* cfrom, float* cto, int* m, int* n, float* a, int* lda, int* info);
void sorm2r_(char* side, char* trans, int* m, int* n, int* k, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* info);
void strevc_(char* side, char* howmny, int* select, int* n, float* t, int* ldt, float* vl, int* ldvl, float* vr, int* ldvr, int* mm, int* m, float* work, int* info);
void strsen_(char* job, char* compq, int* select, int* n, float* t, int* ldt, float* q, int* ldq, float* wr, float* wi, int* m, float* s, float* sep, float* work, int* lwork, int* iwork, int* liwork, int* info);

#endif
