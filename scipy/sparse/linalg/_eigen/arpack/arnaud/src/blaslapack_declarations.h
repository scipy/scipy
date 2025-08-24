#ifndef BLASLAPACK_DECLARATIONS_H
#define BLASLAPACK_DECLARATIONS_H

#include "arnaud/types.h"

// BLAS
void saxpy_(int* n, float* alpha, float* x, int* incx, float* y, int* incy);
void scopy_(int* n, float* x, int* incx, float* y, int* incy);
float sdot_(int* n, float* x, int* incx, float* y, int* incy);
void sgemv_(char* trans, int* m, int* n, float* alpha, float* a, int* lda, float* x, int* incx, float* beta, float* y, int* incy);
void sger_(int* m, int* n, float* alpha, float* x, int* incx, float* y, int* incy, float* a, int* lda);
float snrm2_(int* n, float* x, int* incx);
void srot_(int* n, float* sx, int* incx, float* sy, int* incy, float* c, float* s);
void sscal_(int* n, float* alpha, float* x, int* incx);
void sswap_(int* n, float* x, int* incx, float* y, int* incy);
void strmm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, float* alpha, float* a, int* lda, float* b, int* ldb);


void daxpy_(int* n, double* alpha, double* x, int* incx, double* y, int* incy);
void dcopy_(int* n, double* x, int* incx, double* y, int* incy);
double ddot_(int* n, double* x, int* incx, double* y, int* incy);
void dgemv_(char* trans, int* m, int* n, double* alpha, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
void dger_(int* m, int* n, double* alpha, double* x, int* incx, double* y, int* incy, double* a, int* lda);
double dnrm2_(int* n, double* x, int* incx);
void drot_(int* n, double* sx, int* incx, double* sy, int* incy, double* c, double* s);
void dscal_(int* n, double* alpha, double* x, int* incx);
void dswap_(int* n, double* x, int* incx, double* y, int* incy);
void dtrmm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha, double* a, int* lda, double* b, int* ldb);


void caxpy_(int* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* x, int* incx, ARNAUD_CPLXF_TYPE* y, int* incy);
void ccopy_(int* n, ARNAUD_CPLXF_TYPE* x, int* incx, ARNAUD_CPLXF_TYPE* y, int* incy);
void cgeru_(int* m, int* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* x, int* incx, ARNAUD_CPLXF_TYPE* y, int* incy, ARNAUD_CPLXF_TYPE* a, int* lda);
float scnrm2_(int* n, ARNAUD_CPLXF_TYPE* x, int* incx);
void cscal_(int* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* x, int* incx);
void csscal_(int* n, float* da, ARNAUD_CPLXF_TYPE* zx, int* incx);
void cgemv_(char* trans, int* m, int* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* a, int* lda, ARNAUD_CPLXF_TYPE* x, int* incx, ARNAUD_CPLXF_TYPE* beta, ARNAUD_CPLXF_TYPE* y, int* incy);
void crot_(int* n, ARNAUD_CPLXF_TYPE* cx, int* incx, ARNAUD_CPLXF_TYPE* cy, int* incy, float* c, ARNAUD_CPLXF_TYPE* s);
void ctrmm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* a, int* lda, ARNAUD_CPLXF_TYPE* b, int* ldb);


void zaxpy_(int* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* x, int* incx, ARNAUD_CPLX_TYPE* y, int* incy);
void zcopy_(int* n, ARNAUD_CPLX_TYPE* x, int* incx, ARNAUD_CPLX_TYPE* y, int* incy);
void zgeru_(int* m, int* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* x, int* incx, ARNAUD_CPLX_TYPE* y, int* incy, ARNAUD_CPLX_TYPE* a, int* lda);
double dznrm2_(int* n, ARNAUD_CPLX_TYPE* x, int* incx);
void zscal_(int* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* x, int* incx);
void zdscal_(int* n, double* da, ARNAUD_CPLX_TYPE* zx, int* incx);
void zgemv_(char* trans, int* m, int* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* a, int* lda, ARNAUD_CPLX_TYPE* x, int* incx, ARNAUD_CPLX_TYPE* beta, ARNAUD_CPLX_TYPE* y, int* incy);
void zrot_(int* n, ARNAUD_CPLX_TYPE* cx, int* incx, ARNAUD_CPLX_TYPE* cy, int* incy, double* c, ARNAUD_CPLX_TYPE* s);
void ztrmm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* a, int* lda, ARNAUD_CPLX_TYPE* b, int* ldb);



// LAPACK
void sgeqr2_(int* m, int* n, float* a, int* lda, float* tau, float* work, int* info);
void slacpy_(char* uplo, int* m, int* n, float* a, int* lda, float* b, int* ldb);
void slaev2_(float* a, float* b, float* c, float* rt1, float* rt2, float* cs1, float* sn1);
void slahqr_(int* wantt, int* wantz, int* n, int* ilo, int* ihi, float* h, int* ldh, float* wr, float* wi, int* iloz, int* ihiz, float* z, int* ldz, int* info );
float slanhs_(char* norm, int* n, float* a, int* lda, float* work);
float slanst_(char* norm, int* n, float* d, float* e);
void slarf_(char* side, int* m, int* n, float* v, int* incv, float* tau, float* c, int* ldc, float* work);
void slarfg_(int* n, float* alpha, float* x, int* incx, float* tau);
void slartg_(float* f, float* g, float* c, float* s, float* r);
void slartgp_(float* f, float* g, float* c, float* s, float* r);
void slascl_(char* mtype, int* kl, int* ku, float* cfrom, float* cto, int* m, int* n, float* a, int* lda, int* info);
void slaset_(char* uplo, int* m, int* n, float* alpha, float* beta, float* a, int* lda);
void slasr_(char* side, char* pivot, char* direct, int* m, int* n, float* c, float* s, float* a, int* lda);
void sorm2r_(char* side, char* trans, int* m, int* n, int* k, float* a, int* lda, float* tau, float* c, int* ldc, float* work, int* info);
void ssteqr_(char* compz, int* n, float* d, float* e, float* z, int* ldz, float* work, int* info);
void strevc_(char* side, char* howmny, int* select, int* n, float* t, int* ldt, float* vl, int* ldvl, float* vr, int* ldvr, int* mm, int* m, float* work, int* info);
void strsen_(char* job, char* compq, int* select, int* n, float* t, int* ldt, float* q, int* ldq, float* wr, float* wi, int* m, float* s, float* sep, float* work, int* lwork, int* iwork, int* liwork, int* info);


void dgeqr2_(int* m, int* n, double* a, int* lda, double* tau, double* work, int* info);
void dlacpy_(char* uplo, int* m, int* n, double* a, int* lda, double* b, int* ldb);
void dlaev2_(double* a, double* b, double* c, double* rt1, double* rt2, double* cs1, double* sn1);
void dlahqr_(int* wantt, int* wantz, int* n, int* ilo, int* ihi, double* h, int* ldh, double* wr, double* wi, int* iloz, int* ihiz, double* z, int* ldz, int* info );
double dlanhs_(char* norm, int* n, double* a, int* lda, double* work);
double dlanst_(char* norm, int* n, double* d, double* e);
void dlarf_(char* side, int* m, int* n, double* v, int* incv, double* tau, double* c, int* ldc, double* work);
void dlarfg_(int* n, double* alpha, double* x, int* incx, double* tau);
void dlartg_(double* f, double* g, double* c, double* s, double* r);
void dlartgp_(double* f, double* g, double* c, double* s, double* r);
void dlascl_(char* mtype, int* kl, int* ku, double* cfrom, double* cto, int* m, int* n, double* a, int* lda, int* info);
void dlaset_(char* uplo, int* m, int* n, double* alpha, double* beta, double* a, int* lda);
void dlasr_(char* side, char* pivot, char* direct, int* m, int* n, double* c, double* s, double* a, int* lda);
void dorm2r_(char* side, char* trans, int* m, int* n, int* k, double* a, int* lda, double* tau, double* c, int* ldc, double* work, int* info);
void dsteqr_(char* compz, int* n, double* d, double* e, double* z, int* ldz, double* work, int* info);
void dtrevc_(char* side, char* howmny, int* select, int* n, double* t, int* ldt, double* vl, int* ldvl, double* vr, int* ldvr, int* mm, int* m, double* work, int* info);
void dtrsen_(char* job, char* compq, int* select, int* n, double* t, int* ldt, double* q, int* ldq, double* wr, double* wi, int* m, double* s, double* sep, double* work, int* lwork, int* iwork, int* liwork, int* info);


void cgeqr2_(int* m, int* n, ARNAUD_CPLXF_TYPE* a, int* lda, ARNAUD_CPLXF_TYPE* tau, ARNAUD_CPLXF_TYPE* work, int* info);
void clacpy_(char* uplo, int* m, int* n, ARNAUD_CPLXF_TYPE* a, int* lda, ARNAUD_CPLXF_TYPE* b, int* ldb);
void clahqr_(int* wantt, int* wantz, int* n, int* ilo, int* ihi, ARNAUD_CPLXF_TYPE* h, int* ldh, ARNAUD_CPLXF_TYPE* w, int* iloz, int* ihiz, ARNAUD_CPLXF_TYPE* z, int* ldz, int* info );
float clanhs_(char* norm, int* n, ARNAUD_CPLXF_TYPE* a, int* lda, float* work);
void clarf_(char* side, int* m, int* n, ARNAUD_CPLXF_TYPE* v, int* incv, ARNAUD_CPLXF_TYPE* tau, ARNAUD_CPLXF_TYPE* c, int* ldc, ARNAUD_CPLXF_TYPE* work);
void clarfg_(int* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* x, int* incx, ARNAUD_CPLXF_TYPE* tau);
void clartg_(ARNAUD_CPLXF_TYPE* f, ARNAUD_CPLXF_TYPE* g, float* c, ARNAUD_CPLXF_TYPE* s, ARNAUD_CPLXF_TYPE* r);
void clascl_(char* mtype, int* kl, int* ku, float* cfrom, float* cto, int* m, int* n, ARNAUD_CPLXF_TYPE* a, int* lda, int* info);
void claset_(char* uplo, int* m, int* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* beta, ARNAUD_CPLXF_TYPE* a, int* lda);
void ctrevc_(char* side, char* howmny, int* select, int* n, ARNAUD_CPLXF_TYPE* t, int* ldt, ARNAUD_CPLXF_TYPE* vl, int* ldvl, ARNAUD_CPLXF_TYPE* vr, int* ldvr, int* mm, int* m, ARNAUD_CPLXF_TYPE* work, float* rwork, int* info);
void ctrsen_(char* job, char* compq, int* select, int* n, ARNAUD_CPLXF_TYPE* t, int* ldt, ARNAUD_CPLXF_TYPE* q, int* ldq, ARNAUD_CPLXF_TYPE* w, int* m, float* s, float* sep, ARNAUD_CPLXF_TYPE* work, int* lwork, int* info);
void cunm2r_(char* side, char* trans, int* m, int* n, int* k, ARNAUD_CPLXF_TYPE* a, int* lda, ARNAUD_CPLXF_TYPE* tau, ARNAUD_CPLXF_TYPE* c, int* ldc, ARNAUD_CPLXF_TYPE* work, int* info);


void zgeqr2_(int* m, int* n, ARNAUD_CPLX_TYPE* a, int* lda, ARNAUD_CPLX_TYPE* tau, ARNAUD_CPLX_TYPE* work, int* info);
void zlacpy_(char* uplo, int* m, int* n, ARNAUD_CPLX_TYPE* a, int* lda, ARNAUD_CPLX_TYPE* b, int* ldb);
void zlahqr_(int* wantt, int* wantz, int* n, int* ilo, int* ihi, ARNAUD_CPLX_TYPE* h, int* ldh, ARNAUD_CPLX_TYPE* w, int* iloz, int* ihiz, ARNAUD_CPLX_TYPE* z, int* ldz, int* info );
double zlanhs_(char* norm, int* n, ARNAUD_CPLX_TYPE* a, int* lda, double* work);
void zlarf_(char* side, int* m, int* n, ARNAUD_CPLX_TYPE* v, int* incv, ARNAUD_CPLX_TYPE* tau, ARNAUD_CPLX_TYPE* c, int* ldc, ARNAUD_CPLX_TYPE* work);
void zlarfg_(int* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* x, int* incx, ARNAUD_CPLX_TYPE* tau);
void zlartg_(ARNAUD_CPLX_TYPE* f, ARNAUD_CPLX_TYPE* g, double* c, ARNAUD_CPLX_TYPE* s, ARNAUD_CPLX_TYPE* r);
void zlascl_(char* mtype, int* kl, int* ku, double* cfrom, double* cto, int* m, int* n, ARNAUD_CPLX_TYPE* a, int* lda, int* info);
void zlaset_(char* uplo, int* m, int* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* beta, ARNAUD_CPLX_TYPE* a, int* lda);
void ztrevc_(char* side, char* howmny, int* select, int* n, ARNAUD_CPLX_TYPE* t, int* ldt, ARNAUD_CPLX_TYPE* vl, int* ldvl, ARNAUD_CPLX_TYPE* vr, int* ldvr, int* mm, int* m, ARNAUD_CPLX_TYPE* work, double* rwork, int* info);
void ztrsen_(char* job, char* compq, int* select, int* n, ARNAUD_CPLX_TYPE* t, int* ldt, ARNAUD_CPLX_TYPE* q, int* ldq, ARNAUD_CPLX_TYPE* w, int* m, double* s, double* sep, ARNAUD_CPLX_TYPE* work, int* lwork, int* info);
void zunm2r_(char* side, char* trans, int* m, int* n, int* k, ARNAUD_CPLX_TYPE* a, int* lda, ARNAUD_CPLX_TYPE* tau, ARNAUD_CPLX_TYPE* c, int* ldc, ARNAUD_CPLX_TYPE* work, int* info);


#endif
