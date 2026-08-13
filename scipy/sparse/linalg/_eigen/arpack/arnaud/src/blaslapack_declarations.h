#ifndef BLASLAPACK_DECLARATIONS_H
#define BLASLAPACK_DECLARATIONS_H

#include "arnaud/types.h"

/*
 * BLAS/LAPACK Fortran symbol mangling.
 *
 * When building within SciPy, `ARNAUD_HAS_BLAS_CONFIG` is defined,
 * that sets the `ARNAUD_BLAS` macro to contain the desired symbol
 * mangling functionality (`BLAS_FUNC()` in other submodules).
 *
 * For standalone builds, ARNAUD_BLAS defaults to appending an underscore
 * (standard Fortran name mangling).
 */

#ifdef ARNAUD_HAS_BLAS_CONFIG
#include "arnaud_blas_config.h"
#endif

#ifndef ARNAUD_BLAS
#define ARNAUD_BLAS(name) name ## _
#endif


// BLAS — single real
void ARNAUD_BLAS(saxpy)(ARNAUD_INT* n, float* alpha, float* x, ARNAUD_INT* incx, float* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(scopy)(ARNAUD_INT* n, float* x, ARNAUD_INT* incx, float* y, ARNAUD_INT* incy);
float ARNAUD_BLAS(sdot)(ARNAUD_INT* n, float* x, ARNAUD_INT* incx, float* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(sgemv)(char* trans, ARNAUD_INT* m, ARNAUD_INT* n, float* alpha, float* a, ARNAUD_INT* lda, float* x, ARNAUD_INT* incx, float* beta, float* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(sger)(ARNAUD_INT* m, ARNAUD_INT* n, float* alpha, float* x, ARNAUD_INT* incx, float* y, ARNAUD_INT* incy, float* a, ARNAUD_INT* lda);
float ARNAUD_BLAS(snrm2)(ARNAUD_INT* n, float* x, ARNAUD_INT* incx);
void ARNAUD_BLAS(srot)(ARNAUD_INT* n, float* sx, ARNAUD_INT* incx, float* sy, ARNAUD_INT* incy, float* c, float* s);
void ARNAUD_BLAS(sscal)(ARNAUD_INT* n, float* alpha, float* x, ARNAUD_INT* incx);
void ARNAUD_BLAS(sswap)(ARNAUD_INT* n, float* x, ARNAUD_INT* incx, float* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(strmm)(char* side, char* uplo, char* transa, char* diag, ARNAUD_INT* m, ARNAUD_INT* n, float* alpha, float* a, ARNAUD_INT* lda, float* b, ARNAUD_INT* ldb);


// BLAS — double real
void ARNAUD_BLAS(daxpy)(ARNAUD_INT* n, double* alpha, double* x, ARNAUD_INT* incx, double* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(dcopy)(ARNAUD_INT* n, double* x, ARNAUD_INT* incx, double* y, ARNAUD_INT* incy);
double ARNAUD_BLAS(ddot)(ARNAUD_INT* n, double* x, ARNAUD_INT* incx, double* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(dgemv)(char* trans, ARNAUD_INT* m, ARNAUD_INT* n, double* alpha, double* a, ARNAUD_INT* lda, double* x, ARNAUD_INT* incx, double* beta, double* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(dger)(ARNAUD_INT* m, ARNAUD_INT* n, double* alpha, double* x, ARNAUD_INT* incx, double* y, ARNAUD_INT* incy, double* a, ARNAUD_INT* lda);
double ARNAUD_BLAS(dnrm2)(ARNAUD_INT* n, double* x, ARNAUD_INT* incx);
void ARNAUD_BLAS(drot)(ARNAUD_INT* n, double* sx, ARNAUD_INT* incx, double* sy, ARNAUD_INT* incy, double* c, double* s);
void ARNAUD_BLAS(dscal)(ARNAUD_INT* n, double* alpha, double* x, ARNAUD_INT* incx);
void ARNAUD_BLAS(dswap)(ARNAUD_INT* n, double* x, ARNAUD_INT* incx, double* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(dtrmm)(char* side, char* uplo, char* transa, char* diag, ARNAUD_INT* m, ARNAUD_INT* n, double* alpha, double* a, ARNAUD_INT* lda, double* b, ARNAUD_INT* ldb);


// BLAS — single complex
void ARNAUD_BLAS(caxpy)(ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* x, ARNAUD_INT* incx, ARNAUD_CPLXF_TYPE* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(ccopy)(ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* x, ARNAUD_INT* incx, ARNAUD_CPLXF_TYPE* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(cgeru)(ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* x, ARNAUD_INT* incx, ARNAUD_CPLXF_TYPE* y, ARNAUD_INT* incy, ARNAUD_CPLXF_TYPE* a, ARNAUD_INT* lda);
float ARNAUD_BLAS(scnrm2)(ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* x, ARNAUD_INT* incx);
void ARNAUD_BLAS(cscal)(ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* x, ARNAUD_INT* incx);
void ARNAUD_BLAS(csscal)(ARNAUD_INT* n, float* da, ARNAUD_CPLXF_TYPE* zx, ARNAUD_INT* incx);
void ARNAUD_BLAS(cgemv)(char* trans, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* a, ARNAUD_INT* lda, ARNAUD_CPLXF_TYPE* x, ARNAUD_INT* incx, ARNAUD_CPLXF_TYPE* beta, ARNAUD_CPLXF_TYPE* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(crot)(ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* cx, ARNAUD_INT* incx, ARNAUD_CPLXF_TYPE* cy, ARNAUD_INT* incy, float* c, ARNAUD_CPLXF_TYPE* s);
void ARNAUD_BLAS(ctrmm)(char* side, char* uplo, char* transa, char* diag, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* a, ARNAUD_INT* lda, ARNAUD_CPLXF_TYPE* b, ARNAUD_INT* ldb);


// BLAS — double complex
void ARNAUD_BLAS(zaxpy)(ARNAUD_INT* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* x, ARNAUD_INT* incx, ARNAUD_CPLX_TYPE* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(zcopy)(ARNAUD_INT* n, ARNAUD_CPLX_TYPE* x, ARNAUD_INT* incx, ARNAUD_CPLX_TYPE* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(zgeru)(ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* x, ARNAUD_INT* incx, ARNAUD_CPLX_TYPE* y, ARNAUD_INT* incy, ARNAUD_CPLX_TYPE* a, ARNAUD_INT* lda);
double ARNAUD_BLAS(dznrm2)(ARNAUD_INT* n, ARNAUD_CPLX_TYPE* x, ARNAUD_INT* incx);
void ARNAUD_BLAS(zscal)(ARNAUD_INT* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* x, ARNAUD_INT* incx);
void ARNAUD_BLAS(zdscal)(ARNAUD_INT* n, double* da, ARNAUD_CPLX_TYPE* zx, ARNAUD_INT* incx);
void ARNAUD_BLAS(zgemv)(char* trans, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* a, ARNAUD_INT* lda, ARNAUD_CPLX_TYPE* x, ARNAUD_INT* incx, ARNAUD_CPLX_TYPE* beta, ARNAUD_CPLX_TYPE* y, ARNAUD_INT* incy);
void ARNAUD_BLAS(zrot)(ARNAUD_INT* n, ARNAUD_CPLX_TYPE* cx, ARNAUD_INT* incx, ARNAUD_CPLX_TYPE* cy, ARNAUD_INT* incy, double* c, ARNAUD_CPLX_TYPE* s);
void ARNAUD_BLAS(ztrmm)(char* side, char* uplo, char* transa, char* diag, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* a, ARNAUD_INT* lda, ARNAUD_CPLX_TYPE* b, ARNAUD_INT* ldb);



// LAPACK — single real
void ARNAUD_BLAS(sgeqr2)(ARNAUD_INT* m, ARNAUD_INT* n, float* a, ARNAUD_INT* lda, float* tau, float* work, ARNAUD_INT* info);
void ARNAUD_BLAS(slacpy)(char* uplo, ARNAUD_INT* m, ARNAUD_INT* n, float* a, ARNAUD_INT* lda, float* b, ARNAUD_INT* ldb);
void ARNAUD_BLAS(slaev2)(float* a, float* b, float* c, float* rt1, float* rt2, float* cs1, float* sn1);
void ARNAUD_BLAS(slahqr)(ARNAUD_INT* wantt, ARNAUD_INT* wantz, ARNAUD_INT* n, ARNAUD_INT* ilo, ARNAUD_INT* ihi, float* h, ARNAUD_INT* ldh, float* wr, float* wi, ARNAUD_INT* iloz, ARNAUD_INT* ihiz, float* z, ARNAUD_INT* ldz, ARNAUD_INT* info);
float ARNAUD_BLAS(slanhs)(char* norm, ARNAUD_INT* n, float* a, ARNAUD_INT* lda, float* work);
float ARNAUD_BLAS(slanst)(char* norm, ARNAUD_INT* n, float* d, float* e);
void ARNAUD_BLAS(slarf)(char* side, ARNAUD_INT* m, ARNAUD_INT* n, float* v, ARNAUD_INT* incv, float* tau, float* c, ARNAUD_INT* ldc, float* work);
void ARNAUD_BLAS(slarfg)(ARNAUD_INT* n, float* alpha, float* x, ARNAUD_INT* incx, float* tau);
void ARNAUD_BLAS(slartg)(float* f, float* g, float* c, float* s, float* r);
void ARNAUD_BLAS(slartgp)(float* f, float* g, float* c, float* s, float* r);
void ARNAUD_BLAS(slascl)(char* mtype, ARNAUD_INT* kl, ARNAUD_INT* ku, float* cfrom, float* cto, ARNAUD_INT* m, ARNAUD_INT* n, float* a, ARNAUD_INT* lda, ARNAUD_INT* info);
void ARNAUD_BLAS(slaset)(char* uplo, ARNAUD_INT* m, ARNAUD_INT* n, float* alpha, float* beta, float* a, ARNAUD_INT* lda);
void ARNAUD_BLAS(slasr)(char* side, char* pivot, char* direct, ARNAUD_INT* m, ARNAUD_INT* n, float* c, float* s, float* a, ARNAUD_INT* lda);
void ARNAUD_BLAS(sorm2r)(char* side, char* trans, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_INT* k, float* a, ARNAUD_INT* lda, float* tau, float* c, ARNAUD_INT* ldc, float* work, ARNAUD_INT* info);
void ARNAUD_BLAS(ssteqr)(char* compz, ARNAUD_INT* n, float* d, float* e, float* z, ARNAUD_INT* ldz, float* work, ARNAUD_INT* info);
void ARNAUD_BLAS(strevc)(char* side, char* howmny, ARNAUD_INT* select, ARNAUD_INT* n, float* t, ARNAUD_INT* ldt, float* vl, ARNAUD_INT* ldvl, float* vr, ARNAUD_INT* ldvr, ARNAUD_INT* mm, ARNAUD_INT* m, float* work, ARNAUD_INT* info);
void ARNAUD_BLAS(strsen)(char* job, char* compq, ARNAUD_INT* select, ARNAUD_INT* n, float* t, ARNAUD_INT* ldt, float* q, ARNAUD_INT* ldq, float* wr, float* wi, ARNAUD_INT* m, float* s, float* sep, float* work, ARNAUD_INT* lwork, ARNAUD_INT* iwork, ARNAUD_INT* liwork, ARNAUD_INT* info);


// LAPACK — double real
void ARNAUD_BLAS(dgeqr2)(ARNAUD_INT* m, ARNAUD_INT* n, double* a, ARNAUD_INT* lda, double* tau, double* work, ARNAUD_INT* info);
void ARNAUD_BLAS(dlacpy)(char* uplo, ARNAUD_INT* m, ARNAUD_INT* n, double* a, ARNAUD_INT* lda, double* b, ARNAUD_INT* ldb);
void ARNAUD_BLAS(dlaev2)(double* a, double* b, double* c, double* rt1, double* rt2, double* cs1, double* sn1);
void ARNAUD_BLAS(dlahqr)(ARNAUD_INT* wantt, ARNAUD_INT* wantz, ARNAUD_INT* n, ARNAUD_INT* ilo, ARNAUD_INT* ihi, double* h, ARNAUD_INT* ldh, double* wr, double* wi, ARNAUD_INT* iloz, ARNAUD_INT* ihiz, double* z, ARNAUD_INT* ldz, ARNAUD_INT* info);
double ARNAUD_BLAS(dlanhs)(char* norm, ARNAUD_INT* n, double* a, ARNAUD_INT* lda, double* work);
double ARNAUD_BLAS(dlanst)(char* norm, ARNAUD_INT* n, double* d, double* e);
void ARNAUD_BLAS(dlarf)(char* side, ARNAUD_INT* m, ARNAUD_INT* n, double* v, ARNAUD_INT* incv, double* tau, double* c, ARNAUD_INT* ldc, double* work);
void ARNAUD_BLAS(dlarfg)(ARNAUD_INT* n, double* alpha, double* x, ARNAUD_INT* incx, double* tau);
void ARNAUD_BLAS(dlartg)(double* f, double* g, double* c, double* s, double* r);
void ARNAUD_BLAS(dlartgp)(double* f, double* g, double* c, double* s, double* r);
void ARNAUD_BLAS(dlascl)(char* mtype, ARNAUD_INT* kl, ARNAUD_INT* ku, double* cfrom, double* cto, ARNAUD_INT* m, ARNAUD_INT* n, double* a, ARNAUD_INT* lda, ARNAUD_INT* info);
void ARNAUD_BLAS(dlaset)(char* uplo, ARNAUD_INT* m, ARNAUD_INT* n, double* alpha, double* beta, double* a, ARNAUD_INT* lda);
void ARNAUD_BLAS(dlasr)(char* side, char* pivot, char* direct, ARNAUD_INT* m, ARNAUD_INT* n, double* c, double* s, double* a, ARNAUD_INT* lda);
void ARNAUD_BLAS(dorm2r)(char* side, char* trans, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_INT* k, double* a, ARNAUD_INT* lda, double* tau, double* c, ARNAUD_INT* ldc, double* work, ARNAUD_INT* info);
void ARNAUD_BLAS(dsteqr)(char* compz, ARNAUD_INT* n, double* d, double* e, double* z, ARNAUD_INT* ldz, double* work, ARNAUD_INT* info);
void ARNAUD_BLAS(dtrevc)(char* side, char* howmny, ARNAUD_INT* select, ARNAUD_INT* n, double* t, ARNAUD_INT* ldt, double* vl, ARNAUD_INT* ldvl, double* vr, ARNAUD_INT* ldvr, ARNAUD_INT* mm, ARNAUD_INT* m, double* work, ARNAUD_INT* info);
void ARNAUD_BLAS(dtrsen)(char* job, char* compq, ARNAUD_INT* select, ARNAUD_INT* n, double* t, ARNAUD_INT* ldt, double* q, ARNAUD_INT* ldq, double* wr, double* wi, ARNAUD_INT* m, double* s, double* sep, double* work, ARNAUD_INT* lwork, ARNAUD_INT* iwork, ARNAUD_INT* liwork, ARNAUD_INT* info);


// LAPACK — single complex
void ARNAUD_BLAS(cgeqr2)(ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* a, ARNAUD_INT* lda, ARNAUD_CPLXF_TYPE* tau, ARNAUD_CPLXF_TYPE* work, ARNAUD_INT* info);
void ARNAUD_BLAS(clacpy)(char* uplo, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* a, ARNAUD_INT* lda, ARNAUD_CPLXF_TYPE* b, ARNAUD_INT* ldb);
void ARNAUD_BLAS(clahqr)(ARNAUD_INT* wantt, ARNAUD_INT* wantz, ARNAUD_INT* n, ARNAUD_INT* ilo, ARNAUD_INT* ihi, ARNAUD_CPLXF_TYPE* h, ARNAUD_INT* ldh, ARNAUD_CPLXF_TYPE* w, ARNAUD_INT* iloz, ARNAUD_INT* ihiz, ARNAUD_CPLXF_TYPE* z, ARNAUD_INT* ldz, ARNAUD_INT* info);
float ARNAUD_BLAS(clanhs)(char* norm, ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* a, ARNAUD_INT* lda, float* work);
void ARNAUD_BLAS(clarf)(char* side, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* v, ARNAUD_INT* incv, ARNAUD_CPLXF_TYPE* tau, ARNAUD_CPLXF_TYPE* c, ARNAUD_INT* ldc, ARNAUD_CPLXF_TYPE* work);
void ARNAUD_BLAS(clarfg)(ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* x, ARNAUD_INT* incx, ARNAUD_CPLXF_TYPE* tau);
void ARNAUD_BLAS(clartg)(ARNAUD_CPLXF_TYPE* f, ARNAUD_CPLXF_TYPE* g, float* c, ARNAUD_CPLXF_TYPE* s, ARNAUD_CPLXF_TYPE* r);
void ARNAUD_BLAS(clascl)(char* mtype, ARNAUD_INT* kl, ARNAUD_INT* ku, float* cfrom, float* cto, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* a, ARNAUD_INT* lda, ARNAUD_INT* info);
void ARNAUD_BLAS(claset)(char* uplo, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* alpha, ARNAUD_CPLXF_TYPE* beta, ARNAUD_CPLXF_TYPE* a, ARNAUD_INT* lda);
void ARNAUD_BLAS(ctrevc)(char* side, char* howmny, ARNAUD_INT* select, ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* t, ARNAUD_INT* ldt, ARNAUD_CPLXF_TYPE* vl, ARNAUD_INT* ldvl, ARNAUD_CPLXF_TYPE* vr, ARNAUD_INT* ldvr, ARNAUD_INT* mm, ARNAUD_INT* m, ARNAUD_CPLXF_TYPE* work, float* rwork, ARNAUD_INT* info);
void ARNAUD_BLAS(ctrsen)(char* job, char* compq, ARNAUD_INT* select, ARNAUD_INT* n, ARNAUD_CPLXF_TYPE* t, ARNAUD_INT* ldt, ARNAUD_CPLXF_TYPE* q, ARNAUD_INT* ldq, ARNAUD_CPLXF_TYPE* w, ARNAUD_INT* m, float* s, float* sep, ARNAUD_CPLXF_TYPE* work, ARNAUD_INT* lwork, ARNAUD_INT* info);
void ARNAUD_BLAS(cunm2r)(char* side, char* trans, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_INT* k, ARNAUD_CPLXF_TYPE* a, ARNAUD_INT* lda, ARNAUD_CPLXF_TYPE* tau, ARNAUD_CPLXF_TYPE* c, ARNAUD_INT* ldc, ARNAUD_CPLXF_TYPE* work, ARNAUD_INT* info);


// LAPACK — double complex
void ARNAUD_BLAS(zgeqr2)(ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLX_TYPE* a, ARNAUD_INT* lda, ARNAUD_CPLX_TYPE* tau, ARNAUD_CPLX_TYPE* work, ARNAUD_INT* info);
void ARNAUD_BLAS(zlacpy)(char* uplo, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLX_TYPE* a, ARNAUD_INT* lda, ARNAUD_CPLX_TYPE* b, ARNAUD_INT* ldb);
void ARNAUD_BLAS(zlahqr)(ARNAUD_INT* wantt, ARNAUD_INT* wantz, ARNAUD_INT* n, ARNAUD_INT* ilo, ARNAUD_INT* ihi, ARNAUD_CPLX_TYPE* h, ARNAUD_INT* ldh, ARNAUD_CPLX_TYPE* w, ARNAUD_INT* iloz, ARNAUD_INT* ihiz, ARNAUD_CPLX_TYPE* z, ARNAUD_INT* ldz, ARNAUD_INT* info);
double ARNAUD_BLAS(zlanhs)(char* norm, ARNAUD_INT* n, ARNAUD_CPLX_TYPE* a, ARNAUD_INT* lda, double* work);
void ARNAUD_BLAS(zlarf)(char* side, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLX_TYPE* v, ARNAUD_INT* incv, ARNAUD_CPLX_TYPE* tau, ARNAUD_CPLX_TYPE* c, ARNAUD_INT* ldc, ARNAUD_CPLX_TYPE* work);
void ARNAUD_BLAS(zlarfg)(ARNAUD_INT* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* x, ARNAUD_INT* incx, ARNAUD_CPLX_TYPE* tau);
void ARNAUD_BLAS(zlartg)(ARNAUD_CPLX_TYPE* f, ARNAUD_CPLX_TYPE* g, double* c, ARNAUD_CPLX_TYPE* s, ARNAUD_CPLX_TYPE* r);
void ARNAUD_BLAS(zlascl)(char* mtype, ARNAUD_INT* kl, ARNAUD_INT* ku, double* cfrom, double* cto, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLX_TYPE* a, ARNAUD_INT* lda, ARNAUD_INT* info);
void ARNAUD_BLAS(zlaset)(char* uplo, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_CPLX_TYPE* alpha, ARNAUD_CPLX_TYPE* beta, ARNAUD_CPLX_TYPE* a, ARNAUD_INT* lda);
void ARNAUD_BLAS(ztrevc)(char* side, char* howmny, ARNAUD_INT* select, ARNAUD_INT* n, ARNAUD_CPLX_TYPE* t, ARNAUD_INT* ldt, ARNAUD_CPLX_TYPE* vl, ARNAUD_INT* ldvl, ARNAUD_CPLX_TYPE* vr, ARNAUD_INT* ldvr, ARNAUD_INT* mm, ARNAUD_INT* m, ARNAUD_CPLX_TYPE* work, double* rwork, ARNAUD_INT* info);
void ARNAUD_BLAS(ztrsen)(char* job, char* compq, ARNAUD_INT* select, ARNAUD_INT* n, ARNAUD_CPLX_TYPE* t, ARNAUD_INT* ldt, ARNAUD_CPLX_TYPE* q, ARNAUD_INT* ldq, ARNAUD_CPLX_TYPE* w, ARNAUD_INT* m, double* s, double* sep, ARNAUD_CPLX_TYPE* work, ARNAUD_INT* lwork, ARNAUD_INT* info);
void ARNAUD_BLAS(zunm2r)(char* side, char* trans, ARNAUD_INT* m, ARNAUD_INT* n, ARNAUD_INT* k, ARNAUD_CPLX_TYPE* a, ARNAUD_INT* lda, ARNAUD_CPLX_TYPE* tau, ARNAUD_CPLX_TYPE* c, ARNAUD_INT* ldc, ARNAUD_CPLX_TYPE* work, ARNAUD_INT* info);


#endif
