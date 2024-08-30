#ifndef _ARPACK_N_SINGLE_COMPLEX_H
#define _ARPACK_N_SINGLE_COMPLEX_H

#include "_arpack.h"

// BLAS Routines used
void caxpy_(int* n, ARPACK_CPLXF_TYPE* alpha, ARPACK_CPLXF_TYPE* x, int* incx, ARPACK_CPLXF_TYPE* y, int* incy);
void ccopy_(int* n, ARPACK_CPLXF_TYPE* x, int* incx, ARPACK_CPLXF_TYPE* y, int* incy);
ARPACK_CPLXF_TYPE cdotc_(int* n, ARPACK_CPLXF_TYPE* x, int* incx, ARPACK_CPLXF_TYPE* y, int* incy);
float scnrm2_(int* n, ARPACK_CPLXF_TYPE* x, int* incx);
void cscal_(int* n, ARPACK_CPLXF_TYPE* alpha, ARPACK_CPLXF_TYPE* x, int* incx);
void csscal_(int* n, float* da, ARPACK_CPLXF_TYPE* zx, int* incx);
void cgemv_(char* trans, int* m, int* n, ARPACK_CPLXF_TYPE* alpha, ARPACK_CPLXF_TYPE* a, int* lda, ARPACK_CPLXF_TYPE* x, int* incx, ARPACK_CPLXF_TYPE* beta, ARPACK_CPLXF_TYPE* y, int* incy);


// LAPACK Routines used
void clacpy_(char* uplo, int* m, int* n, ARPACK_CPLXF_TYPE* a, int* lda, ARPACK_CPLXF_TYPE* b, int* ldb);
void clahqr_(int* wantt, int* wantz, int* n, int* ilo, int* ihi, ARPACK_CPLXF_TYPE* h, int* ldh, ARPACK_CPLXF_TYPE* w, int* iloz, int* ihiz, ARPACK_CPLXF_TYPE* z, int* ldz, int* info );
float clanhs_(char* norm, int* n, ARPACK_CPLXF_TYPE* a, int* lda, float* work);
void clarf_(char* side, int* m, int* n, ARPACK_CPLXF_TYPE* v, int* incv, ARPACK_CPLXF_TYPE* tau, ARPACK_CPLXF_TYPE* c, int* ldc, ARPACK_CPLXF_TYPE* work);
void clarfg_(int* n, ARPACK_CPLXF_TYPE* alpha, ARPACK_CPLXF_TYPE* x, int* incx, ARPACK_CPLXF_TYPE* tau);
void clartg_(ARPACK_CPLXF_TYPE* f, ARPACK_CPLXF_TYPE* g, float* c, ARPACK_CPLXF_TYPE* s, ARPACK_CPLXF_TYPE* r);
void clascl_(char* mtype, int* kl, int* ku, float* cfrom, float* cto, int* m, int* n, ARPACK_CPLXF_TYPE* a, int* lda, int* info);
void claset_(char* uplo, int* m, int* n, ARPACK_CPLXF_TYPE* alpha, ARPACK_CPLXF_TYPE* beta, ARPACK_CPLXF_TYPE* a, int* lda);
void ctrevc_(char* side, char* howmny, int* select, int* n, ARPACK_CPLXF_TYPE* t, int* ldt, ARPACK_CPLXF_TYPE* vl, int* ldvl, ARPACK_CPLXF_TYPE* vr, int* ldvr, int* mm, int* m, ARPACK_CPLXF_TYPE* work, float* rwork, int* info);


#endif
