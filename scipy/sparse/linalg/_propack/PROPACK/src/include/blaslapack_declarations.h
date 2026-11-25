#ifndef BLASLAPACK_DECLARATIONS_H
#define BLASLAPACK_DECLARATIONS_H

#include "../include/propack/types.h"


// BLAS
void BLAS_FUNC(saxpy)(CBLAS_INT* n, float* alpha, float* x, CBLAS_INT* incx, float* y, CBLAS_INT* incy);
void BLAS_FUNC(scopy)(CBLAS_INT* n, float* x, CBLAS_INT* incx, float* y, CBLAS_INT* incy);
float BLAS_FUNC(sdot)(CBLAS_INT* n, float* x, CBLAS_INT* incx, float* y, CBLAS_INT* incy);
void BLAS_FUNC(sgemm)(char* transa, char* transb, CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* k, float* alpha, float* a, CBLAS_INT* lda, float* b, CBLAS_INT* ldb, float* beta, float* c, CBLAS_INT* ldc);
void BLAS_FUNC(sgemv)(char* trans, CBLAS_INT* m, CBLAS_INT* n, float* alpha, float* a, CBLAS_INT* lda, float* x, CBLAS_INT* incx, float* beta, float* y, CBLAS_INT* incy);
float BLAS_FUNC(snrm2)(CBLAS_INT* n, float* x, CBLAS_INT* incx);
void BLAS_FUNC(srot)(CBLAS_INT* n, float* sx, CBLAS_INT* incx, float* sy, CBLAS_INT* incy, float* c, float* s);
void BLAS_FUNC(sscal)(CBLAS_INT* n, float* alpha, float* x, CBLAS_INT* incx);

void BLAS_FUNC(daxpy)(CBLAS_INT* n, double* alpha, double* x, CBLAS_INT* incx, double* y, CBLAS_INT* incy);
void BLAS_FUNC(dcopy)(CBLAS_INT* n, double* x, CBLAS_INT* incx, double* y, CBLAS_INT* incy);
double BLAS_FUNC(ddot)(CBLAS_INT* n, double* x, CBLAS_INT* incx, double* y, CBLAS_INT* incy);
void BLAS_FUNC(dgemm)(char* transa, char* transb, CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* k, double* alpha, double* a, CBLAS_INT* lda, double* b, CBLAS_INT* ldb, double* beta, double* c, CBLAS_INT* ldc);
void BLAS_FUNC(dgemv)(char* trans, CBLAS_INT* m, CBLAS_INT* n, double* alpha, double* a, CBLAS_INT* lda, double* x, CBLAS_INT* incx, double* beta, double* y, CBLAS_INT* incy);
double BLAS_FUNC(dnrm2)(CBLAS_INT* n, double* x, CBLAS_INT* incx);
void BLAS_FUNC(drot)(CBLAS_INT* n, double* sx, CBLAS_INT* incx, double* sy, CBLAS_INT* incy, double* c, double* s);
void BLAS_FUNC(dscal)(CBLAS_INT* n, double* alpha, double* x, CBLAS_INT* incx);

void BLAS_FUNC(caxpy)(CBLAS_INT* n, PROPACK_CPLXF_TYPE* alpha, PROPACK_CPLXF_TYPE* x, CBLAS_INT* incx, PROPACK_CPLXF_TYPE* y, CBLAS_INT* incy);
float BLAS_FUNC(scnrm2)(CBLAS_INT* n, PROPACK_CPLXF_TYPE* x, CBLAS_INT* incx);
void BLAS_FUNC(cgemm)(char* transa, char* transb, CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* k, PROPACK_CPLXF_TYPE* alpha, PROPACK_CPLXF_TYPE* a, CBLAS_INT* lda, PROPACK_CPLXF_TYPE* b, CBLAS_INT* ldb, PROPACK_CPLXF_TYPE* beta, PROPACK_CPLXF_TYPE* c, CBLAS_INT* ldc);
void BLAS_FUNC(cgemv)(char* trans, CBLAS_INT* m, CBLAS_INT* n, PROPACK_CPLXF_TYPE* alpha, PROPACK_CPLXF_TYPE* a, CBLAS_INT* lda, PROPACK_CPLXF_TYPE* x, CBLAS_INT* incx, PROPACK_CPLXF_TYPE* beta, PROPACK_CPLXF_TYPE* y, CBLAS_INT* incy);
void BLAS_FUNC(csscal)(CBLAS_INT* n, float* da, PROPACK_CPLXF_TYPE* zx, CBLAS_INT* incx);

void BLAS_FUNC(zaxpy)(CBLAS_INT* n, PROPACK_CPLX_TYPE* alpha, PROPACK_CPLX_TYPE* x, CBLAS_INT* incx, PROPACK_CPLX_TYPE* y, CBLAS_INT* incy);
double BLAS_FUNC(dznrm2)(CBLAS_INT* n, PROPACK_CPLX_TYPE* x, CBLAS_INT* incx);
void BLAS_FUNC(zgemm)(char* transa, char* transb, CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* k, PROPACK_CPLX_TYPE* alpha, PROPACK_CPLX_TYPE* a, CBLAS_INT* lda, PROPACK_CPLX_TYPE* b, CBLAS_INT* ldb, PROPACK_CPLX_TYPE* beta, PROPACK_CPLX_TYPE* c, CBLAS_INT* ldc);
void BLAS_FUNC(zgemv)(char* trans, CBLAS_INT* m, CBLAS_INT* n, PROPACK_CPLX_TYPE* alpha, PROPACK_CPLX_TYPE* a, CBLAS_INT* lda, PROPACK_CPLX_TYPE* x, CBLAS_INT* incx, PROPACK_CPLX_TYPE* beta, PROPACK_CPLX_TYPE* y, CBLAS_INT* incy);
void BLAS_FUNC(zdscal)(CBLAS_INT* n, double* da, PROPACK_CPLX_TYPE* zx, CBLAS_INT* incx);

// LAPACK

void BLAS_FUNC(sbdsdc)(char* uplo, char* compq, CBLAS_INT* n, float* d, float* e, float* u, CBLAS_INT* ldu, float* vt, CBLAS_INT* ldvt, float* q, CBLAS_INT* iq, float* work, CBLAS_INT* iwork, CBLAS_INT* info);
void BLAS_FUNC(sbdsqr)(char* uplo, CBLAS_INT* n, CBLAS_INT* ncvt, CBLAS_INT* nru, CBLAS_INT* ncc, float* d, float* e, float* vt, CBLAS_INT* ldvt, float* u, CBLAS_INT* ldu, float* c, CBLAS_INT* ldc, float* work, CBLAS_INT* info);
void BLAS_FUNC(slartg)(float* f, float* g, float* c, float* s, float* r);
void BLAS_FUNC(slascl)(char* mtype, CBLAS_INT* kl, CBLAS_INT* ku, float* cfrom, float* cto, CBLAS_INT* m, CBLAS_INT* n, float* a, CBLAS_INT* lda, CBLAS_INT* info);
void BLAS_FUNC(slaset)(char* uplo, CBLAS_INT* m, CBLAS_INT* n, float* alpha, float* beta, float* a, CBLAS_INT* lda);

void BLAS_FUNC(dbdsdc)(char* uplo, char* compq, CBLAS_INT* n, double* d, double* e, double* u, CBLAS_INT* ldu, double* vt, CBLAS_INT* ldvt, double* q, CBLAS_INT* iq, double* work, CBLAS_INT* iwork, CBLAS_INT* info);
void BLAS_FUNC(dbdsqr)(char* uplo, CBLAS_INT* n, CBLAS_INT* ncvt, CBLAS_INT* nru, CBLAS_INT* ncc, double* d, double* e, double* vt, CBLAS_INT* ldvt, double* u, CBLAS_INT* ldu, double* c, CBLAS_INT* ldc, double* work, CBLAS_INT* info);
void BLAS_FUNC(dlartg)(double* f, double* g, double* c, double* s, double* r);
void BLAS_FUNC(dlascl)(char* mtype, CBLAS_INT* kl, CBLAS_INT* ku, double* cfrom, double* cto, CBLAS_INT* m, CBLAS_INT* n, double* a, CBLAS_INT* lda, CBLAS_INT* info);
void BLAS_FUNC(dlaset)(char* uplo, CBLAS_INT* m, CBLAS_INT* n, double* alpha, double* beta, double* a, CBLAS_INT* lda);

void BLAS_FUNC(clarfg)(CBLAS_INT* n, PROPACK_CPLXF_TYPE* alpha, PROPACK_CPLXF_TYPE* x, CBLAS_INT* incx, PROPACK_CPLXF_TYPE* tau);
void BLAS_FUNC(clascl)(char* mtype, CBLAS_INT* kl, CBLAS_INT* ku, float* cfrom, float* cto, CBLAS_INT* m, CBLAS_INT* n, PROPACK_CPLXF_TYPE* a, CBLAS_INT* lda, CBLAS_INT* info);

void BLAS_FUNC(zlarfg)(CBLAS_INT* n, PROPACK_CPLX_TYPE* alpha, PROPACK_CPLX_TYPE* x, CBLAS_INT* incx, PROPACK_CPLX_TYPE* tau);
void BLAS_FUNC(zlascl)(char* mtype, CBLAS_INT* kl, CBLAS_INT* ku, double* cfrom, double* cto, CBLAS_INT* m, CBLAS_INT* n, PROPACK_CPLX_TYPE* a, CBLAS_INT* lda, CBLAS_INT* info);


// (c,z)dotc is the complex conjugate dot product of two complex vectors.
// Due some historical reasons, this function can cause segfaults on some
// platforms. Hence implemented here instead of using the BLAS version.
static PROPACK_CPLXF_TYPE
cdotc_(const CBLAS_INT* n, const PROPACK_CPLXF_TYPE* restrict x, const CBLAS_INT* incx, const PROPACK_CPLXF_TYPE* restrict y, const CBLAS_INT* incy)
{
    PROPACK_CPLXF_TYPE result = PROPACK_cplxf(0.0, 0.0);
#ifdef _MSC_VER
    PROPACK_CPLXF_TYPE temp = PROPACK_cplxf(0.0, 0.0);
#endif
    if (*n <= 0) { return result; }
    if ((*incx == 1) && (*incy == 1))
    {
        for (CBLAS_INT i = 0; i < *n; i++)
        {
#ifdef _MSC_VER
            temp = _FCmulcc(conjf(x[i]), y[i]);
            result = PROPACK_cplxf(crealf(result) + crealf(temp), cimagf(result) + cimagf(temp));
#else
            result = result + (conjf(x[i]) * y[i]);
#endif
        }

    } else {

        for (CBLAS_INT i = 0; i < *n; i++)
        {
#ifdef _MSC_VER
            temp = _FCmulcc(conjf(x[i * (*incx)]), y[i * (*incy)]);
            result = PROPACK_cplxf(crealf(result) + crealf(temp), cimagf(result) + cimagf(temp));
#else
            result = result + (conjf(x[i * (*incx)]) * y[i * (*incy)]);
#endif
        }
    }

    return result;
}


static PROPACK_CPLX_TYPE
zdotc_(const CBLAS_INT* n, const PROPACK_CPLX_TYPE* restrict x, const CBLAS_INT* incx, const PROPACK_CPLX_TYPE* restrict y, const CBLAS_INT* incy)
{
    PROPACK_CPLX_TYPE result = PROPACK_cplx(0.0, 0.0);
#ifdef _MSC_VER
    PROPACK_CPLX_TYPE temp = PROPACK_cplx(0.0, 0.0);
#endif
    if (*n <= 0) { return result; }
    if ((*incx == 1) && (*incy == 1))
    {
        for (CBLAS_INT i = 0; i < *n; i++)
        {
#ifdef _MSC_VER
            temp = _Cmulcc(conj(x[i]), y[i]);
            result = PROPACK_cplx(creal(result) + creal(temp), cimag(result) + cimag(temp));
#else
            result = result + (conj(x[i]) * y[i]);
#endif
        }

    } else {

        for (CBLAS_INT i = 0; i < *n; i++)
        {
#ifdef _MSC_VER
            temp = _Cmulcc(conj(x[i * (*incx)]), y[i * (*incy)]);
            result = PROPACK_cplx(creal(result) + creal(temp), cimag(result) + cimag(temp));
#else
            result = result + (conj(x[i * (*incx)]) * y[i * (*incy)]);
#endif
        }
    }

    return result;
}

#endif
