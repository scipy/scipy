#ifndef BLASLAPACK_DECLARATIONS_H
#define BLASLAPACK_DECLARATIONS_H

#include "../include/propack/types.h"


// BLAS
void saxpy_(int* n, float* alpha, float* x, int* incx, float* y, int* incy);
void scopy_(int* n, float* x, int* incx, float* y, int* incy);
float sdot_(int* n, float* x, int* incx, float* y, int* incy);
void sgemm_(char* transa, char* transb, int* m, int* n, int* k, float* alpha, float* a, int* lda, float* b, int* ldb, float* beta, float* c, int* ldc);
void sgemv_(char* trans, int* m, int* n, float* alpha, float* a, int* lda, float* x, int* incx, float* beta, float* y, int* incy);
float snrm2_(int* n, float* x, int* incx);
void srot_(int* n, float* sx, int* incx, float* sy, int* incy, float* c, float* s);
void sscal_(int* n, float* alpha, float* x, int* incx);

void daxpy_(int* n, double* alpha, double* x, int* incx, double* y, int* incy);
void dcopy_(int* n, double* x, int* incx, double* y, int* incy);
double ddot_(int* n, double* x, int* incx, double* y, int* incy);
void dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);
void dgemv_(char* trans, int* m, int* n, double* alpha, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
double dnrm2_(int* n, double* x, int* incx);
void drot_(int* n, double* sx, int* incx, double* sy, int* incy, double* c, double* s);
void dscal_(int* n, double* alpha, double* x, int* incx);

void caxpy_(int* n, PROPACK_CPLXF_TYPE* alpha, PROPACK_CPLXF_TYPE* x, int* incx, PROPACK_CPLXF_TYPE* y, int* incy);
float scnrm2_(int* n, PROPACK_CPLXF_TYPE* x, int* incx);
void cgemm_(char* transa, char* transb, int* m, int* n, int* k, PROPACK_CPLXF_TYPE* alpha, PROPACK_CPLXF_TYPE* a, int* lda, PROPACK_CPLXF_TYPE* b, int* ldb, PROPACK_CPLXF_TYPE* beta, PROPACK_CPLXF_TYPE* c, int* ldc);
void cgemv_(char* trans, int* m, int* n, PROPACK_CPLXF_TYPE* alpha, PROPACK_CPLXF_TYPE* a, int* lda, PROPACK_CPLXF_TYPE* x, int* incx, PROPACK_CPLXF_TYPE* beta, PROPACK_CPLXF_TYPE* y, int* incy);
void csscal_(int* n, float* da, PROPACK_CPLXF_TYPE* zx, int* incx);

void zaxpy_(int* n, PROPACK_CPLX_TYPE* alpha, PROPACK_CPLX_TYPE* x, int* incx, PROPACK_CPLX_TYPE* y, int* incy);
double dznrm2_(int* n, PROPACK_CPLX_TYPE* x, int* incx);
void zgemm_(char* transa, char* transb, int* m, int* n, int* k, PROPACK_CPLX_TYPE* alpha, PROPACK_CPLX_TYPE* a, int* lda, PROPACK_CPLX_TYPE* b, int* ldb, PROPACK_CPLX_TYPE* beta, PROPACK_CPLX_TYPE* c, int* ldc);
void zgemv_(char* trans, int* m, int* n, PROPACK_CPLX_TYPE* alpha, PROPACK_CPLX_TYPE* a, int* lda, PROPACK_CPLX_TYPE* x, int* incx, PROPACK_CPLX_TYPE* beta, PROPACK_CPLX_TYPE* y, int* incy);
void zdscal_(int* n, double* da, PROPACK_CPLX_TYPE* zx, int* incx);

// LAPACK

void sbdsdc_(char* uplo, char* compq, int* n, float* d, float* e, float* u, int* ldu, float* vt, int* ldvt, float* q, int* iq, float* work, int* iwork, int* info);
void sbdsqr_(char* uplo, int* n, int* ncvt, int* nru, int* ncc, float* d, float* e, float* vt, int* ldvt, float* u, int* ldu, float* c, int* ldc, float* work, int* info);
void slartg_(float* f, float* g, float* c, float* s, float* r);
void slascl_(char* mtype, int* kl, int* ku, float* cfrom, float* cto, int* m, int* n, float* a, int* lda, int* info);
void slaset_(char* uplo, int* m, int* n, float* alpha, float* beta, float* a, int* lda);

void dbdsdc_(char* uplo, char* compq, int* n, double* d, double* e, double* u, int* ldu, double* vt, int* ldvt, double* q, int* iq, double* work, int* iwork, int* info);
void dbdsqr_(char* uplo, int* n, int* ncvt, int* nru, int* ncc, double* d, double* e, double* vt, int* ldvt, double* u, int* ldu, double* c, int* ldc, double* work, int* info);
void dlartg_(double* f, double* g, double* c, double* s, double* r);
void dlascl_(char* mtype, int* kl, int* ku, double* cfrom, double* cto, int* m, int* n, double* a, int* lda, int* info);
void dlaset_(char* uplo, int* m, int* n, double* alpha, double* beta, double* a, int* lda);

void clarfg_(int* n, PROPACK_CPLXF_TYPE* alpha, PROPACK_CPLXF_TYPE* x, int* incx, PROPACK_CPLXF_TYPE* tau);
void clascl_(char* mtype, int* kl, int* ku, float* cfrom, float* cto, int* m, int* n, PROPACK_CPLXF_TYPE* a, int* lda, int* info);

void zlarfg_(int* n, PROPACK_CPLX_TYPE* alpha, PROPACK_CPLX_TYPE* x, int* incx, PROPACK_CPLX_TYPE* tau);
void zlascl_(char* mtype, int* kl, int* ku, double* cfrom, double* cto, int* m, int* n, PROPACK_CPLX_TYPE* a, int* lda, int* info);


// (c,z)dotc is the complex conjugate dot product of two complex vectors.
// Due some historical reasons, this function can cause segfaults on some
// platforms. Hence implemented here instead of using the BLAS version.
static PROPACK_CPLXF_TYPE
cdotc_(const int* n, const PROPACK_CPLXF_TYPE* restrict x, const int* incx, const PROPACK_CPLXF_TYPE* restrict y, const int* incy)
{
    PROPACK_CPLXF_TYPE result = PROPACK_cplxf(0.0, 0.0);
#ifdef _MSC_VER
    PROPACK_CPLXF_TYPE temp = PROPACK_cplxf(0.0, 0.0);
#endif
    if (*n <= 0) { return result; }
    if ((*incx == 1) && (*incy == 1))
    {
        for (int i = 0; i < *n; i++)
        {
#ifdef _MSC_VER
            temp = _FCmulcc(conjf(x[i]), y[i]);
            result = PROPACK_cplxf(crealf(result) + crealf(temp), cimagf(result) + cimagf(temp));
#else
            result = result + (conjf(x[i]) * y[i]);
#endif
        }

    } else {

        for (int i = 0; i < *n; i++)
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
zdotc_(const int* n, const PROPACK_CPLX_TYPE* restrict x, const int* incx, const PROPACK_CPLX_TYPE* restrict y, const int* incy)
{
    PROPACK_CPLX_TYPE result = PROPACK_cplx(0.0, 0.0);
#ifdef _MSC_VER
    PROPACK_CPLX_TYPE temp = PROPACK_cplx(0.0, 0.0);
#endif
    if (*n <= 0) { return result; }
    if ((*incx == 1) && (*incy == 1))
    {
        for (int i = 0; i < *n; i++)
        {
#ifdef _MSC_VER
            temp = _Cmulcc(conj(x[i]), y[i]);
            result = PROPACK_cplx(creal(result) + creal(temp), cimag(result) + cimag(temp));
#else
            result = result + (conj(x[i]) * y[i]);
#endif
        }

    } else {

        for (int i = 0; i < *n; i++)
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
