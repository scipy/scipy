#ifndef _SCIPY_COMMON_ARRAY_UTILS_H
#define _SCIPY_COMMON_ARRAY_UTILS_H

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <math.h>
#include "numpy/arrayobject.h"

#if defined(_MSC_VER)
    #include <complex.h>
    #define SCIPY_Z _Dcomplex
    #define SCIPY_C _Fcomplex
    #define CPLX_Z(real, imag) (_Cbuild(real, imag))
    #define CPLX_C(real, imag) (_FCbuild(real, imag))
#else
    #include <complex.h>
    #define SCIPY_Z double complex
    #define SCIPY_C float complex
    #define CPLX_Z(real, imag) (real + imag*I)
    #define CPLX_C(real, imag) (real + imag*I)
#endif

// BLAS and LAPACK functions used
void saxpy_(int* n, float* sa, float* sx, int* incx, float* sy, int* incy);
void scopy_(int* n, float* dx, int* incx, float* dy, int* incy);
void sgees_(char* jobvs, char* sort, int (*select)(float*, float*), int* n, float* a, int* lda, int* sdim, float* wr, float* wi, float* vs, int* ldvs, float* work, int* lwork, int* bwork, int* info);
void sgemm_(char* transa, char* transb, int* m, int* n, int* k, float* alpha, float* a, int* lda, float* b, int* ldb, float* beta, float* c, int* ldc);
void sgemv_(char* trans, int* m, int* n, float* alpha, float* a, int* lda, float* x, int* incx, float* beta, float* y, int* incy);
void sgetrf_(int* m, int* n, float* a, int* lda, int* ipiv, int* info);
void sgetrs_(char* trans, int* n, int* nrhs, float* a, int* lda, int* ipiv, float* b, int* ldb, int* info);
void slacn2_(int* n, float* v, float* x, int* isgn, float* est, int* kase, int* isave);
void slanv2_(float* a, float* b, float* c, float* d, float* rt1r, float* rt1i, float* rt2r, float* rt2i, float* cs, float* sn);
void sscal_(int* n, float* sa, float* sx, int* incx);
void strsyl_(char* trana, char* tranb, int* isgn, int* m, int* n, float* a, int* lda, float* b, int* ldb, float* c, int* ldc, float* scale, int* info);
// void strsyl3_(char* trana, char* tranb, int* isgn, int* m, int* n, float* a, int* lda, float* b, int* ldb, float* c, int* ldc, float* scale, int* iwork, int* liwork, float* swork, int* ldswork, int* info);

void daxpy_(int* n, double* sa, double* sx, int* incx, double* sy, int* incy);
void dcopy_(int* n, double* dx, int* incx, double* dy, int* incy);
void dgees_(char* jobvs, char* sort, int (*select)(double*, double*), int* n, double* a, int* lda, int* sdim, double* wr, double* wi, double* vs, int* ldvs, double* work, int* lwork, int* bwork, int* info);
void dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);
void dgemv_(char* trans, int* m, int* n, double* alpha, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
void dgetrs_(char* trans, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
void dlacn2_(int* n, double* v, double* x, int* isgn, double* est, int* kase, int* isave);
void dlanv2_(double* a, double* b, double* c, double* d, double* rt1r, double* rt1i, double* rt2r, double* rt2i, double* cs, double* sn);
void dscal_(int* n, double* sa, double* sx, int* incx);
void dtrsyl_(char* trana, char* tranb, int* isgn, int* m, int* n, double* a, int* lda, double* b, int* ldb, double* c, int* ldc, double* scale, int* info);
// void dtrsyl3_(char* trana, char* tranb, int* isgn, int* m, int* n, double* a, int* lda, double* b, int* ldb, double* c, int* ldc, double* scale, int* iwork, int* liwork, double* swork, int* ldswork, int* info);

void caxpy_(int* n, SCIPY_C* sa, SCIPY_C* sx, int* incx, SCIPY_C* sy, int* incy);
void ccopy_(int* n, SCIPY_C* dx, int* incx, SCIPY_C* dy, int* incy);
void cgees_(char* jobvs, char* sort, int (*select)(SCIPY_C), int* n, SCIPY_C* a, int* lda, int* sdim, SCIPY_C* w, SCIPY_C* vs, int* ldvs, SCIPY_C* work, int* lwork, float* rwork, int* bwork, int* info);
void cgemm_(char* transa, char* transb, int* m, int* n, int* k, SCIPY_C* alpha, SCIPY_C* a, int* lda, SCIPY_C* b, int* ldb, SCIPY_C* beta, SCIPY_C* c, int* ldc);
void cgemv_(char* trans, int* m, int* n, SCIPY_C* alpha, SCIPY_C* a, int* lda, SCIPY_C* x, int* incx, SCIPY_C* beta, SCIPY_C* y, int* incy);
void cgetrf_(int* m, int* n, SCIPY_C* a, int* lda, int* ipiv, int* info);
void cgetrs_(char* trans, int* n, int* nrhs, SCIPY_C* a, int* lda, int* ipiv, SCIPY_C* b, int* ldb, int* info);
void clacn2_(int* n, SCIPY_C* v, SCIPY_C* x, float* est, int* kase, int* isave);
void crot_(int* n, SCIPY_C* cx, int* incx, SCIPY_C* cy, int* incy, float* c, SCIPY_C* s);
void csscal_(int* n, float* sa, SCIPY_C* sx, int* incx);
void ctrsyl_(char* trana, char* tranb, int* isgn, int* m, int* n, SCIPY_C* a, int* lda, SCIPY_C* b, int* ldb, SCIPY_C* c, int* ldc, float* scale, int* info);
// void ctrsyl3_(char* trana, char* tranb, int* isgn, int* m, int* n, SCIPY_C* a, int* lda, SCIPY_C* b, int* ldb, SCIPY_C* c, int* ldc, float* scale, float* swork, int* ldswork, int* info);

void zaxpy_(int* n, SCIPY_Z* sa, SCIPY_Z* sx, int* incx, SCIPY_Z* sy, int* incy);
void zcopy_(int* n, SCIPY_Z* dx, int* incx, SCIPY_Z* dy, int* incy);
void zgees_(char* jobvs, char* sort, int (*select)(SCIPY_Z), int* n, SCIPY_Z* a, int* lda, int* sdim, SCIPY_Z* w, SCIPY_Z* vs, int* ldvs, SCIPY_Z* work, int* lwork, double* rwork, int* bwork, int* info);
void zgemm_(char* transa, char* transb, int* m, int* n, int* k, SCIPY_Z* alpha, SCIPY_Z* a, int* lda, SCIPY_Z* b, int* ldb, SCIPY_Z* beta, SCIPY_Z* c, int* ldc);
void zgemv_(char* trans, int* m, int* n, SCIPY_Z* alpha, SCIPY_Z* a, int* lda, SCIPY_Z* x, int* incx, SCIPY_Z* beta, SCIPY_Z* y, int* incy);
void zgetrf_(int* m, int* n, SCIPY_Z* a, int* lda, int* ipiv, int* info);
void zgetrs_(char* trans, int* n, int* nrhs, SCIPY_Z* a, int* lda, int* ipiv, SCIPY_Z* b, int* ldb, int* info);
void zlacn2_(int* n, SCIPY_Z* v, SCIPY_Z* x, double* est, int* kase, int* isave);
void zrot_(int* n, SCIPY_Z* cx, int* incx, SCIPY_Z* cy, int* incy, double* c, SCIPY_Z* s);
void zdscal_(int* n, double* sa, SCIPY_Z* sx, int* incx);
void ztrsyl_(char* trana, char* tranb, int* isgn, int* m, int* n, SCIPY_Z* a, int* lda, SCIPY_Z* b, int* ldb, SCIPY_Z* c, int* ldc, double* scale, int* info);
// void ztrsyl3_(char* trana, char* tranb, int* isgn, int* m, int* n, SCIPY_Z* a, int* lda, SCIPY_Z* b, int* ldb, SCIPY_Z* c, int* ldc, double* scale, double* swork, int* ldswork, int* info);

/*
 *  These swap functions are used to convert a C-contiguous n*n array to an F-
 *  contiguous one and back. Recursively halves on the longer dimension each
 *  time until it reaches to a small enough piece such that tries to benefit the
 *  locality of the data. Last letter denotes the LAPACK flavor s,d,c,z.
 */
static inline void
swap_cf_s(float* restrict src, float* restrict dst, const Py_ssize_t r, const Py_ssize_t c, const Py_ssize_t n)
{
    Py_ssize_t i, j, ith_row, r2, c2;
    float *bb = dst;
    float *aa = src;
    if (r < 16) {
        for (j = 0; j < c; j++)
        {
            ith_row = 0;
            for (i = 0; i < r; i++) {
                bb[ith_row] = aa[i];
                ith_row += n;
            }
            aa += n;
            bb += 1;
        }
    } else {
        // If tall
        if (r > c)
        {
            r2 = r/2;
            swap_cf_s(src, dst, r2, c, n);
            swap_cf_s(src + r2, dst+(r2)*n, r-r2, c, n);
        } else {  // Nope
            c2 = c/2;
            swap_cf_s(src, dst, r, c2, n);
            swap_cf_s(src+(c2)*n, dst+c2, r, c-c2, n);
        }
    }
}


static inline void
swap_cf_d(double* restrict src, double* restrict dst, const Py_ssize_t r, const Py_ssize_t c, const Py_ssize_t n)
{
    Py_ssize_t i, j, ith_row, r2, c2;
    double *bb = dst;
    double *aa = src;
    if (r < 16) {
        for (j = 0; j < c; j++)
        {
            ith_row = 0;
            for (i = 0; i < r; i++) {
                bb[ith_row] = aa[i];
                ith_row += n;
            }
            aa += n;
            bb += 1;
        }
    } else {
        // If tall
        if (r > c)
        {
            r2 = r/2;
            swap_cf_d(src, dst, r2, c, n);
            swap_cf_d(src + r2, dst+(r2)*n, r-r2, c, n);
        } else {  // Nope
            c2 = c/2;
            swap_cf_d(src, dst, r, c2, n);
            swap_cf_d(src+(c2)*n, dst+c2, r, c-c2, n);
        }
    }
}


static inline void
swap_cf_c(SCIPY_C* restrict src, SCIPY_C* restrict dst, const Py_ssize_t r, const Py_ssize_t c, const Py_ssize_t n)
{
    Py_ssize_t i, j, ith_row, r2, c2;
    SCIPY_C *bb = dst;
    SCIPY_C *aa = src;
    if (r < 16) {
        for (j = 0; j < c; j++)
        {
            ith_row = 0;
            for (i = 0; i < r; i++) {
                bb[ith_row] = aa[i];
                ith_row += n;
            }
            aa += n;
            bb += 1;
        }
    } else {
        // If tall
        if (r > c)
        {
            r2 = r/2;
            swap_cf_c(src, dst, r2, c, n);
            swap_cf_c(src + r2, dst+(r2)*n, r-r2, c, n);
        } else {  // Nope
            c2 = c/2;
            swap_cf_c(src, dst, r, c2, n);
            swap_cf_c(src+(c2)*n, dst+c2, r, c-c2, n);
        }
    }
}


static inline void
swap_cf_z(SCIPY_Z* restrict src, SCIPY_Z* restrict dst, const Py_ssize_t r, const Py_ssize_t c, const Py_ssize_t n)
{
    Py_ssize_t i, j, ith_row, r2, c2;
    SCIPY_Z *bb = dst;
    SCIPY_Z *aa = src;
    if (r < 16) {
        for (j = 0; j < c; j++)
        {
            ith_row = 0;
            for (i = 0; i < r; i++) {
                bb[ith_row] = aa[i];
                ith_row += n;
            }
            aa += n;
            bb += 1;
        }
    } else {
        // If tall
        if (r > c)
        {
            r2 = r/2;
            swap_cf_z(src, dst, r2, c, n);
            swap_cf_z(src + r2, dst+(r2)*n, r-r2, c, n);
        } else {  // Nope
            c2 = c/2;
            swap_cf_z(src, dst, r, c2, n);
            swap_cf_z(src+(c2)*n, dst+c2, r, c-c2, n);
        }
    }
}

/*
 * Function:  isschur/isschurf
 * ---------------------------
 * checks whether a square array is (quasi)upper triangular and if any 2x2 blocks
 * are present whether these block are in the form of [[a, b], [c, a]] where b*c < 0.
 * This is to check that the 2x2 blocks have actually complex eigenvalues needed
 * for various algorithms.
 *
 *  data: the float/double array to be checked (in Fortran-contiguous layout)
 *  n: the size of the array
 *
 *  returns: a int value indicating it is Schur (1) or not (0).
 *
 * For complex arrays, being Schur is equivalent to being upper triangular.
 */
int
isschurf(const float* restrict data, const Py_ssize_t n)
{
    float prev_a, prev_b;
    int isSchur = 1;
    int expect_zero = 0;
    for (Py_ssize_t j = 0; j < n; j++) {
        if (data[j*n + j+1] == 0.0f) {
            // Subdiagonal is zero
            if (expect_zero) {
                // Previous column had a non-zero entry and this marks the 2x2
                // block. Check the [[a, b], [c, a]] form with b*c < 0.
                if ((prev_a == data[j*n + j]) && (prev_b*data[j*n + j-1] < 0.0))
                {
                    expect_zero = 0; // All good, reset the flag
                } else {
                    // Broken 2x2 form
                    return 0;
                }
            }
            expect_zero = 0;
        } else {
            // Subdiagonal is non-zero. Either first column of 2x2 or not Schur
            if (expect_zero)
            {
                // Already nonzero in the previous col so not Schur.
                return 0;
            } else {
                // First column of 2x2 block
                prev_a = data[j*n + j];
                prev_b = data[j*n + j+1];
                expect_zero = 1; // Next should be zero
            }
        }

        // Check the remaining entries in the col for zero.
        if (j < n - 2)
        {
            for (Py_ssize_t i = j+2; i < n; i++) {
                if (data[j*n + i] != 0.0f) {
                    // Nonzero value below the subdiagonal
                    return 0;
                }
            }
        }
        if (!isSchur) { break; }
    }
    return isSchur;
}


int
isschur(const double* restrict data, const Py_ssize_t n)
{
    double prev_a, prev_b;
    int isSchur = 1;
    int expect_zero = 0;
    for (Py_ssize_t j = 0; j < n; j++) {
        if (data[j*n + j+1] == 0.0) {
            if (expect_zero) {
                if ((prev_a == data[j*n + j]) && (prev_b*data[j*n + j-1] < 0.0))
                {
                    expect_zero = 0;
                } else {
                    return 0;
                }
            }
            expect_zero = 0;
        } else {
            if (expect_zero)
            {
                return 0;
            } else {
                prev_a = data[j*n + j];
                prev_b = data[j*n + j+1];
                expect_zero = 1;
            }
        }
        if (j < n - 2)
        {
            for (Py_ssize_t i = j+2; i < n; i++) {
                if (data[j*n + i] != 0.0) {
                    return 0;
                }
            }
        }
        if (!isSchur) { break; }
    }
    return isSchur;
}

#endif