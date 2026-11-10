#ifndef _SCIPY_COMMON_ARRAY_UTILS_H
#define _SCIPY_COMMON_ARRAY_UTILS_H

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <math.h>
#include "numpy/arrayobject.h"
#include "scipy_blas_defines.h"
#include "scipy_complex_support.h"



// BLAS and LAPACK functions used
void BLAS_FUNC(saxpy)(CBLAS_INT* n, float* sa, float* sx, CBLAS_INT* incx, float* sy, CBLAS_INT* incy);
void BLAS_FUNC(scopy)(CBLAS_INT* n, float* dx, CBLAS_INT* incx, float* dy, CBLAS_INT* incy);
void BLAS_FUNC(sgees)(char* jobvs, char* sort, CBLAS_INT (*select)(float*, float*), CBLAS_INT* n, float* a, CBLAS_INT* lda, CBLAS_INT* sdim, float* wr, float* wi, float* vs, CBLAS_INT* ldvs, float* work, CBLAS_INT* lwork, CBLAS_INT* bwork, CBLAS_INT* info);
void BLAS_FUNC(sgemm)(char* transa, char* transb, CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* k, float* alpha, float* a, CBLAS_INT* lda, float* b, CBLAS_INT* ldb, float* beta, float* c, CBLAS_INT* ldc);
void BLAS_FUNC(sgemv)(char* trans, CBLAS_INT* m, CBLAS_INT* n, float* alpha, float* a, CBLAS_INT* lda, float* x, CBLAS_INT* incx, float* beta, float* y, CBLAS_INT* incy);
void BLAS_FUNC(sgetrf)(CBLAS_INT* m, CBLAS_INT* n, float* a, CBLAS_INT* lda, CBLAS_INT* ipiv, CBLAS_INT* info);
void BLAS_FUNC(sgetrs)(char* trans, CBLAS_INT* n, CBLAS_INT* nrhs, float* a, CBLAS_INT* lda, CBLAS_INT* ipiv, float* b, CBLAS_INT* ldb, CBLAS_INT* info);
void BLAS_FUNC(slacn2)(CBLAS_INT* n, float* v, float* x, CBLAS_INT* isgn, float* est, CBLAS_INT* kase, CBLAS_INT* isave);
void BLAS_FUNC(slanv2)(float* a, float* b, float* c, float* d, float* rt1r, float* rt1i, float* rt2r, float* rt2i, float* cs, float* sn);
void BLAS_FUNC(sscal)(CBLAS_INT* n, float* sa, float* sx, CBLAS_INT* incx);
void BLAS_FUNC(strsyl)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, float* a, CBLAS_INT* lda, float* b, CBLAS_INT* ldb, float* c, CBLAS_INT* ldc, float* scale, CBLAS_INT* info);
// void BLAS_FUNC(strsyl3)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, float* a, CBLAS_INT* lda, float* b, CBLAS_INT* ldb, float* c, CBLAS_INT* ldc, float* scale, CBLAS_INT* iwork, CBLAS_INT* liwork, float* swork, CBLAS_INT* ldswork, CBLAS_INT* info);

void BLAS_FUNC(daxpy)(CBLAS_INT* n, double* sa, double* sx, CBLAS_INT* incx, double* sy, CBLAS_INT* incy);
void BLAS_FUNC(dcopy)(CBLAS_INT* n, double* dx, CBLAS_INT* incx, double* dy, CBLAS_INT* incy);
void BLAS_FUNC(dgees)(char* jobvs, char* sort, CBLAS_INT (*select)(double*, double*), CBLAS_INT* n, double* a, CBLAS_INT* lda, CBLAS_INT* sdim, double* wr, double* wi, double* vs, CBLAS_INT* ldvs, double* work, CBLAS_INT* lwork, CBLAS_INT* bwork, CBLAS_INT* info);
void BLAS_FUNC(dgemm)(char* transa, char* transb, CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* k, double* alpha, double* a, CBLAS_INT* lda, double* b, CBLAS_INT* ldb, double* beta, double* c, CBLAS_INT* ldc);
void BLAS_FUNC(dgemv)(char* trans, CBLAS_INT* m, CBLAS_INT* n, double* alpha, double* a, CBLAS_INT* lda, double* x, CBLAS_INT* incx, double* beta, double* y, CBLAS_INT* incy);
void BLAS_FUNC(dgetrf)(CBLAS_INT* m, CBLAS_INT* n, double* a, CBLAS_INT* lda, CBLAS_INT* ipiv, CBLAS_INT* info);
void BLAS_FUNC(dgetrs)(char* trans, CBLAS_INT* n, CBLAS_INT* nrhs, double* a, CBLAS_INT* lda, CBLAS_INT* ipiv, double* b, CBLAS_INT* ldb, CBLAS_INT* info);
void BLAS_FUNC(dlacn2)(CBLAS_INT* n, double* v, double* x, CBLAS_INT* isgn, double* est, CBLAS_INT* kase, CBLAS_INT* isave);
void BLAS_FUNC(dlanv2)(double* a, double* b, double* c, double* d, double* rt1r, double* rt1i, double* rt2r, double* rt2i, double* cs, double* sn);
void BLAS_FUNC(dscal)(CBLAS_INT* n, double* sa, double* sx, CBLAS_INT* incx);
void BLAS_FUNC(dtrsyl)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, double* a, CBLAS_INT* lda, double* b, CBLAS_INT* ldb, double* c, CBLAS_INT* ldc, double* scale, CBLAS_INT* info);
// void BLAS_FUNC(dtrsyl3)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, double* a, CBLAS_INT* lda, double* b, CBLAS_INT* ldb, double* c, CBLAS_INT* ldc, double* scale, CBLAS_INT* iwork, CBLAS_INT* liwork, double* swork, CBLAS_INT* ldswork, CBLAS_INT* info);

void BLAS_FUNC(caxpy)(CBLAS_INT* n, scipy_c64* sa, scipy_c64* sx, CBLAS_INT* incx, scipy_c64* sy, CBLAS_INT* incy);
void BLAS_FUNC(ccopy)(CBLAS_INT* n, scipy_c64* dx, CBLAS_INT* incx, scipy_c64* dy, CBLAS_INT* incy);
void BLAS_FUNC(cgees)(char* jobvs, char* sort, CBLAS_INT (*select)(scipy_c64*), CBLAS_INT* n, scipy_c64* a, CBLAS_INT* lda, CBLAS_INT* sdim, scipy_c64* w, scipy_c64* vs, CBLAS_INT* ldvs, scipy_c64* work, CBLAS_INT* lwork, float* rwork, CBLAS_INT* bwork, CBLAS_INT* info);
void BLAS_FUNC(cgemm)(char* transa, char* transb, CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* k, scipy_c64* alpha, scipy_c64* a, CBLAS_INT* lda, scipy_c64* b, CBLAS_INT* ldb, scipy_c64* beta, scipy_c64* c, CBLAS_INT* ldc);
void BLAS_FUNC(cgemv)(char* trans, CBLAS_INT* m, CBLAS_INT* n, scipy_c64* alpha, scipy_c64* a, CBLAS_INT* lda, scipy_c64* x, CBLAS_INT* incx, scipy_c64* beta, scipy_c64* y, CBLAS_INT* incy);
void BLAS_FUNC(cgetrf)(CBLAS_INT* m, CBLAS_INT* n, scipy_c64* a, CBLAS_INT* lda, CBLAS_INT* ipiv, CBLAS_INT* info);
void BLAS_FUNC(cgetrs)(char* trans, CBLAS_INT* n, CBLAS_INT* nrhs, scipy_c64* a, CBLAS_INT* lda, CBLAS_INT* ipiv, scipy_c64* b, CBLAS_INT* ldb, CBLAS_INT* info);
void BLAS_FUNC(clacn2)(CBLAS_INT* n, scipy_c64* v, scipy_c64* x, float* est, CBLAS_INT* kase, CBLAS_INT* isave);
void BLAS_FUNC(crot)(CBLAS_INT* n, scipy_c64* cx, CBLAS_INT* incx, scipy_c64* cy, CBLAS_INT* incy, float* c, scipy_c64* s);
void BLAS_FUNC(csscal)(CBLAS_INT* n, float* sa, scipy_c64* sx, CBLAS_INT* incx);
void BLAS_FUNC(ctrsyl)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, scipy_c64* a, CBLAS_INT* lda, scipy_c64* b, CBLAS_INT* ldb, scipy_c64* c, CBLAS_INT* ldc, float* scale, CBLAS_INT* info);
// void BLAS_FUNC(ctrsyl3)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, scipy_c64* a, CBLAS_INT* lda, scipy_c64* b, CBLAS_INT* ldb, scipy_c64* c, CBLAS_INT* ldc, float* scale, float* swork, CBLAS_INT* ldswork, CBLAS_INT* info);

void BLAS_FUNC(zaxpy)(CBLAS_INT* n, scipy_c128* sa, scipy_c128* sx, CBLAS_INT* incx, scipy_c128* sy, CBLAS_INT* incy);
void BLAS_FUNC(zcopy)(CBLAS_INT* n, scipy_c128* dx, CBLAS_INT* incx, scipy_c128* dy, CBLAS_INT* incy);
void BLAS_FUNC(zgees)(char* jobvs, char* sort, CBLAS_INT (*select)(scipy_c128*), CBLAS_INT* n, scipy_c128* a, CBLAS_INT* lda, CBLAS_INT* sdim, scipy_c128* w, scipy_c128* vs, CBLAS_INT* ldvs, scipy_c128* work, CBLAS_INT* lwork, double* rwork, CBLAS_INT* bwork, CBLAS_INT* info);
void BLAS_FUNC(zgemm)(char* transa, char* transb, CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* k, scipy_c128* alpha, scipy_c128* a, CBLAS_INT* lda, scipy_c128* b, CBLAS_INT* ldb, scipy_c128* beta, scipy_c128* c, CBLAS_INT* ldc);
void BLAS_FUNC(zgemv)(char* trans, CBLAS_INT* m, CBLAS_INT* n, scipy_c128* alpha, scipy_c128* a, CBLAS_INT* lda, scipy_c128* x, CBLAS_INT* incx, scipy_c128* beta, scipy_c128* y, CBLAS_INT* incy);
void BLAS_FUNC(zgetrf)(CBLAS_INT* m, CBLAS_INT* n, scipy_c128* a, CBLAS_INT* lda, CBLAS_INT* ipiv, CBLAS_INT* info);
void BLAS_FUNC(zgetrs)(char* trans, CBLAS_INT* n, CBLAS_INT* nrhs, scipy_c128* a, CBLAS_INT* lda, CBLAS_INT* ipiv, scipy_c128* b, CBLAS_INT* ldb, CBLAS_INT* info);
void BLAS_FUNC(zlacn2)(CBLAS_INT* n, scipy_c128* v, scipy_c128* x, double* est, CBLAS_INT* kase, CBLAS_INT* isave);
void BLAS_FUNC(zrot)(CBLAS_INT* n, scipy_c128* cx, CBLAS_INT* incx, scipy_c128* cy, CBLAS_INT* incy, double* c, scipy_c128* s);
void BLAS_FUNC(zdscal)(CBLAS_INT* n, double* sa, scipy_c128* sx, CBLAS_INT* incx);
void BLAS_FUNC(ztrsyl)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, scipy_c128* a, CBLAS_INT* lda, scipy_c128* b, CBLAS_INT* ldb, scipy_c128* c, CBLAS_INT* ldc, double* scale, CBLAS_INT* info);
// void BLAS_FUNC(ztrsyl3)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, scipy_c128* a, CBLAS_INT* lda, scipy_c128* b, CBLAS_INT* ldb, scipy_c128* c, CBLAS_INT* ldc, double* scale, double* swork, CBLAS_INT* ldswork, CBLAS_INT* info);

/*
 *  These swap functions are used to convert a C-contiguous n*n array to an F-
 *  contiguous one and back. Recursively halves on the longer dimension each
 *  time until it reaches to a small enough piece such that tries to benefit the
 *  locality of the data. Last letter denotes the LAPACK flavor s,d,c,z.
 *
 *  There is not much science behind the magical values but empirically they seem
 *  to be at the sweetspot for the tradeoff between recursion overhead and data
 *  locality.
 */
static inline void
swap_cf_s(float* restrict src, float* restrict dst, const Py_ssize_t r, const Py_ssize_t c, const Py_ssize_t n)
{
    Py_ssize_t i, j, ith_row, r2, c2;
    float *bb = dst;
    float *aa = src;
    if (r < 32) {
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
swap_cf_c(scipy_c64* restrict src, scipy_c64* restrict dst, const Py_ssize_t r, const Py_ssize_t c, const Py_ssize_t n)
{
    Py_ssize_t i, j, ith_row, r2, c2;
    scipy_c64 *bb = dst;
    scipy_c64 *aa = src;
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
swap_cf_z(scipy_c128* restrict src, scipy_c128* restrict dst, const Py_ssize_t r, const Py_ssize_t c, const Py_ssize_t n)
{
    Py_ssize_t i, j, ith_row, r2, c2;
    scipy_c128 *bb = dst;
    scipy_c128 *aa = src;
    if (r < 8) {
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
 *  returns: an int value indicating it is Schur (1) or not (0).
 *
 * For complex arrays, being Schur is equivalent to being upper triangular.
 */
static int
isschurf(const float* restrict data, const Py_ssize_t n)
{
    int in_block = 0;
    float block_diag, block_offdiag;

    for (Py_ssize_t k = 0; k < n - 1; k++) {
        // Below-subdiagonal entries in column k must be zero
        for (Py_ssize_t i = k + 2; i < n; i++) {
            if (data[k*n + i] != 0.0f) { return 0; }
        }

        if (in_block) {
            // Second column of 2x2 block: validate using current column only
            if (block_diag != data[k*n + k]) { return 0; }
            if (block_offdiag * data[k*n + k - 1] >= 0.0f) { return 0; }
            if (data[k*n + k + 1] != 0.0f) { return 0; }
            in_block = 0;
        } else if (data[k*n + k + 1] != 0.0f) {
            // First column of 2x2 block: save and mark
            block_diag = data[k*n + k];
            block_offdiag = data[k*n + k + 1];
            in_block = 1;
        }
    }
    // Validate if last block was 2x2.
    if (in_block) {
        Py_ssize_t k = n - 1;
        if (block_diag != data[k*n + k]) { return 0; }
        if (block_offdiag * data[k*n + k - 1] >= 0.0f) { return 0; }
    }
    return 1;
}


static int
isschur(const double* restrict data, const Py_ssize_t n)
{
    int in_block = 0;
    double block_diag, block_offdiag;

    for (Py_ssize_t k = 0; k < n - 1; k++) {
        // Below-subdiagonal entries in column k must be zero
        for (Py_ssize_t i = k + 2; i < n; i++) {
            if (data[k*n + i] != 0.0) { return 0; }
        }

        if (in_block) {
            // Second column of 2x2 block: validate using current column only
            if (block_diag != data[k*n + k]) { return 0; }
            if (block_offdiag * data[k*n + k - 1] >= 0.0) { return 0; }
            if (data[k*n + k + 1] != 0.0) { return 0; }
            in_block = 0;
        } else if (data[k*n + k + 1] != 0.0) {
            // First column of 2x2 block: save and mark
            block_diag = data[k*n + k];
            block_offdiag = data[k*n + k + 1];
            in_block = 1;
        }
    }
    // Validate if last block was 2x2.
    if (in_block) {
        Py_ssize_t k = n - 1;
        if (block_diag != data[k*n + k]) { return 0; }
        if (block_offdiag * data[k*n + k - 1] >= 0.0) { return 0; }
    }
    return 1;
}

#endif
