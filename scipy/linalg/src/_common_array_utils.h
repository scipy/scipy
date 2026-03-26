#ifndef _SCIPY_COMMON_ARRAY_UTILS_H
#define _SCIPY_COMMON_ARRAY_UTILS_H

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <math.h>
#include "numpy/arrayobject.h"
#include "scipy_blas_defines.h"

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

// AVX2 dispatch requires __builtin_cpu_supports and __attribute__((target)),
// available on GCC and Clang on non-Windows x86-64. Clang on Windows uses
// lld-link which cannot resolve __cpu_model from compiler-rt.
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(_WIN32)
#include <immintrin.h>
#define SCIPY_HAVE_AVX2_TARGET 1
#endif

// BLAS and LAPACK functions used
void BLAS_FUNC(saxpy)(CBLAS_INT* n, float* sa, float* sx, CBLAS_INT* incx, float* sy, CBLAS_INT* incy);
void BLAS_FUNC(scopy)(CBLAS_INT* n, float* dx, CBLAS_INT* incx, float* dy, CBLAS_INT* incy);
void BLAS_FUNC(sgees)(char* jobvs, char* sort, int (*select)(float*, float*), CBLAS_INT* n, float* a, CBLAS_INT* lda, CBLAS_INT* sdim, float* wr, float* wi, float* vs, CBLAS_INT* ldvs, float* work, CBLAS_INT* lwork, CBLAS_INT* bwork, CBLAS_INT* info);
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
void BLAS_FUNC(dgees)(char* jobvs, char* sort, int (*select)(double*, double*), CBLAS_INT* n, double* a, CBLAS_INT* lda, CBLAS_INT* sdim, double* wr, double* wi, double* vs, CBLAS_INT* ldvs, double* work, CBLAS_INT* lwork, CBLAS_INT* bwork, CBLAS_INT* info);
void BLAS_FUNC(dgemm)(char* transa, char* transb, CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* k, double* alpha, double* a, CBLAS_INT* lda, double* b, CBLAS_INT* ldb, double* beta, double* c, CBLAS_INT* ldc);
void BLAS_FUNC(dgemv)(char* trans, CBLAS_INT* m, CBLAS_INT* n, double* alpha, double* a, CBLAS_INT* lda, double* x, CBLAS_INT* incx, double* beta, double* y, CBLAS_INT* incy);
void BLAS_FUNC(dgetrf)(CBLAS_INT* m, CBLAS_INT* n, double* a, CBLAS_INT* lda, CBLAS_INT* ipiv, CBLAS_INT* info);
void BLAS_FUNC(dgetrs)(char* trans, CBLAS_INT* n, CBLAS_INT* nrhs, double* a, CBLAS_INT* lda, CBLAS_INT* ipiv, double* b, CBLAS_INT* ldb, CBLAS_INT* info);
void BLAS_FUNC(dlacn2)(CBLAS_INT* n, double* v, double* x, CBLAS_INT* isgn, double* est, CBLAS_INT* kase, CBLAS_INT* isave);
void BLAS_FUNC(dlanv2)(double* a, double* b, double* c, double* d, double* rt1r, double* rt1i, double* rt2r, double* rt2i, double* cs, double* sn);
void BLAS_FUNC(dscal)(CBLAS_INT* n, double* sa, double* sx, CBLAS_INT* incx);
void BLAS_FUNC(dtrsyl)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, double* a, CBLAS_INT* lda, double* b, CBLAS_INT* ldb, double* c, CBLAS_INT* ldc, double* scale, CBLAS_INT* info);
// void BLAS_FUNC(dtrsyl3)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, double* a, CBLAS_INT* lda, double* b, CBLAS_INT* ldb, double* c, CBLAS_INT* ldc, double* scale, CBLAS_INT* iwork, CBLAS_INT* liwork, double* swork, CBLAS_INT* ldswork, CBLAS_INT* info);

void BLAS_FUNC(caxpy)(CBLAS_INT* n, SCIPY_C* sa, SCIPY_C* sx, CBLAS_INT* incx, SCIPY_C* sy, CBLAS_INT* incy);
void BLAS_FUNC(ccopy)(CBLAS_INT* n, SCIPY_C* dx, CBLAS_INT* incx, SCIPY_C* dy, CBLAS_INT* incy);
void BLAS_FUNC(cgees)(char* jobvs, char* sort, int (*select)(SCIPY_C), CBLAS_INT* n, SCIPY_C* a, CBLAS_INT* lda, CBLAS_INT* sdim, SCIPY_C* w, SCIPY_C* vs, CBLAS_INT* ldvs, SCIPY_C* work, CBLAS_INT* lwork, float* rwork, CBLAS_INT* bwork, CBLAS_INT* info);
void BLAS_FUNC(cgemm)(char* transa, char* transb, CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* k, SCIPY_C* alpha, SCIPY_C* a, CBLAS_INT* lda, SCIPY_C* b, CBLAS_INT* ldb, SCIPY_C* beta, SCIPY_C* c, CBLAS_INT* ldc);
void BLAS_FUNC(cgemv)(char* trans, CBLAS_INT* m, CBLAS_INT* n, SCIPY_C* alpha, SCIPY_C* a, CBLAS_INT* lda, SCIPY_C* x, CBLAS_INT* incx, SCIPY_C* beta, SCIPY_C* y, CBLAS_INT* incy);
void BLAS_FUNC(cgetrf)(CBLAS_INT* m, CBLAS_INT* n, SCIPY_C* a, CBLAS_INT* lda, CBLAS_INT* ipiv, CBLAS_INT* info);
void BLAS_FUNC(cgetrs)(char* trans, CBLAS_INT* n, CBLAS_INT* nrhs, SCIPY_C* a, CBLAS_INT* lda, CBLAS_INT* ipiv, SCIPY_C* b, CBLAS_INT* ldb, CBLAS_INT* info);
void BLAS_FUNC(clacn2)(CBLAS_INT* n, SCIPY_C* v, SCIPY_C* x, float* est, CBLAS_INT* kase, CBLAS_INT* isave);
void BLAS_FUNC(crot)(CBLAS_INT* n, SCIPY_C* cx, CBLAS_INT* incx, SCIPY_C* cy, CBLAS_INT* incy, float* c, SCIPY_C* s);
void BLAS_FUNC(csscal)(CBLAS_INT* n, float* sa, SCIPY_C* sx, CBLAS_INT* incx);
void BLAS_FUNC(ctrsyl)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, SCIPY_C* a, CBLAS_INT* lda, SCIPY_C* b, CBLAS_INT* ldb, SCIPY_C* c, CBLAS_INT* ldc, float* scale, CBLAS_INT* info);
// void BLAS_FUNC(ctrsyl3)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, SCIPY_C* a, CBLAS_INT* lda, SCIPY_C* b, CBLAS_INT* ldb, SCIPY_C* c, CBLAS_INT* ldc, float* scale, float* swork, CBLAS_INT* ldswork, CBLAS_INT* info);

void BLAS_FUNC(zaxpy)(CBLAS_INT* n, SCIPY_Z* sa, SCIPY_Z* sx, CBLAS_INT* incx, SCIPY_Z* sy, CBLAS_INT* incy);
void BLAS_FUNC(zcopy)(CBLAS_INT* n, SCIPY_Z* dx, CBLAS_INT* incx, SCIPY_Z* dy, CBLAS_INT* incy);
void BLAS_FUNC(zgees)(char* jobvs, char* sort, int (*select)(SCIPY_Z), CBLAS_INT* n, SCIPY_Z* a, CBLAS_INT* lda, CBLAS_INT* sdim, SCIPY_Z* w, SCIPY_Z* vs, CBLAS_INT* ldvs, SCIPY_Z* work, CBLAS_INT* lwork, double* rwork, CBLAS_INT* bwork, CBLAS_INT* info);
void BLAS_FUNC(zgemm)(char* transa, char* transb, CBLAS_INT* m, CBLAS_INT* n, CBLAS_INT* k, SCIPY_Z* alpha, SCIPY_Z* a, CBLAS_INT* lda, SCIPY_Z* b, CBLAS_INT* ldb, SCIPY_Z* beta, SCIPY_Z* c, CBLAS_INT* ldc);
void BLAS_FUNC(zgemv)(char* trans, CBLAS_INT* m, CBLAS_INT* n, SCIPY_Z* alpha, SCIPY_Z* a, CBLAS_INT* lda, SCIPY_Z* x, CBLAS_INT* incx, SCIPY_Z* beta, SCIPY_Z* y, CBLAS_INT* incy);
void BLAS_FUNC(zgetrf)(CBLAS_INT* m, CBLAS_INT* n, SCIPY_Z* a, CBLAS_INT* lda, CBLAS_INT* ipiv, CBLAS_INT* info);
void BLAS_FUNC(zgetrs)(char* trans, CBLAS_INT* n, CBLAS_INT* nrhs, SCIPY_Z* a, CBLAS_INT* lda, CBLAS_INT* ipiv, SCIPY_Z* b, CBLAS_INT* ldb, CBLAS_INT* info);
void BLAS_FUNC(zlacn2)(CBLAS_INT* n, SCIPY_Z* v, SCIPY_Z* x, double* est, CBLAS_INT* kase, CBLAS_INT* isave);
void BLAS_FUNC(zrot)(CBLAS_INT* n, SCIPY_Z* cx, CBLAS_INT* incx, SCIPY_Z* cy, CBLAS_INT* incy, double* c, SCIPY_Z* s);
void BLAS_FUNC(zdscal)(CBLAS_INT* n, double* sa, SCIPY_Z* sx, CBLAS_INT* incx);
void BLAS_FUNC(ztrsyl)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, SCIPY_Z* a, CBLAS_INT* lda, SCIPY_Z* b, CBLAS_INT* ldb, SCIPY_Z* c, CBLAS_INT* ldc, double* scale, CBLAS_INT* info);
// void BLAS_FUNC(ztrsyl3)(char* trana, char* tranb, CBLAS_INT* isgn, CBLAS_INT* m, CBLAS_INT* n, SCIPY_Z* a, CBLAS_INT* lda, SCIPY_Z* b, CBLAS_INT* ldb, SCIPY_Z* c, CBLAS_INT* ldc, double* scale, double* swork, CBLAS_INT* ldswork, CBLAS_INT* info);


/**
 *
 * Bandwidth detection for general matrices. There are two implementations
 * provided below, a scalar scanning implementation and an AVX2-optimized
 * implementation. The AVX2 implementation is selected at load time when
 * CPU support is detected, otherwise the scalar version is used as a
 * fallback. The AVX2 path is available on GCC and Clang on x86-64; MSVC
 * and other platforms always use the scalar fallback.
 *
 */

// bandwidth_s

static inline void
bandwidth_s_scalar(float* restrict data, npy_intp n, npy_intp m,
                   npy_intp* lower_band, npy_intp* upper_band)
{
    Py_ssize_t lb = 0, ub = 0;
    for (Py_ssize_t r = n-1; r > 0; r--) {
        Py_ssize_t limit = r - lb;
        if (limit > m) { limit = m; }
        for (Py_ssize_t c = 0; c < limit; c++) {
            if (data[r*m + c] != 0.0f) { lb = r - c; break; }
        }
        if (r <= lb) { break; }
    }
    for (Py_ssize_t r = 0; r < n-1; r++) {
        for (Py_ssize_t c = m-1; c > r + ub; c--) {
            if (data[r*m + c] != 0.0f) { ub = c - r; break; }
        }
        if (r + ub + 1 > m) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}

#ifdef SCIPY_HAVE_AVX2_TARGET
__attribute__((target("avx2")))
static inline void
bandwidth_s_avx2(float* restrict data, npy_intp n, npy_intp m,
                 npy_intp* lower_band, npy_intp* upper_band)
{
    const __m256 zero = _mm256_setzero_ps();
    Py_ssize_t lb = 0, ub = 0;

    for (Py_ssize_t r = n-1; r > 0; r--) {
        Py_ssize_t limit = r - lb;
        if (limit > m) { limit = m; }
        if (limit >= 8) {
            Py_ssize_t c = 0;
            for (; c + 7 < limit; c += 8) {
                __m256 v = _mm256_loadu_ps(&data[r*m + c]);
                __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_ps(cmp);
                if (mask != 0xFF) {
                    lb = r - (c + __builtin_ctz(~mask));
                    goto lb_done_s;
                }
            }
            if (c < limit) {
                __m256 v = _mm256_loadu_ps(&data[r*m + limit - 8]);
                __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_ps(cmp);
                if (mask != 0xFF) {
                    lb = r - (limit - 8 + __builtin_ctz(~mask));
                    goto lb_done_s;
                }
            }
        } else {
            for (Py_ssize_t c = 0; c < limit; c++) {
                if (data[r*m + c] != 0.0f) { lb = r - c; goto lb_done_s; }
            }
        }
        lb_done_s:
        if (r <= lb) { break; }
    }

    for (Py_ssize_t r = 0; r < n-1; r++) {
        Py_ssize_t limit = r + ub;
        Py_ssize_t span = m - 1 - limit;
        if (span >= 8) {
            Py_ssize_t c = m - 1;
            for (; c - 7 > limit; c -= 8) {
                __m256 v = _mm256_loadu_ps(&data[r*m + c - 7]);
                __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_ps(cmp);
                if (mask != 0xFF) {
                    int last = 31 - __builtin_clz(~mask & 0xFF);
                    ub = (c - 7 + last) - r;
                    goto ub_done_s;
                }
            }
            if (c > limit) {
                __m256 v = _mm256_loadu_ps(&data[r*m + limit + 1]);
                __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_ps(cmp);
                if (mask != 0xFF) {
                    int last = 31 - __builtin_clz(~mask & 0xFF);
                    ub = (limit + 1 + last) - r;
                    goto ub_done_s;
                }
            }
        } else {
            for (Py_ssize_t c = m - 1; c > limit; c--) {
                if (data[r*m + c] != 0.0f) { ub = c - r; goto ub_done_s; }
            }
        }
        ub_done_s:
        if (r + ub + 1 > m) { break; }
    }

    *lower_band = lb;
    *upper_band = ub;
}

typedef void (*bandwidth_s_fn)(float* restrict, npy_intp, npy_intp, npy_intp*, npy_intp*);
static bandwidth_s_fn bandwidth_s_impl = 0;

__attribute__((constructor))
static void bandwidth_s_resolve(void) {
    bandwidth_s_impl = __builtin_cpu_supports("avx2") ? bandwidth_s_avx2
                                                      : bandwidth_s_scalar;
}
#endif // SCIPY_HAVE_AVX2_TARGET

static inline void
bandwidth_s(float* restrict data, npy_intp n, npy_intp m,
            npy_intp* lower_band, npy_intp* upper_band)
{
#ifdef SCIPY_HAVE_AVX2_TARGET
    bandwidth_s_impl(data, n, m, lower_band, upper_band);
#else
    bandwidth_s_scalar(data, n, m, lower_band, upper_band);
#endif
}


// bandwidth_d

static inline void
bandwidth_d_scalar(double* restrict data, npy_intp n, npy_intp m,
                   npy_intp* lower_band, npy_intp* upper_band)
{
    Py_ssize_t lb = 0, ub = 0;
    for (Py_ssize_t r = n-1; r > 0; r--) {
        Py_ssize_t limit = r - lb;
        if (limit > m) { limit = m; }
        for (Py_ssize_t c = 0; c < limit; c++) {
            if (data[r*m + c] != 0.0) { lb = r - c; break; }
        }
        if (r <= lb) { break; }
    }
    for (Py_ssize_t r = 0; r < n-1; r++) {
        for (Py_ssize_t c = m-1; c > r + ub; c--) {
            if (data[r*m + c] != 0.0) { ub = c - r; break; }
        }
        if (r + ub + 1 > m) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}

#ifdef SCIPY_HAVE_AVX2_TARGET
__attribute__((target("avx2")))
static inline void
bandwidth_d_avx2(double* restrict data, npy_intp n, npy_intp m,
                 npy_intp* lower_band, npy_intp* upper_band)
{
    const __m256d zero = _mm256_setzero_pd();
    Py_ssize_t lb = 0, ub = 0;

    for (Py_ssize_t r = n-1; r > 0; r--) {
        Py_ssize_t limit = r - lb;
        if (limit > m) { limit = m; }
        if (limit >= 4) {
            Py_ssize_t c = 0;
            for (; c + 3 < limit; c += 4) {
                __m256d v = _mm256_loadu_pd(&data[r*m + c]);
                __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_pd(cmp);
                if (mask != 0xF) {
                    lb = r - (c + __builtin_ctz(~mask));
                    goto lb_done_d;
                }
            }
            if (c < limit) {
                __m256d v = _mm256_loadu_pd(&data[r*m + limit - 4]);
                __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_pd(cmp);
                if (mask != 0xF) {
                    lb = r - (limit - 4 + __builtin_ctz(~mask));
                    goto lb_done_d;
                }
            }
        } else {
            for (Py_ssize_t c = 0; c < limit; c++) {
                if (data[r*m + c] != 0.0) { lb = r - c; goto lb_done_d; }
            }
        }
        lb_done_d:
        if (r <= lb) { break; }
    }

    for (Py_ssize_t r = 0; r < n-1; r++) {
        Py_ssize_t limit = r + ub;
        Py_ssize_t span = m - 1 - limit;
        if (span >= 4) {
            Py_ssize_t c = m - 1;
            for (; c - 3 > limit; c -= 4) {
                __m256d v = _mm256_loadu_pd(&data[r*m + c - 3]);
                __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_pd(cmp);
                if (mask != 0xF) {
                    int last = 31 - __builtin_clz(~mask & 0xF);
                    ub = (c - 3 + last) - r;
                    goto ub_done_d;
                }
            }
            if (c > limit) {
                __m256d v = _mm256_loadu_pd(&data[r*m + limit + 1]);
                __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_pd(cmp);
                if (mask != 0xF) {
                    int last = 31 - __builtin_clz(~mask & 0xF);
                    ub = (limit + 1 + last) - r;
                    goto ub_done_d;
                }
            }
        } else {
            for (Py_ssize_t c = m - 1; c > limit; c--) {
                if (data[r*m + c] != 0.0) { ub = c - r; goto ub_done_d; }
            }
        }
        ub_done_d:
        if (r + ub + 1 > m) { break; }
    }

    *lower_band = lb;
    *upper_band = ub;
}

typedef void (*bandwidth_d_fn)(double* restrict, npy_intp, npy_intp, npy_intp*, npy_intp*);
static bandwidth_d_fn bandwidth_d_impl = 0;

__attribute__((constructor))
static void bandwidth_d_resolve(void) {
    bandwidth_d_impl = __builtin_cpu_supports("avx2") ? bandwidth_d_avx2
                                                      : bandwidth_d_scalar;
}
#endif // SCIPY_HAVE_AVX2_TARGET

static inline void
bandwidth_d(double* restrict data, npy_intp n, npy_intp m,
            npy_intp* lower_band, npy_intp* upper_band)
{
#ifdef SCIPY_HAVE_AVX2_TARGET
    bandwidth_d_impl(data, n, m, lower_band, upper_band);
#else
    bandwidth_d_scalar(data, n, m, lower_band, upper_band);
#endif
}


// bandwidth_c
// Complex64: each element is 2 floats. Cast to float*, index as (r*m+c)*2.
// AVX2 loads 8 floats = 4 complex values. Mask bit-pairs: [re,im] per element.

static inline void
bandwidth_c_scalar(SCIPY_C* restrict data, npy_intp n, npy_intp m,
                   npy_intp* lower_band, npy_intp* upper_band)
{
    float* fdata = (float*)data;
    Py_ssize_t lb = 0, ub = 0;
    for (Py_ssize_t r = n-1; r > 0; r--) {
        Py_ssize_t limit = r - lb;
        if (limit > m) { limit = m; }
        for (Py_ssize_t c = 0; c < limit; c++) {
            Py_ssize_t idx = (r*m + c) * 2;
            if (fdata[idx] != 0.0f || fdata[idx + 1] != 0.0f) { lb = r - c; break; }
        }
        if (r <= lb) { break; }
    }
    for (Py_ssize_t r = 0; r < n-1; r++) {
        for (Py_ssize_t c = m-1; c > r + ub; c--) {
            Py_ssize_t idx = (r*m + c) * 2;
            if (fdata[idx] != 0.0f || fdata[idx + 1] != 0.0f) { ub = c - r; break; }
        }
        if (r + ub + 1 > m) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}

#ifdef SCIPY_HAVE_AVX2_TARGET
__attribute__((target("avx2")))
static inline void
bandwidth_c_avx2(SCIPY_C* restrict data, npy_intp n, npy_intp m,
                 npy_intp* lower_band, npy_intp* upper_band)
{
    float* fdata = (float*)data;
    const __m256 zero = _mm256_setzero_ps();
    Py_ssize_t lb = 0, ub = 0;

    // 8 floats = 4 complex64 per tray
    for (Py_ssize_t r = n-1; r > 0; r--) {
        Py_ssize_t limit = r - lb;
        if (limit > m) { limit = m; }
        if (limit >= 4) {
            Py_ssize_t c = 0;
            for (; c + 3 < limit; c += 4) {
                __m256 v = _mm256_loadu_ps(&fdata[(r*m + c) * 2]);
                __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_ps(cmp);
                if (mask != 0xFF) {
                    // 4 complex elements: bits [0,1]=c, [2,3]=c+1, [4,5]=c+2, [6,7]=c+3
                    if ((mask & 0x03) != 0x03) { lb = r - c;       goto lb_done_c; }
                    if ((mask & 0x0C) != 0x0C) { lb = r - (c + 1); goto lb_done_c; }
                    if ((mask & 0x30) != 0x30) { lb = r - (c + 2); goto lb_done_c; }
                    lb = r - (c + 3); goto lb_done_c;
                }
            }
            if (c < limit) {
                __m256 v = _mm256_loadu_ps(&fdata[(r*m + limit - 4) * 2]);
                __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_ps(cmp);
                if (mask != 0xFF) {
                    Py_ssize_t base = limit - 4;
                    if ((mask & 0x03) != 0x03) { lb = r - base;       goto lb_done_c; }
                    if ((mask & 0x0C) != 0x0C) { lb = r - (base + 1); goto lb_done_c; }
                    if ((mask & 0x30) != 0x30) { lb = r - (base + 2); goto lb_done_c; }
                    lb = r - (base + 3); goto lb_done_c;
                }
            }
        } else {
            for (Py_ssize_t c = 0; c < limit; c++) {
                Py_ssize_t idx = (r*m + c) * 2;
                if (fdata[idx] != 0.0f || fdata[idx + 1] != 0.0f) {
                    lb = r - c; goto lb_done_c;
                }
            }
        }
        lb_done_c:
        if (r <= lb) { break; }
    }

    for (Py_ssize_t r = 0; r < n-1; r++) {
        Py_ssize_t limit = r + ub;
        Py_ssize_t span = m - 1 - limit;
        if (span >= 4) {
            Py_ssize_t c = m - 1;
            for (; c - 3 > limit; c -= 4) {
                __m256 v = _mm256_loadu_ps(&fdata[(r*m + c - 3) * 2]);
                __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_ps(cmp);
                if (mask != 0xFF) {
                    // rightmost nonzero: check from high bits down
                    if ((mask & 0xC0) != 0xC0) { ub = c - r;       goto ub_done_c; }
                    if ((mask & 0x30) != 0x30) { ub = (c - 1) - r; goto ub_done_c; }
                    if ((mask & 0x0C) != 0x0C) { ub = (c - 2) - r; goto ub_done_c; }
                    ub = (c - 3) - r; goto ub_done_c;
                }
            }
            if (c > limit) {
                __m256 v = _mm256_loadu_ps(&fdata[(r*m + limit + 1) * 2]);
                __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_ps(cmp);
                if (mask != 0xFF) {
                    Py_ssize_t base = limit + 1;
                    if ((mask & 0xC0) != 0xC0) { ub = (base + 3) - r; goto ub_done_c; }
                    if ((mask & 0x30) != 0x30) { ub = (base + 2) - r; goto ub_done_c; }
                    if ((mask & 0x0C) != 0x0C) { ub = (base + 1) - r; goto ub_done_c; }
                    ub = base - r; goto ub_done_c;
                }
            }
        } else {
            for (Py_ssize_t c = m - 1; c > limit; c--) {
                Py_ssize_t idx = (r*m + c) * 2;
                if (fdata[idx] != 0.0f || fdata[idx + 1] != 0.0f) {
                    ub = c - r; goto ub_done_c;
                }
            }
        }
        ub_done_c:
        if (r + ub + 1 > m) { break; }
    }

    *lower_band = lb;
    *upper_band = ub;
}

typedef void (*bandwidth_c_fn)(SCIPY_C* restrict, npy_intp, npy_intp, npy_intp*, npy_intp*);
static bandwidth_c_fn bandwidth_c_impl = 0;

__attribute__((constructor))
static void bandwidth_c_resolve(void) {
    bandwidth_c_impl = __builtin_cpu_supports("avx2") ? bandwidth_c_avx2
                                                      : bandwidth_c_scalar;
}
#endif // SCIPY_HAVE_AVX2_TARGET

static inline void
bandwidth_c(SCIPY_C* restrict data, npy_intp n, npy_intp m,
            npy_intp* lower_band, npy_intp* upper_band)
{
#ifdef SCIPY_HAVE_AVX2_TARGET
    bandwidth_c_impl(data, n, m, lower_band, upper_band);
#else
    bandwidth_c_scalar(data, n, m, lower_band, upper_band);
#endif
}


// bandwidth_z
// Complex128: each element is 2 doubles. Cast to double*, index as (r*m+c)*2.
// AVX2 loads 4 doubles = 2 complex values. Mask bit-pairs: [re,im] per element.

static inline void
bandwidth_z_scalar(SCIPY_Z* restrict data, npy_intp n, npy_intp m,
                   npy_intp* lower_band, npy_intp* upper_band)
{
    double* ddata = (double*)data;
    Py_ssize_t lb = 0, ub = 0;
    for (Py_ssize_t r = n-1; r > 0; r--) {
        Py_ssize_t limit = r - lb;
        if (limit > m) { limit = m; }
        for (Py_ssize_t c = 0; c < limit; c++) {
            Py_ssize_t idx = (r*m + c) * 2;
            if (ddata[idx] != 0.0 || ddata[idx + 1] != 0.0) { lb = r - c; break; }
        }
        if (r <= lb) { break; }
    }
    for (Py_ssize_t r = 0; r < n-1; r++) {
        for (Py_ssize_t c = m-1; c > r + ub; c--) {
            Py_ssize_t idx = (r*m + c) * 2;
            if (ddata[idx] != 0.0 || ddata[idx + 1] != 0.0) { ub = c - r; break; }
        }
        if (r + ub + 1 > m) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}

#ifdef SCIPY_HAVE_AVX2_TARGET
__attribute__((target("avx2")))
static inline void
bandwidth_z_avx2(SCIPY_Z* restrict data, npy_intp n, npy_intp m,
                 npy_intp* lower_band, npy_intp* upper_band)
{
    double* ddata = (double*)data;
    const __m256d zero = _mm256_setzero_pd();
    Py_ssize_t lb = 0, ub = 0;

    // 4 doubles = 2 complex128 per tray
    for (Py_ssize_t r = n-1; r > 0; r--) {
        Py_ssize_t limit = r - lb;
        if (limit > m) { limit = m; }
        if (limit >= 2) {
            Py_ssize_t c = 0;
            for (; c + 1 < limit; c += 2) {
                __m256d v = _mm256_loadu_pd(&ddata[(r*m + c) * 2]);
                __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_pd(cmp);
                if (mask != 0xF) {
                    if ((mask & 0x3) != 0x3) { lb = r - c;       goto lb_done_z; }
                    else                      { lb = r - (c + 1); goto lb_done_z; }
                }
            }
            if (c < limit) {
                __m256d v = _mm256_loadu_pd(&ddata[(r*m + limit - 2) * 2]);
                __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_pd(cmp);
                if (mask != 0xF) {
                    if ((mask & 0x3) != 0x3) { lb = r - (limit - 2); goto lb_done_z; }
                    else                      { lb = r - (limit - 1); goto lb_done_z; }
                }
            }
        } else {
            for (Py_ssize_t c = 0; c < limit; c++) {
                Py_ssize_t idx = (r*m + c) * 2;
                if (ddata[idx] != 0.0 || ddata[idx + 1] != 0.0) {
                    lb = r - c; goto lb_done_z;
                }
            }
        }
        lb_done_z:
        if (r <= lb) { break; }
    }

    for (Py_ssize_t r = 0; r < n-1; r++) {
        Py_ssize_t limit = r + ub;
        Py_ssize_t span = m - 1 - limit;
        if (span >= 2) {
            Py_ssize_t c = m - 1;
            for (; c - 1 > limit; c -= 2) {
                __m256d v = _mm256_loadu_pd(&ddata[(r*m + c - 1) * 2]);
                __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_pd(cmp);
                if (mask != 0xF) {
                    if ((mask & 0xC) != 0xC) { ub = c - r;       goto ub_done_z; }
                    else                      { ub = (c - 1) - r; goto ub_done_z; }
                }
            }
            if (c > limit) {
                __m256d v = _mm256_loadu_pd(&ddata[(r*m + limit + 1) * 2]);
                __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                int mask = _mm256_movemask_pd(cmp);
                if (mask != 0xF) {
                    if ((mask & 0xC) != 0xC) { ub = (limit + 2) - r; goto ub_done_z; }
                    else                      { ub = (limit + 1) - r; goto ub_done_z; }
                }
            }
        } else {
            for (Py_ssize_t c = m - 1; c > limit; c--) {
                Py_ssize_t idx = (r*m + c) * 2;
                if (ddata[idx] != 0.0 || ddata[idx + 1] != 0.0) {
                    ub = c - r; goto ub_done_z;
                }
            }
        }
        ub_done_z:
        if (r + ub + 1 > m) { break; }
    }

    *lower_band = lb;
    *upper_band = ub;
}

typedef void (*bandwidth_z_fn)(SCIPY_Z* restrict, npy_intp, npy_intp, npy_intp*, npy_intp*);
static bandwidth_z_fn bandwidth_z_impl = 0;

__attribute__((constructor))
static void bandwidth_z_resolve(void) {
    bandwidth_z_impl = __builtin_cpu_supports("avx2") ? bandwidth_z_avx2
                                                      : bandwidth_z_scalar;
}
#endif // SCIPY_HAVE_AVX2_TARGET

static inline void
bandwidth_z(SCIPY_Z* restrict data, npy_intp n, npy_intp m,
            npy_intp* lower_band, npy_intp* upper_band)
{
#ifdef SCIPY_HAVE_AVX2_TARGET
    bandwidth_z_impl(data, n, m, lower_band, upper_band);
#else
    bandwidth_z_scalar(data, n, m, lower_band, upper_band);
#endif
}


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
