#include "_linalg_simd_kernels.h"

// AVX2 dispatch requires __builtin_cpu_supports and __attribute__((target)),
// available on GCC and Clang on non-Windows x86-64. Clang on Windows uses
// lld-link which cannot resolve __cpu_model from compiler-rt.

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(_WIN32)
#include <immintrin.h>
#define SCIPY_HAVE_AVX2_TARGET 1
#endif


/**
 * Bandwidth detection — scalar implementations (always available)
 */

// bandwidth_s_scalar

static void
bandwidth_s_scalar(const float* restrict data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
{
    int64_t lb = 0, ub = 0;
    for (int64_t r = n-1; r > 0; r--) {
        int64_t limit = r - lb;
        if (limit > m) { limit = m; }
        for (int64_t c = 0; c < limit; c++) {
            if (data[r*m + c] != 0.0f) { lb = r - c; break; }
        }
        if (r <= lb) { break; }
    }
    for (int64_t r = 0; r < n-1; r++) {
        for (int64_t c = m-1; c > r + ub; c--) {
            if (data[r*m + c] != 0.0f) { ub = c - r; break; }
        }
        if (r + ub + 1 > m) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}


// bandwidth_d_scalar

static void
bandwidth_d_scalar(const double* restrict data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
{
    int64_t lb = 0, ub = 0;
    for (int64_t r = n-1; r > 0; r--) {
        int64_t limit = r - lb;
        if (limit > m) { limit = m; }
        for (int64_t c = 0; c < limit; c++) {
            if (data[r*m + c] != 0.0) { lb = r - c; break; }
        }
        if (r <= lb) { break; }
    }
    for (int64_t r = 0; r < n-1; r++) {
        for (int64_t c = m-1; c > r + ub; c--) {
            if (data[r*m + c] != 0.0) { ub = c - r; break; }
        }
        if (r + ub + 1 > m) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}


// bandwidth_c_scalar

static void
bandwidth_c_scalar(const SCIPY_C* restrict data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
{
    float* fdata = (float*)data;
    int64_t lb = 0, ub = 0;
    for (int64_t r = n-1; r > 0; r--) {
        int64_t limit = r - lb;
        if (limit > m) { limit = m; }
        for (int64_t c = 0; c < limit; c++) {
            int64_t idx = (r*m + c) * 2;
            if (fdata[idx] != 0.0f || fdata[idx + 1] != 0.0f) { lb = r - c; break; }
        }
        if (r <= lb) { break; }
    }
    for (int64_t r = 0; r < n-1; r++) {
        for (int64_t c = m-1; c > r + ub; c--) {
            int64_t idx = (r*m + c) * 2;
            if (fdata[idx] != 0.0f || fdata[idx + 1] != 0.0f) { ub = c - r; break; }
        }
        if (r + ub + 1 > m) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}


// bandwidth_z_scalar

static void
bandwidth_z_scalar(const SCIPY_Z* restrict data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
{
    double* ddata = (double*)data;
    int64_t lb = 0, ub = 0;
    for (int64_t r = n-1; r > 0; r--) {
        int64_t limit = r - lb;
        if (limit > m) { limit = m; }
        for (int64_t c = 0; c < limit; c++) {
            int64_t idx = (r*m + c) * 2;
            if (ddata[idx] != 0.0 || ddata[idx + 1] != 0.0) { lb = r - c; break; }
        }
        if (r <= lb) { break; }
    }
    for (int64_t r = 0; r < n-1; r++) {
        for (int64_t c = m-1; c > r + ub; c--) {
            int64_t idx = (r*m + c) * 2;
            if (ddata[idx] != 0.0 || ddata[idx + 1] != 0.0) { ub = c - r; break; }
        }
        if (r + ub + 1 > m) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}


/**
 * Bandwidth detection — AVX2 implementations and dynamic dispatch
 */

#ifdef SCIPY_HAVE_AVX2_TARGET

// bandwidth_s_avx2

__attribute__((target("avx2")))
static void
bandwidth_s_avx2(const float* restrict data, int64_t n, int64_t m,
                 int64_t* lower_band, int64_t* upper_band)
{
    const __m256 zero = _mm256_setzero_ps();
    int64_t lb = 0, ub = 0;

    for (int64_t r = n-1; r > 0; r--) {
        int64_t limit = r - lb;
        if (limit > m) { limit = m; }
        if (limit >= 8) {
            int64_t c = 0;
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
            for (int64_t c = 0; c < limit; c++) {
                if (data[r*m + c] != 0.0f) { lb = r - c; goto lb_done_s; }
            }
        }
lb_done_s:
        if (r <= lb) { break; }
    }

    for (int64_t r = 0; r < n-1; r++) {
        int64_t limit = r + ub;
        int64_t span = m - 1 - limit;
        if (span >= 8) {
            int64_t c = m - 1;
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
            for (int64_t c = m - 1; c > limit; c--) {
                if (data[r*m + c] != 0.0f) { ub = c - r; goto ub_done_s; }
            }
        }
ub_done_s:
        if (r + ub + 1 > m) { break; }
    }

    *lower_band = lb;
    *upper_band = ub;
}


// bandwidth_d_avx2

__attribute__((target("avx2")))
static void
bandwidth_d_avx2(const double* restrict data, int64_t n, int64_t m,
                 int64_t* lower_band, int64_t* upper_band)
{
    const __m256d zero = _mm256_setzero_pd();
    int64_t lb = 0, ub = 0;

    for (int64_t r = n-1; r > 0; r--) {
        int64_t limit = r - lb;
        if (limit > m) { limit = m; }
        if (limit >= 4) {
            int64_t c = 0;
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
            for (int64_t c = 0; c < limit; c++) {
                if (data[r*m + c] != 0.0) { lb = r - c; goto lb_done_d; }
            }
        }
lb_done_d:
        if (r <= lb) { break; }
    }

    for (int64_t r = 0; r < n-1; r++) {
        int64_t limit = r + ub;
        int64_t span = m - 1 - limit;
        if (span >= 4) {
            int64_t c = m - 1;
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
            for (int64_t c = m - 1; c > limit; c--) {
                if (data[r*m + c] != 0.0) { ub = c - r; goto ub_done_d; }
            }
        }
ub_done_d:
        if (r + ub + 1 > m) { break; }
    }

    *lower_band = lb;
    *upper_band = ub;
}


// bandwidth_c_avx2

__attribute__((target("avx2")))
static void
bandwidth_c_avx2(const SCIPY_C* restrict data, int64_t n, int64_t m,
                 int64_t* lower_band, int64_t* upper_band)
{
    float* fdata = (float*)data;
    const __m256 zero = _mm256_setzero_ps();
    int64_t lb = 0, ub = 0;

    // 8 floats = 4 complex64 per tray
    for (int64_t r = n-1; r > 0; r--) {
        int64_t limit = r - lb;
        if (limit > m) { limit = m; }
        if (limit >= 4) {
            int64_t c = 0;
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
                    int64_t base = limit - 4;
                    if ((mask & 0x03) != 0x03) { lb = r - base;       goto lb_done_c; }
                    if ((mask & 0x0C) != 0x0C) { lb = r - (base + 1); goto lb_done_c; }
                    if ((mask & 0x30) != 0x30) { lb = r - (base + 2); goto lb_done_c; }
                    lb = r - (base + 3); goto lb_done_c;
                }
            }
        } else {
            for (int64_t c = 0; c < limit; c++) {
                int64_t idx = (r*m + c) * 2;
                if (fdata[idx] != 0.0f || fdata[idx + 1] != 0.0f) {
                    lb = r - c; goto lb_done_c;
                }
            }
        }
lb_done_c:
        if (r <= lb) { break; }
    }

    for (int64_t r = 0; r < n-1; r++) {
        int64_t limit = r + ub;
        int64_t span = m - 1 - limit;
        if (span >= 4) {
            int64_t c = m - 1;
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
                    int64_t base = limit + 1;
                    if ((mask & 0xC0) != 0xC0) { ub = (base + 3) - r; goto ub_done_c; }
                    if ((mask & 0x30) != 0x30) { ub = (base + 2) - r; goto ub_done_c; }
                    if ((mask & 0x0C) != 0x0C) { ub = (base + 1) - r; goto ub_done_c; }
                    ub = base - r; goto ub_done_c;
                }
            }
        } else {
            for (int64_t c = m - 1; c > limit; c--) {
                int64_t idx = (r*m + c) * 2;
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


// bandwidth_z_avx2

__attribute__((target("avx2")))
static void
bandwidth_z_avx2(const SCIPY_Z* restrict data, int64_t n, int64_t m,
                 int64_t* lower_band, int64_t* upper_band)
{
    double* ddata = (double*)data;
    const __m256d zero = _mm256_setzero_pd();
    int64_t lb = 0, ub = 0;

    // 4 doubles = 2 complex128 per tray
    for (int64_t r = n-1; r > 0; r--) {
        int64_t limit = r - lb;
        if (limit > m) { limit = m; }
        if (limit >= 2) {
            int64_t c = 0;
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
            for (int64_t c = 0; c < limit; c++) {
                int64_t idx = (r*m + c) * 2;
                if (ddata[idx] != 0.0 || ddata[idx + 1] != 0.0) {
                    lb = r - c; goto lb_done_z;
                }
            }
        }
lb_done_z:
        if (r <= lb) { break; }
    }

    for (int64_t r = 0; r < n-1; r++) {
        int64_t limit = r + ub;
        int64_t span = m - 1 - limit;
        if (span >= 2) {
            int64_t c = m - 1;
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
            for (int64_t c = m - 1; c > limit; c--) {
                int64_t idx = (r*m + c) * 2;
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


/**
 * Dynamic dispatch — resolved once at load time via constructor
 */

typedef void (*bandwidth_s_fn)(const float* restrict, int64_t, int64_t, int64_t*, int64_t*);
static bandwidth_s_fn bandwidth_s_impl = 0;

typedef void (*bandwidth_d_fn)(const double* restrict, int64_t, int64_t, int64_t*, int64_t*);
static bandwidth_d_fn bandwidth_d_impl = 0;

typedef void (*bandwidth_c_fn)(const SCIPY_C* restrict, int64_t, int64_t, int64_t*, int64_t*);
static bandwidth_c_fn bandwidth_c_impl = 0;

typedef void (*bandwidth_z_fn)(const SCIPY_Z* restrict, int64_t, int64_t, int64_t*, int64_t*);
static bandwidth_z_fn bandwidth_z_impl = 0;

__attribute__((constructor))
static void bandwidth_resolve(void) {
    if (__builtin_cpu_supports("avx2")) {
        bandwidth_s_impl = bandwidth_s_avx2;
        bandwidth_d_impl = bandwidth_d_avx2;
        bandwidth_c_impl = bandwidth_c_avx2;
        bandwidth_z_impl = bandwidth_z_avx2;
    } else {
        bandwidth_s_impl = bandwidth_s_scalar;
        bandwidth_d_impl = bandwidth_d_scalar;
        bandwidth_c_impl = bandwidth_c_scalar;
        bandwidth_z_impl = bandwidth_z_scalar;
    }
}

#endif /* SCIPY_HAVE_AVX2_TARGET */


/**
 * Public entry points
 */

void
bandwidth_s(const float* restrict data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
{
#ifdef SCIPY_HAVE_AVX2_TARGET
    bandwidth_s_impl(data, n, m, lower_band, upper_band);
#else
    bandwidth_s_scalar(data, n, m, lower_band, upper_band);
#endif
}

void
bandwidth_d(const double* restrict data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
{
#ifdef SCIPY_HAVE_AVX2_TARGET
    bandwidth_d_impl(data, n, m, lower_band, upper_band);
#else
    bandwidth_d_scalar(data, n, m, lower_band, upper_band);
#endif
}

void
bandwidth_c(const SCIPY_C* restrict data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
{
#ifdef SCIPY_HAVE_AVX2_TARGET
    bandwidth_c_impl(data, n, m, lower_band, upper_band);
#else
    bandwidth_c_scalar(data, n, m, lower_band, upper_band);
#endif
}

void
bandwidth_z(const SCIPY_Z* restrict data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
{
#ifdef SCIPY_HAVE_AVX2_TARGET
    bandwidth_z_impl(data, n, m, lower_band, upper_band);
#else
    bandwidth_z_scalar(data, n, m, lower_band, upper_band);
#endif
}
