#ifndef SCIPY_LINALG_BANDWIDTH_H
#define SCIPY_LINALG_BANDWIDTH_H

#include <stdint.h>

// Complex definitions depending on whether the header is used with C++ or C and/or MSVC.
#ifndef SCIPY_Z
    #ifdef __cplusplus
        #include <complex>
        #define SCIPY_C std::complex<float>
        #define SCIPY_Z std::complex<double>
    #else
        #include <complex.h>
        #if defined(_MSC_VER)
            #define SCIPY_Z _Dcomplex
            #define SCIPY_C _Fcomplex
        #else
            #define SCIPY_Z double complex
            #define SCIPY_C float complex
        #endif
    #endif
#endif

#ifdef __cplusplus
    #define SCIPY_RESTRICT
#else
    #define SCIPY_RESTRICT restrict
#endif

// AVX2 dispatch requires __builtin_cpu_supports and __attribute__((target)),
// available on GCC and Clang on non-Windows x86-64. Clang on Windows uses
// lld-link which cannot resolve __cpu_model from compiler-rt.

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && !defined(_WIN32)
#include <immintrin.h>
#define SCIPY_HAVE_AVX2_TARGET 1
#endif


/**
 * Bandwidth detection — scalar contiguous loop implementations
 */

static inline void
bandwidth_s_scalar(const float* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
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


static inline void
bandwidth_d_scalar(const double* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
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


static inline void
bandwidth_c_scalar(const SCIPY_C* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
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


static inline void
bandwidth_z_scalar(const SCIPY_Z* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
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
 * Strided scalar bandwidth — C++ only, handles arbitrary strides and all
 * numeric types via templates. Used for non-contiguous arrays and integer types,
 * called from the batched linalg module _bandwidth function.
 */

#ifdef __cplusplus

    template<typename T>
    inline void
    bandwidth_strided_scalar(const void *data, int64_t offset,
                            int64_t n, int64_t m, int64_t s0, int64_t s1,
                            int64_t *lower_band, int64_t *upper_band)
    {
        const T *slice = (const T *)data + offset;
        const T zero = T(0);
        int64_t lb = 0, ub = 0;

        for (int64_t r = n - 1; r > 0; r--) {
            int64_t ncols = r - lb;
            if (ncols > m) { ncols = m; }
            for (int64_t c = 0; c < ncols; c++) {
                if (slice[r * s0 + c * s1] != zero) { lb = r - c; break; }
            }
            if (r <= lb) { break; }
        }

        for (int64_t r = 0; r < n - 1; r++) {
            for (int64_t c = m - 1; c > r + ub; c--) {
                if (slice[r * s0 + c * s1] != zero) { ub = c - r; break; }
            }
            if (r + ub + 1 > m) { break; }
        }

        *lower_band = lb;
        *upper_band = ub;
    }

#endif /* __cplusplus */


/**
 * Bandwidth detection — AVX2 implementations and dynamic dispatch
 */

#ifdef SCIPY_HAVE_AVX2_TARGET

    __attribute__((target("avx2")))
    static void
    bandwidth_s_avx2(const float* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
    {
        int64_t lb = 0, ub = 0;

        // Create a 256-bit AVX vector of zeros
        const __m256 zero = _mm256_setzero_ps();

        // Lower bandwidth:
        // We are scanning the array in C-major order. For the lower bandwidth, we scan
        // from lower left corner going right for each row.
        for (int64_t r = n-1; r > 0; r--) {

            // How many columns in this row to check? If in the previous rows we found
            // a nonzero element we don't need to check beyond past the current detected lb,
            // so we only need to check up to column r-lb. Obviously, we cannot exceed the
            // total number of columns m.
            int64_t ncols = (r - lb > m) ? m : (r - lb);

            // Do we have enough columns to fill a 256-bit vector (8 floats)?
            if (ncols >= 8) {
                int64_t c = 0;

                for (; c + 7 < ncols; c += 8) {

                    // Load 8 floats (8*32bit) from consecutive rows to a vector v
                    __m256 v = _mm256_loadu_ps(&data[r*m + c]);

                    // Compare v to zero, that returns another mask vector cmp,
                    // equality marks 0xFFFFFFFF for true, and 0x00000000 for false.
                    __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);

                    // Extract one bit (sign bit) per float to create an 8-bit mask.
                    // If all 8 floats are zero, the mask will be 0b11111111 = 0xFF.
                    int mask = _mm256_movemask_ps(cmp);

                    // All zero?
                    if (mask != 0xFF) {
                        // No, count the trailing zeros in the inverted mask to find the first nonzero element
                        lb = r - (c + __builtin_ctz(~mask));

                        // We found nonzero element, this row is done.
                        goto row_done_lb_s;
                    }
                }

                // A small overlay trick to avoid falling back to scalar loops.
                // If a tail is left, then reload previous elements to complete the vector:
                // Say, ncols = 10, and we checked elements [0..7] in the loop,
                // then take [2..9] as another 8-vector to check the last 2 elements.
                if (c < ncols) {
                    __m256 v = _mm256_loadu_ps(&data[r*m + ncols - 8]);
                    __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_ps(cmp);
                    if (mask != 0xFF) {
                        lb = r - (ncols - 8 + __builtin_ctz(~mask));
                        goto row_done_lb_s;
                    }
                }
            } else {
                // Not enough elements to fill a vector, hence scalar loop.
                for (int64_t c = 0; c < ncols; c++) {
                    if (data[r*m + c] != 0.0f) { lb = r - c; goto row_done_lb_s; }
                }
            }
    row_done_lb_s:

            // If lb is already, say, 5 then we don't need to check rows above 5th row.
            if (r <= lb) { break; }
        }


        // Upper bandwidth:
        // Scanning from the upper right corner going left for each row.
        for (int64_t r = 0; r < n-1; r++) {

            // We scan right-to-left starting from column m-1. Similar story
            // as lower bandwidth. So we stop the row scan at the current ub
            // which is denoted as stop_col. ncols denotes how many columns until
            // that value from the right edge.
            int64_t stop_col = r + ub;
            int64_t ncols = m - 1 - stop_col;

            if (ncols >= 8) {
                int64_t c = m - 1;
                for (; c - 7 > stop_col; c -= 8) {

                    // Load 8 floats ending at column c (columns [c-7 .. c])
                    __m256 v = _mm256_loadu_ps(&data[r*m + c - 7]);
                    __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_ps(cmp);
                    if (mask != 0xFF) {
                        // First nonzero from the right with (c)ount (l)eading (z)eros in the inverted mask
                        int nz = 31 - __builtin_clz(~mask & 0xFF);
                        ub = (c - 7 + nz) - r;
                        goto row_done_ub_s;
                    }
                }
                if (c > stop_col) {
                    __m256 v = _mm256_loadu_ps(&data[r*m + stop_col + 1]);
                    __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_ps(cmp);
                    if (mask != 0xFF) {
                        int nz = 31 - __builtin_clz(~mask & 0xFF);
                        ub = (stop_col + 1 + nz) - r;
                        goto row_done_ub_s;
                    }
                }
            } else {
                for (int64_t c = m - 1; c > stop_col; c--) {
                    if (data[r*m + c] != 0.0f) { ub = c - r; goto row_done_ub_s; }
                }
            }
    row_done_ub_s:

            // If the known band already reaches the last column, no more rows to check.
            if (r + ub + 1 > m) { break; }
        }

        *lower_band = lb;
        *upper_band = ub;
    }


    __attribute__((target("avx2")))
    static void
    bandwidth_d_avx2(const double* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
    {
        const __m256d zero = _mm256_setzero_pd();
        int64_t lb = 0, ub = 0;

        for (int64_t r = n-1; r > 0; r--) {
            int64_t ncols = (r - lb > m) ? m : (r - lb);
            if (ncols >= 4) {
                int64_t c = 0;
                for (; c + 3 < ncols; c += 4) {
                    __m256d v = _mm256_loadu_pd(&data[r*m + c]);
                    __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_pd(cmp);
                    if (mask != 0xF) {
                        lb = r - (c + __builtin_ctz(~mask));
                        goto row_done_lb_d;
                    }
                }
                if (c < ncols) {
                    __m256d v = _mm256_loadu_pd(&data[r*m + ncols - 4]);
                    __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_pd(cmp);
                    if (mask != 0xF) {
                        lb = r - (ncols - 4 + __builtin_ctz(~mask));
                        goto row_done_lb_d;
                    }
                }
            } else {
                for (int64_t c = 0; c < ncols; c++) {
                    if (data[r*m + c] != 0.0) { lb = r - c; goto row_done_lb_d; }
                }
            }
    row_done_lb_d:
            if (r <= lb) { break; }
        }

        for (int64_t r = 0; r < n-1; r++) {
            int64_t stop_col = r + ub;
            int64_t ncols = m - 1 - stop_col;
            if (ncols >= 4) {
                int64_t c = m - 1;
                for (; c - 3 > stop_col; c -= 4) {
                    __m256d v = _mm256_loadu_pd(&data[r*m + c - 3]);
                    __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_pd(cmp);
                    if (mask != 0xF) {
                        int nz = 31 - __builtin_clz(~mask & 0xF);
                        ub = (c - 3 + nz) - r;
                        goto row_done_ub_d;
                    }
                }
                if (c > stop_col) {
                    __m256d v = _mm256_loadu_pd(&data[r*m + stop_col + 1]);
                    __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_pd(cmp);
                    if (mask != 0xF) {
                        int nz = 31 - __builtin_clz(~mask & 0xF);
                        ub = (stop_col + 1 + nz) - r;
                        goto row_done_ub_d;
                    }
                }
            } else {
                for (int64_t c = m - 1; c > stop_col; c--) {
                    if (data[r*m + c] != 0.0) { ub = c - r; goto row_done_ub_d; }
                }
            }
    row_done_ub_d:
            if (r + ub + 1 > m) { break; }
        }

        *lower_band = lb;
        *upper_band = ub;
    }


    __attribute__((target("avx2")))
    static void
    bandwidth_c_avx2(const SCIPY_C* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
    {
        float* fdata = (float*)data;
        const __m256 zero = _mm256_setzero_ps();
        int64_t lb = 0, ub = 0;

        // 8 floats = 4 complex64 per vector; each complex element is 2 floats
        // so bit positions in the mask map to elements via / 2.
        for (int64_t r = n-1; r > 0; r--) {
            int64_t ncols = (r - lb > m) ? m : (r - lb);
            if (ncols >= 4) {
                int64_t c = 0;
                for (; c + 3 < ncols; c += 4) {
                    __m256 v = _mm256_loadu_ps(&fdata[(r*m + c) * 2]);
                    __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_ps(cmp);
                    if (mask != 0xFF) {
                        lb = r - (c + __builtin_ctz(~mask) / 2);
                        goto row_done_lb_c;
                    }
                }
                if (c < ncols) {
                    __m256 v = _mm256_loadu_ps(&fdata[(r*m + ncols - 4) * 2]);
                    __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_ps(cmp);
                    if (mask != 0xFF) {
                        lb = r - (ncols - 4 + __builtin_ctz(~mask) / 2);
                        goto row_done_lb_c;
                    }
                }
            } else {
                for (int64_t c = 0; c < ncols; c++) {
                    int64_t idx = (r*m + c) * 2;
                    if (fdata[idx] != 0.0f || fdata[idx + 1] != 0.0f) {
                        lb = r - c; goto row_done_lb_c;
                    }
                }
            }
    row_done_lb_c:
            if (r <= lb) { break; }
        }

        for (int64_t r = 0; r < n-1; r++) {
            int64_t stop_col = r + ub;
            int64_t ncols = m - 1 - stop_col;
            if (ncols >= 4) {
                int64_t c = m - 1;
                for (; c - 3 > stop_col; c -= 4) {
                    __m256 v = _mm256_loadu_ps(&fdata[(r*m + c - 3) * 2]);
                    __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_ps(cmp);
                    if (mask != 0xFF) {
                        int nz = (31 - __builtin_clz(~mask & 0xFF)) / 2;
                        ub = (c - 3 + nz) - r;
                        goto row_done_ub_c;
                    }
                }
                if (c > stop_col) {
                    __m256 v = _mm256_loadu_ps(&fdata[(r*m + stop_col + 1) * 2]);
                    __m256 cmp = _mm256_cmp_ps(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_ps(cmp);
                    if (mask != 0xFF) {
                        int nz = (31 - __builtin_clz(~mask & 0xFF)) / 2;
                        ub = (stop_col + 1 + nz) - r;
                        goto row_done_ub_c;
                    }
                }
            } else {
                for (int64_t c = m - 1; c > stop_col; c--) {
                    int64_t idx = (r*m + c) * 2;
                    if (fdata[idx] != 0.0f || fdata[idx + 1] != 0.0f) {
                        ub = c - r; goto row_done_ub_c;
                    }
                }
            }
    row_done_ub_c:
            if (r + ub + 1 > m) { break; }
        }

        *lower_band = lb;
        *upper_band = ub;
    }


    __attribute__((target("avx2")))
    static void
    bandwidth_z_avx2(const SCIPY_Z* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
    {
        double* ddata = (double*)data;
        const __m256d zero = _mm256_setzero_pd();
        int64_t lb = 0, ub = 0;

        // 4 doubles = 2 complex128 per vector
        for (int64_t r = n-1; r > 0; r--) {
            int64_t ncols = (r - lb > m) ? m : (r - lb);
            if (ncols >= 2) {
                int64_t c = 0;
                for (; c + 1 < ncols; c += 2) {
                    __m256d v = _mm256_loadu_pd(&ddata[(r*m + c) * 2]);
                    __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_pd(cmp);
                    if (mask != 0xF) {
                        if ((mask & 0x3) != 0x3) { lb = r - c; }
                        else                     { lb = r - (c + 1); }
                        goto row_done_lb_z;
                    }
                }
                if (c < ncols) {
                    __m256d v = _mm256_loadu_pd(&ddata[(r*m + ncols - 2) * 2]);
                    __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_pd(cmp);
                    if (mask != 0xF) {
                        if ((mask & 0x3) != 0x3) { lb = r - (ncols - 2); }
                        else                     { lb = r - (ncols - 1); }
                        goto row_done_lb_z;
                    }
                }
            } else {
                for (int64_t c = 0; c < ncols; c++) {
                    int64_t idx = (r*m + c) * 2;
                    if (ddata[idx] != 0.0 || ddata[idx + 1] != 0.0) {
                        lb = r - c; goto row_done_lb_z;
                    }
                }
            }
    row_done_lb_z:
            if (r <= lb) { break; }
        }

        for (int64_t r = 0; r < n-1; r++) {
            int64_t stop_col = r + ub;
            int64_t ncols = m - 1 - stop_col;
            if (ncols >= 2) {
                int64_t c = m - 1;
                for (; c - 1 > stop_col; c -= 2) {
                    __m256d v = _mm256_loadu_pd(&ddata[(r*m + c - 1) * 2]);
                    __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_pd(cmp);
                    if (mask != 0xF) {
                        if ((mask & 0xC) != 0xC) { ub = c - r; }
                        else                     { ub = (c - 1) - r; }
                        goto row_done_ub_z;
                    }
                }
                if (c > stop_col) {
                    __m256d v = _mm256_loadu_pd(&ddata[(r*m + stop_col + 1) * 2]);
                    __m256d cmp = _mm256_cmp_pd(v, zero, _CMP_EQ_OQ);
                    int mask = _mm256_movemask_pd(cmp);
                    if (mask != 0xF) {
                        if ((mask & 0xC) != 0xC) { ub = (stop_col + 2) - r; }
                        else                     { ub = (stop_col + 1) - r; }
                        goto row_done_ub_z;
                    }
                }
            } else {
                for (int64_t c = m - 1; c > stop_col; c--) {
                    int64_t idx = (r*m + c) * 2;
                    if (ddata[idx] != 0.0 || ddata[idx + 1] != 0.0) {
                        ub = c - r; goto row_done_ub_z;
                    }
                }
            }
    row_done_ub_z:
            if (r + ub + 1 > m) { break; }
        }

        *lower_band = lb;
        *upper_band = ub;
    }


    /**
    * Dynamic dispatch — resolved once at load time via constructor
    */

    typedef void (*bandwidth_s_fn)(const float* SCIPY_RESTRICT, int64_t, int64_t, int64_t*, int64_t*);
    static bandwidth_s_fn bandwidth_s_impl = 0;

    typedef void (*bandwidth_d_fn)(const double* SCIPY_RESTRICT, int64_t, int64_t, int64_t*, int64_t*);
    static bandwidth_d_fn bandwidth_d_impl = 0;

    typedef void (*bandwidth_c_fn)(const SCIPY_C* SCIPY_RESTRICT, int64_t, int64_t, int64_t*, int64_t*);
    static bandwidth_c_fn bandwidth_c_impl = 0;

    typedef void (*bandwidth_z_fn)(const SCIPY_Z* SCIPY_RESTRICT, int64_t, int64_t, int64_t*, int64_t*);
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
 * Public entry points — use AVX2 when available, scalar otherwise.
 */

static inline void
bandwidth_s(const float* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
{
#ifdef SCIPY_HAVE_AVX2_TARGET
    bandwidth_s_impl(data, n, m, lower_band, upper_band);
#else
    bandwidth_s_scalar(data, n, m, lower_band, upper_band);
#endif
}

static inline void
bandwidth_d(const double* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
{
#ifdef SCIPY_HAVE_AVX2_TARGET
    bandwidth_d_impl(data, n, m, lower_band, upper_band);
#else
    bandwidth_d_scalar(data, n, m, lower_band, upper_band);
#endif
}

static inline void
bandwidth_c(const SCIPY_C* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
{
#ifdef SCIPY_HAVE_AVX2_TARGET
    bandwidth_c_impl(data, n, m, lower_band, upper_band);
#else
    bandwidth_c_scalar(data, n, m, lower_band, upper_band);
#endif
}

static inline void
bandwidth_z(const SCIPY_Z* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band)
{
#ifdef SCIPY_HAVE_AVX2_TARGET
    bandwidth_z_impl(data, n, m, lower_band, upper_band);
#else
    bandwidth_z_scalar(data, n, m, lower_band, upper_band);
#endif
}

#endif /* SCIPY_LINALG_BANDWIDTH_H */
