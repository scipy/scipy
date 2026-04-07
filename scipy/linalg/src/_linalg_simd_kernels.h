#ifndef SCIPY_LINALG_SIMD_KERNELS_H
#define SCIPY_LINALG_SIMD_KERNELS_H

#include <stdint.h>

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

/*
 * SIMD-accelerated kernels for scipy.linalg.
 *
 * This header declares the public entry points. Implementations live in
 * _linalg_simd_kernels.c which is compiled as its own translation unit.
 * On GCC/Clang x86-64 (non-Windows), AVX2 paths are selected at load time
 * via __attribute__((constructor)); otherwise scalar fallbacks are used.
 */

#ifdef __cplusplus
  #define SCIPY_RESTRICT
extern "C" {
#else
  #define SCIPY_RESTRICT restrict
#endif

void bandwidth_s(const float*   SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band);
void bandwidth_d(const double*  SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band);
void bandwidth_c(const SCIPY_C* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band);
void bandwidth_z(const SCIPY_Z* SCIPY_RESTRICT data, int64_t n, int64_t m, int64_t* lower_band, int64_t* upper_band);

#ifdef __cplusplus
}
#endif

#endif /* SCIPY_LINALG_SIMD_KERNELS_H */
