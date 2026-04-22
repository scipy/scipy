/*
 * Portable CMPLX / CMPLXF macros for constructing complex numbers without
 * the NaN-corrupting  real + imag * I  pattern.
 *
 * The C11 standard provides CMPLX/CMPLXF, but not all compilers expose
 * them (notably MSVC).  This header provides fallback definitions so
 * that every translation unit can simply write  CMPLX(re, im)  and get
 * correct behaviour with NaN / Inf arguments on every platform.
 *
 * This header includes <complex.h> itself.
 *
 * Priority order for the fallback:
 *   1. Compiler-provided CMPLX / CMPLXF  (C11, already defined)
 *   2. MSVC intrinsics  (_Cbuild / _FCbuild)
 *   3. __builtin_complex  (GCC ≥ 4.7, Clang)
 *   4. Union type-pun     (portable, avoids  real + imag*I)
 */

#ifndef SCIPY_COMPLEX_SUPPORT_H
#define SCIPY_COMPLEX_SUPPORT_H

#include <complex.h>

/* ---- CMPLX (double complex) ------------------------------------------- */
#ifndef CMPLX
    #if defined(_MSC_VER)
        #define CMPLX(x, y) _Cbuild(x, y)
    #elif defined(__has_builtin)
        #if __has_builtin(__builtin_complex)
            #define CMPLX(x, y) __builtin_complex((double)(x), (double)(y))
        #endif
    #endif
    #ifndef CMPLX
        #define CMPLX(x, y) \
            ((union { double a[2]; double complex z; }){{(x), (y)}}).z
    #endif
#endif

/* ---- CMPLXF (float complex) ------------------------------------------- */
#ifndef CMPLXF
    #if defined(_MSC_VER)
        #define CMPLXF(x, y) _FCbuild(x, y)
    #elif defined(__has_builtin)
        #if __has_builtin(__builtin_complex)
            #define CMPLXF(x, y) __builtin_complex((float)(x), (float)(y))
        #endif
    #endif
    #ifndef CMPLXF
        #define CMPLXF(x, y) \
            ((union { float a[2]; float complex z; }){{(x), (y)}}).z
    #endif
#endif

#endif /* SCIPY_COMPLEX_SUPPORT_H */
