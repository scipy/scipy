/*
 * Complex element types shared by C and C++ code, plus portable CMPLX / CMPLXF
 * macros for C.
 *
 * scipy_c64 and scipy_c128 are the complex64 and complex128 element types:
 *
 *   C++                      std::complex<float>   std::complex<double>
 *   C, MSVC ABI (incl.       _Fcomplex             _Dcomplex
 *   clang-cl)
 *   C, otherwise             float _Complex        double _Complex
 *
 * Each has the layout of float[2] / double[2] (real part first).
 *
 * CMPLX / CMPLXF (C only) construct complex numbers without the NaN-corrupting
 * real + imag * I pattern. The C11 standard provides them, but not all compilers
 * expose them (notably MSVC), so this header provides fallback definitions.
 * Priority order for the fallback:
 *   1. Compiler-provided CMPLX / CMPLXF  (C11, already defined)
 *   2. MSVC intrinsics  (_Cbuild / _FCbuild)
 *   3. __builtin_complex  (GCC >= 4.7, Clang)
 *   4. Union type-pun     (portable, avoids  real + imag*I)
 *
 * In C this header includes <complex.h>; in C++ it includes <complex>.
 */

#ifndef SCIPY_COMPLEX_SUPPORT_H
#define SCIPY_COMPLEX_SUPPORT_H

#ifdef __cplusplus

#include <complex>

typedef std::complex<float>  scipy_c64;
typedef std::complex<double> scipy_c128;

#else

#include <complex.h>

#if defined(_MSC_VER)
typedef _Fcomplex scipy_c64;
typedef _Dcomplex scipy_c128;
#else
typedef float _Complex  scipy_c64;
typedef double _Complex scipy_c128;
#endif

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

#endif /* __cplusplus */

#endif /* SCIPY_COMPLEX_SUPPORT_H */
