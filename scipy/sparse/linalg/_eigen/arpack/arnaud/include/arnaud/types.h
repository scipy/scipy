#ifndef ARNAUD_TYPES_H
#define ARNAUD_TYPES_H

#include <stdint.h>
#include <complex.h>

/*
 * Portable CMPLX / CMPLXF fallbacks.
 *
 * The C11 standard provides CMPLX/CMPLXF, but not all compilers expose
 * them (notably MSVC).  The  real + imag * I  pattern corrupts the real
 * part when imag is NaN.
 *
 * Fallback priority:
 *   1. Compiler-provided CMPLX / CMPLXF  (C11, already defined)
 *   2. MSVC intrinsics  (_Cbuild / _FCbuild)
 *   3. __builtin_complex  (GCC >= 4.7, Clang)
 *   4. Union type-pun     (portable last resort)
 */
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


/* Integer type for BLAS/LAPACK interface */
#ifdef HAVE_BLAS_ILP64
    typedef int64_t ARNAUD_INT;
#else
    typedef int32_t ARNAUD_INT;
#endif


#if defined(_MSC_VER)
    // MSVC definition
    typedef _Fcomplex ARNAUD_CPLXF_TYPE;
    typedef _Dcomplex ARNAUD_CPLX_TYPE;
    #define ARNAUD_cplxf(real, imag) ((_Fcomplex){real, imag})
    #define ARNAUD_cplx(real, imag) ((_Dcomplex){real, imag})
#else
    // C99 compliant compilers
    typedef float complex ARNAUD_CPLXF_TYPE;
    typedef double complex ARNAUD_CPLX_TYPE;
    #define ARNAUD_cplxf(real, imag) CMPLXF(real, imag)
    #define ARNAUD_cplx(real, imag) CMPLX(real, imag)
#endif


#endif
