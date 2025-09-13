#ifndef PROPACK_TYPES_H
#define PROPACK_TYPES_H

#include <stdint.h>
#include <complex.h>


#if defined(_MSC_VER)
    // MSVC definition
    typedef _Fcomplex PROPACK_CPLXF_TYPE;
    typedef _Dcomplex PROPACK_CPLX_TYPE;
    #define PROPACK_cplxf(real, imag) ((_Fcomplex){real, imag})
    #define PROPACK_cplx(real, imag) ((_Dcomplex){real, imag})
#else
    // C99 compliant compilers
    typedef float complex PROPACK_CPLXF_TYPE;
    typedef double complex PROPACK_CPLX_TYPE;
    #define PROPACK_cplxf(real, imag) ((real) + (imag)*I)
    #define PROPACK_cplx(real, imag) ((real) + (imag)*I)
#endif


// Function pointer typedefs for aprod callbacks
typedef void (*PROPACK_aprod_s)(int transa, int m, int n, float* x, float* y, float* dparm, int* iparm);
typedef void (*PROPACK_aprod_d)(int transa, int m, int n, double* x, double* y, double* dparm, int* iparm);
typedef void (*PROPACK_aprod_c)(int transa, int m, int n, PROPACK_CPLXF_TYPE* x, PROPACK_CPLXF_TYPE* y, PROPACK_CPLXF_TYPE* cparm, int* iparm);
typedef void (*PROPACK_aprod_z)(int transa, int m, int n, PROPACK_CPLX_TYPE* x, PROPACK_CPLX_TYPE* y, PROPACK_CPLX_TYPE* zparm, int* iparm);


#endif
