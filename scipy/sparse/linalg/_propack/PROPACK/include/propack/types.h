#ifndef PROPACK_TYPES_H
#define PROPACK_TYPES_H

#include <stdint.h>
#include "scipy_complex_support.h"
#include "scipy_blas_defines.h"


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
    #define PROPACK_cplxf(real, imag) CMPLXF(real, imag)
    #define PROPACK_cplx(real, imag) CMPLX(real, imag)
#endif


// Function pointer typedefs for aprod callbacks
typedef void (*PROPACK_aprod_s)(int transa, CBLAS_INT m, CBLAS_INT n, float* x, float* y, float* dparm, CBLAS_INT* iparm);
typedef void (*PROPACK_aprod_d)(int transa, CBLAS_INT m, CBLAS_INT n, double* x, double* y, double* dparm, CBLAS_INT* iparm);
typedef void (*PROPACK_aprod_c)(int transa, CBLAS_INT m, CBLAS_INT n, PROPACK_CPLXF_TYPE* x, PROPACK_CPLXF_TYPE* y, PROPACK_CPLXF_TYPE* cparm, CBLAS_INT* iparm);
typedef void (*PROPACK_aprod_z)(int transa, CBLAS_INT m, CBLAS_INT n, PROPACK_CPLX_TYPE* x, PROPACK_CPLX_TYPE* y, PROPACK_CPLX_TYPE* zparm, CBLAS_INT* iparm);


#endif
