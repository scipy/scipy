#ifndef ARNAUD_TYPES_H
#define ARNAUD_TYPES_H

#include <stdint.h>
#include <complex.h>


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
    #define ARNAUD_cplxf(real, imag) ((real) + (imag)*I)
    #define ARNAUD_cplx(real, imag) ((real) + (imag)*I)
#endif


#endif
