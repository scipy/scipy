#ifndef ARNAUD_N_DOUBLE_COMPLEX_H
#define ARNAUD_N_DOUBLE_COMPLEX_H

#include "arnaud/arnaud.h"
#include "blaslapack_declarations.h"
#include <complex.h>

#if defined(_MSC_VER)
    // MSVC definitions
    typedef _Dcomplex ARNAUD_CPLX_TYPE;
    #define ARNAUD_cplx(real, imag) ((_Dcomplex){real, imag})

#else
    // C99 compliant compilers
    typedef double complex ARNAUD_CPLX_TYPE;
    #define ARNAUD_cplx(real, imag) ((real) + (imag)*I)

#endif

#endif
