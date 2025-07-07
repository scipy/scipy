#ifndef ARNAUD_N_SINGLE_COMPLEX_H
#define ARNAUD_N_SINGLE_COMPLEX_H

#include "arnaud/arnaud.h"
#include "blaslapack_declarations.h"
#include <complex.h>

#if defined(_MSC_VER)
    // MSVC definitions
    typedef _Fcomplex ARNAUD_CPLXF_TYPE;
    #define ARNAUD_cplxf(real, imag) ((_Fcomplex){real, imag})

#else
    // C99 compliant compilers
    typedef float complex ARNAUD_CPLXF_TYPE;
    #define ARNAUD_cplxf(real, imag) ((real) + (imag)*I)

#endif

#endif
