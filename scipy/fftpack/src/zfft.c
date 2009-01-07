/*
  Interface to various FFT libraries.
  Double complex FFT and IFFT.
  Author: Pearu Peterson, August 2002
 */

#include "fftpack.h"

/* The following macro convert private backend specific function to the public
 * functions exported by the module  */
#define GEN_PUBLIC_API(name) \
void destroy_zfft_cache(void)\
{\
        destroy_z##name##_caches();\
}\
\
void zfft(complex_double *inout, int n, \
        int direction, int howmany, int normalize)\
{\
        zfft_##name(inout, n, direction, howmany, normalize);\
}\
void destroy_cfft_cache(void)\
{\
        destroy_c##name##_caches();\
}\
\
void cfft(complex_float *inout, int n, \
        int direction, int howmany, int normalize)\
{\
        cfft_##name(inout, n, direction, howmany, normalize);\
}

#include "zfft_fftpack.c"
#include "cfft_fftpack.c"
GEN_PUBLIC_API(fftpack)
