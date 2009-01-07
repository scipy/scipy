/*
  Interface to various FFT libraries.
  Double real FFT and IFFT.
  Author: Pearu Peterson, August 2002
 */

#include "fftpack.h"

/* The following macro convert private backend specific function to the public
 * functions exported by the module  */
#define GEN_PUBLIC_API(name) \
void destroy_drfft_cache(void)\
{\
        destroy_dr##name##_caches();\
}\
\
void drfft(double *inout, int n, \
        int direction, int howmany, int normalize)\
{\
        drfft_##name(inout, n, direction, howmany, normalize);\
}\
void destroy_rfft_cache(void)\
{\
        destroy_r##name##_caches();\
}\
\
void rfft(float *inout, int n, \
        int direction, int howmany, int normalize)\
{\
        rfft_##name(inout, n, direction, howmany, normalize);\
}


#include "drfft_fftpack.c"
#include "rfft_fftpack.c"
GEN_PUBLIC_API(fftpack)
