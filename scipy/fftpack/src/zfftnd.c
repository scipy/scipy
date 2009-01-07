/*
  Interface to various FFT libraries.
  Double complex FFT and IFFT, arbitrary dimensions.
  Author: Pearu Peterson, August 2002
 */
#include "fftpack.h"

/* The following macro convert private backend specific function to the public
 * functions exported by the module  */
#define GEN_PUBLIC_API(name) \
void destroy_zfftnd_cache(void)\
{\
        destroy_zfftnd_##name##_caches();\
}\
\
void zfftnd(complex_double * inout, int rank,\
		           int *dims, int direction, int howmany, int normalize)\
{\
        zfftnd_##name(inout, rank, dims, direction, howmany, normalize);\
}\
void destroy_cfftnd_cache(void)\
{\
        destroy_cfftnd_##name##_caches();\
}\
\
void cfftnd(complex_float * inout, int rank,\
		           int *dims, int direction, int howmany, int normalize)\
{\
        cfftnd_##name(inout, rank, dims, direction, howmany, normalize);\
}


#include "zfftnd_fftpack.c"
#include "cfftnd_fftpack.c"
GEN_PUBLIC_API(fftpack)
