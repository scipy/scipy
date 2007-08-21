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
}

#if defined(WITH_FFTW) || defined(WITH_MKL)
static
int equal_dims(int rank,int *dims1,int *dims2) {
  int i;
  for (i=0;i<rank;++i)
    if (dims1[i]!=dims2[i])
      return 0;
  return 1;
}
#endif

#ifdef WITH_FFTW3
    #include "zfftnd_fftw3.c"
    GEN_PUBLIC_API(fftw3)
#elif defined WITH_FFTW
    #include "zfftnd_fftw.c"
    GEN_PUBLIC_API(fftw)
#elif defined WITH_MKL
    #include "zfftnd_mkl.c"
    GEN_PUBLIC_API(mkl)
#else /* Use fftpack by default */
    #include "zfftnd_fftpack.c"
    GEN_PUBLIC_API(fftpack)
#endif
