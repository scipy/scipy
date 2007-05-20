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
}

/* ************** Definition of backend specific functions ********* */

/*
 * To add a backend :
 *  - create a file zfft_name.c, where you define a function zfft_name where
 *  name is the name of your backend. If you do not use the GEN_CACHE macro,
 *  you will need to define a function void destroy_zname_caches(void), 
 *  which can do nothing
 *  - in zfft.c, include the zfft_name.c file, and add the 3 following lines
 *  just after it:
 *  #ifndef WITH_DJBFFT
 *      GEN_PUBLIC_API(name)
 *  #endif
 */

#ifdef WITH_FFTW3
    #include "zfft_fftw3.c"
    #ifndef WITH_DJBFFT
        GEN_PUBLIC_API(fftw3)
    #endif
#elif defined WITH_FFTW
    #include "zfft_fftw.c"
    #ifndef WITH_DJBFFT
        GEN_PUBLIC_API(fftw)
    #endif
#elif defined WITH_MKL
    #include "zfft_mkl.c"
    #ifndef WITH_DJBFFT
        GEN_PUBLIC_API(mkl)
    #endif
#else /* Use fftpack by default */
    #include "zfft_fftpack.c"
    #ifndef WITH_DJBFFT
        GEN_PUBLIC_API(fftpack)
    #endif 
#endif

/* 
 * djbfft must be used at the end, because it needs another backend (defined
 * above) for non 2^n * size 
 */
#ifdef WITH_DJBFFT
    #include "zfft_djbfft.c"
    void destroy_zfft_cache(void)
    {
        destroy_zdjbfft_caches();
        zfft_def_destroy_cache();
    }
    void zfft(complex_double *inout, int n, 
            int direction, int howmany, int normalize)
    {
        zfft_djbfft(inout, n, direction, howmany, normalize);
    }
#endif
