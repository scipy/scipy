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
}

/* ************** Definition of backend specific functions ********* */

/*
 * To add a backend :
 *  - create a file drfft_name.c, where you define a function drfft_name where
 *  name is the name of your backend. If you do not use the GEN_CACHE macro,
 *  you will need to define a function void destroy_drname_caches(void), 
 *  which can do nothing
 *  - in drfft.c, include the drfft_name.c file, and add the 3 following lines
 *  just after it:
 *  #ifndef WITH_DJBFFT
 *      GEN_PUBLIC_API(name)
 *  #endif
 */

#ifdef WITH_FFTW3
    #include "drfft_fftw3.c"
    #ifndef WITH_DJBFFT
        GEN_PUBLIC_API(fftw3)
    #endif
#elif defined WITH_FFTW
    #include "drfft_fftw.c"
    #ifndef WITH_DJBFFT
        GEN_PUBLIC_API(fftw)
    #endif
#else /* Use fftpack by default */
    #include "drfft_fftpack.c"
    #ifndef WITH_DJBFFT
        GEN_PUBLIC_API(fftpack)
    #endif 
#endif

/* 
 * djbfft must be used at the end, because it needs another backend (defined
 * above) for non 2^n * size 
 */
#ifdef WITH_DJBFFT
    #include "drfft_djbfft.c"
    void destroy_drfft_cache(void)
    {
        destroy_drdjbfft_caches();
        drfft_def_destroy_cache();
    }
    void drfft(double *inout, int n, 
            int direction, int howmany, int normalize)
    {
        drfft_djbfft(inout, n, direction, howmany, normalize);
    }
#endif
