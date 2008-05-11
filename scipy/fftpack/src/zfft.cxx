/*
  Interface to various FFT libraries.
  Double complex FFT and IFFT.
  Author: Pearu Peterson, August 2002
 */

#include "fftpack.h"

/* The following macro convert private backend specific function to the public
 * functions exported by the module  */
#define GEN_PUBLIC_API(name) \
extern "C" void zfft(complex_double *inout, int n, \
        int direction, int howmany, int normalize)\
{\
        zfft_##name(inout, n, direction, howmany, normalize);\
}

/* ************** Definition of backend specific functions ********* */

#ifdef WITH_FFTW3
    #include "fftw3/zfft.cxx"
    #ifndef WITH_DJBFFT
        GEN_PUBLIC_API(fftw3)
    #endif
#elif defined WITH_FFTW
    #include "fftw/zfft.cxx"
    #ifndef WITH_DJBFFT
        GEN_PUBLIC_API(fftw)
    #endif
#elif defined WITH_MKL
    #include "mkl/zfft.cxx"
    #ifndef WITH_DJBFFT
        GEN_PUBLIC_API(mkl)
    #endif
#else /* Use fftpack by default */
    #include "fftpack/zfft.cxx"
    #ifndef WITH_DJBFFT
        GEN_PUBLIC_API(fftpack)
    #endif 
#endif

/* 
 * djbfft must be used at the end, because it needs another backend (defined
 * above) for non 2^n * size 
 */
#ifdef WITH_DJBFFT
    #include "djbfft/zfft.cxx"
    extern "C" void zfft(complex_double *inout, int n, 
            int direction, int howmany, int normalize)
    {
        zfft_djbfft(inout, n, direction, howmany, normalize);
    }
#endif
