/*
  Generic functions for computing 1D convolutions of periodic sequences.

  Supported FFT libraries:
    DJBFFT  - optional, used for power-of-two length arrays
    FFTW    - optional
    FFTPACK - used if any of the above libraries is not available

  Author: Pearu Peterson, September 2002
 */

#include "fftpack.h"

#define GEN_CONVOLVE_API(name) \
extern "C" void convolve(int n,double* inout,double* omega,int swap_real_imag)  \
{\
        convolve_##name(n, inout, omega, swap_real_imag);\
}\
extern "C" void convolve_z(int n,double* inout,double* omega_real,double* omega_imag)  \
{\
        convolve_z_##name(n, inout, omega_real, omega_imag);\
}\
extern "C" void init_convolution_kernel(int n,double* omega, int d, \
			     double (*kernel_func)(int), \
			     int zero_nyquist) \
{\
        init_convolution_kernel_##name(n,omega, d, kernel_func, zero_nyquist);\
} \
extern "C" void destroy_convolve_cache(void)  \
{\
}

/**************** FFTW *****************************/
#ifdef WITH_FFTW

#include "fftw/api.h"
#ifndef WITH_DJBFFT
        GEN_CONVOLVE_API(fftw)
#endif

#else
/**************** FFTPACK ZFFT **********************/
#include "fftpack/api.h"

#ifndef WITH_DJBFFT
        GEN_CONVOLVE_API(fftpack)
#endif

#endif

#ifdef WITH_DJBFFT
	#include "djbfft/api.h"
        GEN_CONVOLVE_API(djbfft)
#endif
