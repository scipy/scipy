/*
  Generic functions for computing 1D convolutions of periodic sequences.

  Supported FFT libraries:
    DJBFFT  - optional, used for power-of-two length arrays
    FFTW    - optional
    FFTPACK - used if any of the above libraries is not available

  Author: Pearu Peterson, September 2002
 */

#include "fftpack.h"

/**************** FFTW *****************************/
#ifdef WITH_FFTW
#include "fftw/convolve.cxx"

#ifndef WITH_DJBFFT
extern "C" void destroy_convolve_cache(void) 
{
	destroy_convolve_cache_fftw();
}

extern "C" void convolve(int n,double* inout,double* omega,int swap_real_imag) 
{
	convolve_fftw(n, inout, omega, swap_real_imag); 
}

extern "C" void convolve_z(int n,double* inout,double* omega_real,double* omega_imag) 
{
	convolve_z_fftw(n, inout, omega_real, omega_imag);
}

extern "C" void init_convolution_kernel(int n,double* omega, int d,
			     double (*kernel_func)(int),
			     int zero_nyquist) 
{
	init_convolution_kernel_fftw(n, omega, d, kernel_func, zero_nyquist);
}
#endif

#else
/**************** FFTPACK ZFFT **********************/
#include "fftpack/convolve.cxx"

#ifndef WITH_DJBFFT
extern "C" void destroy_convolve_cache(void) 
{
	destroy_convolve_cache_fftpack();
}

extern "C" void convolve(int n,double* inout,double* omega,int swap_real_imag) 
{
	convolve_fftpack(n, inout, omega, swap_real_imag); 
}

extern "C" void convolve_z(int n,double* inout,double* omega_real,double* omega_imag) 
{
	convolve_z_fftpack(n, inout, omega_real, omega_imag);
}

extern "C" void init_convolution_kernel(int n,double* omega, int d,
			     double (*kernel_func)(int),
			     int zero_nyquist) 
{
	init_convolution_kernel_fftpack(n, omega, d, kernel_func, zero_nyquist);
}
#endif

#endif

#ifdef WITH_DJBFFT
	#include "djbfft/convolve.cxx"
#endif
