#ifndef _SCIPY_FFTPACK_FFTW_API_H_
#define _SCIPY_FFTPACK_FFTW_API_H_

#include "fftpack.h"

/*
 * straight FFT api
 */
void drfft_fftw(double * inout, int n, int direction, int howmany, 
                  int normalize);

void zfft_fftw(complex_double * inout,
		 int n, int direction, int howmany, int normalize);

void zfftnd_fftw(complex_double * inout, int rank,
			 int *dims, int direction, int howmany,
			 int normalize);

/*
 * Convolution api
 */
void convolve_fftw(int n, double *inout, double *omega, int swap_real_imag);
void convolve_z_fftw(int n, double *inout, double *omega_real, 
                       double* omega_imag);

void init_convolution_kernel_fftw(int n, double *omega, int d,
				 double (*kernel_func) (int),
				 int zero_nyquist);

#endif
