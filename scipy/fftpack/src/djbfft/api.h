#ifndef _SCIPY_FFTPACK_DJBFFT_API_H_
#define _SCIPY_FFTPACK_DJBFFT_API_H_

#include "fftpack.h"

/*
 * straight FFT api
 */
void drfft_djbfft(double * inout, int n, int direction, int howmany, 
                  int normalize);

void zfft_djbfft(complex_double * inout,
		 int n, int direction, int howmany, int normalize);

/*
 * Convolution api
 */
void convolve_djbfft(int n, double *inout, double *omega, int swap_real_imag);
void convolve_z_djbfft(int n, double *inout, double *omega_real, 
                       double* omega_imag);

void init_convolution_kernel_djbfft(int n, double *omega, int d,
				 double (*kernel_func) (int),
				 int zero_nyquist);

/*
 * Common headers and def
 */
#define complex8 complex_double
#define COMPLEX8_H

extern "C" {
#include <fftfreq.h>
#include <fftc8.h>
#include <fftr8.h>
};

#endif
