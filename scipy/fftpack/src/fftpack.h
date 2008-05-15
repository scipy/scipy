/*
  Interface to various FFT libraries.
  Author: Pearu Peterson, August 2002
 */

#ifndef FFTPACK_H
#define FFTPACK_H

#include <stdio.h>
#include <stdlib.h>

typedef struct {double r,i;} complex_double;
typedef struct {float r,i;} complex_float;

#ifdef __cplusplus
extern "C" {
#endif
void init_convolution_kernel_fftpack(int n,double* omega, int d, 
			     double (*kernel_func)(int),
			     int zero_nyquist);
void convolve_fftpack(int n,double* inout,double* omega,int swap_real_imag);
void convolve_z_fftpack(int n,double* inout,double* omega_real,double* omega_imag);

void drfft_fftpack(double *inout, int n, int direction, int howmany,
			  int normalize);
void zfft_fftpack(complex_double * inout,
			 int n, int direction, int howmany, int normalize);
void zfftnd_fftpack(complex_double * inout, int rank,
			   int *dims, int direction, int howmany,
			   int normalize);
#ifdef __cplusplus
};
#endif

#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif

#endif
