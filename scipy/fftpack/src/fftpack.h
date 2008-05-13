/*
  Interface to various FFT libraries.
  Author: Pearu Peterson, August 2002
 */

#ifndef FFTPACK_H
#define FFTPACK_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

#ifdef SCIPY_DJBFFT_H
#define WITH_DJBFFT
#endif

#ifdef SCIPY_MKL_H
#define WITH_MKL
#endif

#ifdef SCIPY_FFTW3_H
#define WITH_FFTW3
#endif

#ifdef SCIPY_FFTW_H
#define WITH_FFTW
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

#define COPYRFFTW2STD(SRC,DEST,N) { \
  int j,n2=(N)/2; \
  *(DEST) = *(SRC); \
  for (j=1;j<n2;++j) { \
    *(DEST+2*j-1) = *(SRC+j); \
    *(DEST+2*j) = *(SRC+(N)-j); \
  } \
  if (N>1) { \
    *(DEST+2*n2-1) = *(SRC+n2); \
    if ((N)%2) \
      *(DEST+2*n2) = *(SRC+(N)-n2); \
  } \
}
#define COPYINVRFFTW2STD(SRC,DEST,N) { \
  int j,n2=(N)/2; \
  *(DEST) = *(SRC); \
  for (j=1;j<n2;++j) { \
    *(DEST+j) = *(SRC+2*j-1); \
    *(DEST+(N)-j) = *(SRC+2*j); \
  } \
  if (N>1) {\
    *(DEST+n2) = *(SRC+2*n2-1); \
    if ((N)%2) \
      *(DEST+(N)-n2) = *(SRC+2*n2); \
  } \
}

#endif
