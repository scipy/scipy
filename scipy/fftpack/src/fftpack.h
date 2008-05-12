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
void init_convolution_kernel(int n,double* omega, int d, 
			     double (*kernel_func)(int),
			     int zero_nyquist);
void convolve(int n,double* inout,double* omega,int swap_real_imag);
void convolve_z(int n,double* inout,double* omega_real,double* omega_imag);
#ifdef __cplusplus
};
#endif

extern int ispow2le2e30(int n);
extern int ispow2le2e13(int n);

#ifdef SCIPY_FFTWORK_H
#define WITH_FFTWORK
#include "fftwork/fast_header.h"
#endif

#ifdef SCIPY_DJBFFT_H
#ifdef __cplusplus
extern "C" {
#endif
#define WITH_DJBFFT
#define complex8 complex_double
#define COMPLEX8_H
#include <fftfreq.h>
#include <fftc8.h>
#include <fftr8.h>
#ifdef __cplusplus
}
#endif
#endif

#ifdef SCIPY_MKL_H
#define WITH_MKL
#include <mkl_dfti.h>
#endif

#ifdef SCIPY_FFTW3_H
#define WITH_FFTW3
#include <fftw3.h>
#endif

#ifdef SCIPY_DFFTW_H
#define WITH_FFTW
#include <dfftw.h>
#include <drfftw.h>
#endif

#ifdef SCIPY_FFTW_H
#define WITH_FFTW
#include <fftw.h>
#include <rfftw.h>
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

#define	COPYSTD2DJB(SRC,DEST,N) { \
  int n2 = (N)/2,k,j; \
  *(DEST) = *(SRC); \
  *(DEST+1) = *(SRC+n2); \
  for (j=(N)/2-1,k=2;j>0;--j,k+=2) { \
    *(DEST+k) = *(SRC+n2+j); \
    *(DEST+k+1) = *(SRC+j); \
  } \
}

#define	COPYINVDJB2STD(SRC,DEST,N) { \
  int n2 = (N)/2,k,j; \
  *(DEST) = *(SRC); \
  *(DEST+n2) = *(SRC+1); \
  for (j=(N)/2-1,k=2;j>0;--j,k+=2) { \
    *(DEST+n2+j) = *(SRC+k); \
    *(DEST+j) = *(SRC+k+1); \
  } \
}

#define	COPYINVDJB2STD2(SRC,DEST,N) { \
  int n2 = (N)/2,k,j; \
  *(DEST) = *(SRC); \
  *(DEST+(N)-1) = *(SRC+(N)-1); \
  for (j=1,k=1;j<n2;++j,k+=2) { \
    *(DEST+n2+j-1) = *(SRC+k); \
    *(DEST+j) = *(SRC+k+1); \
  } \
}

#define COPYDJB2STD(SRC,DEST,FRQ,N) { \
  int n2 = (N)/2,k,j; \
  *(DEST) = *(SRC); \
  *(DEST+N-1) = *(SRC+1); \
  for (k=2;k<N-1;k+=2) { \
    j = FRQ[k]; \
    if (j>n2) { \
      j = 2*(N-j); \
      *(DEST+j-1) = *(SRC+k); \
      *(DEST+j) = -*(SRC+k+1); \
    } else { \
      j *= 2; \
      *(DEST+j-1) = *(SRC+k); \
      *(DEST+j) = *(SRC+k+1); \
    } \
  } \
}
#define COPYINVSTD2DJB(SRC,DEST,NORMALIZE,FRQ,N) { \
  int n2 = (N)/2,k,j; \
  if (NORMALIZE) { \
    *(DEST) = *(SRC); \
    *(DEST+1) = *(SRC+N-1); \
  } else { \
    *(DEST) = (*(SRC))*0.5; \
    *(DEST+1) = (*(SRC+N-1))*0.5; \
  } \
  for (k=2;k<N-1;k+=2) { \
    j = FRQ[k]; \
    if (j>n2) { \
      j = 2*(N-j); \
      *(DEST+k) = *(SRC+j-1); \
      *(DEST+k+1) = -*(SRC+j); \
    } else { \
      j *= 2; \
      *(DEST+k) = *(SRC+j-1); \
      *(DEST+k+1) = *(SRC+j); \
    } \
  } \
}
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
