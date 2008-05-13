#ifndef _SCIPY_FFTPACK_FFTW3_API_H_
#define _SCIPY_FFTPACK_FFTW3_API_H_

#include "fftpack.h"

/*
 * straight FFT api
 */
void drfft_fftw3(double * inout, int n, int direction, int howmany, 
                  int normalize);

void zfft_fftw3(complex_double * inout,
		 int n, int direction, int howmany, int normalize);

void zfftnd_fftw3(complex_double * inout, int rank,
			 int *dims, int direction, int howmany,
			 int normalize);

#endif
