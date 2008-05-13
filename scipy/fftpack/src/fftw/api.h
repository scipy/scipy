#ifndef _SCIPY_FFTPACK_FFTW_API_H_
#define _SCIPY_FFTPACK_FFTW_API_H_

#include "fftpack.h"

void drfft_fftw(double * inout, int n, int direction, int howmany, 
                  int normalize);

void zfft_fftw(complex_double * inout,
		 int n, int direction, int howmany, int normalize);

void zfftnd_fftw(complex_double * inout, int rank,
			 int *dims, int direction, int howmany,
			 int normalize);

#endif
