#ifndef _SCIPY_FFTPACK_FFTPACK_API_H_
#define _SCIPY_FFTPACK_FFTPACK_API_H_

#include "fftpack.h"

void drfft_fftpack(double * inout, int n, int direction, int howmany, 
                  int normalize);

void zfft_fftpack(complex_double * inout,
		 int n, int direction, int howmany, int normalize);

void zfftnd_fftpack(complex_double * inout, int rank,
			 int *dims, int direction, int howmany,
			 int normalize);

#endif
