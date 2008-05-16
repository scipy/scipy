#ifndef _SCIPY_FFTPACK_MKL_API_H_
#define _SCIPY_FFTPACK_MKL_API_H_

#include "misc.h"

/*
 * straight FFT api
 */
extern "C" {
void zfft_mkl(complex_double * inout,
		 int n, int direction, int howmany, int normalize);

void zfftnd_mkl(complex_double * inout, int rank,
			 int *dims, int direction, int howmany,
			 int normalize);
};

#endif
