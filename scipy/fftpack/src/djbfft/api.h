#ifndef _SCIPY_FFTPACK_DJBFFT_API_H_
#define _SCIPY_FFTPACK_DJBFFT_API_H_

#include "fftpack.h"

void drfft_djbfft(double * inout, int n, int direction, int howmany, 
                  int normalize);

void zfft_djbfft(complex_double * inout,
		 int n, int direction, int howmany, int normalize);

#define complex8 complex_double
#define COMPLEX8_H

extern "C" {
#include <fftfreq.h>
#include <fftc8.h>
#include <fftr8.h>
};

#endif
