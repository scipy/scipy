#ifndef _SCIPY_FFTPACK_DJBFFT_API_H_
#define _SCIPY_FFTPACK_DJBFFT_API_H_

void drfft_djbfft(double * inout, int n, int direction, int howmany, 
                  int normalize);

void zfft_djbfft(complex_double * inout,
		 int n, int direction, int howmany, int normalize);

#endif
