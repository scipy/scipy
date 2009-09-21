/*
  Interface to various FFT libraries.
  Double complex FFT and IFFT.
  Author: Pearu Peterson, August 2002
 */

#include "fftpack.h"

extern void F_FUNC(zfftf,ZFFTF)(int*,double*,double*);
extern void F_FUNC(zfftb,ZFFTB)(int*,double*,double*);
extern void F_FUNC(zffti,ZFFTI)(int*,double*);
extern void F_FUNC(cfftf,CFFTF)(int*,float*,float*);
extern void F_FUNC(cfftb,CFFTB)(int*,float*,float*);
extern void F_FUNC(cffti,CFFTI)(int*,float*);

GEN_CACHE(zfft,(int n)
	  ,double* wsave;
	  ,(caches_zfft[i].n==n)
	  ,caches_zfft[id].wsave = (double*)malloc(sizeof(double)*(4*n+15));
	   F_FUNC(zffti,ZFFTI)(&n,caches_zfft[id].wsave);
	  ,free(caches_zfft[id].wsave);
	  ,10)

GEN_CACHE(cfft,(int n)
	  ,float* wsave;
	  ,(caches_cfft[i].n==n)
	  ,caches_cfft[id].wsave = (float*)malloc(sizeof(float)*(4*n+15));
	   F_FUNC(cffti,CFFTI)(&n,caches_cfft[id].wsave);
	  ,free(caches_cfft[id].wsave);
	  ,10)

void zfft(complex_double * inout, int n, int direction, int howmany,
		int normalize)
{
	int i;
	complex_double *ptr = inout;
	double *wsave = NULL;

	wsave = caches_zfft[get_cache_id_zfft(n)].wsave;

	switch (direction) {
	case 1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			F_FUNC(zfftf,ZFFTF)(&n, (double *) (ptr), wsave);
		}
		break;

	case -1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			F_FUNC(zfftb,ZFFTB)(&n, (double *) (ptr), wsave);
		}
		break;
	default:
		fprintf(stderr, "zfft: invalid direction=%d\n", direction);
	}

	if (normalize) {
		ptr = inout;
		for (i = n * howmany - 1; i >= 0; --i) {
			*((double *) (ptr)) /= n;
			*((double *) (ptr++) + 1) /= n;
		}
	}
}

void cfft(complex_float * inout, int n, int direction, int howmany,
	int normalize)
{
	int i;
	complex_float *ptr = inout;
	float *wsave = NULL;

	wsave = caches_cfft[get_cache_id_cfft(n)].wsave;

	switch (direction) {
	case 1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			F_FUNC(cfftf, CFFTF)(&n, (float *) (ptr), wsave);

		}
		break;

	case -1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			F_FUNC(cfftb, CFFTB)(&n, (float *) (ptr), wsave);
		}
		break;
	default:
		fprintf(stderr, "cfft: invalid direction=%d\n", direction);
	}

	if (normalize) {
		ptr = inout;
		for (i = n * howmany - 1; i >= 0; --i) {
			*((float *) (ptr)) /= n;
			*((float *) (ptr++) + 1) /= n;
		}
	}
}
