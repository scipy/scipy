/*
  Interface to various FFT libraries.
  Double real FFT and IFFT.
  Author: Pearu Peterson, August 2002
 */

#include "fftpack.h"

extern void F_FUNC(dfftf, DFFTF) (int *, double *, double *);
extern void F_FUNC(dfftb, DFFTB) (int *, double *, double *);
extern void F_FUNC(dffti, DFFTI) (int *, double *);
extern void F_FUNC(rfftf, RFFTF) (int *, float *, float *);
extern void F_FUNC(rfftb, RFFTB) (int *, float *, float *);
extern void F_FUNC(rffti, RFFTI) (int *, float *);


GEN_CACHE(drfft, (int n)
	  , double *wsave;
	  , (caches_drfft[i].n == n)
	  , caches_drfft[id].wsave =
	  (double *) malloc(sizeof(double) * (2 * n + 15));
	  F_FUNC(dffti, DFFTI) (&n, caches_drfft[id].wsave);
	  , free(caches_drfft[id].wsave);
	  , 10)

GEN_CACHE(rfft, (int n)
	  , float *wsave;
	  , (caches_rfft[i].n == n)
	  , caches_rfft[id].wsave =
	  (float *) malloc(sizeof(float) * (2 * n + 15));
	  F_FUNC(rffti, RFFTI) (&n, caches_rfft[id].wsave);
	  , free(caches_rfft[id].wsave);
	  , 10)

void drfft(double *inout, int n, int direction, int howmany,
			  int normalize)
{
    int i;
    double *ptr = inout;
    double *wsave = NULL;
    wsave = caches_drfft[get_cache_id_drfft(n)].wsave;


    switch (direction) {
        case 1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            F_FUNC(dfftf,DFFTF)(&n, ptr, wsave);
        }
        break;

    case -1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            F_FUNC(dfftb,DFFTB)(&n, ptr, wsave);
        }
        break;

    default:
        fprintf(stderr, "drfft: invalid direction=%d\n", direction);
    }

    if (normalize) {
        double d = 1.0 / n;
        ptr = inout;
        for (i = n * howmany - 1; i >= 0; --i) {
            (*(ptr++)) *= d;
        }
    }
}

void rfft(float *inout, int n, int direction, int howmany,
			 int normalize)
{
    int i;
    float *ptr = inout;
    float *wsave = NULL;
    wsave = caches_rfft[get_cache_id_rfft(n)].wsave;


    switch (direction) {
        case 1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            F_FUNC(rfftf,RFFTF)(&n, ptr, wsave);
        }
        break;

    case -1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            F_FUNC(rfftb,RFFTB)(&n, ptr, wsave);
        }
        break;

    default:
        fprintf(stderr, "rfft: invalid direction=%d\n", direction);
    }

    if (normalize) {
        float d = 1.0 / n;
        ptr = inout;
        for (i = n * howmany - 1; i >= 0; --i) {
            (*(ptr++)) *= d;
        }
    }
}
