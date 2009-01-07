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


GEN_CACHE(drfftpack, (int n)
	  , double *wsave;
	  , (caches_drfftpack[i].n == n)
	  , caches_drfftpack[id].wsave =
	  (double *) malloc(sizeof(double) * (2 * n + 15));
	  F_FUNC(dffti, DFFTI) (&n, caches_drfftpack[id].wsave);
	  , free(caches_drfftpack[id].wsave);
	  , 10)

GEN_CACHE(rfftpack, (int n)
	  , float *wsave;
	  , (caches_rfftpack[i].n == n)
	  , caches_rfftpack[id].wsave =
	  (float *) malloc(sizeof(float) * (2 * n + 15));
	  F_FUNC(rffti, RFFTI) (&n, caches_rfftpack[id].wsave);
	  , free(caches_rfftpack[id].wsave);
	  , 10)

void drfft_fftpack(double *inout, int n, int direction, int howmany,
			  int normalize)
{
    int i;
    double *ptr = inout;
    double *wsave = NULL;
    wsave = caches_drfftpack[get_cache_id_drfftpack(n)].wsave;


    switch (direction) {
        case 1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            dfftf_(&n, ptr, wsave);
        }
        break;

    case -1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            dfftb_(&n, ptr, wsave);
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

void rfft_fftpack(float *inout, int n, int direction, int howmany,
			 int normalize)
{
    int i;
    float *ptr = inout;
    float *wsave = NULL;
    wsave = caches_rfftpack[get_cache_id_rfftpack(n)].wsave;


    switch (direction) {
        case 1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            rfftf_(&n, ptr, wsave);
        }
        break;

    case -1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            rfftb_(&n, ptr, wsave);
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
