/*
 * Last Change: Wed Aug 01 07:00 PM 2007 J
 *
 * FFTPACK implementation
 *
 * Original code by Pearu Peterson.
 */

extern void F_FUNC(dfftf, DFFTF) (int *, double *, double *);
extern void F_FUNC(dfftb, DFFTB) (int *, double *, double *);
extern void F_FUNC(dffti, DFFTI) (int *, double *);
GEN_CACHE(drfftpack, (int n)
	  , double *wsave;
	  , (caches_drfftpack[i].n == n)
	  , caches_drfftpack[id].wsave =
	  (double *) malloc(sizeof(double) * (2 * n + 15));
	  F_FUNC(dffti, DFFTI) (&n, caches_drfftpack[id].wsave);
	  , free(caches_drfftpack[id].wsave);
	  , 10)

static void drfft_fftpack(double *inout, int n, int direction, int howmany,
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
