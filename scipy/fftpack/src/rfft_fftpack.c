/*
 * Last Change: Wed Aug 01 07:00 PM 2007 J
 *
 * FFTPACK implementation
 *
 * Original code by Pearu Peterson.
 */

extern void F_FUNC(rfftf, RFFTF) (int *, float *, float *);
extern void F_FUNC(rfftb, RFFTB) (int *, float *, float *);
extern void F_FUNC(rffti, RFFTI) (int *, float *);
GEN_CACHE(rfftpack, (int n)
	  , float *wsave;
	  , (caches_rfftpack[i].n == n)
	  , caches_rfftpack[id].wsave =
	  (float *) malloc(sizeof(float) * (2 * n + 15));
	  F_FUNC(rffti, RFFTI) (&n, caches_rfftpack[id].wsave);
	  , free(caches_rfftpack[id].wsave);
	  , 10)

static void rfft_fftpack(float *inout, int n, int direction, int howmany,
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
