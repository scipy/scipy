/*
 * Last Change: Wed Aug 01 07:00 PM 2007 J
 *
 * FFTW2 implementation
 *
 * Original code by Pearu Peterson.
 */

GEN_CACHE(drfftw, (int n, int d, int flags)
	  , int direction;
	  int flags;
	  rfftw_plan plan;
	  double *ptr;, ((caches_drfftw[i].n == n) &&
			 (caches_drfftw[i].direction == d) &&
			 (caches_drfftw[i].flags == flags))
	  , caches_drfftw[id].direction = d;
	  caches_drfftw[id].flags = flags;
	  caches_drfftw[id].plan = rfftw_create_plan(n,
						     (d >
						      0 ?
						      FFTW_REAL_TO_COMPLEX
						      :
						      FFTW_COMPLEX_TO_REAL),
						     flags);
	  caches_drfftw[id].ptr =
	  (double *) malloc(sizeof(double) * (n));,
	  rfftw_destroy_plan(caches_drfftw[id].plan);
	  free(caches_drfftw[id].ptr);, 10)

static void drfft_fftw(double *inout, int n, int dir, int
			howmany, int normalize)
{
    int i;
    double *ptr = inout;
    double *ptrc = NULL;
    rfftw_plan plan = NULL;

    i = get_cache_id_drfftw(n, dir, FFTW_IN_PLACE | FFTW_ESTIMATE);
    plan = caches_drfftw[i].plan;
    ptrc = caches_drfftw[i].ptr;

    switch (dir) {
    case 1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            memcpy(ptrc, ptr, sizeof(double) * n);
            rfftw(plan, 1, (fftw_real *) ptrc, 1, 1, NULL, 1, 1);
            COPYRFFTW2STD(ptrc, ptr, n);
        }
        break;

    case -1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            COPYINVRFFTW2STD(ptr, ptrc, n);
            rfftw(plan, 1, (fftw_real *) ptrc, 1, 1, NULL, 1, 1);
            memcpy(ptr, ptrc, sizeof(double) * n);
        }
        break;

    default:
        fprintf(stderr, "drfft: invalid direction=%d\n", dir);
    }

    if (normalize) {
        double d = 1.0 / n;
        ptr = inout;
        for (i = n * howmany - 1; i >= 0; --i) {
            (*(ptr++)) *= d;
        }
    }
}
