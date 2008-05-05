/*
 * Last Change: Wed Aug 01 07:00 PM 2007 J
 *
 * FFTW3 implementation
 *
 * Original code by Pearu Peterson.
 */

GEN_CACHE(drfftw3, (int n, int d, int flags)
	  , int direction;
	  int flags;
	  fftw_plan plan;
	  double *ptr;, ((caches_drfftw3[i].n == n) &&
			 (caches_drfftw3[i].direction == d) &&
			 (caches_drfftw3[i].flags == flags))
	  , caches_drfftw3[id].direction = d;
	  caches_drfftw3[id].flags = flags;
	  caches_drfftw3[id].ptr =
	  (double *) fftw_malloc(sizeof(double) * (n));
	  caches_drfftw3[id].plan =
	  fftw_plan_r2r_1d(n, caches_drfftw3[id].ptr, caches_drfftw3[id].ptr,
			   (d > 0 ? FFTW_R2HC : FFTW_HC2R), flags);,
	  fftw_destroy_plan(caches_drfftw3[id].plan);
	  fftw_free(caches_drfftw3[id].ptr);, 10)

static void drfft_fftw3(double *inout, int n, int direction, int
			howmany, int normalize)
{
    int i;
    double *ptr = inout;

    double *ptrc = NULL;
    fftw_plan plan = NULL;

    i = get_cache_id_drfftw3(n, direction, FFTW_ESTIMATE);
    plan = caches_drfftw3[i].plan;
    ptrc = caches_drfftw3[i].ptr;
    switch (direction) {
    case 1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            memcpy(ptrc, ptr, sizeof(double) * n);
            fftw_execute(plan);
            COPYRFFTW2STD(ptrc, ptr, n);
        }
        break;

    case -1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            COPYINVRFFTW2STD(ptr, ptrc, n);
            fftw_execute(plan);
            memcpy(ptr, ptrc, sizeof(double) * n);
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
