/*
 * fftw2 backend for multi dimensional fft
 *
 * Original code by Pearu Peaterson
 *
 * Last Change: Thu Sep 06 05:00 PM 2007 J
 */

GEN_CACHE(zfftnd_fftw, (int n, int *dims, int d, int flags)
	  , int direction;
	  int *dims;
	  fftwnd_plan plan;, ((caches_zfftnd_fftw[i].n == n) &&
			      (caches_zfftnd_fftw[i].direction == d) &&
			      (equal_dims
			       (n, caches_zfftnd_fftw[i].dims, dims)))
	  , caches_zfftnd_fftw[id].direction = d;
	  caches_zfftnd_fftw[id].n = n;
	  caches_zfftnd_fftw[id].dims = (int *) malloc(sizeof(int) * n);
	  memcpy(caches_zfftnd_fftw[id].dims, dims, sizeof(int) * n);
	  caches_zfftnd_fftw[id].plan =
	  fftwnd_create_plan(n, dims,
			     (d > 0 ? FFTW_FORWARD : FFTW_BACKWARD),
			     flags);,
	  fftwnd_destroy_plan(caches_zfftnd_fftw[id].plan);
	  free(caches_zfftnd_fftw[id].dims);, 10)


extern void zfftnd_fftw(complex_double * inout, int rank,
		       int *dims, int direction, int howmany,
		       int normalize)
{
    int i, sz;
    complex_double *ptr = inout;
    fftwnd_plan plan = NULL;

    sz = 1;
    for (i = 0; i < rank; ++i) {
	sz *= dims[i];
    }
    i = get_cache_id_zfftnd_fftw(rank, dims, direction,
				 FFTW_IN_PLACE | FFTW_ESTIMATE);
    plan = caches_zfftnd_fftw[i].plan;
    for (i = 0; i < howmany; ++i, ptr += sz) {
	fftwnd_one(plan, (fftw_complex *) ptr, NULL);
    }
    if (normalize) {
	ptr = inout;
	for (i = sz * howmany - 1; i >= 0; --i) {
	    *((double *) (ptr)) /= sz;
	    *((double *) (ptr++) + 1) /= sz;
	}
    }
}
