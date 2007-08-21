GEN_CACHE(zfftw,(int n,int d)
	  ,int direction;
	   fftw_plan plan;
	  ,((caches_zfftw[i].n==n) && 
	    (caches_zfftw[i].direction==d))
	  ,caches_zfftw[id].direction = d;
	   caches_zfftw[id].plan = fftw_create_plan(n,
		(d>0?FFTW_FORWARD:FFTW_BACKWARD),
		FFTW_IN_PLACE|FFTW_ESTIMATE);
	  ,fftw_destroy_plan(caches_zfftw[id].plan);
	  ,10)

extern void zfft_fftw(complex_double * inout, int n,
		      int dir, int howmany, int normalize)
{
	int i;
	complex_double *ptr = inout;
	fftw_plan plan = NULL;
	plan = caches_zfftw[get_cache_id_zfftw(n, dir)].plan;

	switch (dir) {
	case 1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			fftw_one(plan, (fftw_complex *) ptr, NULL);
		}
		break;
	case -1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			fftw_one(plan, (fftw_complex *) ptr, NULL);
		}
		break;
	default:
		fprintf(stderr, "zfft: invalid dir=%d\n", dir);
	}

	if (normalize) {
		ptr = inout;
		for (i = n * howmany - 1; i >= 0; --i) {
			*((double *) (ptr)) /= n;
			*((double *) (ptr++) + 1) /= n;
		}
	}
}
