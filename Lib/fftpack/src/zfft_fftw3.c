GEN_CACHE(zfftw3,(int n,int d)
	,int direction;
	fftw_plan plan;
	fftw_complex* ptr;
	,((caches_zfftw3[i].n==n) &&
	    (caches_zfftw3[i].direction==d))
	,caches_zfftw3[id].direction = d;
	caches_zfftw3[id].ptr = fftw_malloc(sizeof(fftw_complex)*(n));
	    caches_zfftw3[id].plan = fftw_plan_dft_1d(n, caches_zfftw3[id].ptr,
	caches_zfftw3[id].ptr,
		(d>0?FFTW_FORWARD:FFTW_BACKWARD),
		FFTW_ESTIMATE);
	,fftw_destroy_plan(caches_zfftw3[id].plan);
	fftw_free(caches_zfftw3[id].ptr);
	,10)

static void zfft_fftw3(complex_double * inout, int n, int dir, int
howmany, int normalize)
{
	complex_double *ptr = inout;
	fftw_complex *ptrm = NULL;
	fftw_plan plan = NULL;

	int i;

	plan = caches_zfftw3[get_cache_id_zfftw3(n, dir)].plan;

	switch (dir) {
	case 1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			ptrm =
			    caches_zfftw3[get_cache_id_zfftw3(n, dir)].ptr;
			memcpy(ptrm, ptr, sizeof(double) * 2 * n);
			fftw_execute(plan);
			memcpy(ptr, ptrm, sizeof(double) * 2 * n);
		}
		break;

	case -1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			ptrm =
			    caches_zfftw3[get_cache_id_zfftw3(n, dir)].ptr;
			memcpy(ptrm, ptr, sizeof(double) * 2 * n);
			fftw_execute(plan);
			memcpy(ptr, ptrm, sizeof(double) * 2 * n);
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
