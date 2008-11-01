extern void F_FUNC(zfftf,ZFFTF)(int*,double*,double*);
extern void F_FUNC(zfftb,ZFFTB)(int*,double*,double*);
extern void F_FUNC(zffti,ZFFTI)(int*,double*);
GEN_CACHE(zfftpack,(int n)
	  ,double* wsave;
	  ,(caches_zfftpack[i].n==n)
	  ,caches_zfftpack[id].wsave = (double*)malloc(sizeof(double)*(4*n+15));
	   F_FUNC(zffti,ZFFTI)(&n,caches_zfftpack[id].wsave);
	  ,free(caches_zfftpack[id].wsave);
	  ,10)

static void zfft_fftpack(complex_double * inout,
			 int n, int direction, int howmany, int normalize)
{
	int i;
	complex_double *ptr = inout;
	double *wsave = NULL;

	wsave = caches_zfftpack[get_cache_id_zfftpack(n)].wsave;

	switch (direction) {
	case 1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			zfftf_(&n, (double *) (ptr), wsave);

		}
		break;

	case -1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			zfftb_(&n, (double *) (ptr), wsave);
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
