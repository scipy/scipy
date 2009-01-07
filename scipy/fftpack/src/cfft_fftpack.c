/* XXX: use the .src mechanism: zfft_cfft are the same */
extern void F_FUNC(cfftf,CFFTF)(int*,float*,float*);
extern void F_FUNC(cfftb,CFFTB)(int*,float*,float*);
extern void F_FUNC(cffti,CFFTI)(int*,float*);
GEN_CACHE(cfftpack,(int n)
	  ,float* wsave;
	  ,(caches_cfftpack[i].n==n)
	  ,caches_cfftpack[id].wsave = (float*)malloc(sizeof(float)*(4*n+15));
	   F_FUNC(cffti,CFFTI)(&n,caches_cfftpack[id].wsave);
	  ,free(caches_cfftpack[id].wsave);
	  ,10)

static void cfft_fftpack(complex_float * inout,
			 int n, int direction, int howmany, int normalize)
{
	int i;
	complex_float *ptr = inout;
	float *wsave = NULL;

	wsave = caches_cfftpack[get_cache_id_cfftpack(n)].wsave;

	switch (direction) {
	case 1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			cfftf_(&n, (float *) (ptr), wsave);

		}
		break;

	case -1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			cfftb_(&n, (float *) (ptr), wsave);
		}
		break;
	default:
		fprintf(stderr, "cfft: invalid direction=%d\n", direction);
	}

	if (normalize) {
		ptr = inout;
		for (i = n * howmany - 1; i >= 0; --i) {
			*((float *) (ptr)) /= n;
			*((float *) (ptr++) + 1) /= n;
		}
	}
}
