GEN_CACHE(zmkl,(int n)
	  ,DFTI_DESCRIPTOR_HANDLE desc_handle;
	  ,(caches_zmkl[i].n==n)
      ,DftiCreateDescriptor(&caches_zmkl[id].desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, (long)n); 
       DftiCommitDescriptor(caches_zmkl[id].desc_handle);
	  ,DftiFreeDescriptor(&caches_zmkl[id].desc_handle);
	  ,10)

static void zfft_mkl(complex_double * inout,
		 int n, int direction, int howmany, int normalize)
{
	int i;
	complex_double *ptr = inout;
	DFTI_DESCRIPTOR_HANDLE desc_handle;
	desc_handle = caches_zmkl[get_cache_id_zmkl(n)].desc_handle;

	switch (direction) {

	case 1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			DftiComputeForward(desc_handle, (double *) ptr);
		}
		break;

	case -1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			DftiComputeBackward(desc_handle, (double *) ptr);
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
