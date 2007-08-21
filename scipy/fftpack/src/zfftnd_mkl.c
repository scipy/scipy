/*
 * MKL backend for multi dimensional fft
 *
 * Original code by David M. Cooke
 *
 * Last Change: Wed Aug 08 03:00 PM 2007 J
 */

static long *convert_dims(int n, int *dims)
{
    long *ndim;
    int i;
    ndim = (long *) malloc(sizeof(long) * n);
    for (i = 0; i < n; i++) {
        ndim[i] = (long) dims[i];
    }
    return ndim;
}

GEN_CACHE(zfftnd_mkl, (int n, int *dims)
	  , DFTI_DESCRIPTOR_HANDLE desc_handle;
	  int *dims;
	  long *ndims;, ((caches_zfftnd_mkl[i].n == n) &&
			 (equal_dims(n, caches_zfftnd_mkl[i].dims, dims)))
	  , caches_zfftnd_mkl[id].ndims = convert_dims(n, dims);
	  caches_zfftnd_mkl[id].n = n;
	  caches_zfftnd_mkl[id].dims = (int *) malloc(sizeof(int) * n);
	  memcpy(caches_zfftnd_mkl[id].dims, dims, sizeof(int) * n);
	  DftiCreateDescriptor(&caches_zfftnd_mkl[id].desc_handle,
			       DFTI_DOUBLE, DFTI_COMPLEX, (long) n,
			       caches_zfftnd_mkl[id].ndims);
	  DftiCommitDescriptor(caches_zfftnd_mkl[id].desc_handle);,
	  DftiFreeDescriptor(&caches_zfftnd_mkl[id].desc_handle);
	  free(caches_zfftnd_mkl[id].dims);
	  free(caches_zfftnd_mkl[id].ndims);, 10)

extern void zfftnd_mkl(complex_double * inout, int rank,
		       int *dims, int direction, int howmany,
		       int normalize)
{
    int i, sz;
    complex_double *ptr = inout;

    DFTI_DESCRIPTOR_HANDLE desc_handle;
    sz = 1;
    for (i = 0; i < rank; ++i) {
        sz *= dims[i];
    }

    desc_handle =
	caches_zfftnd_mkl[get_cache_id_zfftnd_mkl(rank, dims)].desc_handle;
    for (i = 0; i < howmany; ++i, ptr += sz) {
        if (direction == 1) {
            DftiComputeForward(desc_handle, (double *) ptr);
        } else if (direction == -1) {
            DftiComputeBackward(desc_handle, (double *) ptr);
        }
    }
    if (normalize) {
        ptr = inout;
        for (i = sz * howmany - 1; i >= 0; --i) {
            *((double *) (ptr)) /= sz;
            *((double *) (ptr++) + 1) /= sz;
        }
    }
}
