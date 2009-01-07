/*
 * fftpack backend for multi dimensional fft
 *
 * Original code by Pearu Peaterson
 *
 * Last Change: Wed Aug 08 02:00 PM 2007 J
 */

GEN_CACHE(cfftnd_fftpack, (int n, int rank)
	  , complex_float * ptr; int *iptr; int rank;
	  , ((caches_cfftnd_fftpack[i].n == n)
	     && (caches_cfftnd_fftpack[i].rank == rank))
	  , caches_cfftnd_fftpack[id].n = n;
	  caches_cfftnd_fftpack[id].ptr =
	  (complex_float *) malloc(2 * sizeof(float) * n);
	  caches_cfftnd_fftpack[id].iptr =
	  (int *) malloc(4 * rank * sizeof(int));
	  ,
	  free(caches_cfftnd_fftpack[id].ptr);
	  free(caches_cfftnd_fftpack[id].iptr);
	  , 10)

static
void sflatten(complex_float * dest, complex_float * src,
	     int rank, int strides_axis, int dims_axis, int unflat,
	     int *tmp)
{
    int *new_strides = tmp + rank;
    int *new_dims = tmp + 2 * rank;
    int *ia = tmp + 3 * rank;
    int rm1 = rank - 1, rm2 = rank - 2;
    int i, j, k;
    for (i = 0; i < rm2; ++i)
	ia[i] = 0;
    ia[rm2] = -1;
    j = 0;
    if (unflat) {
        while (next_comb(ia, new_dims, rm2)) {
            k = 0;
            for (i = 0; i < rm1; ++i) {
                k += ia[i] * new_strides[i];
            }
            for (i = 0; i < dims_axis; ++i) {
                *(dest + k + i * strides_axis) = *(src + j++);
            }
        }
    } else {
        while (next_comb(ia, new_dims, rm2)) {
            k = 0;
            for (i = 0; i < rm1; ++i) {
                k += ia[i] * new_strides[i];
            }
            for (i = 0; i < dims_axis; ++i) {
                *(dest + j++) = *(src + k + i * strides_axis);
            }
        }
    }
}

extern void cfft(complex_float * inout,
		 int n, int direction, int howmany, int normalize);

extern void cfftnd_fftpack(complex_float * inout, int rank,
			   int *dims, int direction, int howmany,
			   int normalize)
{
    int i, sz;
    complex_float *ptr = inout;
    int axis;
    complex_float *tmp;
    int *itmp;
    int k, j;

    sz = 1;
    for (i = 0; i < rank; ++i) {
        sz *= dims[i];
    }
    cfft(ptr, dims[rank - 1], direction, howmany * sz / dims[rank - 1],
	 normalize);

    i = get_cache_id_cfftnd_fftpack(sz, rank);
    tmp = caches_cfftnd_fftpack[i].ptr;
    itmp = caches_cfftnd_fftpack[i].iptr;

    itmp[rank - 1] = 1;
    for (i = 2; i <= rank; ++i) {
        itmp[rank - i] = itmp[rank - i + 1] * dims[rank - i + 1];
    }

    for (i = 0; i < howmany; ++i, ptr += sz) {
        for (axis = 0; axis < rank - 1; ++axis) {
            for (k = j = 0; k < rank; ++k) {
                if (k != axis) {
                    *(itmp + rank + j) = itmp[k];
                    *(itmp + 2 * rank + j++) = dims[k] - 1;
                }
            }
            sflatten(tmp, ptr, rank, itmp[axis], dims[axis], 0, itmp);
            cfft(tmp, dims[axis], direction, sz / dims[axis], normalize);
            sflatten(ptr, tmp, rank, itmp[axis], dims[axis], 1, itmp);
        }
    }

}
