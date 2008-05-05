/*
 * fftpack backend for multi dimensional fft
 *
 * Original code by Pearu Peaterson
 *
 * Last Change: Wed Aug 08 02:00 PM 2007 J
 */

GEN_CACHE(zfftnd_fftpack, (int n, int rank)
	  , complex_double * ptr; int *iptr; int rank;
	  , ((caches_zfftnd_fftpack[i].n == n)
	     && (caches_zfftnd_fftpack[i].rank == rank))
	  , caches_zfftnd_fftpack[id].n = n;
	  caches_zfftnd_fftpack[id].ptr =
	  (complex_double *) malloc(2 * sizeof(double) * n);
	  caches_zfftnd_fftpack[id].iptr =
	  (int *) malloc(4 * rank * sizeof(int));
	  ,
	  free(caches_zfftnd_fftpack[id].ptr);
	  free(caches_zfftnd_fftpack[id].iptr);
	  , 10)

static
/*inline : disabled because MSVC6.0 fails to compile it. */
int next_comb(int *ia, int *da, int m)
{
    while (m >= 0 && ia[m] == da[m]) {
        ia[m--] = 0;
    }
    if (m < 0) {
        return 0;
    }
    ia[m]++;
    return 1;
}

static
void flatten(complex_double * dest, complex_double * src,
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

extern void zfft(complex_double * inout,
		 int n, int direction, int howmany, int normalize);

extern void zfftnd_fftpack(complex_double * inout, int rank,
			   int *dims, int direction, int howmany,
			   int normalize)
{
    int i, sz;
    complex_double *ptr = inout;
    int axis;
    complex_double *tmp;
    int *itmp;
    int k, j;

    sz = 1;
    for (i = 0; i < rank; ++i) {
        sz *= dims[i];
    }
    zfft(ptr, dims[rank - 1], direction, howmany * sz / dims[rank - 1],
	 normalize);

    i = get_cache_id_zfftnd_fftpack(sz, rank);
    tmp = caches_zfftnd_fftpack[i].ptr;
    itmp = caches_zfftnd_fftpack[i].iptr;

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
            flatten(tmp, ptr, rank, itmp[axis], dims[axis], 0, itmp);
            zfft(tmp, dims[axis], direction, sz / dims[axis], normalize);
            flatten(ptr, tmp, rank, itmp[axis], dims[axis], 1, itmp);
        }
    }

}
