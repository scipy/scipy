/*
 * fftpack backend for multi dimensional fft
 *
 * Original code by Pearu Peaterson
 *
 * Last Change: Sun May 11 09:00 PM 2008 J
 */
#include <new>

#include "cycliccache.h"

extern "C" {
extern void zfft(complex_double * inout,
		 int n, int direction, int howmany, int normalize);
};

static int next_comb(int *ia, int *da, int m);
static void flatten(complex_double * dest, complex_double * src,
                int rank, int strides_axis, int dims_axis, int unflat,
                int *tmp);


using namespace fft;

class NDFFTPackCacheId : public CacheId {
        public:
                NDFFTPackCacheId(int n, int rank) : 
                        CacheId(n),
                        m_rank(rank)  
                {
                };

                virtual bool operator==(const NDFFTPackCacheId &other) const 
                {
                        return is_equal(other);
                };

                virtual bool is_equal(const NDFFTPackCacheId &other) const;

        public:
                int m_rank;
};

bool NDFFTPackCacheId::is_equal(const NDFFTPackCacheId & other) const
{
        return m_n == other.m_n && m_rank == other.m_rank;
}

class NDFFTPackCache: public Cache<NDFFTPackCacheId> {
        public:
                NDFFTPackCache(const NDFFTPackCacheId& id);
                virtual ~NDFFTPackCache();

                int compute(complex_double * inout, int sz, int* dims, 
                            int direction, int normalize, int howmany) const;

        protected:
                complex_double* m_wsave;
                int* m_iptr;

        private:
                int prepare(int *dims) const;
};

NDFFTPackCache::NDFFTPackCache(const NDFFTPackCacheId& id)
:	Cache<NDFFTPackCacheId>(id)
{
        int n = id.m_n;
        int rank = id.m_rank;

        m_wsave = (complex_double *)malloc(sizeof(*m_wsave) * (2 * n));
        if (m_wsave == NULL) {
                goto fail;
        }

        m_iptr = (int*)malloc(4 * rank * sizeof(*m_iptr));
        if (m_iptr == NULL) {
                goto clean_wsave;
        }

        return;

clean_wsave:
        free(m_wsave);
fail:
	throw std::bad_alloc();
}

NDFFTPackCache::~NDFFTPackCache()
{
        free(m_iptr);
        free(m_wsave);
}

int NDFFTPackCache::compute(complex_double *inout, int sz, int *dims, 
                int direction, int normalize, 
                int howmany) const 
{
        int rank = m_id.m_rank;
        int i, axis, k, j;
        complex_double *tmp = m_wsave;
        complex_double *ptr = inout;

        zfft(inout, dims[rank - 1], direction, howmany * sz / dims[rank - 1],
             normalize);
        prepare(dims);

        for (i = 0; i < howmany; ++i, ptr += sz) {
                for (axis = 0; axis < rank - 1; ++axis) {
                        for (k = j = 0; k < rank; ++k) {
                                if (k != axis) {
                                        *(m_iptr + rank + j) = m_iptr[k];
                                        *(m_iptr + 2 * rank + j++) = dims[k] - 1;
                                }
                        }
                        flatten(tmp, ptr, rank, m_iptr[axis], dims[axis], 0, m_iptr);
                        zfft(tmp, dims[axis], direction, sz / dims[axis], normalize);
                        flatten(ptr, tmp, rank, m_iptr[axis], dims[axis], 1, m_iptr);
                }
        }
        return 0;
}

int NDFFTPackCache::prepare(int *dims) const 
{
        int rank = m_id.m_rank;
        int i;

        m_iptr[rank - 1] = 1;
        for (i = 2; i <= rank; ++i) {
                m_iptr[rank - i] = m_iptr[rank - i + 1] * dims[rank - i + 1];
        }

        return 0;
}

static CacheManager<NDFFTPackCacheId, NDFFTPackCache> ndfftpack_cmgr(10);

#if 0
/* stub to make PUBLIC_GEN_API happy */
static void destroy_zfftnd_fftpack_caches()
{
}
#endif

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

        for (i = 0; i < rm2; ++i) {
                ia[i] = 0;
        }

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

extern void zfftnd_fftpack(complex_double * inout, int rank,
			   int *dims, int direction, int howmany,
			   int normalize)
{
        int i, sz;
        complex_double *ptr = inout;
        NDFFTPackCache* cache;

        sz = 1;
        for (i = 0; i < rank; ++i) {
                sz *= dims[i];
        }

        cache = ndfftpack_cmgr.get_cache(NDFFTPackCacheId(sz, rank));
        cache->compute(ptr, sz, dims, direction, normalize, howmany);

}
