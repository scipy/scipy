/*
 * MKL backend for multi dimensional fft
 *
 * Original code by David M. Cooke
 *
 * Last Change: Tue May 13 12:00 PM 2008 J
 */
#include <new>

#include "cycliccache.h"

using namespace fft;

class NDMKLCacheId {
        public:
                NDMKLCacheId(int rank, int *dim);
                virtual ~NDMKLCacheId();

                NDMKLCacheId(const NDMKLCacheId &);

                virtual bool operator==(const NDMKLCacheId & other) const {
                        return is_equal(other);
                };

                virtual bool is_equal(const NDMKLCacheId & other) const;

        public:
                int m_rank;
                int *m_dims;

        private:
                int init(int rank, int *dims);
};

int NDMKLCacheId::init(int rank, int *dims)
{
	m_dims = (int *) malloc(sizeof(int) * rank);
	if (m_dims == NULL) {
		return -1;
	}
	memcpy(m_dims, dims, rank * sizeof(*m_dims));

	return 0;

}

NDMKLCacheId::NDMKLCacheId(int rank, int *dims) :
        m_rank(rank)
{
        if (init(rank, dims)) {
                goto fail;
        }

fail:
        std::bad_alloc();
}

NDMKLCacheId::NDMKLCacheId(const NDMKLCacheId & copy) :
        m_rank(copy.m_rank)
{
	if (init(copy.m_rank, copy.m_dims)) {
		goto fail;
	}

fail:
	std::bad_alloc();
}

NDMKLCacheId::~NDMKLCacheId()
{
	free(m_dims);
}

bool NDMKLCacheId::is_equal(const NDMKLCacheId & other) const
{
	bool res;

	if (m_rank == other.m_rank) {
                res = equal_dims(m_rank, m_dims, other.m_dims);
	} else {
		return false;
	}

	return res;
}

/*
 * Cache class for nd-MKL
 */
class NDMKLCache:public Cache < NDMKLCacheId > {
        public:
                NDMKLCache(const NDMKLCacheId & id);
                virtual ~ NDMKLCache();

                int compute_forward(double * inout) const
                {
                        DftiComputeForward(m_hdl, inout);
                        return 0;
                };

                int compute_backward(double * inout) const
                {
                        DftiComputeBackward(m_hdl, inout);
                        return 0;
                };

        protected:
                int m_rank;
                int *m_dims;
                long *m_ndims;
                DFTI_DESCRIPTOR_HANDLE m_hdl;

        private:
                long *convert_dims(int n, int *dims) const;

};

NDMKLCache::NDMKLCache(const NDMKLCacheId & id)
:  Cache < NDMKLCacheId > (id)
{
        m_rank = id.m_rank;
        m_ndims = convert_dims(id.m_rank, id.m_dims);
        m_dims = (int *) malloc(sizeof(int) * m_rank);
        if (m_dims == NULL) {
                goto fail;
        }

	memcpy(m_dims, id.m_dims, sizeof(int) * m_rank);
        DftiCreateDescriptor(&m_hdl, DFTI_DOUBLE, DFTI_COMPLEX, (long) m_rank,
                             m_ndims);
        DftiCommitDescriptor(m_hdl);

        return;

fail:
        throw std::bad_alloc();
}

NDMKLCache::~NDMKLCache()
{
        DftiFreeDescriptor(&m_hdl);
        free(m_dims);
        free(m_ndims);
}

long* NDMKLCache::convert_dims(int n, int *dims) const
{
        long *ndim;
        int i;

        ndim = (long *) malloc(sizeof(*ndim) * n);
        for (i = 0; i < n; i++) {
                ndim[i] = (long) dims[i];
        }
        return ndim;
}

static CacheManager < NDMKLCacheId, NDMKLCache > ndmkl_cmgr(10);

extern void zfftnd_mkl(complex_double * inout, int rank,
		       int *dims, int direction, int howmany,
		       int normalize)
{
        int i, sz;
        complex_double *ptr = inout;
	NDMKLCache *cache;

        sz = 1;
        for (i = 0; i < rank; ++i) {
                sz *= dims[i];
        }

        cache = ndmkl_cmgr.get_cache(NDMKLCacheId(rank, dims));
        switch(direction) {
                case 1:
                        for (i = 0; i < howmany; ++i, ptr += sz) {
                                cache->compute_forward((double*)ptr);
                        }
                        break;
                case -1:
                        for (i = 0; i < howmany; ++i, ptr += sz) {
                                cache->compute_backward((double*)ptr);
                        }
                        break;
                default:
                        fprintf(stderr,
                                "nd mkl:Wrong direction (this is a bug)\n");
                        return;
        }
        if (normalize) {
                ptr = inout;
                for (i = sz * howmany - 1; i >= 0; --i) {
                        *((double *) (ptr)) /= sz;
                        *((double *) (ptr++) + 1) /= sz;
                }
        }
}
