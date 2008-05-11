#include <new>

#include "cycliccache.h"

using namespace fft;

class MKLCacheId : public CacheId {
        public:
                MKLCacheId(int n) : CacheId(n) {};
};


class MKLCache: public Cache<MKLCacheId> {
        public:
                MKLCache(const MKLCacheId& id);
                virtual ~MKLCache();

                int compute_forward(complex_double * inout) const;
                int compute_backward(complex_double * inout) const;

        protected:
                DFTI_DESCRIPTOR_HANDLE m_hdl;
};

MKLCache::MKLCache(const MKLCacheId& id)
:	Cache<MKLCacheId>(id)
{
        int n = id.m_n;

        DftiCreateDescriptor(&m_hdl, DFTI_DOUBLE, DFTI_COMPLEX, 1, (long)n); 
        DftiCommitDescriptor(m_hdl);

        return;
}

MKLCache::~MKLCache()
{
        DftiFreeDescriptor(&m_hdl);
}

int MKLCache::compute_forward(complex_double *inout) const 
{
        DftiComputeForward(m_hdl, (double *) inout);
        return 0;
}

int MKLCache::compute_backward(complex_double *inout) const 
{
        DftiComputeBackward(m_hdl, (double *) inout);
        return 0;
}

CacheManager<MKLCacheId, MKLCache> mkl_cmgr(10);

static void zfft_mkl(complex_double * inout,
		 int n, int direction, int howmany, int normalize)
{
	int i;
	complex_double *ptr = inout;
        MKLCache *cache;

        cache = mkl_cmgr.get_cache(MKLCacheId(n));

	switch (direction) {

	case 1:
		for (i = 0; i < howmany; ++i, ptr += n) {
                        cache->compute_forward(ptr);
		}
		break;

	case -1:
		for (i = 0; i < howmany; ++i, ptr += n) {
                        cache->compute_backward(ptr);
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
