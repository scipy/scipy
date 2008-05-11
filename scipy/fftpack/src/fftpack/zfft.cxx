#include <new>

#include "cycliccache.h"

extern "C" {
        extern void F_FUNC(zfftf,ZFFTF)(int*,double*,double*);
        extern void F_FUNC(zfftb,ZFFTB)(int*,double*,double*);
        extern void F_FUNC(zffti,ZFFTI)(int*,double*);
};

using namespace fft;

class FFTPackCacheId : public CacheId {
        public:
                FFTPackCacheId(int n) : CacheId(n) {};
};

class FFTPackCache: public Cache<FFTPackCacheId> {
        public:
                FFTPackCache(const FFTPackCacheId& id);
                virtual ~FFTPackCache();

                int compute_forward(complex_double * inout) const;
                int compute_backward(complex_double * inout) const;

        protected:
                double* m_wsave;
};

FFTPackCache::FFTPackCache(const FFTPackCacheId& id)
:	Cache<FFTPackCacheId>(id)
{
        int n = id.m_n;

        m_wsave = (double *)malloc(sizeof(*m_wsave) * (4 * n + 15));
        if (m_wsave == NULL) {
                goto fail;
        }

        F_FUNC(zffti,ZFFTI)(&n, m_wsave);

        return;

fail:
	throw std::bad_alloc();
}

FFTPackCache::~FFTPackCache()
{
        free(m_wsave);
}

int FFTPackCache::compute_forward(complex_double *inout) const 
{
        int n = m_id.m_n;

        zfftf_(&n, (double *)(inout), m_wsave);
        return 0;
}

int FFTPackCache::compute_backward(complex_double *inout) const 
{
        int n = m_id.m_n;

        zfftb_(&n, (double *)(inout), m_wsave);
        return 0;
}

static CacheManager<FFTPackCacheId, FFTPackCache> fftpack_cmgr(10);

static void zfft_fftpack(complex_double * inout,
			 int n, int direction, int howmany, int normalize)
{
	int i;
	complex_double *ptr = inout;
        FFTPackCache *cache;

	cache = fftpack_cmgr.get_cache(FFTPackCacheId(n));

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
