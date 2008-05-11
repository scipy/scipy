/*
 * Last Change: Sun May 11 09:00 PM 2008 J
 *
 * FFTPACK implementation
 *
 * Original code by Pearu Peterson.
 */
#include <new>

#include "cycliccache.h"

extern "C" {
extern void F_FUNC(dfftf, DFFTF) (int *, double *, double *);
extern void F_FUNC(dfftb, DFFTB) (int *, double *, double *);
extern void F_FUNC(dffti, DFFTI) (int *, double *);
};

using namespace fft;

class RFFTPackCacheId : public CacheId {
        public:
                RFFTPackCacheId(int n) : CacheId(n) {};
};

class RFFTPackCache: public Cache<RFFTPackCacheId> {
        public:
                RFFTPackCache(const RFFTPackCacheId& id);
                virtual ~RFFTPackCache();

                int compute_forward(double * inout) const;
                int compute_backward(double * inout) const;

        protected:
                double* m_wsave;
};

RFFTPackCache::RFFTPackCache(const RFFTPackCacheId& id)
:	Cache<RFFTPackCacheId>(id)
{
        int n = id.m_n;

        m_wsave = (double *)malloc(sizeof(*m_wsave) * (2 * n + 15));
        if (m_wsave == NULL) {
                goto fail;
        }

        F_FUNC(dffti, DFFTI)(&n, m_wsave);

        return;

fail:
	throw std::bad_alloc();
}

RFFTPackCache::~RFFTPackCache()
{
        free(m_wsave);
}

int RFFTPackCache::compute_forward(double *inout) const 
{
        int n = m_id.m_n;

        F_FUNC(dfftf, DFFTF)(&n, inout, m_wsave);
        return 0;
}

int RFFTPackCache::compute_backward(double *inout) const 
{
        int n = m_id.m_n;

        F_FUNC(dfftb, DFFTB)(&n, inout, m_wsave);
        return 0;
}

static CacheManager<RFFTPackCacheId, RFFTPackCache> rfftpack_cmgr(10);

static void drfft_fftpack(double *inout, int n, int direction, int howmany,
			  int normalize)
{
        int i;
        double *ptr = inout;
        RFFTPackCache* cache;

        cache = rfftpack_cmgr.get_cache(RFFTPackCacheId(n));

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
                        fprintf(stderr, "drfft: invalid direction=%d\n", direction);
                        return;
        }

        if (normalize) {
                double d = 1.0 / n;
                ptr = inout;
                for (i = n * howmany - 1; i >= 0; --i) {
                        (*(ptr++)) *= d;
                }
        }
}
