#ifndef _SCIPY_FFTPACK_FFTPACK_COMMON_H
#define _SCIPY_FFTPACK_FFTPACK_COMMON_H

#include <new>

#include "cycliccache.h"

namespace fft {

extern "C" {
extern void F_FUNC(dfftf, DFFTF) (int *, double *, double *);
extern void F_FUNC(dfftb, DFFTB) (int *, double *, double *);
extern void F_FUNC(dffti, DFFTI) (int *, double *);
};

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

};

#endif
