/*
 * Last Change: Sun May 11 09:00 PM 2008 J
 *
 * FFTW2 implementation
 *
 * Original code by Pearu Peterson.
 */
#include <new>

#include "cycliccache.h"

using namespace fft;

class RFFTWCacheId : public CacheId {
	public:
		RFFTWCacheId(int n, int dir, int flags);

		virtual bool operator==(const RFFTWCacheId& other) const
		{
			return is_equal(other);
		}

		virtual bool is_equal(const RFFTWCacheId& other) const
		{
			const CacheId *ot = &other;
			const CacheId *th = this;

                        return m_dir == other.m_dir &&  
                               m_flags == other.m_flags && 
                               th->is_equal(*ot);
		}

	public:
		int m_dir;
		int m_flags;
};

RFFTWCacheId::RFFTWCacheId(int n, int dir, int flags):
        CacheId(n),
        m_dir(dir),
        m_flags(flags)
{
}

class RFFTWCache : public Cache<RFFTWCacheId> {
	public:	
		RFFTWCache(const RFFTWCacheId& id);
		virtual ~RFFTWCache();

		int compute(double* inout) const;

	protected:
		rfftw_plan m_plan;
                double *m_wrk;
};

RFFTWCache::RFFTWCache(const RFFTWCacheId& id)
:	Cache<RFFTWCacheId>(id)
{
        m_wrk = (double *) malloc(sizeof(double) * id.m_n);
        if(m_wrk == NULL) {
                goto fail;
        }

	m_plan = rfftw_create_plan(id.m_n, 
				  (id.m_dir > 0 ?  FFTW_REAL_TO_COMPLEX :
                                   FFTW_COMPLEX_TO_REAL), 
				  id.m_flags);

	if (m_plan == NULL) {
		goto clean_wrk;
	}

        return;

clean_wrk:
        free(m_wrk);
fail:
	throw std::bad_alloc();
}

RFFTWCache::~RFFTWCache()
{
        free(m_wrk);
	rfftw_destroy_plan(m_plan);
}

int RFFTWCache::compute(double* inout) const
{
        if(m_id.m_dir == 1) {
                memcpy(m_wrk, inout, sizeof(double) * m_id.m_n);
                rfftw(m_plan, 1, (fftw_real *) m_wrk, 1, 1, NULL, 1, 1);
                COPYRFFTW2STD(m_wrk, inout, m_id.m_n);
        } else {
                COPYINVRFFTW2STD(inout, m_wrk, m_id.m_n);
                rfftw(m_plan, 1, (fftw_real *) m_wrk, 1, 1, NULL, 1, 1);
                memcpy(inout, m_wrk, sizeof(double) * m_id.m_n);
        }
        return 0;
};

CacheManager<RFFTWCacheId, RFFTWCache> rfftw_cmgr(10);

static void drfft_fftw(double *inout, int n, int dir, 
                       int howmany, int normalize)
{
        int i;
        double *ptr = inout;
        RFFTWCache *cache;

        cache = rfftw_cmgr.get_cache(
                        RFFTWCacheId(n, dir, FFTW_IN_PLACE | FFTW_ESTIMATE));

        if (dir != -1 && dir != 1) {
		fprintf(stderr, "drfft: invalid dir=%d\n", dir);
                return;
        } else {
                for (i = 0; i < howmany; ++i, ptr += n) {
                        cache->compute(ptr);
                }
        }

        if (normalize) {
                double d = 1.0 / n;
                ptr = inout;
                for (i = n * howmany - 1; i >= 0; --i) {
                        (*(ptr++)) *= d;
                }
        }
}
