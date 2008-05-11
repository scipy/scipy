#include <new>

#include "cycliccache.h"

using namespace fft;

class FFTWCacheId : public CacheId {
	public:
		FFTWCacheId(int n, int dir) : CacheId(n), m_dir(dir) {};

		virtual bool operator==(const FFTWCacheId& other) const
		{
			return is_equal(other);
		}

		virtual bool is_equal(const FFTWCacheId& other) const
		{
			const CacheId *ot = &other;
			const CacheId *th = this;

			return m_dir == other.m_dir &&  th->is_equal(*ot);
		}

	public:
		int m_dir;
};

class FFTWCache : public Cache<FFTWCacheId> {
	public:	
		FFTWCache(const FFTWCacheId& id);
		virtual ~FFTWCache();

		int compute(fftw_complex* inout) 
		{
			fftw_one(m_plan, inout, NULL);
			return 0;
		};

	protected:
		fftw_plan m_plan;	
};

FFTWCache::FFTWCache(const FFTWCacheId& id)
:	Cache<FFTWCacheId>(id)
{
	m_plan = fftw_create_plan(id.m_n, 
				  (id.m_dir > 0 ?  FFTW_FORWARD:FFTW_BACKWARD), 
				  FFTW_ESTIMATE | FFTW_IN_PLACE);

	if (m_plan == NULL) {
		goto fail;
	}

	return ;

fail:
	throw std::bad_alloc();
}

FFTWCache::~FFTWCache()
{
	fftw_destroy_plan(m_plan);
}

CacheManager<FFTWCacheId, FFTWCache> fftw_cmgr(10);

/* stub to make GEN_PUBLIC_API happy */
static void destroy_zfftw_caches()
{
}

extern void zfft_fftw(complex_double * inout, int n,
		      int dir, int howmany, int normalize)
{
	int i;
	complex_double *ptr = inout;
        FFTWCache* cache;

	cache = fftw_cmgr.get_cache(FFTWCacheId(n, dir));

        if (dir != -1 && dir != 1) {
		fprintf(stderr, "zfft: invalid dir=%d\n", dir);
        } else {
		for (i = 0; i < howmany; ++i, ptr += n) {
			cache->compute((fftw_complex *) ptr);
		}
        }

	if (normalize) {
		ptr = inout;
		for (i = n * howmany - 1; i >= 0; --i) {
			*((double *) (ptr)) /= n;
			*((double *) (ptr++) + 1) /= n;
		}
	}
}
