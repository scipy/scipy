#include <new>
#include <cassert>

#include "common.h"

using namespace fft;

class FFTW3Cache : public Cache<FFTW3CacheId> {
	public:	
		FFTW3Cache(const FFTW3CacheId& id);
		virtual ~FFTW3Cache();

		int compute(fftw_complex* inout) const
		{
                        assert (m_id.m_isalign ? is_simd_aligned(inout) :
                                                 true);
			fftw_execute_dft(m_plan, inout, inout);
			return 0;
		};

	protected:
		fftw_plan m_plan;	
		fftw_complex *m_wrk;	
};

FFTW3Cache::FFTW3Cache(const FFTW3CacheId& id)
:	Cache<FFTW3CacheId>(id)
{
        int flags = FFTW_ESTIMATE;

	m_wrk = (fftw_complex*)fftw_malloc(id.m_n * sizeof(double) * 2);
	if (m_wrk == NULL) {
		goto fail_wrk;
	}

        if (!m_id.m_isalign) {
                flags |= FFTW_UNALIGNED;
        } 

	m_plan = fftw_plan_dft_1d(id.m_n, m_wrk, m_wrk, 
				  (id.m_dir > 0 ?  FFTW_FORWARD :
                                                   FFTW_BACKWARD),
				  flags);

	if (m_plan == NULL) {
		goto clean_wrk;
	}

	return ;

clean_wrk:
	fftw_free(m_wrk);
fail_wrk:
	throw std::bad_alloc();
}

FFTW3Cache::~FFTW3Cache()
{
	fftw_destroy_plan(m_plan);
	fftw_free(m_wrk);
}

static CacheManager<FFTW3CacheId, FFTW3Cache> fftw3_cmgr(10);

static void zfft_fftw3(complex_double * inout, int n, int dir, int howmany, 
                       int normalize)
{
	fftw_complex *ptr = (fftw_complex*)inout;
        double factor = 1./n;
        FFTW3Cache *cache;
        bool isaligned;

	int i;

        isaligned = is_simd_aligned(ptr); 

        if (howmany > 1) {
                /* 
                 * If executing for several consecutive buffers, we have to
                 * check that the shifting one buffer does not make it
                 * unaligned 
                 */
                isaligned = isaligned && is_simd_aligned(ptr + n);
        }
	cache = fftw3_cmgr.get_cache(FFTW3CacheId(n, dir, isaligned));

	switch (dir) {
	case 1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			cache->compute(ptr);
		}
		break;

	case -1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			cache->compute(ptr);
		}
		break;

	default:
		fprintf(stderr, "zfft: invalid dir=%d\n", dir);
	}

	if (normalize) {
                ptr =(fftw_complex*)inout;
		for (i = n * howmany - 1; i >= 0; --i) {
			*((double *) (ptr)) *= factor;
			*((double *) (ptr++) + 1) *= factor;
		}
	}
}
