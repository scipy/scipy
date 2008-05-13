/*
 * Last Change: Tue May 13 12:00 PM 2008 J
 *
 * Original code by Pearu Peterson.
 */

/*
 * DJBFFT only implements size 2^N !
 *
 * drfft_def and drfft_def_destroy_cache are the functions used for size different
 * than 2^N
 */
#include <new>
#include <cassert>

#include "common.h"

#ifdef WITH_FFTW3
#define drfft_def drfft_fftw3
#elif defined WITH_FFTW
#define drfft_def drfft_fftw
#else
#define drfft_def drfft_fftpack
#endif

using namespace fft;

class RDJBFFTCache: public Cache<DJBFFTCacheId> {
        public:
                RDJBFFTCache(const DJBFFTCacheId& id);
                virtual ~RDJBFFTCache();

                int compute_forward(double * inout) const;
                int compute_backward(double * inout, int normalize) const;

        protected:
                unsigned int* m_f;
                double* m_ptr;
};

RDJBFFTCache::RDJBFFTCache(const DJBFFTCacheId& id)
:	Cache<DJBFFTCacheId>(id)
{
        int n = id.m_n;

        m_f = (unsigned int*)malloc(sizeof(*m_f) * n);
        if (m_f == NULL) {
                goto fail_f;
        }

        m_ptr = (double *)malloc(sizeof(*m_ptr) * 2 * n);
        if (m_ptr == NULL) {
                goto clean_f;
        }

        fftfreq_rtable(m_f, id.m_n);
        return;

clean_f:
        free(m_f);
fail_f:
	throw std::bad_alloc();
}

RDJBFFTCache::~RDJBFFTCache()
{
        free(m_ptr);
        free(m_f);
}

int RDJBFFTCache::compute_forward(double *inout) const
{
	double *ptr = inout;
	double *ptrc = m_ptr;

	int n =  m_id.m_n;
	COPYSTD2DJB(ptr, ptrc, n);
	switch (n) {
#define TMPCASE(N) case N: fftr8_##N(ptrc); break
		TMPCASE(2);
		TMPCASE(4);
		TMPCASE(8);
		TMPCASE(16);
		TMPCASE(32);
		TMPCASE(64);
		TMPCASE(128);
		TMPCASE(256);
		TMPCASE(512);
		TMPCASE(1024);
		TMPCASE(2048);
		TMPCASE(4096);
		TMPCASE(8192);
#undef TMPCASE
	}
	COPYDJB2STD(ptrc, ptr, m_f, n);
	return 0;
}

int RDJBFFTCache::compute_backward(double *inout, int normalize) const
{
	double *ptr = inout;
	double *ptrc = m_ptr;

	int n =  m_id.m_n;

	COPYINVSTD2DJB(ptr, ptrc, normalize, m_f, n);
	switch (n) {

#define TMPCASE(N)case N:if(normalize)fftr8_scale##N(ptrc);fftr8_un##N(ptrc);break
		TMPCASE(2);
		TMPCASE(4);
		TMPCASE(8);
		TMPCASE(16);
		TMPCASE(32);
		TMPCASE(64);
		TMPCASE(128);
		TMPCASE(256);
		TMPCASE(512);
		TMPCASE(1024);
		TMPCASE(2048);
		TMPCASE(4096);
		TMPCASE(8192);
#undef TMPCASE
	}
	COPYINVDJB2STD(ptrc, ptr, n);
	return 0;
}

static CacheManager<DJBFFTCacheId, RDJBFFTCache> rdjbfft_cmgr(10);

/**************** ZFFT function **********************/
void drfft_djbfft(double * inout,
			 int n, int direction, int howmany, int normalize)
{
	int i;
	double *ptr = inout;
	RDJBFFTCache *cache;

	switch (n) {
		case 2:;
		case 4:;
		case 8:;
		case 16:;
		case 32:;
		case 64:;
		case 128:;
		case 256:;
		case 512:;
		case 1024:;
		case 2048:;
		case 4096:;
		case 8192:
			  cache = rdjbfft_cmgr.get_cache(DJBFFTCacheId(n));
			  break;
		default:
			  drfft_def(inout, n, direction, howmany, normalize);
			  return;
	}

	switch (direction) {
		case 1:
			for (i = 0; i < howmany; ++i, ptr += n) {
				cache->compute_forward(ptr);
			}
			break;

		case -1:
			for (i = 0; i < howmany; ++i, ptr += n) {
				cache->compute_backward(ptr, normalize);
			}
			break;

		default:
			fprintf(stderr, "drfft: invalid direction=%d\n", direction);
	}

	if (normalize && direction == 1) {
		double d = 1.0 / n;
		ptr = inout;
		for (i = n * howmany - 1; i >= 0; --i) {
			(*(ptr++)) *= d;
		}
	}
}
