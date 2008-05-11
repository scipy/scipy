/*
 * fftw2 backend for multi dimensional fft
 *
 * Original code by Pearu Peaterson
 *
 * Last Change: Sun May 11 09:00 PM 2008 J
 */
#include <new>
#include <cassert>

#include <cycliccache.h>

using namespace fft;

class NDFFTWCacheId {
        public:
                NDFFTWCacheId(int rank, int *dims, int dir, int flags);
                virtual ~NDFFTWCacheId();

                NDFFTWCacheId(const NDFFTWCacheId &);

                virtual bool operator==(const NDFFTWCacheId & other) const {
                        return is_equal(other);
                };

                virtual bool is_equal(const NDFFTWCacheId & other) const;

        public:
                int m_dir;
                int m_rank;
                int *m_dims;
                int m_flags;

        private:
                int init(int rank, int *dims);
};

int NDFFTWCacheId::init(int rank, int *dims)
{
	m_dims = (int *) malloc(sizeof(int) * rank);
	if (m_dims == NULL) {
		return -1;
	}
	memcpy(m_dims, dims, rank * sizeof(*m_dims));

	return 0;

}

NDFFTWCacheId::NDFFTWCacheId(int rank, int *dims, int dir, int flags) :
        m_dir(dir),
        m_rank(rank),
        m_flags(flags)
{
        if (init(rank, dims)) {
                goto fail;
        }

fail:
        std::bad_alloc();
}

NDFFTWCacheId::NDFFTWCacheId(const NDFFTWCacheId & copy) :
        m_dir(copy.m_dir),
        m_rank(copy.m_rank),
        m_flags(copy.m_flags)
{
	if (init(copy.m_rank, copy.m_dims)) {
		goto fail;
	}

fail:
	std::bad_alloc();
}

NDFFTWCacheId::~NDFFTWCacheId()
{
	free(m_dims);
}

bool NDFFTWCacheId::is_equal(const NDFFTWCacheId & other) const
{
	bool res;

	res = (other.m_dir == m_dir);
	res = res && (other.m_flags == m_flags);

	if (m_rank == other.m_rank) {
                res = res && equal_dims(m_rank, m_dims, other.m_dims);
	} else {
		return false;
	}

	return res;
}

class NDFFTWCache : public Cache < NDFFTWCacheId > {
        public:
                NDFFTWCache(const NDFFTWCacheId & id);
                virtual ~ NDFFTWCache();

                int compute(fftw_complex * inout) const
                {
                        fftwnd_one(m_plan, inout, NULL);
                        return 0;
                };

        protected:
                fftwnd_plan m_plan;
};

NDFFTWCache::NDFFTWCache(const NDFFTWCacheId & id)
:  Cache < NDFFTWCacheId > (id)
{
	int flags = FFTW_ESTIMATE | FFTW_IN_PLACE;
#if 0
	int sz;
	int i;

	sz = 1;
	for (i = 0; i < m_id.m_rank; ++i) {
		sz *= m_id.m_dims[i];
	}
#endif

	m_plan = fftwnd_create_plan(m_id.m_rank,
				    m_id.m_dims,
				    (id.m_dir > 0 ?
                                     FFTW_FORWARD : FFTW_BACKWARD),
				    flags);

	if (m_plan == NULL) {
		goto fail;
	}

	return;

fail:
        throw std::bad_alloc();
}

NDFFTWCache::~NDFFTWCache()
{
	fftwnd_destroy_plan(m_plan);
}

static CacheManager < NDFFTWCacheId, NDFFTWCache > fftwnd_cmgr(10);

extern void zfftnd_fftw(complex_double * inout, int rank,
		       int *dims, int direction, int howmany,
		       int normalize)
{
        int i, sz;
        complex_double *ptr = inout;
        NDFFTWCache *cache;

        sz = 1;
        for (i = 0; i < rank; ++i) {
                sz *= dims[i];
        }

        cache = fftwnd_cmgr.get_cache(
                        NDFFTWCacheId(rank, dims, direction,
                                      FFTW_IN_PLACE | FFTW_ESTIMATE));

        for (i = 0; i < howmany; ++i, ptr += sz) {
                cache->compute((fftw_complex*)ptr);
        }

        if (normalize) {
                ptr = inout;
                for (i = sz * howmany - 1; i >= 0; --i) {
                        *((double *) (ptr)) /= sz;
                        *((double *) (ptr++) + 1) /= sz;
                }
        }
}
