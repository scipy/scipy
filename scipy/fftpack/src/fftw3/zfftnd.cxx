/*
 * fftw3 backend for multi dimensional fft
 *
 * Original code by Pearu Peaterson
 *
 * Last Change: Wed Aug 08 02:00 PM 2007 J
 */
#include <new>
#include <cassert>

#include "common.h"

using namespace fft;

class NDFFTW3CacheId {
	public:
		NDFFTW3CacheId(int rank, int* dims, int howmany, int dir, bool isalign);
                virtual ~NDFFTW3CacheId();

		virtual bool operator==(const NDFFTW3CacheId& other) const
		{
			return is_equal(other);
		};

		virtual bool is_equal(const NDFFTW3CacheId& other) const;

	public:
		int m_rank;
		int *m_dims;
                int m_howmany;
                int m_dir;
                bool m_isalign;
		
};

NDFFTW3CacheId::NDFFTW3CacheId(int rank, int* dims, int howmany, int dir, bool isalign) :
        m_rank(rank),
        m_howmany(howmany),
        m_dir(dir),
        m_isalign(isalign)
{
        m_dims = (int*)fftw_malloc(m_rank * sizeof(*m_dims));
        if (m_dims == NULL) {
                goto fail;
        }
        memcpy(m_dims, dims, m_rank * sizeof(*m_dims));

        return;

fail:
        std::bad_alloc();
}

NDFFTW3CacheId::~NDFFTW3CacheId()
{
        fftw_free(m_dims);
}

bool NDFFTW3CacheId::is_equal(const NDFFTW3CacheId& other) const
{
        bool res;
        int i;

        res = (other.m_dir == m_dir);
        res = res && (other.m_isalign == m_isalign);
        res = res && (m_howmany == other.m_howmany);
        
        if (m_rank == other.m_rank) {
                fprintf(stderr, "rank is %d, %d\n", m_rank, other.m_rank);
                for (i = 0; i < m_rank; ++i) {
                        res = res && (m_dims[i] == other.m_dims[i]);
                }
        } else {
                return false;
        }

        return res;
}

/* stub because fftw3 has no cache mechanism (yet) */
static void destroy_zfftnd_fftw3_caches(void) {}

class NDFFTW3Cache : public Cache<NDFFTW3CacheId> {
	public:	
		NDFFTW3Cache(const NDFFTW3CacheId& id);
		virtual ~NDFFTW3Cache();

		int compute(fftw_complex* inout) const
		{
                        assert (m_id.m_isalign ? is_simd_aligned(inout) : true);
			fftw_execute_dft(m_plan, inout, inout);
			return 0;
		};

	protected:
		fftw_plan m_plan;
		fftw_complex *m_wrk;	
};

NDFFTW3Cache::NDFFTW3Cache(const NDFFTW3CacheId& id)
:	Cache<NDFFTW3CacheId>(id)
{
        int flags = FFTW_ESTIMATE;
        int sz;
        int i;

        sz = 1;
        for (i = 0; i < m_id.m_rank; ++i) {
                sz *= m_id.m_dims[i];
        }

	m_wrk = (fftw_complex*)fftw_malloc(sz * sizeof(fftw_complex) );
	if (m_wrk == NULL) {
		goto fail_wrk;
	}

        if (!m_id.m_isalign) {
                flags |= FFTW_UNALIGNED;
        } 

	m_plan = fftw_plan_many_dft(m_id.m_rank, 
                                m_id.m_dims, m_id.m_howmany,
                                m_wrk, NULL, 1, sz, 
                                m_wrk, NULL, 1, sz, 
				(id.m_dir > 0 ?  FFTW_FORWARD:FFTW_BACKWARD), 
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

NDFFTW3Cache::~NDFFTW3Cache()
{
	fftw_destroy_plan(m_plan);
	fftw_free(m_wrk);
}

static CacheManager<NDFFTW3CacheId, NDFFTW3Cache> fftw3_cmgr(10);

extern void zfftnd_fftw3(complex_double * inout, int rank,
			   int *dims, int direction, int howmany,
			   int normalize)
{
    int i, sz;
    fftw_complex *ptr = (fftw_complex*)inout;
    NDFFTW3Cache *cache;

    sz = 1;
    for (i = 0; i < rank; ++i) {
        sz *= dims[i];
    }


    cache = fftw3_cmgr.get_cache(NDFFTW3CacheId(rank, dims, howmany, direction, is_simd_aligned(inout)));

    cache->compute(ptr);

    if (normalize) {
        ptr = (fftw_complex*)inout;
        for (i = sz * howmany - 1; i >= 0; --i) {
            *((double *) (ptr)) /= sz;
            *((double *) (ptr++) + 1) /= sz;
        }
    }
}
