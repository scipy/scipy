/*
* DJBFFT only implements size 2^N !
*
* zfft_def and zfft_def_destroy_cache are the functions
* used for size different than 2^N
*/
#ifdef WITH_FFTWORK
#define zfft_def zfft_fftwork
#define zfft_def_destroy_cache destroy_zfftwork_cache
#elif defined WITH_FFTW3
#define zfft_def zfft_fftw3
#define zfft_def_destroy_cache destroy_zfftw3_caches
#elif defined WITH_FFTW
#define zfft_def zfft_fftw
#define zfft_def_destroy_cache destroy_zfftw_caches
#else
#define zfft_def zfft_fftpack
#define zfft_def_destroy_cache destroy_zfftpack_caches
#endif

class DJBFFTCacheId : public CacheId {
        public:
                DJBFFTCacheId(int n) : CacheId(n) {};
};

class DJBFFTCache: public Cache<DJBFFTCacheId> {
        public:
                DJBFFTCache(const DJBFFTCacheId& id);
                virtual ~DJBFFTCache();

                int compute_forward(complex_double * inout) const;
                int compute_backward(complex_double * inout) const;
                int normalize(complex_double * inout) const;

        protected:
                unsigned int* m_f;
                double* m_ptr;
};

DJBFFTCache::DJBFFTCache(const DJBFFTCacheId& id)
:	Cache<DJBFFTCacheId>(id)
{
        int i;
        int n = id.m_n;

        m_f = (unsigned int*)malloc(sizeof(*m_f) * n);
        if (m_f == NULL) {
                goto fail_f;
        }

        m_ptr = (double *)malloc(sizeof(*m_ptr) * 2 * n);
        if (m_ptr == NULL) {
                goto clean_f;
        }

        fftfreq_ctable(m_f, id.m_n);
        for(i = 0; i < n; ++i) {
                m_f[i] = (id.m_n - m_f[i]) % id.m_n;
        }
        return;

clean_f:
        free(m_f);
fail_f:
	throw std::bad_alloc();
}

DJBFFTCache::~DJBFFTCache()
{
        free(m_ptr);
        free(m_f);
}

int DJBFFTCache::compute_forward(complex_double *inout) const 
{
        const int n = m_id.m_n;
        int j;

        complex_double *ptrc = NULL;
        complex_double *ptr = inout;

        ptrc = (complex_double*)m_ptr;

        memcpy(ptrc, ptr, 2 * n * sizeof(double));
        switch (n) {
#define TMPCASE(N) case N:  fftc8_##N(ptrc); break
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
        for (j = 0; j < n; ++j) {
                *(ptr + m_f[j]) = *(ptrc + j);
        }

        return 0;
}

int DJBFFTCache::compute_backward(complex_double *inout) const 
{
        const int n = m_id.m_n;
        int j;

        complex_double *ptrc = NULL;
        complex_double *ptr = inout;

        ptrc = (complex_double*)m_ptr;

        for (j = 0; j < n; ++j) {
                *(ptrc + j) = *(ptr + m_f[j]);
        }
        switch (n) {
#define TMPCASE(N) case N:  fftc8_un##N(ptrc); break
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
        memcpy(ptr, ptrc, 2 * n * sizeof(double));

        return 0;
}

int DJBFFTCache::normalize(complex_double *ptr) const 
{
        int n = m_id.m_n;

        switch (n) {
#define TMPCASE(N) case N:  fftc8_scale##N(ptr); break
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

        return 0;
}

static CacheManager<DJBFFTCacheId, DJBFFTCache> djbfft_cmgr(10);

/* stub to make GEN_PUBLIC_API happy */
static void destroy_zdjbfft_caches()
{
}

/**************** ZFFT function **********************/
static void zfft_djbfft(complex_double * inout,
		 int n, int direction, int howmany, int normalize)
{
	int i;
	complex_double *ptr = inout;
        DJBFFTCache *cache;

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
		cache = djbfft_cmgr.get_cache(DJBFFTCacheId(n));
                break;
        default:
                /* For sizes not handled by djbfft, use default implementation
                 * and returns */
                zfft_def(inout, n, direction, howmany, normalize);
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
                        cache->compute_backward(ptr);
		}
		break;
	default:
		fprintf(stderr, "zfft: invalid direction=%d\n", direction);
	}

	if (normalize) {
		ptr = inout;
                for (i = 0; i < howmany; ++i, ptr += n) {
                        cache->normalize(ptr);
                }
	}
}
