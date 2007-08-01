/*
 * Last Change: Wed Aug 01 08:00 PM 2007 J
 *
 * Original code by Pearu Peterson.
 */

/*
 * DJBFFT only implements size 2^N !
 *
 * drfft_def and drfft_def_destroy_cache are the functions used for size different
 * than 2^N
 */
#ifdef WITH_FFTW3
#define drfft_def drfft_fftw3
#define drfft_def_destroy_cache destroy_drfftw3_caches
#elif defined WITH_FFTW
#define drfft_def drfft_fftw
#define drfft_def_destroy_cache destroy_drfftw_caches
#else
#define drfft_def drfft_fftpack
#define drfft_def_destroy_cache destroy_drfftpack_caches
#endif

GEN_CACHE(drdjbfft, (int n)
	    , unsigned int *f;
	    double *ptr;, 
        caches_drdjbfft[i].n == n, 
        caches_drdjbfft[id].f = (unsigned int *) malloc(sizeof(unsigned int) * (n));
	    caches_drdjbfft[id].ptr = (double *) malloc(sizeof(double) * n);
	    fftfreq_rtable(caches_drdjbfft[id].f, n);,
	    free(caches_drdjbfft[id].f); 
        free(caches_drdjbfft[id].ptr);, 
        10)

/**************** ZFFT function **********************/
static void drfft_djbfft(double * inout,
			 int n, int direction, int howmany, int normalize)
{
    int i;
    double *ptr = inout;
    double *ptrc = NULL;
    unsigned int *f = NULL;

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
        i = get_cache_id_drdjbfft(n);
        f = caches_drdjbfft[i].f;
        ptrc = caches_drdjbfft[i].ptr;
    }
    if (f == NULL) {
        drfft_def(inout, n, direction, howmany, normalize);
    }

    switch (direction) {
    case 1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            if (f != NULL) {
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
                COPYDJB2STD(ptrc, ptr, f, n);
            } 
        }
        break;

    case -1:
        for (i = 0; i < howmany; ++i, ptr += n) {
            if (f != NULL) {
                COPYINVSTD2DJB(ptr, ptrc, normalize, f, n);
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
            } 
        }
        break;

    default:
        fprintf(stderr, "drfft: invalid direction=%d\n", direction);
    }

    if (normalize && f != NULL && direction == 1) {
        double d = 1.0 / n;
        ptr = inout;
        for (i = n * howmany - 1; i >= 0; --i) {
            (*(ptr++)) *= d;
        }
    }
}
