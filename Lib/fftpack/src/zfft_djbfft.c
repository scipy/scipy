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

GEN_CACHE(zdjbfft,(int n)
	  ,unsigned int* f;
	   double* ptr;
	  ,caches_zdjbfft[i].n==n
	  ,caches_zdjbfft[id].f = (unsigned int*)malloc(sizeof(unsigned int)*(n));
	   caches_zdjbfft[id].ptr = (double*)malloc(sizeof(double)*(2*n));
	   fftfreq_ctable(caches_zdjbfft[id].f,n);
           for(i=0;i<n;++i) 
	     caches_zdjbfft[id].f[i] = (n-caches_zdjbfft[id].f[i])%n;
	  ,free(caches_zdjbfft[id].f);
	   free(caches_zdjbfft[id].ptr);
	  ,10)

/**************** ZFFT function **********************/
static void zfft_djbfft(complex_double * inout,
		 int n, int direction, int howmany, int normalize)
{
	int i;
	complex_double *ptr = inout;
	int j;
	complex_double *ptrc = NULL;
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
		i = get_cache_id_zdjbfft(n);
		f = caches_zdjbfft[i].f;
		ptrc = (complex_double *) caches_zdjbfft[i].ptr;
	}
	if (f == 0) {
                zfft_def(inout, n, direction, howmany, normalize);
	}

	switch (direction) {
	case 1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			if (f != NULL) {
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
					*(ptr + f[j]) = *(ptrc + j);
				}
                        }

		}
		break;

	case -1:
		for (i = 0; i < howmany; ++i, ptr += n) {
			if (f != NULL) {
				for (j = 0; j < n; ++j) {
					*(ptrc + j) = *(ptr + f[j]);
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
			}
		}
		break;
	default:
		fprintf(stderr, "zfft: invalid direction=%d\n", direction);
	}

	if (normalize) {
		ptr = inout;
		if (f != NULL) {
			for (i = 0; i < howmany; ++i, ptr += n) {
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
			}
		}
	}
}
