/*
  Interface to various FFT libraries.
  Double complex FFT and IFFT.
  Author: Pearu Peterson, August 2002
 */

#include "fftpack.h"

/**************** FFTWORK *****************************/

#ifdef WITH_FFTWORK
GEN_CACHE(zfftwork,(int n)
	  ,coef_dbl* coef;
	  ,caches_zfftwork[i].n==n
	  ,caches_zfftwork[id].coef = (coef_dbl*)malloc(sizeof(coef_dbl)*(n));
	   fft_coef_dbl(caches_zfftwork[id].coef,n);
	  ,free(caches_zfftwork[id].coef);
	  ,10)
#endif

/**************** DJBFFT *****************************/
#ifdef WITH_DJBFFT
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
#endif

/**************** FFTW *****************************/
#ifdef WITH_FFTW
GEN_CACHE(zfftw,(int n,int d)
	  ,int direction;
	   fftw_plan plan;
	  ,((caches_zfftw[i].n==n) && 
	    (caches_zfftw[i].direction==d))
	  ,caches_zfftw[id].direction = d;
	   caches_zfftw[id].plan = fftw_create_plan(n,
		(d>0?FFTW_FORWARD:FFTW_BACKWARD),
		FFTW_IN_PLACE|FFTW_ESTIMATE);
	  ,fftw_destroy_plan(caches_zfftw[id].plan);
	  ,10)
#else
/**************** FFTPACK ZFFT **********************/
extern void F_FUNC(zfftf,ZFFTF)(int*,double*,double*);
extern void F_FUNC(zfftb,ZFFTB)(int*,double*,double*);
extern void F_FUNC(zffti,ZFFTI)(int*,double*);
GEN_CACHE(zfftpack,(int n)
	  ,double* wsave;
	  ,(caches_zfftpack[i].n==n)
	  ,caches_zfftpack[id].wsave = (double*)malloc(sizeof(double)*(4*n+15));
	   F_FUNC(zffti,ZFFTI)(&n,caches_zfftpack[id].wsave);
	  ,free(caches_zfftpack[id].wsave);
	  ,10)
#endif

extern void destroy_zfft_cache(void) {
#ifdef WITH_FFTWORK
  destroy_zfftwork_caches();
#endif
#ifdef WITH_DJBFFT
  destroy_zdjbfft_caches();
#endif
#ifdef WITH_FFTW
  destroy_zfftw_caches();
#else
  destroy_zfftpack_caches();
#endif
}

/**************** ZFFT function **********************/
extern void zfft(complex_double *inout,
		 int n,int direction,int howmany,int normalize) {
  int i;
  complex_double *ptr = inout;
#ifdef WITH_FFTW
  fftw_plan plan = NULL;
#else
  double* wsave = NULL;
#endif
#ifdef WITH_FFTWORK
  coef_dbl* coef = NULL;
#endif
#ifdef WITH_DJBFFT
  int j;
  complex_double *ptrc = NULL;
  unsigned int *f = NULL;
#endif
#ifdef WITH_FFTWORK
  if (ispow2le2e30(n)) {
    i = get_cache_id_zfftwork(n);
    coef = caches_zfftwork[i].coef;
  } else
#endif
#ifdef WITH_DJBFFT
  switch (n) {
  case 2:;case 4:;case 8:;case 16:;case 32:;case 64:;case 128:;case 256:;
  case 512:;case 1024:;case 2048:;case 4096:;case 8192:
    i = get_cache_id_zdjbfft(n);
    f = caches_zdjbfft[i].f;
    ptrc = (complex_double*)caches_zdjbfft[i].ptr;
  }
  if (f==0)
#endif
#ifdef WITH_FFTW
    plan = caches_zfftw[get_cache_id_zfftw(n,direction)].plan;
#else
    wsave = caches_zfftpack[get_cache_id_zfftpack(n)].wsave;
#endif

  switch (direction) {

  case 1:
    for (i=0;i<howmany;++i,ptr+=n) {
#ifdef WITH_FFTWORK
      if (coef!=NULL) {
	fft_for_cplx_flt((cplx_dbl*)ptr,coef,n);
      } else
#endif
#ifdef WITH_DJBFFT
      if (f!=NULL) {
	memcpy(ptrc,ptr,2*n*sizeof(double));
	switch (n) {
#define TMPCASE(N) case N:  fftc8_##N(ptrc); break
	  TMPCASE(2);TMPCASE(4);TMPCASE(8);TMPCASE(16);TMPCASE(32);
	  TMPCASE(64);TMPCASE(128);TMPCASE(256);TMPCASE(512);
	  TMPCASE(1024);TMPCASE(2048);TMPCASE(4096);TMPCASE(8192);
#undef TMPCASE
	}
	for (j=0;j<n;++j) *(ptr+f[j]) = *(ptrc+j);
      } else
#endif
#ifdef WITH_FFTW
        fftw_one(plan,(fftw_complex*)ptr,NULL);
#else
	F_FUNC(zfftf,ZFFTF)(&n,(double*)(ptr),wsave);
#endif
    }
    break;

  case -1:
    for (i=0;i<howmany;++i,ptr+=n) {
#ifdef WITH_FFTWORK
      if (coef!=NULL) {
	fft_bak_cplx_flt((cplx_dbl*)ptr,coef,n);
      } else
#endif
#ifdef WITH_DJBFFT
      if (f!=NULL) {
	for (j=0;j<n;++j) *(ptrc+j) = *(ptr+f[j]);
	switch (n) {
#define TMPCASE(N) case N:  fftc8_un##N(ptrc); break
	  TMPCASE(2);TMPCASE(4);TMPCASE(8);TMPCASE(16);TMPCASE(32);
	  TMPCASE(64);TMPCASE(128);TMPCASE(256);TMPCASE(512);
	  TMPCASE(1024);TMPCASE(2048);TMPCASE(4096);TMPCASE(8192);
#undef TMPCASE
	}
	memcpy(ptr,ptrc,2*n*sizeof(double));
      } else
#endif
#ifdef WITH_FFTW
        fftw_one(plan,(fftw_complex*)ptr,NULL);
#else
	F_FUNC(zfftb,ZFFTB)(&n,(double*)(ptr),wsave);
#endif
    }
    break;

  default:
    fprintf(stderr,"zfft: invalid direction=%d\n",direction);
  }

  if (normalize) {
    ptr = inout;
#ifdef WITH_DJBFFT
    if (f!=NULL) {
      for (i=0;i<howmany;++i,ptr+=n) {
	switch (n) {
#define TMPCASE(N) case N:  fftc8_scale##N(ptr); break
	  TMPCASE(2);TMPCASE(4);TMPCASE(8);TMPCASE(16);TMPCASE(32);
	  TMPCASE(64);TMPCASE(128);TMPCASE(256);TMPCASE(512);
	  TMPCASE(1024);TMPCASE(2048);TMPCASE(4096);TMPCASE(8192);
#undef TMPCASE
	}
      }
    } else
#endif
      for (i=n*howmany-1;i>=0;--i) {
	*((double*)(ptr)) /= n;
	*((double*)(ptr++)+1) /= n;
      }
  }
}
