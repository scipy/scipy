/*
  Interface to various FFT libraries.
  Double real FFT and IFFT.
  Author: Pearu Peterson, August 2002
 */

#include "fftpack.h"

/**************** DJBFFT *****************************/
#ifdef WITH_DJBFFT
GEN_CACHE(ddjbfft,(int n)
	  ,unsigned int* f;
	   double* ptr;
	  ,caches_ddjbfft[i].n==n
	  ,caches_ddjbfft[id].f = (unsigned int*)malloc(sizeof(unsigned int)*(n));
	   caches_ddjbfft[id].ptr = (double*)malloc(sizeof(double)*n);
	   fftfreq_rtable(caches_ddjbfft[id].f,n);
	  ,free(caches_ddjbfft[id].f);
	   free(caches_ddjbfft[id].ptr);
	  ,10)
#endif

/**************** FFTW *****************************/
#ifdef WITH_FFTW
GEN_CACHE(drfftw,(int n,int d,int flags)
	  ,int direction;
	   int flags;
	   rfftw_plan plan;
	   double *ptr;
	  ,((caches_drfftw[i].n==n) && 
	   (caches_drfftw[i].direction==d) &&
	   (caches_drfftw[i].flags==flags))
	  ,caches_drfftw[id].direction = d;
	   caches_drfftw[id].flags = flags;
	   caches_drfftw[id].plan = rfftw_create_plan(n,
		(d>0?FFTW_REAL_TO_COMPLEX:FFTW_COMPLEX_TO_REAL),flags);
	   caches_drfftw[id].ptr = (double*)malloc(sizeof(double)*(n));
	  ,rfftw_destroy_plan(caches_drfftw[id].plan);
	  ,10)
#else
/**************** FFTPACK ZFFT **********************/
extern void F_FUNC(dfftf,DFFTF)(int*,double*,double*);
extern void F_FUNC(dfftb,DFFTB)(int*,double*,double*);
extern void F_FUNC(dffti,DFFTI)(int*,double*);
GEN_CACHE(dfftpack,(int n)
	  ,double* wsave;
	  ,(caches_dfftpack[i].n==n)
	  ,caches_dfftpack[id].wsave = (double*)malloc(sizeof(double)*(2*n+15));
	   F_FUNC(dffti,DFFTI)(&n,caches_dfftpack[id].wsave);
	  ,free(caches_dfftpack[id].wsave);
	  ,10)
#endif

extern void destroy_drfft_cache(void) {
#ifdef WITH_DJBFFT
  destroy_ddjbfft_caches();
#endif
#ifdef WITH_FFTW
  destroy_drfftw_caches();
#else
  destroy_dfftpack_caches();
#endif
}

/**************** DRFFT function **********************/


extern void drfft(double *inout,
		  int n,int direction,int howmany,int normalize) {
  int i;
  double *ptr = inout;
#if defined(WITH_FFTW) || defined(WITH_DJBFFT)
  double *ptrc = NULL;
#endif
#ifdef WITH_FFTW
  rfftw_plan plan = NULL;
#else
  double* wsave = NULL;
#endif
#ifdef WITH_DJBFFT
  unsigned int *f = NULL;
#endif

#ifdef WITH_DJBFFT
  switch (n) {
  case 2:;case 4:;case 8:;case 16:;case 32:;case 64:;case 128:;case 256:;
  case 512:;case 1024:;case 2048:;case 4096:;case 8192:
    i = get_cache_id_ddjbfft(n);
    f = caches_ddjbfft[i].f;
    ptrc = caches_ddjbfft[i].ptr;
  }
  if (f==NULL)
#endif
#ifdef WITH_FFTW
    {
      i = get_cache_id_drfftw(n,direction,FFTW_IN_PLACE|FFTW_ESTIMATE);
      plan = caches_drfftw[i].plan;
      ptrc = caches_drfftw[i].ptr;
    }
#else
    wsave = caches_dfftpack[get_cache_id_dfftpack(n)].wsave;
#endif

  switch (direction) {

  case 1:
    for (i=0;i<howmany;++i,ptr+=n) {
#ifdef WITH_DJBFFT
      if (f!=NULL) {
	COPYSTD2DJB(ptr,ptrc,n);
	switch (n) {
#define TMPCASE(N) case N: fftr8_##N(ptrc); break
	  TMPCASE(2);TMPCASE(4);TMPCASE(8);TMPCASE(16);TMPCASE(32);
	  TMPCASE(64);TMPCASE(128);TMPCASE(256);TMPCASE(512);
	  TMPCASE(1024);TMPCASE(2048);TMPCASE(4096);TMPCASE(8192);
#undef TMPCASE
	}
	COPYDJB2STD(ptrc,ptr,f,n);
      } else
#endif
#ifdef WITH_FFTW
	{
	  memcpy(ptrc,ptr,sizeof(double)*n);
	  rfftw(plan,1,(fftw_real*)ptrc,1,1,NULL,1,1);
	  COPYRFFTW2STD(ptrc,ptr,n);
	}
#else
	F_FUNC(dfftf,DFFTF)(&n,ptr,wsave);
#endif
    }
    break;

  case -1:
    for (i=0;i<howmany;++i,ptr+=n) {
#ifdef WITH_DJBFFT
      if (f!=NULL) {
	COPYINVSTD2DJB(ptr,ptrc,normalize,f,n);
	switch (n) {
#define TMPCASE(N)case N:if(normalize)fftr8_scale##N(ptrc);fftr8_un##N(ptrc);break
	  TMPCASE(2);TMPCASE(4);TMPCASE(8);TMPCASE(16);TMPCASE(32);
	  TMPCASE(64);TMPCASE(128);TMPCASE(256);TMPCASE(512);
	  TMPCASE(1024);TMPCASE(2048);TMPCASE(4096);TMPCASE(8192);
#undef TMPCASE
	}
	COPYINVDJB2STD(ptrc,ptr,n);
      } else
#endif
#ifdef WITH_FFTW
	{
	  COPYINVRFFTW2STD(ptr,ptrc,n);
	  rfftw(plan,1,(fftw_real*)ptrc,1,1,NULL,1,1);
	  memcpy(ptr,ptrc,sizeof(double)*n);
	}
#else
	F_FUNC(dfftb,DFFTB)(&n,ptr,wsave);
#endif
    }
    break;

  default:
    fprintf(stderr,"drfft: invalid direction=%d\n",direction);
  }

  if (normalize 
#ifdef WITH_DJBFFT
      && (f==NULL||direction==1)
#endif
      ) {
    double d = 1.0/n;
    ptr = inout;
    for (i=n*howmany-1;i>=0;--i)
      (*(ptr++)) *= d;
  }
}



