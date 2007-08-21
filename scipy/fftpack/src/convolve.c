/*
  Generic functions for computing 1D convolutions of periodic sequences.

  Supported FFT libraries:
    DJBFFT  - optional, used for power-of-two length arrays
    FFTW    - optional
    FFTPACK - used if any of the above libraries is not available

  Author: Pearu Peterson, September 2002
 */

#include "fftpack.h"

/**************** DJBFFT *****************************/
#ifdef WITH_DJBFFT
GEN_CACHE(ddjbfft,(int n)
	  ,double* ptr;
	  ,(caches_ddjbfft[i].n==n)
	  ,caches_ddjbfft[id].ptr = (double*)malloc(sizeof(double)*n);
	  ,free(caches_ddjbfft[id].ptr);
	  ,20)
#endif

/**************** FFTW *****************************/
#ifdef WITH_FFTW
GEN_CACHE(drfftw,(int n)
	  ,rfftw_plan plan1;
	   rfftw_plan plan2;
	  ,(caches_drfftw[i].n==n)
	  ,caches_drfftw[id].plan1 = rfftw_create_plan(n,
		FFTW_REAL_TO_COMPLEX,
		FFTW_IN_PLACE|FFTW_ESTIMATE);
	   caches_drfftw[id].plan2 = rfftw_create_plan(n,
		FFTW_COMPLEX_TO_REAL,
		FFTW_IN_PLACE|FFTW_ESTIMATE);
	  ,rfftw_destroy_plan(caches_drfftw[id].plan1);
  	   rfftw_destroy_plan(caches_drfftw[id].plan2);
	  ,20)
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
	  ,20)
#endif
extern void destroy_convolve_cache(void) {
#ifdef WITH_DJBFFT
  destroy_ddjbfft_caches();
#endif
#ifdef WITH_FFTW
  destroy_drfftw_caches();
#else
  destroy_dfftpack_caches();
#endif
}

/**************** convolve **********************/
extern
void convolve(int n,double* inout,double* omega,int swap_real_imag) {
  int i;
#ifdef WITH_DJBFFT
  double* ptr = NULL;
#endif
#ifdef WITH_FFTW
  rfftw_plan plan1 = NULL;
  rfftw_plan plan2 = NULL;
#else
  double* wsave = NULL;
#endif
#ifdef WITH_DJBFFT
  switch (n) {
  case 2:;case 4:;case 8:;case 16:;case 32:;case 64:;case 128:;case 256:;
  case 512:;case 1024:;case 2048:;case 4096:;case 8192:
    i = get_cache_id_ddjbfft(n);
    ptr = caches_ddjbfft[i].ptr;
    COPYSTD2DJB(inout,ptr,n);
    switch (n) {
#define TMPCASE(N) case N: fftr8_##N(ptr); break
      TMPCASE(2);TMPCASE(4);TMPCASE(8);TMPCASE(16);TMPCASE(32);
      TMPCASE(64);TMPCASE(128);TMPCASE(256);TMPCASE(512);
      TMPCASE(1024);TMPCASE(2048);TMPCASE(4096);TMPCASE(8192);
#undef TMPCASE
    }
    if (swap_real_imag) {
      int n1 = n-1;
      double c;
      ptr[0] *= omega[0];
      ptr[1] *= omega[1];
      for(i=2;i<n1;i+=2) {
	c = ptr[i] * omega[i];
	ptr[i] = ptr[i+1] * omega[i+1];
	ptr[i+1] = c;
      }
    }
    else
      for(i=0;i<n;++i)
	ptr[i] *= omega[i];
    switch (n) {
#define TMPCASE(N)case N:fftr8_un##N(ptr);break
      TMPCASE(2);TMPCASE(4);TMPCASE(8);TMPCASE(16);TMPCASE(32);
      TMPCASE(64);TMPCASE(128);TMPCASE(256);TMPCASE(512);
      TMPCASE(1024);TMPCASE(2048);TMPCASE(4096);TMPCASE(8192);
#undef TMPCASE
    }
    COPYINVDJB2STD2(ptr,inout,n);
    return;
  }
#endif
  {
#ifdef WITH_FFTW
    int l = (n-1)/2+1;
    i = get_cache_id_drfftw(n);
    plan1 = caches_drfftw[i].plan1;
    plan2 = caches_drfftw[i].plan2;
    rfftw_one(plan1, (fftw_real *)inout, NULL);
    if (swap_real_imag) {
      double c;
      inout[0] *= omega[0];
      if (!(n%2))
	inout[n/2] *= omega[n/2];
      for(i=1;i<l;++i) {
	c = inout[i] * omega[i];
	inout[i] = omega[n-i] * inout[n-i];
	inout[n-i] = c;
      }
    }
    else
      for(i=0;i<n;++i)
	inout[i] *= omega[i];
    rfftw_one(plan2, (fftw_real *)inout, NULL);
#else
    i = get_cache_id_dfftpack(n);
    wsave = caches_dfftpack[i].wsave;
    F_FUNC(dfftf,DFFTF)(&n,inout,wsave);
    if (swap_real_imag) {
      double c;
      int n1 = n-1;
      inout[0] *= omega[0];
      if (!(n%2))
	inout[n-1] *= omega[n-1];
      for(i=1;i<n1;i+=2) {
	c = inout[i] * omega[i];
	inout[i] = inout[i+1] * omega[i+1];
	inout[i+1] = c;
      }
    }
    else
      for(i=0;i<n;++i)
	inout[i] *= omega[i];
    F_FUNC(dfftb,DFFTB)(&n,inout,wsave);
#endif
  }
}

/**************** convolve **********************/
extern
void convolve_z(int n,double* inout,double* omega_real,double* omega_imag) {
  int i;
#ifdef WITH_DJBFFT
  double* ptr = NULL;
#endif
#ifdef WITH_FFTW
  rfftw_plan plan1 = NULL;
  rfftw_plan plan2 = NULL;
#else
  double* wsave = NULL;
#endif
#ifdef WITH_DJBFFT
  switch (n) {
  case 2:;case 4:;case 8:;case 16:;case 32:;case 64:;case 128:;case 256:;
  case 512:;case 1024:;case 2048:;case 4096:;case 8192:
    i = get_cache_id_ddjbfft(n);
    ptr = caches_ddjbfft[i].ptr;
    COPYSTD2DJB(inout,ptr,n);
    switch (n) {
#define TMPCASE(N) case N: fftr8_##N(ptr); break
      TMPCASE(2);TMPCASE(4);TMPCASE(8);TMPCASE(16);TMPCASE(32);
      TMPCASE(64);TMPCASE(128);TMPCASE(256);TMPCASE(512);
      TMPCASE(1024);TMPCASE(2048);TMPCASE(4096);TMPCASE(8192);
#undef TMPCASE
    }
    {
      int n1 = n-1;
      double c;
      ptr[0] *= (omega_real[0]+omega_imag[0]);
      ptr[1] *= (omega_real[1]+omega_imag[1]);
      for(i=2;i<n1;i+=2) {
	c = ptr[i] * omega_imag[i];
	ptr[i] *= omega_real[i];
	ptr[i] += ptr[i+1] * omega_imag[i+1];
	ptr[i+1] *= omega_real[i+1];
	ptr[i+1] += c;
      }
    }
    switch (n) {
#define TMPCASE(N)case N:fftr8_un##N(ptr);break
      TMPCASE(2);TMPCASE(4);TMPCASE(8);TMPCASE(16);TMPCASE(32);
      TMPCASE(64);TMPCASE(128);TMPCASE(256);TMPCASE(512);
      TMPCASE(1024);TMPCASE(2048);TMPCASE(4096);TMPCASE(8192);
#undef TMPCASE
    }
    COPYINVDJB2STD2(ptr,inout,n);
    return;
  }
#endif
  {
#ifdef WITH_FFTW
    int l = (n-1)/2+1;
    i = get_cache_id_drfftw(n);
    plan1 = caches_drfftw[i].plan1;
    plan2 = caches_drfftw[i].plan2;
    rfftw_one(plan1, (fftw_real *)inout, NULL);
    {
      double c;
      inout[0] *= (omega_real[0]+omega_imag[0]);
      if (!(n%2))
	inout[n/2] *= (omega_real[n/2]+omega_imag[n/2]);
      for(i=1;i<l;++i) {
	c = inout[i] * omega_imag[i];
	inout[i] *= omega_real[i];
	inout[i] += omega_imag[n-i] * inout[n-i];
	inout[n-i] *= omega_real[n-i];
	inout[n-i] += c;
      }
    }
    rfftw_one(plan2, (fftw_real *)inout, NULL);
#else
    i = get_cache_id_dfftpack(n);
    wsave = caches_dfftpack[i].wsave;
    F_FUNC(dfftf,DFFTF)(&n,inout,wsave);
    {
      double c;
      int n1 = n-1;
      inout[0] *= (omega_real[0]+omega_imag[0]);
      if (!(n%2))
	inout[n-1] *= (omega_real[n-1]+omega_imag[n-1]);
      for(i=1;i<n1;i+=2) {
	c = inout[i] * omega_imag[i];
	inout[i] *= omega_real[i];
	inout[i] += inout[i+1] * omega_imag[i+1];
	inout[i+1] *= omega_real[i+1];
	inout[i+1] += c;
      }
    }
    F_FUNC(dfftb,DFFTB)(&n,inout,wsave);
#endif
  }
}

extern
void init_convolution_kernel(int n,double* omega, int d,
			     double (*kernel_func)(int),
			     int zero_nyquist) {
  /*
    omega[k] = pow(sqrt(-1),d) * kernel_func(k)
    omega[0] = kernel_func(0)
    conjugate(omega[-k]) == omega[k]
   */
#ifdef WITH_DJBFFT
  switch (n) {
  case 2:;case 4:;case 8:;case 16:;case 32:;case 64:;case 128:;case 256:;
  case 512:;case 1024:;case 2048:;case 4096:;case 8192:
    {
      int k,n2=n/2, *f = (int*)malloc(sizeof(int)*(n));
      fftfreq_rtable(f,n);
      for (k=1;k<n;++k)
	if (f[k]>n2) f[k] -= n;
      omega[0] = (*kernel_func)(0)/n;
      switch (d%4) {
      case 0:
	for (k=2;k<n-1;k+=2) {
	  omega[k] = (*kernel_func)(f[k])/n2;
	  omega[k+1] = -omega[k];
	}
	omega[1] = (zero_nyquist?0.0:(*kernel_func)(n2)/n);
	break;
      case 1:;case -3:
	for (k=2;k<n-1;k+=2)
	  omega[k] = omega[k+1] = -(*kernel_func)(f[k])/n2;
	omega[1] = (zero_nyquist?0.0:(*kernel_func)(n2)/n);
	break;
      case 2:;case -2:
	for (k=2;k<n-1;k+=2) {
	  omega[k] = -(*kernel_func)(f[k])/n2;
	  omega[k+1] = -omega[k];
	}
	omega[1] = (zero_nyquist?0.0:-(*kernel_func)(n2)/n);
	break;
      case 3:;case -1:
	for (k=2;k<n-1;k+=2)
	  omega[k] = omega[k+1] = (*kernel_func)(f[k])/n2;
	omega[1] = (zero_nyquist?0.0:-(*kernel_func)(n2)/n);
	break;
      }
      free(f);
    }
    return;
  }
#endif
#ifdef WITH_FFTW
  {
    int k,l=(n-1)/2+1;
    omega[0] = (*kernel_func)(0)/n;;
    switch (d%4) {
    case 0:
      for (k=1;k<l;++k)
	omega[k] = omega[n-k] = (*kernel_func)(k)/n;
      if (!(n%2)) 
	omega[n/2] = (zero_nyquist?0.0:(*kernel_func)(n/2)/n);
      break;
    case 1:;case -3:
      for (k=1;k<l;++k) {
	omega[k] = (*kernel_func)(k)/n;
	omega[n-k] = -omega[k];
      }
      if (!(n%2))
	omega[n/2] = (zero_nyquist?0.0:(*kernel_func)(n/2)/n);
      break;
    case 2:;case -2:
      for (k=1;k<l;++k)
	omega[k] = omega[n-k] = -(*kernel_func)(k)/n;
      if (!(n%2))
	omega[n/2] = (zero_nyquist?0.0:-(*kernel_func)(n/2)/n);
      break;
    case 3:;case -1:
      for (k=1;k<l;++k) {
	omega[k] = -(*kernel_func)(k)/n;
	omega[n-k] = -omega[k];
      }
      if (!(n%2))
        omega[n/2] = (zero_nyquist?0.0:-(*kernel_func)(n/2)/n);
      break;
    }
  }
#else
  {
    int j,k,l=(n%2?n:n-1);
    omega[0] = (*kernel_func)(0)/n;
    switch (d%4) {
    case 0:
      for (k=j=1;j<l;j+=2,++k)
	omega[j] = omega[j+1] = (*kernel_func)(k)/n;
      if (!(n%2))
	omega[n-1] = (zero_nyquist?0.0:(*kernel_func)(k)/n);
      break;
    case 1:;case -3:
      for (k=j=1;j<l;j+=2,++k) {
	omega[j] = (*kernel_func)(k)/n;
	omega[j+1] = -omega[j];
      }
      if (!(n%2))
	omega[n-1] = (zero_nyquist?0.0:(*kernel_func)(k)/n);
      break;
    case 2:;case -2:
      for (k=j=1;j<l;j+=2,++k)
	omega[j] = omega[j+1] = -(*kernel_func)(k)/n;
      if (!(n%2))
	omega[n-1] = (zero_nyquist?0.0:-(*kernel_func)(k)/n);
      break;
    case 3:;case -1:
      for (k=j=1;j<l;j+=2,++k) {
	omega[j] = -(*kernel_func)(k)/n;
	omega[j+1] = -omega[j];
      }
      if (!(n%2))
	omega[n-1] = (zero_nyquist?0.0:-(*kernel_func)(k)/n);
      break;
    }
  }
#endif
}
