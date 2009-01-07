/*
  Interface to various FFT libraries.
  Double complex FFT and IFFT with zero imaginary input.
  Author: Pearu Peterson, August 2002
 */

#include "fftpack.h"

extern void drfft(double *inout,int n,int direction,int howmany,int normalize);
extern void rfft(float *inout,int n,int direction,int howmany,int normalize);

extern void zrfft(complex_double *inout,
		  int n,int direction,int howmany,int normalize) {
  int i,j,k;
  double* ptr = (double *)inout;
  switch (direction) {
    case 1:
      for (i=0;i<howmany;++i,ptr+=2*n) {
	*(ptr+1) = *ptr;
	for(j=2,k=3;j<n;++j,++k)
	  *(ptr+k) = *(ptr+2*j);
	drfft(ptr+1,n,1,1,normalize);
	*ptr = *(ptr+1);
	*(ptr+1) = 0.0;
	if (!(n%2))
	  *(ptr+n+1) = 0.0;
	for(j=2,k=2*n-2;j<n;j+=2,k-=2) {
	  *(ptr+k) = *(ptr+j);
	  *(ptr+k+1) = -(*(ptr+j+1));
	}
      }
      break;
  case -1:
    for (i=0;i<howmany;++i,ptr+=2*n) {
      *(ptr+1) = (*ptr);
      for(j=1,k=2;j<n;++j,++k)
	*(ptr+k) = (*(ptr+2*j));
      drfft(ptr+1,n,1,1,normalize);
      *ptr = *(ptr+1);
      *(ptr+1) = 0.0;
      if (!(n%2))
	*(ptr+n+1) = 0.0;
      for(j=2,k=2*n-2;j<n;j+=2,k-=2) {
	double d;
	*(ptr+k) = *(ptr+j);
	d = *(ptr+j+1);
	*(ptr+k+1) = d; 
	*(ptr+j+1) = -d;
      }
    }
    break;
  default:
    fprintf(stderr,"zrfft: invalid direction=%d\n",direction);
  }
}

extern void crfft(complex_float *inout,
		  int n,int direction,int howmany,int normalize) {
  int i,j,k;
  float* ptr = (float *)inout;
  switch (direction) {
    case 1:
      for (i=0;i<howmany;++i,ptr+=2*n) {
	*(ptr+1) = *ptr;
	for(j=2,k=3;j<n;++j,++k)
	  *(ptr+k) = *(ptr+2*j);
	rfft(ptr+1,n,1,1,normalize);
	*ptr = *(ptr+1);
	*(ptr+1) = 0.0;
	if (!(n%2))
	  *(ptr+n+1) = 0.0;
	for(j=2,k=2*n-2;j<n;j+=2,k-=2) {
	  *(ptr+k) = *(ptr+j);
	  *(ptr+k+1) = -(*(ptr+j+1));
	}
      }
      break;
  case -1:
    for (i=0;i<howmany;++i,ptr+=2*n) {
      *(ptr+1) = (*ptr);
      for(j=1,k=2;j<n;++j,++k)
	*(ptr+k) = (*(ptr+2*j));
      rfft(ptr+1,n,1,1,normalize);
      *ptr = *(ptr+1);
      *(ptr+1) = 0.0;
      if (!(n%2))
	*(ptr+n+1) = 0.0;
      for(j=2,k=2*n-2;j<n;j+=2,k-=2) {
	float d;
	*(ptr+k) = *(ptr+j);
	d = *(ptr+j+1);
	*(ptr+k+1) = d; 
	*(ptr+j+1) = -d;
      }
    }
    break;
  default:
    fprintf(stderr,"crfft: invalid direction=%d\n",direction);
  }
}
