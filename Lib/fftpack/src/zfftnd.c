/*
  Interface to various FFT libraries.
  Double complex FFT and IFFT, arbitrary dimensions.
  Author: Pearu Peterson, August 2002
 */
#include "fftpack.h"

/**************** FFTW *****************************/
#ifdef WITH_FFTW
static
int equal_dims(int rank,int *dims1,int *dims2) {
  int i;
  for (i=0;i<rank;++i)
    if (dims1[i]!=dims2[i])
      return 0;
  return 1;
}
GEN_CACHE(zfftwnd,(int n,int *dims,int d,int flags)
	  ,int direction;
	   int *dims;
	   fftwnd_plan plan;
	  ,((caches_zfftwnd[i].n==n) && 
	   (caches_zfftwnd[i].direction==d) &&
	    (equal_dims(n,caches_zfftwnd[i].dims,dims)))
	  ,caches_zfftwnd[id].direction = d;
	   caches_zfftwnd[id].n = n;
           caches_zfftwnd[id].dims = (int*)malloc(sizeof(int)*n);
           memcpy(caches_zfftwnd[id].dims,dims,sizeof(int)*n);
	   caches_zfftwnd[id].plan = fftwnd_create_plan(n,dims,
		(d>0?FFTW_FORWARD:FFTW_BACKWARD),flags);
	  ,fftwnd_destroy_plan(caches_zfftwnd[id].plan);
	   free(caches_zfftwnd[id].dims);
	  ,10)
#else
GEN_CACHE(zfftnd,(int n,int rank)
	  ,complex_double *ptr;
	  int *iptr;
	  int rank;
	  ,((caches_zfftnd[i].n==n)&&(caches_zfftnd[i].rank==rank))
	  ,caches_zfftnd[id].n = n;
           caches_zfftnd[id].ptr = (complex_double*)malloc(2*sizeof(double)*n);
	   caches_zfftnd[id].iptr = (int*)malloc(4*rank*sizeof(int));
	  ,free(caches_zfftnd[id].ptr);
	  free(caches_zfftnd[id].iptr);
	  ,10)
#endif

extern void destroy_zfftnd_cache(void) {
#ifdef WITH_FFTW
  destroy_zfftwnd_caches();
#else
#endif
}
#ifdef WITH_FFTW
#else
static
/*inline : disabled because MSVC6.0 fails to compile it. */
int next_comb(int *ia,int *da,int m) {
  while (m>=0 && ia[m]==da[m]) ia[m--] = 0;
  if (m<0) return 0;
  ia[m]++;
  return 1;
}
static
void flatten(complex_double *dest,complex_double *src,
	     int rank,int strides_axis,int dims_axis,int unflat,int *tmp) {
  int *new_strides = tmp+rank;
  int *new_dims = tmp+2*rank;
  int *ia = tmp+3*rank;
  int rm1=rank-1,rm2=rank-2;
  int i,j,k;
  for (i=0;i<rm2;++i) ia[i]=0;ia[rm2] = -1;
  j = 0;
  if (unflat)
    while (next_comb(ia,new_dims,rm2)) {
      k = 0;
      for(i=0;i<rm1;++i)
	k += ia[i] * new_strides[i];
      for(i=0;i<dims_axis;++i)
	*(dest+k+i*strides_axis) = *(src+j++);
    }
  else
    while (next_comb(ia,new_dims,rm2)) {
      k = 0;
      for(i=0;i<rm1;++i)
	k += ia[i] * new_strides[i];
      for(i=0;i<dims_axis;++i)
	*(dest+j++) = *(src+k+i*strides_axis);
    }
}
extern void zfft(complex_double *inout,
		 int n,int direction,int howmany,int normalize);
#endif
/**************** ZFFTND function **********************/
extern void zfftnd(complex_double *inout,int rank,
		   int *dims,int direction,int howmany,int normalize) {
  int i,sz;
  complex_double *ptr = inout;
#ifdef WITH_FFTW
  fftwnd_plan plan = NULL;
#else
  int axis;
  complex_double *tmp;
  int *itmp;
  int k,j;
#endif
  sz = 1;
  for(i=0;i<rank;++i)
    sz *= dims[i];
#ifdef WITH_FFTW
  i = get_cache_id_zfftwnd(rank,dims,direction,FFTW_IN_PLACE|FFTW_ESTIMATE);
  plan = caches_zfftwnd[i].plan;
  for (i=0;i<howmany;++i,ptr+=sz)
    fftwnd_one(plan,(fftw_complex*)ptr,NULL);
  if (normalize) {
    ptr = inout;
    for (i=sz*howmany-1;i>=0;--i) {
      *((double*)(ptr)) /= sz;
      *((double*)(ptr++)+1) /= sz;
    }
  }
#else
  zfft(ptr,dims[rank-1],direction,howmany*sz/dims[rank-1],normalize);

  i = get_cache_id_zfftnd(sz,rank); /*Get cache*/
  tmp = caches_zfftnd[i].ptr;
  itmp = caches_zfftnd[i].iptr;

  itmp[rank-1] = 1; /*Calculate strides*/
  for (i=2;i<=rank;++i)
    itmp[rank-i] = itmp[rank-i+1]*dims[rank-i+1];

  for (i=0;i<howmany;++i,ptr+=sz) {
    for(axis=0;axis<rank-1;++axis) {
      for (k=j=0;k<rank;++k)
	if (k!=axis) {
	  *(itmp+rank+j) = itmp[k];
	  *(itmp+2*rank+j++)= dims[k]-1;
	}
      flatten(tmp,ptr,rank,itmp[axis],dims[axis],0,itmp);
      zfft(tmp,dims[axis],direction,sz/dims[axis],normalize);
      flatten(ptr,tmp,rank,itmp[axis],dims[axis],1,itmp);
    }
  }
#endif
}
