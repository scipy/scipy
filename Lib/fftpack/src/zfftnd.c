/*
  Interface to various FFT libraries.
  Double complex FFT and IFFT, arbitrary dimensions.
  Author: Pearu Peterson, August 2002
 */
#include "fftpack.h"

/* The following macro convert private backend specific function to the public
 * functions exported by the module  */
#define GEN_PUBLIC_API(name) \
void destroy_zfftnd_cache(void)\
{\
        destroy_zfftnd_##name##_caches();\
}\
\
void zfftnd(complex_double * inout, int rank,\
		           int *dims, int direction, int howmany, int normalize)\
{\
        zfftnd_##name(inout, rank, dims, direction, howmany, normalize);\
}

#if defined(WITH_FFTW) || defined(WITH_MKL)
static
int equal_dims(int rank,int *dims1,int *dims2) {
  int i;
  for (i=0;i<rank;++i)
    if (dims1[i]!=dims2[i])
      return 0;
  return 1;
}
#endif

#ifdef WITH_FFTW3
    #include "zfftnd_fftw3.c"
    GEN_PUBLIC_API(fftw3)
#elif defined WITH_FFTW
    #include "zfftnd_fftw.c"
    GEN_PUBLIC_API(fftw)
#elif defined WITH_MKL
    #include "zfftnd_mkl.c"
    GEN_PUBLIC_API(mkl)
#else /* Use fftpack by default */
    #include "zfftnd_fftpack.c"
    GEN_PUBLIC_API(fftpack)
#endif

#if 0
/**************** INTEL MKL **************************/
#ifdef WITH_MKL
long* convert_dims(int n, int *dims)
{
    long * ndim;
    int i;
    ndim = (long*)malloc(sizeof(long)*n);
    for(i=0;i<n;i++){ ndim[i] = (long)dims[i];}
    return ndim;
}

GEN_CACHE(zmklfftnd,(int n,int *dims)
	  ,DFTI_DESCRIPTOR_HANDLE desc_handle;
	   int *dims;
       long *ndims;
      ,((caches_zmklfftnd[i].n==n) &&
	   (equal_dims(n,caches_zmklfftnd[i].dims,dims)))
	  ,caches_zmklfftnd[id].ndims = convert_dims(n, dims);
       caches_zmklfftnd[id].n = n;
       caches_zmklfftnd[id].dims = (int*)malloc(sizeof(int)*n);
       memcpy(caches_zmklfftnd[id].dims,dims,sizeof(int)*n);
       DftiCreateDescriptor(&caches_zmklfftnd[id].desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, (long)n, caches_zmklfftnd[id].ndims); 
       DftiCommitDescriptor(caches_zmklfftnd[id].desc_handle);
	  ,DftiFreeDescriptor(&caches_zmklfftnd[id].desc_handle);
	   free(caches_zmklfftnd[id].dims);
       free(caches_zmklfftnd[id].ndims);
	  ,10)

/**************** FFTW3 *****************************/
#elif defined WITH_FFTW3
/* Don't worry about caching for fftw3 - plans take specific arrays and
 * keeping around a lot of memory for such a small speed up isn't
 * worth it.
 */
/**************** FFTW2 *****************************/
#elif defined WITH_FFTW
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
#ifdef WITH_MKL
  destroy_zmklfftnd_caches();
#elif defined WITH_FFTW3
#elif defined WITH_FFTW
  destroy_zfftwnd_caches();
#else
  destroy_zfftnd_caches();
#endif
}
#if defined(WITH_FFTW) || defined(WITH_FFTW3) || defined(WITH_MKL)
#else
#endif
/**************** ZFFTND function **********************/
extern void zfftnd(complex_double *inout,int rank,
		   int *dims,int direction,int howmany,int normalize) {
  int i,sz;
  complex_double *ptr = inout;
#if defined WITH_MKL
  DFTI_DESCRIPTOR_HANDLE desc_handle;
#elif defined WITH_FFTW3
  fftw_plan plan = NULL;
#elif defined WITH_FFTW
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
#ifdef WITH_MKL
  desc_handle = caches_zmklfftnd[get_cache_id_zmklfftnd(rank, dims)].desc_handle;
  for (i=0;i<howmany;++i,ptr+=sz){
    if (direction == 1){
      DftiComputeForward(desc_handle, (double *)ptr);
    }else if (direction == -1){
      DftiComputeBackward(desc_handle, (double *)ptr);
    }
  }
  if (normalize) {
    ptr = inout;
    for (i=sz*howmany-1;i>=0;--i) {
      *((double*)(ptr)) /= sz;
      *((double*)(ptr++)+1) /= sz;
    }
  }
#elif defined WITH_FFTW3
  plan = fftw_plan_many_dft(rank,dims,howmany,
			    (fftw_complex*)ptr,NULL,1,sz,
			    (fftw_complex*)ptr,NULL,1,sz,
			    (direction>0?FFTW_FORWARD:FFTW_BACKWARD),
			    FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
    /* note that fftw_malloc of array *could* lead
     * to faster fft here for processors with SIMD acceleration,
     * but would require more memory and an array memcpy
     */
  if (normalize) {
    ptr = inout;
    for (i=sz*howmany-1;i>=0;--i) {
      *((double*)(ptr)) /= sz;
      *((double*)(ptr++)+1) /= sz;
    }
  }
#elif defined WITH_FFTW
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
#endif
