/*
  Interface to various FFT libraries.
  Author: Pearu Peterson, August 2002
 */

#ifndef FFTPACK_H
#define FFTPACK_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef struct {double r,i;} complex_double;
typedef struct {float r,i;} complex_float;

extern
void init_convolution_kernel(int n,double* omega, int d, 
			     double (*kernel_func)(int),
			     int zero_nyquist);
extern
void convolve(int n,double* inout,double* omega,int swap_real_imag);
extern
void convolve_z(int n,double* inout,double* omega_real,double* omega_imag);

extern int ispow2le2e30(int n);
extern int ispow2le2e13(int n);

#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif

/*
  Simple cyclic cache.
 */
#define GEN_CACHE(name,CACHEARG,CACHETYPE,CHECK,MALLOC,FREE,CACHESIZE) \
typedef struct {\
  int n;\
  CACHETYPE \
} cache_type_##name;\
static cache_type_##name caches_##name[CACHESIZE];\
static int nof_in_cache_##name = 0;\
static int last_cache_id_##name = 0;\
static int get_cache_id_##name CACHEARG { \
  int i,id = -1; \
  for (i=0;i<nof_in_cache_##name;i++) \
    if (CHECK) { \
      id=i; \
      break; \
    } \
  if (id>=0) goto exit;\
  if (nof_in_cache_##name<CACHESIZE) {\
    id = nof_in_cache_##name++;\
  } else {\
    id = (last_cache_id_##name<CACHESIZE-1)?last_cache_id_##name+1:0;\
    /*fprintf(stderr,"Removing cache item n=%d\n",caches_##name[id].n);*/\
    FREE \
    caches_##name[id].n = 0;\
  }\
  /*fprintf(stderr,"New cache item n=%d\n",n);*/\
  caches_##name[id].n = n;\
  MALLOC \
 exit:\
  last_cache_id_##name = id;\
  return id;\
}\
void destroy_##name##_cache(void) {\
  int id;\
  for (id=0;id<nof_in_cache_##name;++id) {\
    FREE \
    caches_##name[id].n = 0;\
  }\
  nof_in_cache_##name = last_cache_id_##name = 0;\
}

#endif
