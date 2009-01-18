/*
  Interface to various FFT libraries.
  Double complex FFT and IFFT.
  Author: Pearu Peterson, August 2002
 */
#include <math.h>

#include "fftpack.h"

enum normalize {
    DCT_NORMALIZE_NO = 0,
    DCT_NORMALIZE_ORTHONORMAL = 1
};

extern void F_FUNC(dcosti,DCOSTI)(int*,double*);
extern void F_FUNC(dcost,DCOST)(int*,double*,double*);
extern void F_FUNC(dcosqi,DCOSQI)(int*,double*);
extern void F_FUNC(dcosqb,DCOSQB)(int*,double*,double*);
extern void F_FUNC(dcosqf,DCOSQF)(int*,double*,double*);

extern void F_FUNC(costi,DCOSTI)(int*,float*);
extern void F_FUNC(cost,COST)(int*,float*,float*);
extern void F_FUNC(cosqi,COSQI)(int*,float*);
extern void F_FUNC(cosqb,COSQB)(int*,float*,float*);
extern void F_FUNC(cosqf,COSQF)(int*,float*,float*);

GEN_CACHE(dct1,(int n)
      ,double* wsave;
      ,(caches_dct1[i].n==n)
      ,caches_dct1[id].wsave = (double*)malloc(sizeof(double)*(3*n+15));
       F_FUNC(dcosti,DCOSTI)(&n,caches_dct1[id].wsave);
      ,free(caches_dct1[id].wsave);
      ,10)

GEN_CACHE(dct2,(int n)
      ,double* wsave;
      ,(caches_dct2[i].n==n)
      ,caches_dct2[id].wsave = (double*)malloc(sizeof(double)*(3*n+15));
       F_FUNC(dcosqi,DCOSQI)(&n,caches_dct2[id].wsave);
      ,free(caches_dct2[id].wsave);
      ,10)

GEN_CACHE(fdct1,(int n)
      ,float* wsave;
      ,(caches_fdct1[i].n==n)
      ,caches_fdct1[id].wsave = (float*)malloc(sizeof(float)*(3*n+15));
       F_FUNC(costi,COSTI)(&n,caches_fdct1[id].wsave);
      ,free(caches_fdct1[id].wsave);
      ,10)

GEN_CACHE(fdct2,(int n)
      ,float* wsave;
      ,(caches_fdct2[i].n==n)
      ,caches_fdct2[id].wsave = (float*)malloc(sizeof(float)*(3*n+15));
       F_FUNC(cosqi,DCOSQI)(&n,caches_fdct2[id].wsave);
      ,free(caches_fdct2[id].wsave);
      ,10)

void dct1(double * inout, int n, int howmany, int normalize)
{
    int i;
    double *ptr = inout;
    double *wsave = NULL;

    wsave = caches_dct1[get_cache_id_dct1(n)].wsave;

    for (i = 0; i < howmany; ++i, ptr += n) {
        dcost_(&n, (double*)(ptr), wsave);
    }

    if (normalize) {
        fprintf(stderr, "dct1: normalize not yet supported=%d\n",
                normalize);
    } else {
        ptr = inout;
        /* 0.5 coeff comes from fftpack defining DCT as
         * 4 * sum(cos(something)), whereas most definition
         * use 2 */
        for (i = n * howmany - 1; i >= 0; --i, ++ptr) {
            *((double *) (ptr)) *= 0.5;
        }
    }
}

void dct2(double * inout, int n, int howmany, int normalize)
{
    int i, j;
    double *ptr = inout;
    double *wsave = NULL;
    double n1, n2;

    wsave = caches_dct2[get_cache_id_dct2(n)].wsave;

    for (i = 0; i < howmany; ++i, ptr += n) {
        dcosqb_(&n, (double *) (ptr), wsave);

    }

    switch (normalize) {
        case DCT_NORMALIZE_NO:
            ptr = inout;
            /* 0.5 coeff comes from fftpack defining DCT as
             * 4 * sum(cos(something)), whereas most definition
             * use 2 */
            for (i = n * howmany - 1; i >= 0; --i, ++ptr) {
                *((double *) (ptr)) *= 0.5;
            }
            break;
        case DCT_NORMALIZE_ORTHONORMAL:
            ptr = inout;
            /* 0.5 coeff comes from fftpack defining DCT as
             * 4 * sum(cos(something)), whereas most definition
             * use 2 */
            n1 = 0.25 * sqrt(1./n);
            n2 = 0.25 * sqrt(2./n);
            for (i = 0; i < howmany; ++i, ptr+=n) {
                ptr[0] *= n1;
                for (j = 1; j < n; ++j) {
                    ptr[j] *= n2;
                }
            }
            break;
        default:
            fprintf(stderr, "dct2: normalize not yet supported=%d\n",
                    normalize);
            break;
    }
}

void fdct1(float * inout, int n, int howmany, int normalize)
{
    int i;
    float *ptr = inout;
    float *wsave = NULL;

    wsave = caches_fdct1[get_cache_id_fdct1(n)].wsave;

    for (i = 0; i < howmany; ++i, ptr += n) {
        cost_(&n, (float*)(ptr), wsave);
    }

    if (normalize) {
        fprintf(stderr, "dct1: normalize not yet supported=%d\n",
                normalize);
    } else {
        ptr = inout;
        /* 0.5 coeff comes from fftpack defining DCT as
         * 4 * sum(cos(something)), whereas most definition
         * use 2 */
        for (i = n * howmany - 1; i >= 0; --i, ++ptr) {
            *((float *) (ptr)) *= 0.5;
        }
    }
}

void fdct2(float * inout, int n, int howmany, int normalize)
{
    int i, j;
    float *ptr = inout;
    float *wsave = NULL;
    float n1, n2;

    wsave = caches_fdct2[get_cache_id_fdct2(n)].wsave;

    for (i = 0; i < howmany; ++i, ptr += n) {
        cosqb_(&n, (float *) (ptr), wsave);

    }

    switch (normalize) {
        case DCT_NORMALIZE_NO:
            ptr = inout;
            /* 0.5 coeff comes from fftpack defining DCT as
             * 4 * sum(cos(something)), whereas most definition
             * use 2 */
            for (i = n * howmany - 1; i >= 0; --i, ++ptr) {
                *((float *) (ptr)) *= 0.5;
            }
            break;
        case DCT_NORMALIZE_ORTHONORMAL:
            ptr = inout;
            /* 0.5 coeff comes from fftpack defining DCT as
             * 4 * sum(cos(something)), whereas most definition
             * use 2 */
            n1 = 0.25 * sqrt(1./n);
            n2 = 0.25 * sqrt(2./n);
            for (i = 0; i < howmany; ++i, ptr+=n) {
                ptr[0] *= n1;
                for (j = 1; j < n; ++j) {
                    ptr[j] *= n2;
                }
            }
            break;
        default:
            fprintf(stderr, "dct2: normalize not yet supported=%d\n",
                    normalize);
            break;
    }
}
