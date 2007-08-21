/*
 * Last Change: Tue Nov 28 03:00 PM 2006 J
 * 
 * vim:syntax=c
 */
#include <stdlib.h>     /* for malloc and co */
#include <stdio.h>   

#include "levinson.h"
#include "autocorr_nofft.h"

#include "lpc.h"


/* 
 * (float version) Compute lpc coeff (order x coefficients) of a 
 * contiguous array 
 *
 * err is a double, coeff must be able to contain order+1 elements, and kcoeff
 * order elements
 */
int flt_lpc(const float* signal, 
        size_t size, size_t order, float* coeff, 
        float* kcoeff, float* err)
{
    size_t  i, nbuff, ncache;
    float *buff, *cache, biasnorm;
    int     status;

    biasnorm    = 1.0/size;
    nbuff       = order + 1;
    ncache      = order;

    buff    = malloc(sizeof(*buff) * nbuff);
    if (buff == NULL) {
        status  = -2;
        goto fail_buff_malloc;
    }

    cache   = malloc(sizeof(*cache) * ncache);
    if (cache == NULL) {
        status  = -2;
        goto fail_cache_malloc;
    }

    /*
     * Compute the autocorreleation up to lag order, normalized by the 
     * size of the signal
     */
    flt_xcorr_nofft_1d(signal, size, buff, order);
    for(i = 0; i < nbuff; ++i) {
        buff[i] *= biasnorm;
    }

    /*
     * Compute the inverse coefficients using (simple) levinson recursive algo
     */
    status  = flt_levinson1d(buff, order, 
                coeff, err, kcoeff, cache);
    if (status) {
        status  = -1;
        goto fail_levinson;
    }
              
    free(cache);
    free(buff);

    return 0;

fail_levinson:
    free(cache);
fail_cache_malloc:
    free(buff);
fail_buff_malloc:
    fprintf(stderr, "Failure\n");
    return status;
}


/* 
 * (double version) Compute lpc coeff (order x coefficients) of a 
 * contiguous array 
 *
 * err is a double, coeff must be able to contain order+1 elements, and kcoeff
 * order elements
 */
int dbl_lpc(const double* signal, 
        size_t size, size_t order, double* coeff, 
        double* kcoeff, double* err)
{
    size_t  i, nbuff, ncache;
    double *buff, *cache, biasnorm;
    int     status;

    biasnorm    = 1.0/size;
    nbuff       = order + 1;
    ncache      = order;

    buff    = malloc(sizeof(*buff) * nbuff);
    if (buff == NULL) {
        status  = -2;
        goto fail_buff_malloc;
    }

    cache   = malloc(sizeof(*cache) * ncache);
    if (cache == NULL) {
        status  = -2;
        goto fail_cache_malloc;
    }

    /*
     * Compute the autocorreleation up to lag order, normalized by the 
     * size of the signal
     */
    dbl_xcorr_nofft_1d(signal, size, buff, order);
    for(i = 0; i < nbuff; ++i) {
        buff[i] *= biasnorm;
    }

    /*
     * Compute the inverse coefficients using (simple) levinson recursive algo
     */
    status  = dbl_levinson1d(buff, order, 
                coeff, err, kcoeff, cache);
    if (status) {
        status  = -1;
        goto fail_levinson;
    }
              
    free(cache);
    free(buff);

    return 0;

fail_levinson:
    free(cache);
fail_cache_malloc:
    free(buff);
fail_buff_malloc:
    fprintf(stderr, "Failure\n");
    return status;
}


