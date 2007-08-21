[+ AutoGen5 template c +]
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

[+ For float_type +]
/* 
 * ([+ (get "type_name") +] version) Compute lpc coeff (order x coefficients) of a 
 * contiguous array 
 *
 * err is a double, coeff must be able to contain order+1 elements, and kcoeff
 * order elements
 */
int [+ (get "short_name") +]_lpc(const [+ (get "type_name") +]* signal, 
        size_t size, size_t order, [+ (get "type_name") +]* coeff, 
        [+ (get "type_name") +]* kcoeff, [+ (get "type_name") +]* err)
{
    size_t  i, nbuff, ncache;
    [+ (get "type_name") +] *buff, *cache, biasnorm;
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
    [+ (get "short_name") +]_xcorr_nofft_1d(signal, size, buff, order);
    for(i = 0; i < nbuff; ++i) {
        buff[i] *= biasnorm;
    }

    /*
     * Compute the inverse coefficients using (simple) levinson recursive algo
     */
    status  = [+ (get "short_name") +]_levinson1d(buff, order, 
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

[+ ENDFOR float_type +]
