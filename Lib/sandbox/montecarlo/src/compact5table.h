/* Include file for discrete random variable generation with Marsaglia's compact
 * 5-table algorithm. 
 *
 * Copyright: Ed Schofield, 2005
 * License: BSD-style (see LICENSE.txt at root of scipy tree)
 */

#include <stdint.h>
#include "mtrand/randomkit.h"

typedef struct sampler_t
{
    int32_t t1,t2,t3,t4;           /* limits for table lookups */
    int32_t *AA,*BB,*CC,*DD,*EE;   /* tables for condensed table-lookup */
    int32_t sizeEE;
    long prob1event;     /* If this is >= 0, this sample point x=prob1event has
                          * p(x) = 1, so the lookup table sampling method won't
                          * work.  (We just let Dran() spit out this value
                          * instead.)  */
    rk_state state;                /* state variable for the RandomKit RNG. */
} Sampler;

/* Represent probabilities as 30-bit integers and create the 5 tables.
 * The prob mass function is specified as n values p(j) = weights[j]. */
Sampler* init_sampler5tbl(double* weights, long n, unsigned long seed);

/* Deallocate it */
void destroy_sampler5tbl(Sampler*);

/* Discrete random variable generating functions using 5 compact tables */
uint32_t Dran(Sampler* sampler);  /* for a single variate */
void Dran_array(Sampler*, uint32_t* output, long samplesize);  /* array version */



