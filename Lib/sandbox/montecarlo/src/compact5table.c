/* Source file for fast variate generation from a discrete distribution.  Uses
 * Marsaglia's compact 5-table method, as described in his paper 'Fast
 * generation of discrete random variables' in the Journal of Statistical
 * Software, July 2004, vol 11, issue 3.  This code is based upon the C
 * implementation that accompanies that paper, but is simpler, and uses a
 * different random number generator, the Mersenne Twister in Jean-Sebastien
 * Roy's RandomKit.
 *
 * Copyright: Ed Schofield, 2005
 * License: BSD-style (see LICENSE.txt at root of scipy tree)
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "compact5table.h"
#include "mtrand/randomkit.h"

#define dg(m,k) ((m>>(30-6*k))&63)  /* gets kth digit of m (base 64) */
#define MAX(a,b)  ((a)>(b)?(a):(b))

/* Represent probabilities as 30-bit integers and create the 5 tables.
   We assume the n weights are all non-negative and that their sum is strictly
   positive. */
Sampler* init_sampler5tbl(double* weights, long n, unsigned long seed)
{
    Sampler* sampler = NULL;
    int32_t *AA,*BB,*CC,*DD,*EE;      /* Tables for condensed table-lookup */

    int32_t *P;         /* probabilities as an array of 30-bit integers */
    int32_t sizeAA=0, sizeBB=0, sizeCC=0, sizeDD=0, sizeEE=0;
    long i;             /* the number of table elements */
    int32_t m, k, j;    /* the probs are stored as 30-bit integers in a 32-bit
                           data type */
    double sum = 0.0;
    int seedhigh, seedlow;

    uint32_t random_value;
    double random_double;

    if(sum <= 0)
    {
        fprintf(stderr, "Error: invalid arguments to init_sampler5tbl()." \
                "The sum of the probabilities is non-positive.  Aborting!\n");
        return NULL;
    }

    /* Allocate space */
    P = malloc(n*sizeof(int32_t));
    
    /* fill P array, (30-bit integers) */
    for (i=0; i<n; i++)
        P[i] = weights[i]/sum * (1<<30);
    
    sampler = malloc(sizeof(Sampler));
    /*  Initialize prob1event to -1, indicating that there's no event with
     *  prob=1 (in which case the lookup table sampling method wouldn't work)
     */
    sampler->prob1event = -1;

    /* Now seed the RandomKit RNG */
    rk_seed(seed, &sampler->state);
    
    /* Normalize weights */
    for (i=0; i<n; i++)
        sum += weights[i];
    
    /* Now create the 5 tables */ 
    /* Get table sizes, then malloc */
    for(i=0; i<n; i++)
    {
        m=P[i];
        sizeAA+=dg(m,1);
        sizeBB+=dg(m,2);
        sizeCC+=dg(m,3);
        sizeDD+=dg(m,4);
        sizeEE+=dg(m,5);
    }
    
    if (sizeAA+sizeBB+sizeCC+sizeDD+sizeEE <= 0)
    {
        /* Either an event has probability one, or all events have probability
         * zero.  If an event has prob 1, find it! */
        for (i=0; i<n; i++)
        {
            m = P[i];
            if (((m>>30) & 63) == 1)
            {
                /* Found event with probability 1 */
                sampler->prob1event = i;
                break;
            }
        }
        // if (sampler->prob1event == -1)
        // {
        //     /* This has been tested above, so I think this should never
        //        occur...  */
        //     fprintf(stderr, "Error: invalid arguments to init_sampler5tbl()." \
        //             "The sum of the probabilities is zero.  Aborting!\n");
        //     free(P);
        //     free(sampler);
        //     return NULL;
        // }
    }
    else
    {
        sampler->AA = AA = malloc(sizeAA*sizeof(int32_t)); 
        sampler->BB = BB = malloc(sizeBB*sizeof(int32_t));
        sampler->CC = CC = malloc(sizeCC*sizeof(int32_t));
        sampler->DD = DD = malloc(sizeDD*sizeof(int32_t));
        sampler->EE = EE = malloc(sizeEE*sizeof(int32_t));
        sampler->sizeEE = sizeEE;

        sampler->t1 = sizeAA<<24;
        sampler->t2 = sampler->t1+(sizeBB<<18);
        sampler->t3 = sampler->t2+(sizeCC<<12);
        sampler->t4 = sampler->t3+(sizeDD<<6);
    
        int32_t na=0, nb=0, nc=0, nd=0, ne=0;
        
        /* Fill tables AA,BB,CC,DD,EE */
        for(i=0; i<n; i++)
        {
            m=P[i]; k=i; // +offset;
            for(j=0;j<dg(m,1);j++) 
                AA[na+j]=k;
            na+=dg(m,1);

            for(j=0;j<dg(m,2);j++)
                BB[nb+j]=k;
            nb+=dg(m,2);

            for(j=0;j<dg(m,3);j++)
                CC[nc+j]=k;
            nc+=dg(m,3);

            for(j=0;j<dg(m,4);j++)
                DD[nd+j]=k;
            nd+=dg(m,4);

            for(j=0;j<dg(m,5);j++)
                EE[ne+j]=k;
            ne+=dg(m,5);
        }
    }
    free(P);
    return sampler;
}


/* Destroys the sampler created by start() and frees memory */
void destroy_sampler5tbl(Sampler* sampler)
{
    if (sampler->prob1event == -1  /* if there is no prob1event */ )
    {
        free(sampler->AA);  sampler->AA = NULL;
        free(sampler->BB);  sampler->BB = NULL;
        free(sampler->CC);  sampler->CC = NULL;
        free(sampler->DD);  sampler->DD = NULL;
        free(sampler->EE);  sampler->EE = NULL;
    }
    free(sampler);
}


/* Discrete random variable generating function. Uses the provided compact
 * 5-table sampler. */
uint32_t Dran(Sampler* sampler)
{
    uint32_t j;
    if (sampler->prob1event > -1)
    {
        /* One event x has probability p(x) = 1.  Return it. */
        return sampler->prob1event;
    }

    /* Generate a 30-bit random number: */
    j = rk_random(&sampler->state) >> 2;
    
    if(j < sampler->t1) return sampler->AA[j>>24];
    if(j < sampler->t2) return sampler->BB[(j-sampler->t1)>>18];
    if(j < sampler->t3) return sampler->CC[(j-sampler->t2)>>12];
    if(j < sampler->t4) return sampler->DD[(j-sampler->t3)>>6];
    
    /* It seems we need to deal specially with the boundary cases where
     *     j >= 2^30-1 == 1073741823
     */
    if (j - sampler->t4 >= sampler->sizeEE)
    {
        /* The random number generated is larger than the sizes of all tables.
         * This should happen only very rarely. For now, just generate another
         * random number. */
        fprintf(stderr, 
            "Debug: random number is larger than the sizes of all tables!");
        return Dran(sampler);
    }
    else
        return sampler->EE[j - sampler->t4];
}


/* Array version of discrete random variable generating function using the
 * compact 5-table representation. */
void Dran_array(Sampler* sampler, uint32_t* output, long samplesize)
{
    uint32_t j;
    long i;
    
    if (sampler->prob1event > -1)
    {
        /* One event x has probability p(x) = 1.  Spit it out! */
        for (i=0; i<samplesize; i++)
            output[i] = sampler->prob1event;
        return;
    }
    for (i=0; i<samplesize; i++)
    {
        /* Generate a 30-bit random number: */
        j = rk_random(&sampler->state) >> 2;
    
        if(j < sampler->t1) {
            output[i] = sampler->AA[j>>24];
        }
        else if(j < sampler->t2) {
            output[i] = sampler->BB[(j-sampler->t1)>>18];
        }
        else if(j < sampler->t3) {
            output[i] = sampler->CC[(j-sampler->t2)>>12];
        }
        else if(j < sampler->t4) {
            output[i] = sampler->DD[(j-sampler->t3)>>6];
        }
        /* It seems we need to deal specially with the boundary cases where
         *     j >= 2^30-1 == 1073741823
         */
        else if (j - sampler->t4 >= sampler->sizeEE)
        {
            /* The random number generated is larger than the sizes of all
             * tables.  This should happen only very rarely. For now, just
             * generate another random number.
             */
            fprintf(stderr,
                "Debug: random number is larger than the sizes of all tables!");
            i--;
        }
        else
            output[i] = sampler->EE[j - sampler->t4];
    }
}


