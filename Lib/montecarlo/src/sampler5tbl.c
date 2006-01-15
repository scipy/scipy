/* Source file for fast variate generation from a discrete distribution.  Uses
 * Marsaglia's compact 5-table method.  Adapted from Marsaglia's paper 'Fast
 * generation of discrete random variables' in the Journal of Statistical
 * Software, July 2004, vol 11, issue 3.
 *
 * Test this by linking with sampler5tbltest.c
 *
 * Copyright: Ed Schofield, 2005
 * License: BSD-style (see LICENSE.txt at root of scipy tree)
 */


#include <stdio.h>
//#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include "sampler5tbl.h"

#define dg(m,k) ((m>>(30-6*k))&63)  /* gets kth digit of m (base 64) */
#define MAX(a,b)  ((a)>(b)?(a):(b))

// static int t1,t2,t3,t4,offset=0,last;  /* limits for table lookups */
static unsigned long jxr=182736531; /* Xorshift RNG */
//static unsigned long maxj=0;
// int prob1event = -1;                 /* If this is >= 0, this sample point x=prob1event has p(x) = 1, so the lookup table sampling method won't work.  (We just spit out this value in Dran() instead.) */

/* Represent probabilities as 30-bit integers and create the 5 tables.
   Assumes the weights are all non-negative and that their sum is strictly
   positive. */
Sampler* init_sampler5tbl(double* weights, long n)
{
    Sampler* sampler = NULL;
    int32_t *AA,*BB,*CC,*DD,*EE;      /* Tables for condensed table-lookup */

    int32_t *P;   /* Probabilities as an array of 30-bit integers*/
    int32_t sizeAA=0, sizeBB=0, sizeCC=0, sizeDD=0, sizeEE=0;
    long i;    /* the number of table elements has type 'long' */ 
    int32_t m, k, j;   /* the probs are stored as 30-bit integers in a 32-bit
                        data type */
    double sum = 0.0;
    
    /* Seed our faster XORshift RNG */
    srand48( (unsigned)time( NULL ) );
    jxr = lrand48()>>2;
    //printf("jxr initially set to %u\n",jxr);
    //jxr = 183025040;  this initial value leads to a SEGV with jxr==2^30-1
    //printf("Random seed is: %ld\n", jxr);
    
    /* Normalize weights */
    for (i=0; i<n; i++)
        sum += weights[i];
    
    assert(sum > 0);

    /* Allocate space */
    P = malloc(n*sizeof(int32_t));
    
    /* fill P array, (30-bit integers) */
    for (i=0; i<n; i++)
        P[i] = weights[i]/sum * (1<<30);
    
    /*P[0] = 1./1800. * (1<<30);
    P[1] = 1699./1800. * (1<<30);
    P[2] = 100./1800. * (1<<30); */

    sampler = malloc(sizeof(Sampler));
    //  Initialize prob1event to -1, indicating that there's no event with prob=1
    //  (in which case the lookup table sampling method wouldn't work)
    sampler->prob1event = -1;

    /* Now create the 5 tables */ 
    /* get table sizes, malloc */
    for(i=0; i<n; i++)
    {
        m=P[i];
        sizeAA+=dg(m,1);sizeBB+=dg(m,2);sizeCC+=dg(m,3);sizeDD+=dg(m,4);sizeEE+=dg(m,5);
    }
    
    if (sizeAA+sizeBB+sizeCC+sizeDD+sizeEE <= 0)
    {
        // Either an event has probability one, or all events have probability zero
        // If an event has prob 1, find it!
        for (i=0; i<n; i++)
        {
            //printf("in loop\n");
            m = P[i];
            //printf("((m>>30) & 63) == %d\n", (m>>30) & 63);
            if (((m>>30) & 63) == 1)
            {
                // printf("(found event with probability 1 ...)\n");
                sampler->prob1event = i;
                break;
            }
        }
        if (sampler->prob1event == -1)
        {
            fprintf(stderr, "Error: invalid arguments to init_sampler5tbl()." \
                    "The sum of the probabilities is zero.  Aborting!\n");
            free(P);
            free(sampler);
            return NULL;
        }
    }
    else
    {
        sampler->AA = AA = malloc(sizeAA*sizeof(int32_t)); // should be 32 !!
        sampler->BB = BB = malloc(sizeBB*sizeof(int32_t));
        sampler->CC = CC = malloc(sizeCC*sizeof(int32_t));
        sampler->DD = DD = malloc(sizeDD*sizeof(int32_t));
        sampler->EE = EE = malloc(sizeEE*sizeof(int32_t));
        sampler->sizeEE = sizeEE;
        //printf(" Table sizes: %d, %d, %d, %d, %d, total=%d\n",sizeAA,sizeBB,sizeCC,sizeDD,sizeEE,sizeAA+sizeBB+sizeCC+sizeDD+sizeEE);
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
        free(sampler->AA);  sampler->AA = NULL;  //sizeAA = 0;
        free(sampler->BB);  sampler->BB = NULL;  //sizeBB = 0;
        free(sampler->CC);  sampler->CC = NULL;  //sizeCC = 0;
        free(sampler->DD);  sampler->DD = NULL;  //sizeDD = 0;
        free(sampler->EE);  sampler->EE = NULL;  //sizeEE = 0;
    }
    free(sampler);
}

/* Discrete random variable generating function */
/* Uses 5 compact tables */

long Dran(Sampler* sampler)
{
    unsigned long j;
    if (sampler->prob1event > -1)
    {
        /* One event x has probability p(x) = 1.  Spit it out! */
        return sampler->prob1event;
    }
    jxr^=jxr<<13;
    jxr^=jxr>>17;
    jxr^=jxr<<5;
    j=(jxr>>2);
    
    // The C way would be this:
    // j = lrand48()>>1;
    
    if(j < sampler->t1) return sampler->AA[j>>24];
    if(j < sampler->t2) return sampler->BB[(j-sampler->t1)>>18];
    if(j < sampler->t3) return sampler->CC[(j-sampler->t2)>>12];
    if(j < sampler->t4) return sampler->DD[(j-sampler->t3)>>6];
    // It seems we need to deal specially with the boundary case where
    //     jxr == 2^30-1 == 1073741823
    //if (jxr >= 1073741823)
    if (j - sampler->t4 >= sampler->sizeEE)
    {
        // The random number generated is larger than the sizes of all tables.
        // This should happen only very rarely. For now, just generate another
        // random number.
        
        return Dran(sampler);
    }
    else
        return sampler->EE[j - sampler->t4];
}

/* Array version of discrete random variable generating function */
/* Uses 5 compact tables */
void Dran_array(Sampler* sampler, long* output, long samplesize)
{
    unsigned long j;
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
        jxr^=jxr<<13;
        jxr^=jxr>>17;
        jxr^=jxr<<5;
        j=(jxr>>2);
        
        // The C way would be this:
        // j = lrand48()>>1;
        
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
        // It seems we need to deal specially with the boundary case where
        //     jxr == 2^30-1 == 1073741823
        //if (jxr >= 1073741823)
        else if (j - sampler->t4 >= sampler->sizeEE)
        {
            // The random number generated is larger than the sizes of all tables.
            // This should happen only very rarely. For now, just generate another
            // random number.
            i--;
        }
        else
            output[i] = sampler->EE[j - sampler->t4];
    }
}


