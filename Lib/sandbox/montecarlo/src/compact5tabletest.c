/* Copyright: Ed Schofield (2005-6)
 *
 * A test of using Marsaglia's compact 5-table method to generate N=10^8 variates
 * from the distribution:
 *     x        0       1         2         3         4 
 *    p(x)   12/100   58/100      0       15/100    15/100
 */



#include <stdio.h>
#include <malloc.h>
#include "compact5table.h"

int main()
{
    Sampler *mysampler;
    unsigned long* sample_array;
    
    long i, j=0,nsmpls=10000000;
    long size = 5;
    double weights[] = {12.0, 58.0, 0.0, 15.0, 15.0};
    long counts[] = {0,0,0,0,0};
    long counts2[] = {0,0,0,0,0};
    unsigned long seed = 100;

    /* Init */
    mysampler = init_sampler5tbl(weights, size, seed);
    
    /* Sample and count # occurrences */
    for (i=0; i<nsmpls; i++)
    {
        j=Dran(mysampler);
        counts[j]++;
        //printf("Sample %d is: %d\n", i, j);
    }
    for (j=0; j<size; j++)
        printf("%d occurrences of %d\n", counts[j], j);
    
    printf("***********************************************\n");
    /* Now repeat, with the array version */
    /* Sample */
    sample_array = (unsigned long*) malloc(nsmpls*sizeof(long));
    Dran_array(mysampler, &sample_array[0], nsmpls);

    /* Count # occurrences */
    for (i=0; i<nsmpls; i++)
    {
        j=sample_array[i];
        counts2[j]++;
    }
    for (j=0; j<size; j++)
        printf("%d occurrences of %d\n", counts2[j], j);
    
    //printf("Max rand # is: %d\n", maxj);
    
    /* Destroy */
    free(sample_array);
    destroy_sampler5tbl(mysampler);
    
    /* Done */
    return 0;
}

