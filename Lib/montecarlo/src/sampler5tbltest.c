/*
 * We use the 5 tables method of Marsaglia (2004) to generate N=10^8 variates
 * from the distribution:
 *     x        0       1         2         3         4 
 *    p(x)   10/180   140/180     0       15/180    15/180
 */



#include <stdio.h>
#include <malloc.h>
#include "sampler5tbl.h"

int main()
{
    Sampler *mysampler;
    int* sample_array;
    
    int i, j=0,nsmpls=10000000;
    int size = 5;
    double weights[] = {10.0, 140.0, 0.0, 15.0, 15.0};
    int counts[] = {0,0,0,0,0};
    int counts2[] = {0,0,0,0,0};
    
    /* Init */
    mysampler = init(weights, size);
    
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
    sample_array = (int*) malloc(nsmpls*sizeof(int));
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
    destroy(mysampler);
    
    /* Die */
    return 0;
}

