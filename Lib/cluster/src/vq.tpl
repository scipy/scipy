[+ AutoGen5 template c +]
/*
 * vim:syntax=c
 *
 * This file implements vq for float and double in C. It is a direct
 * translation from the swig interface which could not be generated anymore
 * with recent swig
 */
#include <stddef.h>
#include <math.h>

#include "vq.h"
/*
 * results is put into code, which contains initially the initial code
 *
 * mdist and code should have at least n elements
 */
const static double rbig = 1e100;

[+ FOR data_type +]
#if 0
static int [+ (get "type_name") +]_vq_1d(const [+ (get "type_name") +] *in, int n, 
    const [+ (get "type_name") +] *init, int ncode, 
    npy_intp *code, [+ (get "type_name") +] *mdist)
{
    int i, j;
    [+ (get "data_type") +] m, d;

    for (i = 0; i < n; ++i) {
        m = ([+ (get "data_type") +])rbig;
        /* Compute the minimal distance for obsvervation i */
        for (j = 0; j < ncode; ++j) {
            d = (in[i] - init[j]);
            d *= d;
            if ( d < m) {
                m = d;
            }
        }
        mdist[i] = m;
        code[i] = j;
    }
    return 0;
}
#endif

static int [+ (get "type_name") +]_vq_obs(const [+ (get "data_type") +] *obs,
    [+ (get "data_type") +] *code_book, int Ncodes, int Nfeatures,
       npy_intp* code, [+ (get "data_type") +] *lowest_dist)
{
	int i,j,k=0;
	[+ (get "data_type") +] dist, diff;

	*lowest_dist = ([+ (get "data_type") +]) rbig;
	for(i = 0; i < Ncodes; i++) {
		dist = 0;
		for(j=0; j < Nfeatures; j++) {
			diff = code_book[k] - obs[j];
			dist += diff*diff;
			k++;
		}
		dist = ([+ (get "data_type") +])sqrt(dist);
		if (dist < *lowest_dist) {
			*code = i;
			*lowest_dist = dist;
		}
	}

    return 0;
}

int [+ (get "type_name") +]_tvq(
    [+ (get "data_type") +]* obs,
    [+ (get "data_type") +]* code_book, 
    int Nobs, int Ncodes, int Nfeatures,
    npy_intp* codes, [+ (get "data_type") +]* lowest_dist)
{
    int i;
	for( i = 0; i < Nobs; i++) {		
		[+ (get "type_name") +]_vq_obs(
                    &(obs[i*Nfeatures]),
                    code_book,Ncodes, Nfeatures,
                    &(codes[i]), &(lowest_dist[i]));
	}
    return 0;
}
[+ ENDFOR data_type +]
