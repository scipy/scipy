[+ AutoGen5 template c +]
/*
 * Last Change: Tue Nov 28 04:00 PM 2006 J
 * 
 * vim:syntax=c
 */
#include <float.h>
#include <math.h> /* for isfinite */
#include <stdio.h>
#include <stdlib.h>

#include "common.h"

#include "levinson.h"

[+ For float_type +]
/* The actual computation :
 *	- in	: the input vector which defines the toeplitz matrix
 *	- size	: size of in (ie number of elements)
 *	- order	: size of the system to solve. order must be < size -1
 *	- acoeff: solution (ie ar coefficients). Size must be at last order+1
 *	- err	: *prediction* error (scalar)
 *	- kcoeff: reflexion coefficients. Size must be at last equal to equal to order.
 *	- tmp   : cache, mnust have at least order elements  
 */
 
/* 
 * this function assume all arrays are allocated with the 
 * right size, and that the parameters make sense. No checking
 * is done, must be done before calling this function 
 *
 * Returns 0 on success, -1 if a compuation error happened (overflow, underflow
 * for error calculation)
 */

int [+ (get "short_name") +]_levinson1d(const [+ (get "type_name") +]* in, 
        size_t order, [+ (get "type_name") +]* acoeff, 
        [+ (get "type_name") +]* err, [+ (get "type_name") +]* kcoeff, 
        [+ (get "type_name") +]* tmp)
{
	/* TODO: to check if first element of corr is 0*/
	
	size_t	i, j;
	[+ (get "type_name") +] acc;
    int     ret = 0;

	/* 
	 * order 0 
	 */
	acoeff[0]	= ([+ (get "type_name") +])1.0;
	*err		= in[0];

	/* 
	 * order >= 1
	 */
	for ( i = 1; i <= order; ++i) {
		acc	= in[i];
		for ( j = 1; j <= i-1; ++j) {
			acc	+= acoeff[j]*in[i-j];
		}
		kcoeff[i-1]	= -acc/(*err);
		acoeff[i]	= kcoeff[i-1];

		for ( j = 0; j < order; ++j) {
			tmp[j]	= acoeff[j];
		}

		for (j = 1; j < i;  ++j) {
			acoeff[j]	+= kcoeff[i-1]*tmp[i-j];
		}
		*err	*= (1-kcoeff[i-1]*kcoeff[i-1]); 
	}
	
	return ret;
}

int [+ (get "short_name") +]_levinson1d_check(const [+ (get "type_name") +]* in, 
        size_t order, [+ (get "type_name") +]* acoeff, 
        [+ (get "type_name") +]* err, [+ (get "type_name") +]* kcoeff, 
        [+ (get "type_name") +]* tmp)
{
	/* TODO: to check if first element of corr is 0*/
	
	size_t	i, j;
	[+ (get "type_name") +] acc;
    int     ret = 0;

	/* 
	 * order 0 
	 */
	acoeff[0]	= ([+ (get "type_name") +])1.0;
	*err		= in[0];

	/* 
	 * order >= 1
	 */
	for ( i = 1; i <= order; ++i) {
		acc	= in[i];
		for ( j = 1; j <= i-1; ++j) {
			acc	+= acoeff[j]*in[i-j];
		}
		kcoeff[i-1]	= -acc/(*err);
        if (!isfinite(kcoeff[i-1])) {
            fprintf(stderr, "%s:%s, kcoeff is not finite, err is %e\n", 
                __FILE__, __func__, *err);
            ret = -1;
        }
		acoeff[i]	= kcoeff[i-1];

		for ( j = 0; j < order; ++j) {
			tmp[j]	= acoeff[j];
		}

		for (j = 1; j < i;  ++j) {
			acoeff[j]	+= kcoeff[i-1]*tmp[i-j];
		}
		*err	*= (1-kcoeff[i-1]*kcoeff[i-1]); 
	}
	
	return ret;
}

/* 
 * For rank 2 arrays, contiguous cases
 * [+ (get "type_name") +] version; out must have a size of dim0 * (lag + 1), 
 * already pre allocated 
 *
 * order should be < dim1
 *
 * returns 0 is succesfull, other value otherwise.
 */
int [+ (get "short_name") +]_levinson2d(const [+ (get "type_name") +] *in, 
    size_t dim0, size_t dim1, 
    [+ (get "type_name") +] *acoeff, 
    [+ (get "type_name") +] *err, 
    [+ (get "type_name") +] *kcoeff)
{
    size_t  i;
    size_t  order = dim1 - 1;
    [+ (get "type_name") +]  *caaxis, *ceaxis, *ckaxis, *buff;

    buff    = malloc(sizeof(*buff) * order);
    if (buff == NULL) {
        goto fail_malloc;
    }

#if 0
    for(i = 0; i < dim0; ++i) {
        fprintf(stdout, "%d 1d levinson, first element is %f\n", i, in[i * dim1]);
    }
#endif
    for(i = 0; i < dim0; ++i) {
        caaxis  = acoeff + i * (order + 1);
        ckaxis  = kcoeff + i * (order);
        ceaxis  = err + i;
        [+ (get "short_name") +]_levinson1d(in + i * dim1, order, caaxis, ceaxis, ckaxis, buff);
    }

    free(buff);
	return 0;

fail_malloc:
    return -1;
}
[+ ENDFOR float_type +]
