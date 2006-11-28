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

int flt_levinson1d(const float* in, 
        size_t order, float* acoeff, 
        float* err, float* kcoeff, 
        float* tmp)
{
	/* TODO: to check if first element of corr is 0*/
	
	size_t	i, j;
	float acc;
    int     ret = 0;

	/* 
	 * order 0 
	 */
	acoeff[0]	= (float)1.0;
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

int flt_levinson1d_check(const float* in, 
        size_t order, float* acoeff, 
        float* err, float* kcoeff, 
        float* tmp)
{
	/* TODO: to check if first element of corr is 0*/
	
	size_t	i, j;
	float acc;
    int     ret = 0;

	/* 
	 * order 0 
	 */
	acoeff[0]	= (float)1.0;
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
 * float version; out must have a size of dim0 * (lag + 1), 
 * already pre allocated 
 *
 * order should be < dim1
 *
 * returns 0 is succesfull, other value otherwise.
 */
int flt_levinson2d(const float *in, 
    size_t dim0, size_t dim1, 
    float *acoeff, 
    float *err, 
    float *kcoeff)
{
    size_t  i;
    size_t  order = dim1 - 1;
    float  *caaxis, *ceaxis, *ckaxis, *buff;

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
        flt_levinson1d(in + i * dim1, order, caaxis, ceaxis, ckaxis, buff);
    }

    free(buff);
	return 0;

fail_malloc:
    return -1;
}

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

int dbl_levinson1d(const double* in, 
        size_t order, double* acoeff, 
        double* err, double* kcoeff, 
        double* tmp)
{
	/* TODO: to check if first element of corr is 0*/
	
	size_t	i, j;
	double acc;
    int     ret = 0;

	/* 
	 * order 0 
	 */
	acoeff[0]	= (double)1.0;
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

int dbl_levinson1d_check(const double* in, 
        size_t order, double* acoeff, 
        double* err, double* kcoeff, 
        double* tmp)
{
	/* TODO: to check if first element of corr is 0*/
	
	size_t	i, j;
	double acc;
    int     ret = 0;

	/* 
	 * order 0 
	 */
	acoeff[0]	= (double)1.0;
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
 * double version; out must have a size of dim0 * (lag + 1), 
 * already pre allocated 
 *
 * order should be < dim1
 *
 * returns 0 is succesfull, other value otherwise.
 */
int dbl_levinson2d(const double *in, 
    size_t dim0, size_t dim1, 
    double *acoeff, 
    double *err, 
    double *kcoeff)
{
    size_t  i;
    size_t  order = dim1 - 1;
    double  *caaxis, *ceaxis, *ckaxis, *buff;

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
        dbl_levinson1d(in + i * dim1, order, caaxis, ceaxis, ckaxis, buff);
    }

    free(buff);
	return 0;

fail_malloc:
    return -1;
}

