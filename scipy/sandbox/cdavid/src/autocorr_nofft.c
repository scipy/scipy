/*
 * Last Change: Tue Nov 28 03:00 PM 2006 J
 * vim:syntax=c
 *
 * TODO: is size_t 64 bits long on 64 bits machines ?
 */
#include <stddef.h> /* for size_t */
#include <stdio.h> /* for size_t */

#include "autocorr_nofft.h"

/*
 * NOFFT auto correlation
 */


/* 
 * float version; out must have a size of lag+1, already pre allocated 
 *
 * lag should be < size
 *
 * returns 0 is succesfull, other value otherwise.
 */
int flt_xcorr_nofft_1d(const float *in, 
    const size_t size, float *out, const size_t lag)
{
	size_t	i, j;
	float acc;
	
	/* lag 0 */
	acc	= 0;
	for (i = 0; i <  size; ++i) {
		acc	+= in[i]*in[i];
	}
	out[0] = acc;

	/* lag : 1 -> lag */
	for (i = 1; i <= lag; ++i) {
		acc	= 0;
		for (j = i; j < size; ++j) {
			acc	+= 	in[j-i]*in[j];
		}
		out[i]	= acc;
	}

	return 0;
}

/* 
 * double version; out must have a size of lag+1, already pre allocated 
 *
 * lag should be < size
 *
 * returns 0 is succesfull, other value otherwise.
 */
int dbl_xcorr_nofft_1d(const double *in, 
    const size_t size, double *out, const size_t lag)
{
	size_t	i, j;
	double acc;
	
	/* lag 0 */
	acc	= 0;
	for (i = 0; i <  size; ++i) {
		acc	+= in[i]*in[i];
	}
	out[0] = acc;

	/* lag : 1 -> lag */
	for (i = 1; i <= lag; ++i) {
		acc	= 0;
		for (j = i; j < size; ++j) {
			acc	+= 	in[j-i]*in[j];
		}
		out[i]	= acc;
	}

	return 0;
}



/* 
 * float version for non contiguous arrays; the corresponding 
 * array should have at least in_size elements. 
 *
 * Constraints:
 *  - lag should be < in_size
 *  - strides in bytes
 *  - TODO: check if should be aligned ?
 * 
 * returns 0 is succesfull, other value otherwise.
 */
int flt_xcorr_nofft_1d_noncontiguous(const float *in, size_t in_size, 
    size_t in_stride, float *out, size_t out_stride, size_t lag)
{
	size_t	i, j, clag;
    size_t  istride = in_stride / sizeof(float);
    size_t  ostride = out_stride / sizeof(float);
	float acc;
	
	/* lag 0 */
	acc	= 0;
	for (i = 0; i < in_size * istride; i+= istride) {
		acc	+= in[i]*in[i];
	}
	out[0] = acc;

	/* lag : 1 -> lag */
	for (i = 1; i <= lag ; ++i) {
		acc	    = 0;
        clag    = i * istride;
        for (j = clag; j < in_size * istride; j += istride) {
			acc	+= 	in[j-clag]*in[j];
		}
		out[i * ostride]	= acc;
	}

    return 0;
}

/* 
 * double version for non contiguous arrays; the corresponding 
 * array should have at least in_size elements. 
 *
 * Constraints:
 *  - lag should be < in_size
 *  - strides in bytes
 *  - TODO: check if should be aligned ?
 * 
 * returns 0 is succesfull, other value otherwise.
 */
int dbl_xcorr_nofft_1d_noncontiguous(const double *in, size_t in_size, 
    size_t in_stride, double *out, size_t out_stride, size_t lag)
{
	size_t	i, j, clag;
    size_t  istride = in_stride / sizeof(double);
    size_t  ostride = out_stride / sizeof(double);
	double acc;
	
	/* lag 0 */
	acc	= 0;
	for (i = 0; i < in_size * istride; i+= istride) {
		acc	+= in[i]*in[i];
	}
	out[0] = acc;

	/* lag : 1 -> lag */
	for (i = 1; i <= lag ; ++i) {
		acc	    = 0;
        clag    = i * istride;
        for (j = clag; j < in_size * istride; j += istride) {
			acc	+= 	in[j-clag]*in[j];
		}
		out[i * ostride]	= acc;
	}

    return 0;
}



/* 
 * For rank 2 arrays, contiguous cases
 * float version; out must have a size of lag+1, already pre allocated 
 *
 * lag should be < size
 *
 * returns 0 is succesfull, other value otherwise.
 */
int flt_xcorr_nofft_2d(const float *in, 
    size_t dim0, size_t dim1, float *out, const size_t lag)
{
    size_t  i;
    float  *coaxis;

#if 0
    for(i = 0; i < dim0; ++i) {
        fprintf(stdout, "%d 1d autocorr, first element is %f\n", i, in[i * dim1]);
    }
#endif
    for(i = 0; i < dim0; ++i) {
        coaxis  = out + i * (lag + 1);
        flt_xcorr_nofft_1d(in + i * dim1, dim1, coaxis, lag);
    }

	return 0;
}

/* 
 * For rank 2 arrays, contiguous cases
 * double version; out must have a size of lag+1, already pre allocated 
 *
 * lag should be < size
 *
 * returns 0 is succesfull, other value otherwise.
 */
int dbl_xcorr_nofft_2d(const double *in, 
    size_t dim0, size_t dim1, double *out, const size_t lag)
{
    size_t  i;
    double  *coaxis;

#if 0
    for(i = 0; i < dim0; ++i) {
        fprintf(stdout, "%d 1d autocorr, first element is %f\n", i, in[i * dim1]);
    }
#endif
    for(i = 0; i < dim0; ++i) {
        coaxis  = out + i * (lag + 1);
        dbl_xcorr_nofft_1d(in + i * dim1, dim1, coaxis, lag);
    }

	return 0;
}



/* 
 * For rank 2 arrays, non contiguous cases
 * float version; out must have a size of lag+1, already pre allocated 
 *
 * lag should be < size
 *
 * returns 0 is succesfull, other value otherwise.
 */
int flt_xcorr_nofft_2d_noncontiguous(const float *in, 
    size_t dim0, size_t dim1, size_t in_stride0, size_t in_stride1, 
    float *out, size_t out_stride0, size_t out_stride1,
    const size_t lag)
{
    size_t  i;

    size_t  istride0    = in_stride0 / sizeof(float);
    size_t  ostride0    = out_stride0 / sizeof(float);

    float  *coaxis;
#if 0
    fprintf(stdout, "%s: shape is (%d, %d)\n", __func__, dim0, dim1);
    fprintf(stdout, "%s: istrides are (%d, %d)\n", __func__, istride0, istride1);
    
    fprintf(stdout, "%s: ostrides are (%d, %d)\n", __func__, ostride0, ostride1);
    for(i = 0; i < dim0; ++i) {
        ciaxis  = in + i * istride0;
        coaxis  = out + i * istride0;
        fprintf(stdout, "%d 1d autocorr, first element is %f, last is %f (%d el)\n", 
            i, ciaxis[0], ciaxis[(dim1-1) * istride1], dim1);
    }
#endif

    for(i = 0; i < dim0; ++i) {
        coaxis  = out + i * ostride0;
        flt_xcorr_nofft_1d_noncontiguous(in + i * istride0, dim1, in_stride1,
            coaxis, out_stride1, lag);
    }
	return 0;
}

/* 
 * For rank 2 arrays, non contiguous cases
 * double version; out must have a size of lag+1, already pre allocated 
 *
 * lag should be < size
 *
 * returns 0 is succesfull, other value otherwise.
 */
int dbl_xcorr_nofft_2d_noncontiguous(const double *in, 
    size_t dim0, size_t dim1, size_t in_stride0, size_t in_stride1, 
    double *out, size_t out_stride0, size_t out_stride1,
    const size_t lag)
{
    size_t  i;

    size_t  istride0    = in_stride0 / sizeof(double);
    size_t  ostride0    = out_stride0 / sizeof(double);

    double  *coaxis;
#if 0
    fprintf(stdout, "%s: shape is (%d, %d)\n", __func__, dim0, dim1);
    fprintf(stdout, "%s: istrides are (%d, %d)\n", __func__, istride0, istride1);
    
    fprintf(stdout, "%s: ostrides are (%d, %d)\n", __func__, ostride0, ostride1);
    for(i = 0; i < dim0; ++i) {
        ciaxis  = in + i * istride0;
        coaxis  = out + i * istride0;
        fprintf(stdout, "%d 1d autocorr, first element is %f, last is %f (%d el)\n", 
            i, ciaxis[0], ciaxis[(dim1-1) * istride1], dim1);
    }
#endif

    for(i = 0; i < dim0; ++i) {
        coaxis  = out + i * ostride0;
        dbl_xcorr_nofft_1d_noncontiguous(in + i * istride0, dim1, in_stride1,
            coaxis, out_stride1, lag);
    }
	return 0;
}


