[+ AutoGen5 template c +]
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

[+ For float_type +]
/* 
 * [+ (get "type_name") +] version; out must have a size of lag+1, already pre allocated 
 *
 * lag should be < size
 *
 * returns 0 is succesfull, other value otherwise.
 */
int [+ (get "short_name") +]_xcorr_nofft_1d(const [+ (get "type_name") +] *in, 
    const size_t size, [+ (get "type_name") +] *out, const size_t lag)
{
	size_t	i, j;
	[+ (get "type_name") +] acc;
	
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
[+ ENDFOR float_type +]

[+ For float_type +]
/* 
 * [+ (get "type_name") +] version for non contiguous arrays; the corresponding 
 * array should have at least in_size elements. 
 *
 * Constraints:
 *  - lag should be < in_size
 *  - strides in bytes
 *  - TODO: check if should be aligned ?
 * 
 * returns 0 is succesfull, other value otherwise.
 */
int [+ (get "short_name") +]_xcorr_nofft_1d_noncontiguous(const [+ (get "type_name") +] *in, size_t in_size, 
    size_t in_stride, [+ (get "type_name") +] *out, size_t out_stride, size_t lag)
{
	size_t	i, j, clag;
    size_t  istride = in_stride / sizeof([+ (get "type_name") +]);
    size_t  ostride = out_stride / sizeof([+ (get "type_name") +]);
	[+ (get "type_name") +] acc;
	
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
[+ ENDFOR float_type +]

[+ For float_type +]
/* 
 * For rank 2 arrays, contiguous cases
 * [+ (get "type_name") +] version; out must have a size of lag+1, already pre allocated 
 *
 * lag should be < size
 *
 * returns 0 is succesfull, other value otherwise.
 */
int [+ (get "short_name") +]_xcorr_nofft_2d(const [+ (get "type_name") +] *in, 
    size_t dim0, size_t dim1, [+ (get "type_name") +] *out, const size_t lag)
{
    size_t  i;
    [+ (get "type_name") +]  *coaxis;

#if 0
    for(i = 0; i < dim0; ++i) {
        fprintf(stdout, "%d 1d autocorr, first element is %f\n", i, in[i * dim1]);
    }
#endif
    for(i = 0; i < dim0; ++i) {
        coaxis  = out + i * (lag + 1);
        [+ (get "short_name") +]_xcorr_nofft_1d(in + i * dim1, dim1, coaxis, lag);
    }

	return 0;
}
[+ ENDFOR float_type +]

[+ For float_type +]
/* 
 * For rank 2 arrays, non contiguous cases
 * [+ (get "type_name") +] version; out must have a size of lag+1, already pre allocated 
 *
 * lag should be < size
 *
 * returns 0 is succesfull, other value otherwise.
 */
int [+ (get "short_name") +]_xcorr_nofft_2d_noncontiguous(const [+ (get "type_name") +] *in, 
    size_t dim0, size_t dim1, size_t in_stride0, size_t in_stride1, 
    [+ (get "type_name") +] *out, size_t out_stride0, size_t out_stride1,
    const size_t lag)
{
    size_t  i;

    size_t  istride0    = in_stride0 / sizeof([+ (get "type_name") +]);
    size_t  ostride0    = out_stride0 / sizeof([+ (get "type_name") +]);

    [+ (get "type_name") +]  *coaxis;
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
        [+ (get "short_name") +]_xcorr_nofft_1d_noncontiguous(in + i * istride0, dim1, in_stride1,
            coaxis, out_stride1, lag);
    }
	return 0;
}
[+ ENDFOR float_type +]

