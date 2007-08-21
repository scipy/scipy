/*
 * Last Change: Mon Nov 27 07:00 PM 2006 J
 */

#ifndef _C_AUTOCORR_NOFFT_H_
#define _C_AUTOCORR_NOFFT_H_

/*
 * Direct implementation of auto correlation (faster when a few lags only are
 * necessary). One side only, out should have at least lag+1 elements allocated
 *
 * Expect in and out to be contiguous
 */
int flt_xcorr_nofft_1d(const float *in, size_t size, float *out, size_t lag);
int dbl_xcorr_nofft_1d(const double *in, size_t size, double *out, size_t lag);

/*
 * Direct implementation of auto correlation (faster when a few lags only are
 * necessary). One side only, out should have at least lag+1 elements allocated
 *
 * Expect in and out need not to be contiguous
 */
int flt_xcorr_nofft_1d_noncontiguous(const float *in, size_t in_size, size_t in_stride, 
        float *out, size_t out_stride, size_t lag);
int dbl_xcorr_nofft_1d_noncontiguous(const double *in, size_t in_size, size_t in_stride, 
        double *out, size_t out_stride, size_t lag);

/*
 * 1d autocorrelation for rank 2 arrays
 */
int flt_xcorr_nofft_2d(const float *in, size_t dim0, size_t dim1, 
        float *out, size_t lag);
int dbl_xcorr_nofft_2d(const double *in, size_t dim0, size_t dim1, 
        double *out, size_t lag);

/*
 * 1d autocorrelation for rank 2 arrays, non contiguous cases
 */
int dbl_xcorr_nofft_2d_noncontiguous(const double *in,
    size_t dim0, size_t dim1, size_t in_stride0, size_t in_stride1, 
    double *out, size_t out_stride0, size_t out_stride1,
    const size_t lag);
int flt_xcorr_nofft_2d_noncontiguous(const float *in,
    size_t dim0, size_t dim1, size_t in_stride0, size_t in_stride1, 
    float *out, size_t out_stride0, size_t out_stride1,
    const size_t lag);
#endif
