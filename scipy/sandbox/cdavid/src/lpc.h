/*
 * Last Change: Tue Nov 28 12:00 PM 2006 J
 */

/*
 * Functions to compute lpc coefficients with contiguous arrays.
 *
 * input is signal, output are coeff, kcoeff and err.
 *
 * requirements:
 *  - signal must have size elements
 *  - order < size
 *  - coeff must have order + 1 elements at least
 *  - kcoeff must have order elements at least
 *  - err must have at least one element
 */
int dbl_lpc(const double* signal, size_t size, size_t order, double* coeff, 
        double* kcoeff, double* err);
int flt_lpc(const float* signal, size_t size, size_t order, float* coeff, 
        float* kcoeff, float* err);
