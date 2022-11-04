#ifndef NI_SPLINES_H
#define NI_SPLINES_H


/*
 * For context on the theory behind this code, see:
 *
 * Unser, Michael. "Splines: A perfect fit for signal and image
 * processing." IEEE Signal processing magazine 16.6 (1999): 22-38.
 */

/*
 * Fills `weights` with the interpolation weights for the B-splines at
 * the `origin + 1` knots closest to the point with coordinate `x`.
 * Because B-splines are symmetrical, the above turns out to be the
 * same as computing the values of a B-spline centered at `x` at the
 * knots.
 */
int
get_spline_interpolation_weights(double x, int order, double *weights);


/* The maximum number of poles for supported spline orders. */
#define MAX_SPLINE_FILTER_POLES 2


/*
 * Stores in `poles` the values of the poles of a spline filter of the
 * selected order, and in `npoles` the number of poles written. `poles`
 * should be at least `MAX_SPLINE_FILTER_POLES` long. Returns 1 if the
 * order is not supported, 0 on success.
 */
int
get_filter_poles(int order, int *npoles, double *poles);


/*
 * Applies the causal and anticausal filters for all poles to the array
 * of coefficients. This is equivalent to solving the banded linear
 * system of equations that computes the spline coefficients for an
 * input array.
 */
void
apply_filter(double *coefficients, const npy_intp len, const double *poles,
             int npoles, NI_ExtendMode mode);


#endif
