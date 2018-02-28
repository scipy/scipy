#ifndef NI_SPLINES_H
#define NI_SPLINES_H


/*
 * For context in the theory behind this code, see:
 *
 * M. Unser, A. Aldroubi and M. Eden, "Fast B-spline transforms for
 * continuous image representation and interpolation," in IEEE
 * Transactions on Pattern Analysis and Machine Intelligence, vol. 13,
 * no. 3, pp. 277-285, Mar 1991.
 *
 * M. Unser, "Splines: A Perfect Fit for Signal and Image Processing,"
 * in IEEE Signal Processing Magazine, vol. 16, no. 6, pp. 22-38, 1999.
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
 * Computes the filter gain for the given poles.
 */
double
filter_gain(const double *poles, int npoles);


/*
 * Multiplies all coefficients in-place by the gain.
 */
void
apply_gain(double gain, double *coefficeints, npy_intp len);


/*
 * Sets the first coefficient to the right value to begin applying a
 * causal filter for the given pole to the sequence. If tolerance is
 * zero the calculation will be exact, if non-zero, the series will be
 * truncated for items  smaller than the tolerance.
 *
 * -----
 *
 * If s[i] are the values to filter, and z is the pole being considered,
 * the initial coefficient for the causal filter can be computed as:
 *
 *   c[0] = Sum(i=0..inf, s[i] * z**i)
 *
 * Since |z| < 1, this can be computed approximately by adding terms
 * until z**i is below a user defined tolerance.
 *
 * It can also be computed exactly by taking into account the boundary
 * conditions coming from the extend mode. Currently NI_EXTEND_MIRROR
 * is the only supported mode.
 *
 * NI_EXTEND_MIRROR : a b c d | c b a b c d c b ...
 * The s[i], when extended beyond the given n items, are periodic with
 * period T, so we can rewrite:
 *
 *   c[0] = Sum(k=0..inf, z**(T * k)) * Sum(i=0..T, s[i] * z**i) =
 *        = 1 / (1 - z**T) * S
 *
 * In this case T = 2*n - 2, and c[-k] = c[k], so S can be computed as:
 *
 *   S = s[0] + s[n-1] * z**(n-1) +
 *       Sum(i=1..n-1, s[i]* (z**i + z**(2*n - 2 - i)))
 *
 */
void
set_initial_causal_coefficient(double *coefficients, npy_intp len,
                               double pole, double tolerance);


/*
 * Sets the last coefficient to the right value to begin applying an
 * anticausal filter for the given pole to the sequence.
 *
 * -----
 *
 * The causal, c+, and anticausal, c-, filter coefficients, and the
 * filtered sequence, c, are computed from the sequence, s, using the
 * following recursive relations:
 *
 *   c+[i] = s[i] + z * c+[i-1]
 *   c-[i] = s[i] + z * c-[i+1]
 *   c[i] = z / (1 - z**2) * (c+[i] + c-[i] - s[i])
 *
 * For the anticausal filter, if the extension mode is NI_EXTEND_MIRROR,
 * we also know that c+[n-1] = c-[n-1], so:
 *
 *   c[n-1] = z / (1 - z**2) * (z * c+[n-2] + c+[n-1])
 *
 */
void
set_initial_anticausal_coefficient(double *coefficients, npy_intp len,
                                   double pole);

#endif

