#include "ni_support.h"
#include "ni_splines.h"
#include <math.h>


int
get_spline_interpolation_weights(double x, int order, double *weights)
{
    int i;
    double y, z, t;

    /* Convert x to the delta to the middle knot. */
    x -= floor(order & 1 ? x : x + 0.5);
    y = x;
    z = 1.0 - x;

    switch (order) {
        case 1:
            /* 0 <= x < 1*/
            weights[0] = 1.0 - x;
            break;
        case 2:
            /* -0.5 < x <= 0.5 */
            weights[1] = 0.75 - x * x;
            /* For weights[0] we'd normally do:
             *
             *   y = 1 + x  # 0.5 < y <= 1.5
             *   yy = 1.5 - y  # yy = 0.5 - x
             *   weights[0] = 0.5 * yy * yy
             *
             * So we set y = 0.5 - x directly instead.
             */
            y = 0.5 - x;
            weights[0] = 0.5 * y * y;
            break;
        case 3:
            /* y = x, 0 <= y < 1 */
            weights[1] = (y * y * (y - 2.0) * 3.0 + 4.0) / 6.0;
            /* z = 1 - x, 0 < z <= 1 */
            weights[2] = (z * z * (z - 2.0) * 3.0 + 4.0) / 6.0;
            /*
             * For weights[0] we would normally do:
             *
             *   y += 1.0  # y = 1.0 + x, 1 <= y < 2
             *   yy = 2.0 - y  # yy = 1 - x
             *   weights[0] = yy * yy * yy / 6.0
             *
             * But we already have yy in z.
             */
            weights[0] = z * z * z / 6.0;
            break;
        case 4:
            /* -0.5 < x <= 0.5 */
            t = x * x;
            weights[2] = t * (t * 0.25 - 0.625) + 115.0 / 192.0;
            /* y = 1 + x, 0.5 < y <= 1.5 */
            y = 1.0 + x;
            weights[1] = y * (y * (y * (5.0 - y) / 6.0 - 1.25) + 5.0 / 24.0) +
                         55.0 / 96.0;
            /* z = 1 - x, 0.5 <= z < 1.5 */
            weights[3] = z * (z * (z * (5.0 - z) / 6.0 - 1.25) + 5.0 / 24.0) +
                         55.0 / 96.0;
            /*
             * For weights[0] we would normally do:
             *
             *   y += 1.0  # y = 2.0 + x, 1.5 <= y < 2.5
             *   yy = 2.5 - y  # yy = 0.5 - x
             *  weights[0] = yy**4 / 24.0
             *
             * So we set y = 0.5 - x directly instead.
             */
            y = 0.5 - x;
            t = y * y;
            weights[0] = t * t / 24.0;
            break;
        case 5:
            /* y = x, 0 <= y < 1 */
            t = y * y;
            weights[2] = t * (t * (0.25 - y / 12.0) - 0.5) + 0.55;
            /* z = 1 - x, 0 < z <= 1 */
            t = z * z;
            weights[3] = t * (t * (0.25 - z / 12.0) - 0.5) + 0.55;
            /* y = 1 + x, 1 <= y < 2 */
            y += 1.0;
            weights[1] = y * (y * (y * (y * (y / 24.0 - 0.375) + 1.25) - 1.75)
                              + 0.625) + 0.425;
            /* z = 2 - x, 1 < z <= 2 */
            z += 1.0;
            weights[4] = z * (z * (z * (z * (z / 24.0 - 0.375) + 1.25) - 1.75)
                              + 0.625) + 0.425;
            /* For weights[0] we would normally do:
             *
             *   y += 1.0  # y = 2.0 + x, 2 <= y < 3
             *   yy = 3.0 - y  # yy = 1.0 - x
             *   weights[0] = yy**5 / 120.0
             *
             * So we set y = 2.0 - y = 1.0 - x directly instead.
             */
            y = 1.0 - x;
            t = y * y;
            weights[0] = y * t * t / 120.0;
            break;
        default:
            return 1; /* Unsupported spline order. */
    }

    /* All interpolation weights add to 1.0, so use it for the last one. */
    weights[order] = 1.0;
    for (i = 0; i < order; ++i) {
        weights[order] -= weights[i];
    }

    return 0;
}


int
get_filter_poles(int order, int *npoles, double *poles)
{
    *npoles = order / 2;
    /*
     * If this assert is triggered, someone added more orders here but
     * MAX_SPLIINE_FILTER_POLES was not kept in sync.
     */
    assert(*npoles <= MAX_SPLINE_FILTER_POLES);

    switch (order) {
        case 2:
            /* sqrt(8.0) - 3.0 */
            poles[0] = -0.171572875253809902396622551580603843;
            break;
        case 3:
            /* sqrt(3.0) - 2.0 */
            poles[0] = -0.267949192431122706472553658494127633;
            break;
        case 4:
            /* sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0 */
            poles[0] = -0.361341225900220177092212841325675255;
            /* sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0 */
            poles[1] = -0.013725429297339121360331226939128204;
            break;
        case 5:
            /* sqrt(67.5 - sqrt(4436.25)) + sqrt(26.25) - 6.5 */
            poles[0] = -0.430575347099973791851434783493520110;
            /* sqrt(67.5 + sqrt(4436.25)) - sqrt(26.25) - 6.5 */
            poles[1] = -0.043096288203264653822712376822550182;
            break;
        default:
            return 1; /* Unsupported order. */
    };

    return 0;
}


typedef void (init_fn)(npy_double*, const npy_intp, const double);


static void
_init_causal_mirror(double *c, const npy_intp n, const double z)
{
    npy_intp i;
    double z_i = z;
    const double z_n_1 = pow(z, n - 1);

    c[0] = c[0] + z_n_1 * c[n - 1];
    for (i = 1; i < n - 1; ++i) {
        c[0] += z_i * (c[i] + z_n_1 * c[n - 1 - i]);
        z_i *= z;
    }
    c[0] /= 1 - z_n_1 * z_n_1;
}


static void
_init_anticausal_mirror(double *c, const npy_intp n, const double z)
{
    c[n - 1] = (z * c[n - 2] + c[n - 1]) * z / (z * z - 1);
}


static void
_init_causal_wrap(double *c, const npy_intp n, const double z)
{
    npy_intp i;
    double z_i = z;

    for (i = 1; i < n; ++i) {
        c[0] += z_i * c[n - i];
        z_i *= z;
    }
    c[0] /= 1 - z_i; /* z_i = pow(z, n) */
}


static void
_init_anticausal_wrap(double *c, const npy_intp n, const double z)
{
    npy_intp i;
    double z_i = z;

    for (i = 0; i < n - 1; ++i) {
        c[n - 1] += z_i * c[i];
        z_i *= z;
    }
    c[n - 1] *= z / (z_i - 1); /* z_i = pow(z, n) */
}


static void
_init_causal_reflect(double *c, const npy_intp n, const double z)
{
    npy_intp i;
    double z_i = z;
    const double z_n = pow(z, n);
    const double c0 = c[0];

    c[0] = c[0] + z_n * c[n - 1];
    for (i = 1; i < n; ++i) {
        c[0] += z_i * (c[i] + z_n * c[n - 1 - i]);
        z_i *= z;
    }
    c[0] *= z / (1 - z_n * z_n);
    c[0] += c0;
}


void
_init_anticausal_reflect(double *c, const npy_intp n, const double z)
{
    c[n - 1] *= z / (z - 1);
}


/*
 * The application of the filter is performed in two passes over the
 * coefficient array, one forward and the other backward. For a
 * detailed discussion of the method see e.g.:
 *
 * Unser, Michael, Akram Aldroubi, and Murray Eden. "Fast B-spline
 * transforms for continuous image representation and interpolation."
 * IEEE Transactions on pattern analysis and machine intelligence 13.3
 * (1991): 277-285.
 *
 * A key part of the process is initializing the first coefficient for
 * each pass, which depends on the boundary conditions chosen to extend
 * the input image beyond its boundaries. The method to initialize these
 * values for the NI_EXTEND_MIRROR mode is discussed in the above paper.
 * For NI_EXTEND_WRAP and NI_EXTEND_REFLECT, the unpublished method was
 * obtained from a private communication with Dr. Philippe Thévenaz.
 */
static void
_apply_filter(double *c, npy_intp n, double z, init_fn *causal_init,
              init_fn *anticausal_init)
{
    npy_intp i;

    causal_init(c, n, z);
    for (i = 1; i < n; ++i) {
        c[i] += z * c[i - 1];
    }
    anticausal_init(c, n, z);
    for (i = n - 2; i >= 0; --i) {
        c[i] = z * (c[i + 1] - c[i]);
    }
}


static void
_apply_filter_gain(double *c, npy_intp n, const double *zs, int nz)
{
    double gain = 1.0;

    while (nz--) {
        const double z = *zs++;
        gain *= (1.0 - z) * (1.0 - 1.0 / z);
    }

    while (n--) {
        *c++ *= gain;
    }
}


void
apply_filter(double *coefficients, const npy_intp len, const double *poles,
             int npoles, NI_ExtendMode mode)
{
    init_fn *causal = NULL;
    init_fn *anticausal = NULL;

    //Note: This switch statement should match the settings used for
	//      the spline_mode variable in NI_GeometricTransform
    switch(mode) {
        case NI_EXTEND_GRID_CONSTANT:
        case NI_EXTEND_CONSTANT:
        case NI_EXTEND_MIRROR:
        case NI_EXTEND_WRAP:
            causal = &_init_causal_mirror;
            anticausal = &_init_anticausal_mirror;
            break;
        case NI_EXTEND_GRID_WRAP:
            causal = &_init_causal_wrap;
            anticausal = &_init_anticausal_wrap;
            break;
        case NI_EXTEND_NEAREST:
        case NI_EXTEND_REFLECT:
            causal = &_init_causal_reflect;
            anticausal = &_init_anticausal_reflect;
            break;
        default:
            assert(0); /* We should never get here. */
    }

    _apply_filter_gain(coefficients, len, poles, npoles);

    while (npoles--) {
        _apply_filter(coefficients, len, *poles++, causal, anticausal);
    }
}
