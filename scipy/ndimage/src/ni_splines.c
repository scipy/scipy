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


double
filter_gain(const double *poles, int npoles)
{
    double gain = 1.0;

    while (npoles--) {
        const double pole = *poles++;
        gain *= (1.0 - pole) * (1.0 - 1.0 / pole);
    }

    return gain;
}


void
apply_gain(double gain, double *coefficients, npy_intp len)
{
    while (len--) {
        *coefficients++ *= gain;
    }
}


void
set_initial_causal_coefficient(double *coefficients, npy_intp len,
                               double pole, double tolerance)
{
    int i;
    int last_coeff = len;
    double sum = 0.0;

    if (tolerance > 0.0) {
        last_coeff = ceil(log(tolerance)) / log(fabs(pole));
    }
    if (last_coeff < len) {
        double z_i = pole;

        sum = coefficients[0];
        for (i = 1; i < last_coeff; ++i) {
            sum += z_i * coefficients[i];
            z_i *= pole;
        }
    }
    else {
        double z_i = pole;
        const double inv_z = 1.0 / pole;
        const double z_n_1 = pow(pole, len - 1);
        double z_2n_2_i = z_n_1 * z_n_1 * inv_z;

        sum = coefficients[0] + coefficients[len - 1] * z_n_1;
        for (i = 1; i < len - 1; ++i) {
            sum += (z_i + z_2n_2_i) * coefficients[i];
            z_i *= pole;
            z_2n_2_i *= inv_z;
        }
        sum /= (1 - z_n_1 * z_n_1);
    }

    coefficients[0] = sum;
}


void
set_initial_anticausal_coefficient(double *coefficients, npy_intp len,
                                   double pole)
{
    coefficients[len - 1] = pole / (pole * pole - 1.0) *
                            (pole * coefficients[len - 2] +
                             coefficients[len - 1]) ;
}
