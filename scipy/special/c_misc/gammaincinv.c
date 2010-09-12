#include <Python.h>
#include <numpy/npy_math.h>

#include <stdio.h>
#include <math.h>

#include "../cephes.h"
#undef fabs
#include "misc.h"

/* Limits after which to issue warnings about non-convergence */
#define ALLOWED_ATOL (1e-306)
#define ALLOWED_RTOL (1e-6)

void scipy_special_raise_warning(char *fmt, ...);

/*
  Inverse of the (regularised) incomplete Gamma integral.

  Given a, find x such that igam(a, x) = y.
  For y not small, we just use igami(a, 1-y) (inverse of the complemented
  incomplete Gamma integral). For y small, however, 1-y is about 1, and we
  lose digits.

*/

extern double MACHEP, MAXNUM;

static double
gammainc(double x, double params[2])
{
    return cephes_igam(params[0], x) - params[1];
}

double
gammaincinv(double a, double y)
{
    double lo = 0.0, hi;
    double flo = -y, fhi = 0.25 - y;
    double params[2];
    double best_x, best_f, errest;
    fsolve_result_t r;

    if (a <= 0.0 || y <= 0.0 || y >= 0.25) {
        return cephes_igami(a, 1-y);
    }

    /* Note: flo and fhi must have different signs (and be != 0),
     *       otherwise fsolve terminates with an error.
     */

    params[0] = a;
    params[1] = y;
    hi = cephes_igami(a, 0.75);
    /* I found Newton to be unreliable. Also, after we generate a small
       interval by bisection above, false position will do a large step
       from an interval of width ~1e-4 to ~1e-14 in one step (a=10, x=0.05,
       but similiar for other values).
     */

    r = false_position(&lo, &flo, &hi, &fhi,
                       (objective_function)gammainc, params,
                       2*MACHEP, 2*MACHEP, 1e-2*a,
                       &best_x, &best_f, &errest);
    if (!(r == FSOLVE_CONVERGED || r == FSOLVE_EXACT) &&
            errest > ALLOWED_ATOL + ALLOWED_RTOL*fabs(best_x)) {
        scipy_special_raise_warning(
            "gammaincinv: failed to converge at (a, y) = (%.20g, %.20g): got %g +- %g, code %d\n",
            a, y, best_x, errest, r);
        best_x = NPY_NAN;
    }
    return best_x;
}
