#include <stdio.h>
#include <math.h>
#include "../cephes.h"
#undef fabs
#include "misc.h"

/*
  Inverse of the (regularised) incomplete Gamma integral.

  Given a, find x such that igam(a, x) = y.
  For y not small, we just use igami(a, 1-y) (inverse of the complemented
  incomplete Gamma integral). For y small, however, 1-y is about 1, and we
  lose digits.

*/

extern double MACHEP;

static double
gammainc(double x, double params[2])
{
    return cephes_igam(params[0], x) - params[1];
}

double
gammaincinv(double a, double y)
{
    if (a <= 0.0 || y <= 0.0 || y > 0.25) {
        return cephes_igami(a, 1-y);
    }

    /* I found Newton to be unreliable. Also, after we generate a small
       interval by bisection above, false position will do a large step
       from an interval of width ~1e-4 to ~1e-14 in one step (a=10, x=0.05,
       but similiar for other values).
     */

    double lo = 0.0, hi = cephes_igami(a, 0.75);
    double flo = -y, fhi = 0.25 - y;
    double params[2] = {a, y};
    double best_x, best_f;
    fsolve_result_t r;

    r = false_position(&lo, &flo, &hi, &fhi,
                       (objective_function)gammainc, params,
                       MACHEP, MACHEP, 1e-2*a,
                       &best_x, &best_f);
    if (r == FSOLVE_NOT_BRACKET) {
        best_x = 0.0;
    }
    return best_x;
}
