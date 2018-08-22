#include <math.h>
#include "newton.h"
#include "../../Zeros/zeros.h"

double
halley(callback_type func, double p0, callback_type fprime, default_parameters *params, double tol, int maxiter, callback_type fprime2)
{
    int i;
    double fder, fval, p, newton_step, fder2;

    if (maxiter < 0) {
        params->error_num = SIGNERR;
        return p0;
    }
    if (tol < 0) {
        params->error_num = SIGNERR;
        return p0;
    }

    params->funcalls = 0;
    params->iterations = 0;
    for (i=0; i<maxiter; i++) {
        params->iterations++;
        fval = (*func)(p0, params);
        params->funcalls++;
        if (fval == 0) {
            return p0;
        }
        fder = (*fprime)(p0, params);
        params->funcalls++;
        if (fder == 0) {
            params->error_num = CONVERR;
            return p0;
        }
        newton_step = fval / fder;
        fder2 = (*fprime2)(p0, params);
        params->funcalls++;
        p = p0 + newton_step / (1.0 - 0.5 * newton_step * fder2 / fder);
        if (fabs(p - p0) < tol) {
            return p;
        }
        p0 = p;
    }
    params->error_num = CONVERR;
    return p;
}
