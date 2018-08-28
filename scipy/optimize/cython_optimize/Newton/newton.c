#include <math.h>
#include "newton.h"
#include "../../Zeros/zeros.h"

double
newton(callback_type func, double p0, callback_type fprime, default_parameters *params, double tol, int maxiter)
{
    int i;
    double fder, fval, p;

    if (maxiter <= 0) {
        params->error_num = SIGNERR;
        return p0;
    }
    if (tol <= 0) {
        params->error_num = SIGNERR;
        return p0;
    }

    params->funcalls = 0;
    params->iterations = 0;
    params->error_num = 0;
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
        p = p0 - fval / fder;
        if (fabs(p - p0) < tol) {
            return p;
        }
        p0 = p;
    }
    params->error_num = CONVERR;
    return p;
}
