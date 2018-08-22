#include <math.h>
#include "newton.h"
#include "../../Zeros/zeros.h"

#define DX 1e-4

double
secant(callback_type func, double p0, default_parameters *params, double tol, int maxiter)
{
    int i;
    double p, p1, q0, q1;

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
    p1 = p0 * (1 + DX);
    if (p0 >= 0) {
        p1 += DX;
    } else {
        p1 -= DX;
    }
    q0 = (*func)(p0, params);
    params->funcalls++;
    q1 = (*func)(p1, params);
    params->funcalls++;
    for (i=0; i<maxiter; i++) {
        params->iterations++;
        if (q1 == q0) {
            if (p1 != p0) {
                params->error_num = CONVERR;
                return (p0 + p1) / 2.0;
            }
        } else {
            p = p1 - q1 * (p1 - p0) / (q1 - q0);
        }
        if (fabs(p - p1) < tol) {
            return p;
        }
        p0 = p1;
        q0 = q1;
        p1 = p;
        q1 = (*func)(p1, params);
        params->funcalls++;
    }
    params->error_num = CONVERR;
    return p;
}
