/*
 * Originally written by Charles Harris charles.harris@sdl.usu.edu.
 * Modified by Travis Oliphant to not depend on Python.
 */

#include <math.h>
#include "zeros.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define SIGN(a) ((a) > 0. ? 1. : -1.)

/* Sets solver_stats->error_num
    SIGNERR for sign_error;
    CONVERR for convergence_error;
*/

double
ridder(callback_type f, double xa, double xb, double xtol, double rtol,
       int iter, void *func_data, scipy_zeros_info *solver_stats)
{
    int i;
    double dm,dn,xm,xn=0.0,fn,fm,fa,fb,tol;
    solver_stats->error_num = INPROGRESS;

    tol = xtol + rtol*MIN(fabs(xa), fabs(xb));
    fa = (*f)(xa, func_data);
    fb = (*f)(xb, func_data);
    solver_stats->funcalls = 2;
    if (fa*fb > 0) {
        solver_stats->error_num = SIGNERR;
        return 0.;
    }
    if (fa == 0) {
        solver_stats->error_num = CONVERGED;
        return xa;
    }
    if (fb == 0) {
        solver_stats->error_num = CONVERGED;
        return xb;
    }

    solver_stats->iterations=0;
    for (i=0; i<iter; i++) {
        solver_stats->iterations++;
        dm = 0.5*(xb - xa);
        xm = xa + dm;
        fm = (*f)(xm, func_data);
        dn = SIGN(fb - fa)*dm*fm/sqrt(fm*fm - fa*fb);
        xn = xm - SIGN(dn) * MIN(fabs(dn), fabs(dm) - .5*tol);
        fn = (*f)(xn, func_data);
        solver_stats->funcalls += 2;
        if (fn*fm < 0.0) {
            xa = xn; fa = fn; xb = xm; fb = fm;
        }
        else if (fn*fa < 0.0) {
            xb = xn; fb = fn;
        }
        else {
            xa = xn; fa = fn;
        }
        tol = xtol + rtol*xn;
        if (fn == 0.0 || fabs(xb - xa) < tol) {
            solver_stats->error_num = CONVERGED;
            return xn;
        }
    }
    solver_stats->error_num = CONVERR;
    return xn;
}
