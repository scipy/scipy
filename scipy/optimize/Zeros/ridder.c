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
       int iter, void *func_data_param, scipy_zeros_info *solver_stats)
{
    int i;
    double dm,dn,xm,xn=0.0,fn,fm,fa,fb,tol;
    solver_stats->error_num = INPROGRESS;

    tol = xtol + rtol*MIN(fabs(xa), fabs(xb));
    fa = (*f)(xa, func_data_param);
    fb = (*f)(xb, func_data_param);
    solver_stats->funcalls = 2;
    if (fa == 0) {
        solver_stats->error_num = CONVERGED;
        return xa;
    }
    if (fb == 0) {
        solver_stats->error_num = CONVERGED;
        return xb;
    }
    if (signbit(fa)==signbit(fb)) {
        solver_stats->error_num = SIGNERR;
        return 0.;
    }

    solver_stats->iterations=0;
    for (i=0; i<iter; i++) {
        solver_stats->iterations++;
        
        dm = 0.5*(xb - xa);
        xm = xa + dm;
        fm = (*f)(xm, func_data_param);
        
        /* * Exiting if the midpoint is the root, to avoid zero division errors
         */
        if (fm == 0.0) {
            solver_stats->error_num = CONVERGED;
            return xm;
        }

        /* * FIX: Implement Equation 6 from Ridders' paper to avoid underflow.
         * Normalizing by fa to avoid underflow in intermediate terms (fm*fm).
         * Using -dm because the sign of (fm/fa) is opposite to 
         * the logic used in the faulty SIGN(fb-fa) implementation.
         */
        double ratio = fm / fa;
        dn = -dm * ratio / sqrt(ratio * ratio - fb / fa);

        xn = xm - SIGN(dn) * MIN(fabs(dn), fabs(dm) - .5*tol);
        fn = (*f)(xn, func_data_param);
        solver_stats->funcalls += 2;
        if (signbit(fn) != signbit(fm)) {
            xa = xn; fa = fn; xb = xm; fb = fm;
        }
        else if (signbit(fn) != signbit(fa)) {
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