/*
 * Originally written by Charles Harris charles.harris@sdl.usu.edu.
 * Modified by Travis Oliphant to not depend on Python.
 */

#include <math.h>
#include "zeros.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define SIGN(a) ((a) > 0. ? 1. : -1.)

/* Sets params->error_num SIGNERR for sign_error;
                         CONVERR for convergence_error;
*/

double
ridder(callback_type f, double xa, double xb, double xtol, double rtol,
       int iter, default_parameters *params)
{
    int i;
    double dm,dn,xm,xn=0.0,fn,fm,fa,fb,tol;

    tol = xtol + rtol*(fabs(xa) + fabs(xb));
    fa = (*f)(xa,params);
    fb = (*f)(xb,params);
    params->funcalls = 2;
    if (fa*fb > 0) {
        params->error_num = SIGNERR;
        return 0.;
    }
    if (fa == 0) {
        return xa;
    }
    if (fb == 0) {
        return xb;
    }

    params->iterations=0;
    for (i=0; i<iter; i++) {
        params->iterations++;
        dm = 0.5*(xb - xa);
        xm = xa + dm;
        fm = (*f)(xm,params);
        dn = SIGN(fb - fa)*dm*fm/sqrt(fm*fm - fa*fb);
        xn = xm - SIGN(dn) * MIN(fabs(dn), fabs(dm) - .5*tol);
        fn = (*f)(xn,params);
        params->funcalls += 2;
        if (fn*fm < 0.0) {
            xa = xn; fa = fn; xb = xm; fb = fm;
        }
        else if (fn*fa < 0.0) {
            xb = xn; fb = fn;
        }
        else {
            xa = xn; fa = fn;
        }
        if (fn == 0.0 || fabs(xb - xa) < tol) {
            return xn;
        }
    }
    params->error_num = CONVERR;
    return xn;
}
