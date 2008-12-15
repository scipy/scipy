/* Written by Charles Harris charles.harris@sdl.usu.edu */

#include "zeros.h"

double
bisect(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)
{
    int i;
    double dm,xm,fm,fa,fb,tol;

    tol = xtol + rtol*(fabs(xa) + fabs(xb));

    fa = (*f)(xa,params);
    fb = (*f)(xb,params);
    params->funcalls = 2;
    if (fa*fb > 0) {ERROR(params,SIGNERR,0.0);}
    if (fa == 0) return xa;
    if (fb == 0) return xb;
    dm = xb - xa;
    params->iterations = 0;
    for(i=0; i<iter; i++) {
        params->iterations++;
        dm *= .5;
        xm = xa + dm;
        fm = (*f)(xm,params);
        params->funcalls++;
        if (fm*fa >= 0) {
            xa = xm;
        }
        if (fm == 0 || fabs(dm) < tol)
            return xm;
    }
    ERROR(params,CONVERR,xa);
}
