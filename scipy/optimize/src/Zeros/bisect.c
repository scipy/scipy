/* Written by Charles Harris charles.harris@sdl.usu.edu */

#include <math.h>
#include "zeros.h"

double
bisect(callback_type f, double xa, double xb, double xtol, double rtol,
       int iter, void *func_data_param, scipy_zeros_info *solver_stats)
{
    int i;
    double dm,xm,fm,fa,fb;
    solver_stats->error_num = INPROGRESS;

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
    dm = xb - xa;
    solver_stats->iterations = 0;
    for (i=0; i<iter; i++) {
        solver_stats->iterations++;
        dm *= .5;
        xm = xa + dm;
        fm = (*f)(xm, func_data_param);
        solver_stats->funcalls++;
        if (signbit(fm)==signbit(fa)) {
            xa = xm;
        }
        if (fm == 0 || fabs(dm) < xtol + rtol*fabs(xm)) {
            solver_stats->error_num = CONVERGED;
            return xm;
        }
    }
    solver_stats->error_num = CONVERR;
    return xa;
}
