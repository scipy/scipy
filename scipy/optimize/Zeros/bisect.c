/* Written by Charles Harris charles.harris@sdl.usu.edu */

#include <math.h>
#include "zeros.h"
#include <stdio.h>

#define OFILE stdout

double
bisect(callback_type f, double xa, double xb, double xtol, double rtol,
       int iter, void *func_data, scipy_zeros_info *solver_stats)
{
    int i;
    double dm,xm,fm,fa,fb;
    fprintf(OFILE, "BISECT:Start\n");
    fprintf(OFILE, "a=%f b=%f\n", xa, xb);
    fprintf(OFILE, "f=%p func_data=%p\n", (void *)f, func_data);
    fflush(OFILE);

    fa = (*f)(xa, func_data);
    fb = (*f)(xb, func_data);
    solver_stats->funcalls = 2;
    if (fa*fb > 0) {
        solver_stats->error_num = SIGNERR;
        return 0.;
    }
    solver_stats->error_num = CONVERGED;
    if (fa == 0) {
        return xa;
    }
    if (fb == 0) {
        return xb;
    }
    dm = xb - xa;
    solver_stats->iterations = 0;
    for (i=0; i<iter; i++) {
        solver_stats->iterations++;
        dm *= .5;
        xm = xa + dm;
        fm = (*f)(xm, func_data);
        solver_stats->funcalls++;
        if (fm*fa >= 0) {
            xa = xm;
        }
        if (fm == 0 || fabs(dm) < xtol + rtol*fabs(xm)) {
            fprintf(OFILE, "Converged!\n"); fflush(OFILE);
            return xm;
        }
        fprintf(OFILE, "Still working: i=%2d: xa=%f, xm=%f fm=%f\n", i, xa, xm, fm); fflush(OFILE);
    }
    fprintf(OFILE, "FAIL: Not converged!\n"); fflush(OFILE);
    solver_stats->error_num = CONVERR;
    return xa;
}
