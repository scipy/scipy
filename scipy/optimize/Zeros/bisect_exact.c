/* Written by Charles Harris charles.harris@sdl.usu.edu */

#include <math.h>
#include <stdint.h>
#include <string.h>
#include "zeros.h"

double
bisect_exact(callback_type f, double xa, double xb, double xtol, double rtol,
       int iter, void *func_data_param, scipy_zeros_info *solver_stats)
{
    int i;
    double xm,fm,fa,fb,fz;
    uint64_t xa_int,xb_int,xm_int;
    solver_stats->error_num = INPROGRESS;

    fa = (*f)(xa, func_data_param);
    fb = (*f)(xb, func_data_param);
    fz = (*f)(0, func_data_param);
    solver_stats->funcalls = 3;
    if (fa == 0) {
        solver_stats->error_num = CONVERGED;
        return xa;
    }
    if (fb == 0) {
        solver_stats->error_num = CONVERGED;
        return xb;
    }
    if (fz == 0) {
        solver_stats->error_num = CONVERGED;
        return 0.;
    }
    if (signbit(fa)==signbit(fb)) {
        solver_stats->error_num = SIGNERR;
        return 0.;
    }
    if(signbit(fz)==signbit(fb)){
        xb = copysign(0.0, xa);
    }else{
        xa = copysign(0.0, xb);
    }
    memcpy(&xa_int,&xa,sizeof xa_int);
    memcpy(&xb_int,&xb,sizeof xb_int);
    solver_stats->iterations = 0;
    for (i=0; i<iter; i++) {
        solver_stats->iterations++;
        xm_int = (xb_int-xa_int)/2+xa_int;
        memcpy(&xm,&xm_int,sizeof xm);
        fm = (*f)(xm, func_data_param);
        solver_stats->funcalls++;
        if (signbit(fm)==signbit(fa)) {
            xa_int = xm_int;
        }else{
            xb_int = xm_int;
        }
        if(fm==0){
            return xm;
        }
    }
    solver_stats->error_num = CONVERGED;
    memcpy(&xa,&xa_int,sizeof xa);
    return xa;
}
