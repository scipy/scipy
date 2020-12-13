/* Written by Charles Harris charles.harris@sdl.usu.edu */

#include <math.h>
#include "zeros.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

/*
 At the top of the loop the situation is the following:

    1. the root is bracketed between xa and xb
    2. xa is the most recent estimate
    3. xp is the previous estimate
    4. |fp| < |fb|

  The order of xa and xp doesn't matter, but assume xp < xb. Then xa lies to
  the right of xp and the assumption is that xa is increasing towards the root.
  In this situation we will attempt quadratic extrapolation as long as the
  condition

  *  |fa| < |fp| < |fb|

  is satisfied. That is, the function value is decreasing as we go along.
  Note the 4 above implies that the right inequlity already holds.

  The first check is that xa is still to the left of the root. If not, xb is
  replaced by xp and the interval reverses, with xb < xa. In this situation
  we will try linear interpolation. That this has happened is signaled by the
  equality xb == xp;


  The second check is that |fa| < |fb|. If this is not the case, we swap
  xa and xb and resort to bisection.

*/

double
brenth(callback_type f, double xa, double xb, double xtol, double rtol,
       int iter, void *func_data, scipy_zeros_info *solver_stats)
{
    double xpre = xa, xcur = xb;
    double xblk = 0., fpre, fcur, fblk = 0., spre = 0., scur = 0., sbis;
    /* the tolerance is 2*delta */
    double delta;
    double stry, dpre, dblk;
    int i;
    solver_stats->error_num = INPROGRESS;

    fpre = (*f)(xpre,func_data);
    fcur = (*f)(xcur,func_data);
    solver_stats->funcalls = 2;
    if (fpre*fcur > 0) {
        solver_stats->error_num = SIGNERR;
        return 0.;
    }
    if (fpre == 0) {
        solver_stats->error_num = CONVERGED;
        return xpre;
    }
    if (fcur == 0) {
        solver_stats->error_num = CONVERGED;
        return xcur;
    }
    solver_stats->iterations = 0;
    for (i = 0; i < iter; i++) {
        solver_stats->iterations++;
        if (fpre*fcur < 0) {
            xblk = xpre;
            fblk = fpre;
            spre = scur = xcur - xpre;
        }
        if (fabs(fblk) < fabs(fcur)) {
            xpre = xcur;
            xcur = xblk;
            xblk = xpre;

            fpre = fcur;
            fcur = fblk;
            fblk = fpre;
        }

        delta = (xtol + rtol*fabs(xcur))/2;
        sbis = (xblk - xcur)/2;
        if (fcur == 0 || fabs(sbis) < delta) {
            solver_stats->error_num = CONVERGED;
            return xcur;
        }

        if (fabs(spre) > delta && fabs(fcur) < fabs(fpre)) {
            if (xpre == xblk) {
                /* interpolate */
                stry = -fcur*(xcur - xpre)/(fcur - fpre);
            }
            else {
                /* extrapolate */
                dpre = (fpre - fcur)/(xpre - xcur);
                dblk = (fblk - fcur)/(xblk - xcur);
                stry = -fcur*(fblk - fpre)/(fblk*dpre - fpre*dblk);
            }

            if (2*fabs(stry) < MIN(fabs(spre), 3*fabs(sbis) - delta)) {
                /* accept step */
                spre = scur;
                scur = stry;
            }
            else {
                /* bisect */
                spre = sbis;
                scur = sbis;
            }
        }
        else {
            /* bisect */
            spre = sbis;
            scur = sbis;
        }

        xpre = xcur;
        fpre = fcur;
        if (fabs(scur) > delta) {
            xcur += scur;
        }
        else {
            xcur += (sbis > 0 ? delta : -delta);
        }

        fcur = (*f)(xcur, func_data);
        solver_stats->funcalls++;
    }
    solver_stats->error_num = CONVERR;
    return xcur;
}
