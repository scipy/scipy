#include "misc.h"
#include <math.h>

#define MAX_ITERATIONS              100
#define FP_CMP_WITH_BISECT_NITER    4
#define FP_CMP_WITH_BISECT_WIDTH    4.0

static double
max(double a, double b)
{
    return (a > b ? a : b);
}

/*
   Use a combination of bisection and false position to find a root
   of a function within a given interval. This is guaranteed to converge,
   and always keeps a bounding interval, unlike Newton's method.

   The false position steps are either unmodified, or modified with
   the Anderson-Bjorck method as appropiate. Theoretically, this has
   a "speed of convergence" of 1.7 (bisection is 1, Newton is 2).

   Input
   -----
    a, b:   initial bounding interval
    fa, fb: value of f() at a and b
    f, f_extra: function to find root of is f(x, f_extra)
    abserr, relerr: absolute and relative errors on the bounding interval
    bisect_til: If > 0.0, perform bisection until the width of the
                bounding interval is less than this.

   Output
   ------
    a, b, fa, fb: Final bounding interval and function values
    best_x, best_f: Best root approximation and the function value there

   Returns
   -------
    FSOLVE_CONVERGED: Bounding interval is smaller than required error.
    FSOLVE_NOT_BRACKET: Initial interval is not a bounding interval.
    FSOLVE_EXACT: An exact root was found (best_f = 0)


   Note that this routine was designed initially to work with gammaincinv, so
   it may not be tuned right for other problems. Don't use it blindly.
 */
fsolve_result_t
false_position(double *a, double *fa, double *b, double *fb,
               objective_function f, void *f_extra,
               double abserr, double relerr, double bisect_til,
               double *best_x, double *best_f, double *errest)
{
    double x1=*a, f1=*fa, x2=*b, f2=*fb;
    fsolve_result_t r = FSOLVE_CONVERGED;
    double gamma = 1.0;
    enum {bisect, falsep} state = bisect;
    int n_falsep = 0;
    double x3, f3;
    double w, last_bisect_width;
    double tol;
    int niter;

    if (f1*f2 >= 0.0) {
        return FSOLVE_NOT_BRACKET;
    }
    if (bisect_til > 0.0) {
        state = bisect;
    } else {
        state = falsep;
    }
    w = fabs(x2 - x1);
    last_bisect_width = w;
    for (niter=0; niter < MAX_ITERATIONS; niter++) {
        switch (state) {
        case bisect: {
            x3 = 0.5 * (x1 + x2);
            if (x3 == x1 || x3 == x2) {
                /* i.e., x1 and x2 are successive floating-point numbers. */
                *best_x = x3;
                *best_f = (x3==x1) ? f1 : f2;
                goto finish;
            }
            f3 = f(x3, f_extra);
            if (f3 == 0.0) {
                goto exact_soln;
            }
            if (f3*f2 < 0.0) {
                x1 = x2; f1 = f2;
            }
            x2 = x3; f2 = f3;
            w = fabs(x2 - x1);
            last_bisect_width = w;
            if (bisect_til > 0.0) {
                if (w < bisect_til) {
                    bisect_til = -1.0;
                    gamma = 1.0;
                    n_falsep = 0;
                    state = falsep;
                }
            } else {
                gamma = 1.0;
                n_falsep = 0;
                state = falsep;
            }
            break;
        }
        case falsep: {
            double s12 = (f2 - gamma*f1) / (x2 - x1);
            x3 = x2 - f2/s12;
            f3 = f(x3, f_extra);
            if (f3 == 0.0) {
                goto exact_soln;
            }
            n_falsep += 1;
            if (f3*f2 < 0.0) {
                gamma = 1.0;
                x1 = x2; f1 = f2;
            } else {
                /* Anderson-Bjorck method */
                double g = 1.0 - f3 / f2;
                if (g <= 0.0) { g = 0.5; }
                /* It's not really clear from the sources I've looked at,
                   but I believe this is *= instead of =. */
                gamma *= g;
            }
            x2 = x3; f2 = f3;
            w = fabs(x2 - x1);
            /* Sanity check. For every 4 false position checks, see if we
               really are decreasing the interval by comparing to what
               bisection would have achieved (or, rather, a bit more lenient
               than that -- interval decreased by 4 instead of by 16, as
               the fp could be decreasing gamma for a bit).

               Note that this should guarantee convergence, as it makes
               sure that we always end up decreasing the interval width
               with a bisection.
             */
            if (n_falsep > FP_CMP_WITH_BISECT_NITER) {
                if (w*FP_CMP_WITH_BISECT_WIDTH > last_bisect_width) {
                    state = bisect;
                }
                n_falsep = 0;
                last_bisect_width = w;
            }
            break;
        }
        }
        tol = abserr + relerr*max(max(fabs(x1), fabs(x2)), 1.0);
        if (w <= tol) {
            if (fabs(f1) < fabs(f2)) {
                *best_x = x1; *best_f = f1;
            } else {
                *best_x = x2; *best_f = f2;
            }
            goto finish;
        }
    }
    r = FSOLVE_MAX_ITERATIONS;
    *best_x = x3; *best_f = f3;
    goto finish;
exact_soln:
    *best_x = x3; *best_f = 0.0;
    r = FSOLVE_EXACT;
finish:
    *a = x1; *fa = f1; *b = x2; *fb = f2;
    *errest = w;
    return r;
}
