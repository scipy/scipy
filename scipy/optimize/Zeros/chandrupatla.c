/*
 * copied from ridder.c and adapted for Chandrupatla method
 * doi:10.1016/S0965-9978(96)00051-8
 */

#include "zeros.h"
#include <math.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define SIGN(a) ((a) > 0. ? 1. : -1.)

/* Sets solver_stats->error_num
    SIGNERR for sign_error;
    CONVERR for convergence_error;
*/

double chandrupatla(callback_type f, double xa, double xb, double xtol,
                    double rtol, int iter, void *func_data_param,
                    scipy_zeros_info *solver_stats) {
  int i;
  double dm, xm, xc, fm, fa, fb, fc, tol, t, phi, xi;
  solver_stats->error_num = INPROGRESS;

  tol = xtol + rtol * MIN(fabs(xa), fabs(xb));
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
  if (signbit(fa) == signbit(fb)) {
    solver_stats->error_num = SIGNERR;
    return 0.;
  }

  t = 0.5;
  xm = xa + t * (xb - xa); // to avoid maybe uninitialized warning
  solver_stats->iterations = 0;
  for (i = 0; i < iter; i++) {
    solver_stats->iterations++;
    dm = t * (xb - xa);
    xm = xa + dm;
    fm = (*f)(xm, func_data_param);
    solver_stats->funcalls += 1;

    if (signbit(fa) != signbit(fm)) { // warning: is noted wrongly in the flow
                                      // chart of original paper
      xc = xb;
      xb = xa;
      xa = xm;
      fc = fb;
      fb = fa;
      fa = fm;
    } else {
      xc = xa;
      xa = xm;
      fc = fa;
      fa = fm;
    }

    tol = xtol + rtol * xm;
    if (fm == 0.0 || fabs(xb - xa) < tol) {
      solver_stats->error_num = CONVERGED;
      return xm;
    }

    // test if inverse quadratic interpolation is applicable
    phi = (fa - fb) / (fc - fb);
    xi = (xa - xb) / (xc - xb);
    if (1 - sqrt(1 - xi) < phi && phi < sqrt(xi)) {
      // do inverse quadratic interpolation
      t = (fa / (fa - fb)) * (fc / (fc - fb)) -
          ((xc - xa) / (xb - xa)) * (fa / (fc - fa)) * (fb / (fb - fc));
    } else {
      // do bisection
      t = 0.5;
    }
  }
  solver_stats->error_num = CONVERR;
  return xm;
}
