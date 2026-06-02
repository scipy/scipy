/* TNC : Minimization example */
/* $Jeannot: example.c,v 1.19 2005/01/28 18:27:31 js Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tnc.h"

static tnc_function function;

static int function(double x[], double *f, double g[], void *state)
{
  *f = pow(x[0],2.0)+pow(fabs(x[1]),3.0);
  g[0] = 2.0*x[0];
  g[1] = 3.0*pow(fabs(x[1]),2.0);
  if(x[1]<0) g[1] = -g[1];
  return 0;
}

int main(int argc, char **argv)
{
  int i, rc, maxCGit = 2, maxnfeval = 20, nfeval;
  double fopt = 1.0, f, g[2],
    x[2] = {-7.0, 3.0},
    xopt[2] = {0.0, 1.0},
    low[2], up[2],
    eta = -1.0, stepmx = -1.0,
    accuracy = -1.0, fmin = 0.0, ftol = -1.0, xtol = -1.0, pgtol = -1.0,
    rescale = -1.0;

  low[0] = - HUGE_VAL; low[1] = 1.0;
  up[0] = HUGE_VAL; up[1] = HUGE_VAL;

  rc = tnc(2, x, &f, g, function, NULL, low, up, NULL, NULL, TNC_MSG_ALL,
    maxCGit, maxnfeval, eta, stepmx, accuracy, fmin, ftol, xtol, pgtol,
    rescale, &nfeval, NULL);

  printf("After %d function evaluations, TNC returned:\n%s\n", nfeval,
    tnc_rc_string[rc - TNC_MINRC]);

  for (i = 0; i < 2; i++)
    printf("x[%d] = %.15f / xopt[%d] = %.15f\n", i, x[i], i, xopt[i]);

  printf("\n");
  printf("f    = %.15f / fopt    = %.15f\n", f, fopt);

  return 0;
}
