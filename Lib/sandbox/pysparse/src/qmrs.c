#include <math.h>
#include <stdio.h>
#include "Python.h"
#include "pysparse/blas.h"
#include "pysparse/fortran.h"
#define SPMATRIX_UNIQUE_SYMBOL itsolvers_spmatrix
#include "pysparse/spmatrix.h"
#include "pysparse/qmrs.h"

#define SpMatrix_PRECON(prec_obj, n, x, y) \
        {if (SpMatrix_Precon((prec_obj),(n),(x),(y))) return -1;}
#define SpMatrix_MATVEC(mat_obj, n, x, m, y) \
        {if (SpMatrix_Matvec((mat_obj), (n), (x), (m), (y))) return -1;}

/*******************************************************************************
  parameters:
  -----------
  n         order of linear system
  b         right hand side of linear system
  x         solution vector
  matvec    matrix vector product function, y := A*x
  precon    preconditioning function,  y := M^-1*x
            my be a NULL pointer
  work      work array of size 6*n
  tol       requested error tolerance
  maxit     maximum allowed iteration steps
  iter      number of iterations 
  err       error of solution

  simplified decoupled QMR method for J-symmetric matrices (Freund & Nachtigal)

********************************************************************************/

int Itsolvers_qmrs_kernel(int n, 
			  double *b, 
			  double *x, 
			  double *work, 
			  double tol, 
			  int maxit, 
			  int *iter, 
			  double *err,
			  PyObject *mat_obj,
			  PyObject *prec_obj)
{
  /* Local variables */
  int ONE = 1;
  double beta;
  double res_init;
  int i;
  double delta, theta;
  double c0, c1, theta0, cc, xi1, rho1inv, tau, eta0, eps0, rho0, rho1, d__1;
  /* work arrays */
  double *wrk1;
  double *p;
  double *d;
  double *v1;
  double *t;
  double *g;

  /* setup pointers into work */
  wrk1 = work;
  p    = work +   n;
  d    = work + 2*n;
  v1   = work + 3*n;
  t    = work + 4*n;
  g    = work + 5*n;

  F77(dcopy)(&n, b, &ONE, v1, &ONE);
  rho0 = F77(dnrm2)(&n, v1, &ONE);
  tau = rho0;
  for (i = 0; i < n; ++i) {
    v1[i] /= rho0;
    p[i] = 0.0;
    g[i] = 0.0;
    d[i] = 0.0;
    x[i] = 0.0;
  }
  c0 = 1.0;
  eps0 = 1.0;
  xi1 = 1.0;
  theta0 = 0.0;
  eta0 = -1.0;
  res_init = rho0;
  *err = 1.0;
  *iter = 0;
  while(*err > tol && *iter < maxit) {
    ++(*iter);
    if (eps0 == 0.0) {
      return -6;
    }
    if (prec_obj != NULL)
      SpMatrix_PRECON(prec_obj, n, v1, wrk1)
    else
      F77(dcopy)(&n, v1, &ONE, wrk1, &ONE);
    delta = F77(ddot)(&n, wrk1, &ONE, v1, &ONE);
    if (delta == 0.0) {
      return -2;
    }
    cc = xi1 * (delta / eps0);
    /* it is: g = P\p */
    for (i = 0; i < n; ++ i) {
      p[i] = v1[i] - p[i] * cc;
      g[i] = wrk1[i] - g[i] * cc;
    }
    SpMatrix_MATVEC(mat_obj, n, g, n, t);
    eps0 = F77(ddot)(&n, g, &ONE, t, &ONE);
    beta = eps0 / delta;
    for (i = 0; i < n; ++ i) {
      v1[i] = t[i] - v1[i] * beta;
    }
    rho1 = F77(dnrm2)(&n, v1, &ONE);
    xi1 = rho1;
    if (c0 * fabs(beta) == 0.0) {
      return -6;
    }
    theta = rho1 / (c0 * fabs(beta));
    c1 = 1.0 / sqrt(theta * theta + 1.0);
    if (beta * (c0 * c0) == 0.0) {
      return -6;
    }
    eta0 = -eta0 * rho0 * (c1 * c1) / (beta * (c0 * c0));
    tau = tau * theta * c1;
    if (rho1 == 0.0) {
      return -6;
    }
    /* Computing 2nd power */
    d__1 = theta0 * c1;
    cc = d__1 * d__1;
    rho1inv = 1.0 / rho1;
    for (i = 0; i < n; ++i) {
      d[i] = p[i] * eta0 + d[i] * cc;
      x[i] += d[i];
      v1[i] *= rho1inv;
    }
    if (xi1 == 0.0) {
      return -6;
    }
    rho0 = rho1;
    *err = tau / res_init;
    c0 = c1;
    theta0 = theta;
    /*       print *,iter,err,tau/res_init */
  }
  /*       print *,iter,err */
  /*       if (err.gt.tol) print *,'convergence not reached' */
  if (prec_obj != NULL) {
    SpMatrix_PRECON(prec_obj, n, x, wrk1);
    F77(dcopy)(&n, wrk1, &ONE, x, &ONE);
  }
  if (*err < tol)
    return 0;
  else
    return -1;
}
