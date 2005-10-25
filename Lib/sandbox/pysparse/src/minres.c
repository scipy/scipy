/**************************************************************************
*                                                                         *
*               Swiss Federal Institute of Technology (ETH),              *
*                       CH-8092 Zuerich, Switzerland                      *
*                                                                         *
*                       (C) 1999 All Rights Reserved                      *
*                                                                         *
*                                NOTICE                                   *
*                                                                         *
*  Permission to use, copy, modify, and distribute this software and      *
*  its documentation for any purpose and without fee is hereby granted    *
*  provided that the above copyright notice appear in all copies and      *
*  that both the copyright notice and this permission notice appear in    *
*  supporting documentation.                                              *
*                                                                         *
*  Neither the Swiss Federal Institute of Technology nor the author make  *
*  any representations about the suitability of this software for any     *
*  purpose.  This software is provided ``as is'' without express or       *
*  implied warranty.                                                      *
*                                                                         *
***************************************************************************
*
*  Implementation of the MINRES iterative linear solver
*
*  $Id: minres.c,v 1.1.1.1 2003/11/09 12:35:47 geus Exp $
*
*
**************************************************************************/
#include <math.h>
#include <stdio.h>
#include "Python.h"
#include "pysparse/fortran.h"
#include "pysparse/blas.h"
#define SPMATRIX_UNIQUE_SYMBOL itsolvers_spmatrix
#include "pysparse/spmatrix.h"
#include "pysparse/minres.h"

#define SpMatrix_PRECON(prec_obj, n, x, y) \
        {if (SpMatrix_Precon((prec_obj),(n),(x),(y))) return -1;}
#define SpMatrix_MATVEC(mat_obj, n, x, m, y) \
        {if (SpMatrix_Matvec((mat_obj), (n), (x), (m), (y))) return -1;}

int Itsolvers_minres_kernel(int n, double errtol, int it_max,
			    int *it, double *nrm_res, int clvl,
			    double *x, double *b, double *work,
			    PyObject *mat_obj,
			    PyObject *prec_obj) {
  int ONE = 1;

  double norm_r0, beta, beta_old, c, c_old, c_oold, s, s_old, s_oold, eta, norm_rmr, alpha, dconst1, dconst2, r1, r1_hat, r2, r3, tmp;
  int i;
  double *v, *v_hat, *v_hat_old, *av, *y, *w, *w_old;
  
  v_hat_old = work;
  v_hat = work + n;
  y = work + 2*n;
  w = work + 3*n;
  w_old = work + 4*n;
  v = work + 5*n;
  av = work + 6*n;
  
  *it = 0;
  /* v_hat_old=zeros(N,1); */
  for (i = 0; i < n; i ++)
    v_hat_old[i] = 0.0;
  /* v_hat = b - A*x0; */
  SpMatrix_MATVEC(mat_obj, n, x, n, v_hat);
  for (i = 0; i < n; i ++)
    v_hat[i] = b[i] - v_hat[i];
  /* norm_r0=norm(v_hat); */
  norm_r0 = F77(dnrm2)(&n, v_hat, &ONE);
  /* y = M\v_hat; */
  if (prec_obj)
    SpMatrix_PRECON(prec_obj, n, v_hat, y)
  else
    F77(dcopy)(&n, v_hat, &ONE, y, &ONE);
  /* beta=sqrt(v_hat'*y); beta_old=1; */
  beta = F77(ddot)(&n, v_hat, &ONE, y, &ONE);
  if (beta < 0.0)
    return -3;			/* preconditioner is not SPD */
  beta = sqrt(beta);
  beta_old = 1.0;
  
  /* c=1; c_old=1; s_old=0; s=0; */
  c = 1.0; c_old = 1.0; s = 0.0; s_old = 0.0;
  /* w=zeros(N,1); w_old=w; eta=beta; */
  for (i = 0; i < n; i ++)
    w[i] = 0.0;
  for (i = 0; i < n; i ++)
    w_old[i] = 0.0;
  eta = beta;
  
  /* xMR=x0; norm_rMR=norm_r0; */
  norm_rmr = norm_r0;
  
  while (1) {

    if (clvl >= 1) {
      if (*it == 0) {
	printf("MINRES.            Solution of symmetric  Ax = b\n"
	       "N      =%7d\n"
	       "IT_MAX =%7d     R_TOL =%11.2e\n\n", 
	       n, it_max, errtol*norm_r0);
	printf("      ITN             NORM(R)\n");
      }
      printf("    %5d %19.10e\n", *it, norm_rmr);
      if (*it % 10 == 0)
	printf("\n");
    }
    
    /*
     * Lanczos
     */
    if (*it >= it_max || norm_rmr < errtol*norm_r0)
      break;
    *it = *it + 1;

    /*
     * Lanczos
     */

    /* v=y/beta; y=v_hat; */
    for (i = 0; i < n; i ++)
      v[i] = y[i] / beta;
    F77(dcopy)(&n, v_hat, &ONE, y, &ONE);
    /* Av = A*v */
    SpMatrix_MATVEC(mat_obj, n, v, n, av);
    /* alpha=v'*Av; */
    alpha = F77(ddot)(&n, v, &ONE, av, &ONE);
    /* v_hat=Av-(alpha/beta)*v_hat-(beta/beta_old)*v_hat_old; */
    dconst1 = alpha/beta; dconst2 = beta/beta_old;
    for (i = 0; i < n; i ++)
      v_hat[i] = av[i] - dconst1*v_hat[i] - dconst2*v_hat_old[i];
    /* v_hat_old=y; */
    F77(dcopy)(&n, y, &ONE, v_hat_old, &ONE);
    /* y = M\v_hat; */
    if (prec_obj)
      SpMatrix_PRECON(prec_obj, n, v_hat, y)
    else
      F77(dcopy)(&n, v_hat, &ONE, y, &ONE);
    /* beta_old=beta; beta=sqrt(v_hat'*y); */
    beta_old=beta;
    beta = F77(ddot)(&n, v_hat, &ONE, y, &ONE);
    if (beta < 0.0)
      return -3;	       	/* preconditioner is not SPD */
    beta = sqrt(beta);
    
    /*
     * QR factorization
     */
    c_oold = c_old; c_old = c; s_oold = s_old; s_old = s;
    
    r1_hat = c_old*alpha - c_oold*s_old*beta_old;
    r1 = sqrt(r1_hat*r1_hat + beta*beta);
    r2 = s_old*alpha + c_oold*c_old*beta_old;
    r3 = s_oold*beta_old;
    
    /*
     * Givens rotation
     */
    if (r1 == 0.0)
      return -6;
    c = r1_hat/r1;
    s = beta/r1;

    /*
     * Update
     */
    /* w_oold=w_old; w_old=w; */
    /* w=(v-r3*w_oold-r2*w_old)/r1; */
    /* Vector w_oold eliminated */
    for (i = 0; i < n; i ++) {
      tmp = w[i];
      w[i] = (v[i] - r3*w_old[i] - r2*tmp)/r1;
      w_old[i] = tmp;
    }
    /* xMR=xMR+c*eta*w; eta=-s*eta; */
    dconst1 = c*eta;
    for (i = 0; i < n; i ++)
      x[i] += dconst1*w[i];
    eta = -s*eta;
    
    /*
     * Norm
     * on illconditioned problems the estimate of norm_rMR may
     * depart from the true norm_rMR and thus cause premature
     * termination.
     * To overcome this, norm_rMR can be computed from xMR:
     *     if isstr(A), norm_rMR = norm(b - feval(A,xMR,varargin{:})); 
     *     else norm_rMR=norm(b - A*xMR); end
     */
    norm_rmr *= fabs(s);  /* Updated norm w.r.t. preconditioned system */
  }
  
  *nrm_res = norm_rmr / norm_r0;
  if (norm_rmr < errtol*norm_r0)
    return 0;
  else
    return -1;
}
