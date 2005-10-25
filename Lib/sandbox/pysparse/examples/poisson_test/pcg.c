/*-------------------------------------------------------------------------

  PCG    Preconditioned Conjugate Gradients Method

  X = PCG(A,B) attempts to solve the system of linear equations A*X=B
  for X.  The coefficient matrix A must be symmetric and positive
  definite and the right hand side (column) vector B must have length
  N, where A is N-by-N.  PCG will start iterating from an initial
  guess. Iterates are produced until the method either converges,
  fails, or has computed the maximum number of iterations.
  Convergence is achieved when an iterate X has relative residual
  NORM(B-A*X)/NORM(B) less than or equal to the tolerance of the
  method. If PCG converged, then a message to that effect is
  displayed.  If PCG failed to converge after the maximum number of
  iterations or halted for any reason, then a message is printed
  displaying the relative residual and the iteration number at which
  the method stopped or failed.

  Parameters:

  N        problem size, A is an N-by-N matrix, B has length N
  X        on entry: the initial guess
           on exit:  the computed approximate solution
  B        the right hand side of the system
  TOL      error tolerance
  MAXIT    maximum number of iterations
  CLVL     verbosity of output (0: no output, 1: some output, 2: more output)
  ITER     number of iterations of iterations
  RELRES   norm of residual upon exit
  FLAG     exit code (see below)
  WORK     work array. must be of size 4*N
  MATVEC   function that performs multiplication with matrix A
  PRECON   function that applies preconditioner (can be NULL)

  FLAG describes the convergence of PCG.  If FLAG is

   0 then PCG converged to the desired tolerance TOL within MAXIT
     iterations without failing for any reason.

   1 then PCG iterated MAXIT times but did not converge.

   2 then a system of equations of the form M*Y = R was ill-conditioned.

   3 then PCG stagnated (two consecutive iterates were the same).

   4 then one of the scalar quantities calculated during PCG became
     too small or too large to continue computing.

  $Id: pcg.c,v 1.1.1.1 2003/11/09 12:35:47 geus Exp $
  Roman Geus, ETHZ

--------------------------------------------------------------------------*/

#include <assert.h>
#include "blas.h"
#include "fortran.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pcg.h"

/* function prototypes */
static void itermsg(double tol,
		    int maxit,
		    int flag,
		    int iter,
		    double relres);

/* PCG - Conjugate Gradients Algorithm
 */
void pcg(int n, 
	 double *x, 
	 double *b,
	 double tol, 
	 int maxit,
	 int clvl,
	 int *iter, 
	 double *relres, 
	 int *flag,
	 double *work,
	 void (*matvec)(double *, double *),
	 void (*precon)(double *, double *))
{
  double ALPHA;			/* used for passing parameters */
  int ONE = 1;			/* to BLAS routines */

  double n2b;			/* norm of rhs vector */
  double tolb;			/* requested tolerance for residual */
  double normr;			/* residual norm */
  double alpha, beta;
  double rho, rho1;
  double pq;
  double dmax, ddum;		/* used to detect stagnation */
  int stag;			/* flag to indicate stagnation */
  int it;			/* current iteration number */
  int i;			/* index variable */
  double *r, *z, *p, *q;	/* pointers to vectors in PCG algorithm */
  
  /* setup pointers into work */
  r = work;
  z = work + n;
  p = work + 2*n;
  q = work + 3*n;

  /* Check for all zero right hand side vector => all zero solution */
  n2b = F77(dnrm2)(&n, b, &ONE);/* Norm of rhs vector, b */
  if (n2b == 0.0) {		/* if rhs vector is all zeros */
    for (i = 0; i < n; i ++)	/* then  solution is all zeros */
      x[i] = 0.0;
    *flag = 0;			/* a valid solution has been obtained */
    *relres = 0.0;		/* the relative residual is actually 0/0 */
    *iter = 0;			/* no iterations need be performed */
    if (clvl)
      itermsg(tol,maxit,*flag,*iter,*relres);
    return;
  }
  
  /* Set up for the method */
  *flag = 1;
  tolb = tol * n2b;		/* Relative tolerance */
  matvec(x, r);			/* Zero-th residual: r = b - A * x*/
  for (i = 0; i < n; i ++)	/* then  solution is all zeros */
    r[i] = b[i] - r[i];
  normr = F77(dnrm2)(&n, r, &ONE); /* Norm of residual */
  
  if (normr <= tolb) {		/* Initial guess is a good enough solution */
    *flag = 0;
    *relres = normr / n2b;
    *iter = 0;
    if (clvl)
      itermsg(tol,maxit,*flag,*iter,*relres);
    return;
  }

  rho = 1.0;
  stag = 0;			/* stagnation of the method */

  /* loop over maxit iterations (unless convergence or failure) */
  
  for (it = 1; it <= maxit; it ++) {
    
    if (precon) {
      precon(r, z);
      /*
	if isinf(norm(y,inf))
	flag = 2;
	break
	end
      */
    } else {
      F77(dcopy)(&n, r, &ONE, z, &ONE);
    }
   
    rho1 = rho;
    rho = F77(ddot)(&n, r, &ONE, z, &ONE);
    if (rho == 0.0) {		/* or isinf(rho) */
      *flag = 4;
      break;
    }
    if (it == 1) {
      F77(dcopy)(&n, z, &ONE, p, &ONE);
    } else {
      beta = rho / rho1;
      if (beta == 0.0) {	/* | isinf(beta) */
	*flag = 4;
	break;
      }
      for (i = 0; i < n; i ++)	/* p = z + beta * p; */
	p[i] = z[i] + beta * p[i];
    }
    matvec(p, q);		/* q = A * p */
    pq = F77(ddot)(&n, p, &ONE, q, &ONE); /* pq = p' * q */
    if (pq == 0.0) {		/* | isinf(pq) */
      *flag = 4;
      break;
    } else {
      alpha = rho / pq;
    }
    /* 
       if isinf(alpha)
       flag = 4;
       break
       end
    */
    if (alpha == 0.0)		/* stagnation of the method */
      stag = 1;
   
    /* Check for stagnation of the method */
    if (stag == 0) {
      dmax = 0.0;
      for (i = 0; i < n; i ++)
	if (x[i] != 0.0) {
	  ddum = fabs(alpha * p[i]/x[i]);
	  if (ddum > dmax)
	    dmax = ddum;
	} else
	  if (p[i] != 0.0)
	    dmax = 1.0;
      stag = (1.0 + dmax == 1.0);
    }
    
    F77(daxpy)(&n, &alpha, p, &ONE, x, &ONE); /* form new iterate */
    ALPHA = -alpha;
    F77(daxpy)(&n, &ALPHA, q, &ONE, r, &ONE); /* r = r - alpha * q */
    
    /* check for convergence */
#ifdef EXPENSIVE_CRIT
    matvec(x, z);		/* normr = norm(b - A * x) */
    for (i = 0; i < n; i ++)
      z[i] = b[i] - z[i];
    normr = F77(dnrm2)(&n, z, &ONE);
#else
    normr = F77(dnrm2)(&n, r, &ONE); /* normr = norm(r) */
#endif
    if (normr <= tolb) {
      *flag = 0;
      break;
    }
    
    if (stag == 1) {
      *flag = 3;
      break;
    }
  } /* for it = 1 : maxit */
  
  *iter = it;
  *relres = normr / n2b;

  if (clvl)
    itermsg(tol,maxit,*flag,*iter,*relres);
}

/* PCG_F77 - Fortran interface to PCG
 */
void F77(pcg_f77)(int *n, 
		  double *x, 
		  double *b,
		  double *tol, 
		  int *maxit,
		  int *clvl,
		  int *iter, 
		  double *relres, 
		  int *flag,
		  double *work,
		  void (*matvec)(double *, double *),
		  void (*precon)(double *, double *)) {
  pcg(*n,
      x,
      b,
      *tol,
      *maxit,
      *clvl,
      iter,
      relres,
      flag,
      work,
      matvec,
      precon);
}

  
/* ITERMSG - Displays the final message for PCG method
 */
static void itermsg(double tol,
		    int maxit,
		    int flag,
		    int iter,
		    double relres) {
  if (flag != 0) {
    printf("PCG stopped at iteration %d without converging to the desired tolerance %0.2g", iter, tol);
  }
  
  switch(flag) {
  case 0:
    if (iter == 0)
      printf("The initial guess has relative residual %0.2g which is within\nthe desired tolerance %0.2g so PCG returned it without iterating.",
	     relres, tol);
    else
      printf("PCG converged at iteration %d to a solution with relative residual %0.2g", iter, relres);
    break;
  case 1:
    printf("\nbecause the maximum number of iterations was reached.");
    break;
  case 2:
    printf("\nbecause the system involving the preconditioner was ill conditioned.");
    break;
  case 3:
    printf("\nbecause the method stagnated.");
    break;
  case 4:
    printf("\nbecause a scalar quantity became too small or too large to continue computing.");
    break;
  }
  
  if (flag != 0)
    printf("\nThe iterate returned (number %d) has relative residual %0.2g",iter,relres);

  printf("\n");
}
