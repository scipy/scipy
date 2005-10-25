#include "Python.h"
#include "pysparse/blas.h"
#include "pysparse/fortran.h"
#define SPMATRIX_UNIQUE_SYMBOL itsolvers_spmatrix
#include "pysparse/spmatrix.h"
#include "pysparse/pcg.h"

#define SpMatrix_PRECON(prec_obj, n, x, y) \
        {if (SpMatrix_Precon((prec_obj),(n),(x),(y))) return -1;}
#define SpMatrix_MATVEC(mat_obj, n, x, m, y) \
        {if (SpMatrix_Matvec((mat_obj), (n), (x), (m), (y))) return -1;}

/* function prototypes */
static void itermsg(double tol,
		    int maxit,
		    int flag,
		    int iter,
		    double relres);

/* PCG - Conjugate Gradients Algorithm
 */
int Itsolvers_pcg_kernel(int n, 
			 double *x, 
			 double *b,
			 double tol, 
			 int maxit,
			 int clvl,
			 int *iter, 
			 double *relres, 
			 int *flag,
			 double *work,
			 PyObject *mat_obj,
			 PyObject *prec_obj)
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
    return 0;
  }
  
  /* Set up for the method */
  *flag = -1;
  tolb = tol * n2b;		/* Relative tolerance */
  SpMatrix_MATVEC(mat_obj, n, x, n, r); /* Zero-th residual: r = b - A * x*/
  for (i = 0; i < n; i ++)	/* then  solution is all zeros */
    r[i] = b[i] - r[i];
  normr = F77(dnrm2)(&n, r, &ONE); /* Norm of residual */
  
  if (normr <= tolb) {		/* Initial guess is a good enough solution */
    *flag = 0;
    *relres = normr / n2b;
    *iter = 0;
    if (clvl)
      itermsg(tol,maxit,*flag,*iter,*relres);
    return 0;
  }

  rho = 1.0;
  stag = 0;			/* stagnation of the method */

  /* loop over maxit iterations (unless convergence or failure) */
  
  for (it = 1; it <= maxit; it ++) {
    
    if (prec_obj) {
      SpMatrix_PRECON(prec_obj, n, r, z);
    } else {
      F77(dcopy)(&n, r, &ONE, z, &ONE);
    }
   
    rho1 = rho;
    rho = F77(ddot)(&n, r, &ONE, z, &ONE);
    if (rho == 0.0) {		/* or isinf(rho) */
      *flag = -2;
      break;
    }
    if (it == 1) {
      F77(dcopy)(&n, z, &ONE, p, &ONE);
    } else {
      beta = rho / rho1;
      if (beta == 0.0) {	/* | isinf(beta) */
	*flag = -6;
	break;
      }
      for (i = 0; i < n; i ++)	/* p = z + beta * p; */
	p[i] = z[i] + beta * p[i];
    }
    SpMatrix_MATVEC(mat_obj, n, p, n, q); /* q = A * p */
    pq = F77(ddot)(&n, p, &ONE, q, &ONE); /* pq = p' * q */
    if (pq == 0.0) {		/* | isinf(pq) */
      *flag = -6;
      break;
    } else {
      alpha = rho / pq;
    }
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
    SpMatrix_MATVEC(mat_obj, n, x, n, z); /* normr = norm(b - A * x) */
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
      *flag = -5;
      break;
    }
  } /* for it = 1 : maxit */
  
  *iter = it;
  *relres = normr / n2b;

  if (clvl)
    itermsg(tol,maxit,*flag,*iter,*relres);
  return 0;
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
  case -1:
    printf("\nbecause the maximum number of iterations was reached.");
    break;
  case -2:
    printf("\nbecause the system involving the preconditioner was ill conditioned.");
    break;
  case -5:
    printf("\nbecause the method stagnated.");
    break;
  case -6:
    printf("\nbecause a scalar quantity became too small or too large to continue computing.");
    break;
  }
  
  if (flag != 0)
    printf("\nThe iterate returned (number %d) has relative residual %0.2g",iter,relres);

  printf("\n");
}
