#include "Python.h"
#include "pysparse/blas.h"
#include "pysparse/fortran.h"
#define SPMATRIX_UNIQUE_SYMBOL itsolvers_spmatrix
#include "pysparse/spmatrix.h"
#include "pysparse/bicgstab.h"

/* function prototypes */
static void itermsg(double tol,
		    int maxit,
		    int info,
		    int iter,
		    double relres);

void v_plus_cw(int n,
	       double *v,
	       double *w,
	       double c,
	       double *z);


#define SpMatrix_PRECON(prec_obj, n, x, y) \
        {if (SpMatrix_Precon((prec_obj),(n),(x),(y))) return -1;}
#define SpMatrix_MATVEC(mat_obj, n, x, m, y) \
        {if (SpMatrix_Matvec((mat_obj), (n), (x), (m), (y))) return -1;}


/* BICGSTAB - stabilized BiConjugate Gradients algorithm
 */
int Itsolvers_bicgstab_kernel2(int n, 
			      double *x, 
			      double *b,
			      double tol, 
			      int maxit,
			      int clvl,
			      int *iter, 
			      double *relres, 
			      int *info,
			      double *work,
			      PyObject *mat_obj,
			      PyObject *prec_obj)
{
  int ONE = 1;			/* to BLAS routines */
  int i;

  double n2b; /* 2 nrm of residual */
  double *r,*rhat,*p,*v,*w,*z,*bb; 
  double c,rho,old_rho,alpha,beta,omega,old_omega;
  double res,res0;

  *info = -6;

  /* Check for all zero right hand side vector => all zero solution */
  n2b = F77(dnrm2)(&n, b, &ONE);/* Norm of rhs vector, b */
  if (n2b == 0.0) {		/* if rhs vector is all zeros */
    for (i = 0; i < n; i ++)	/* then  solution is all zeros */
      x[i] = 0.0;
    *info = 0;			/* a valid solution has been obtained */
    *relres = 0.0;		/* the relative residual is actually 0/0 */
    *iter = 0;			/* no iterations need be performed */
    if (clvl)
      itermsg(tol,maxit,*info,*iter,*relres);
    return(0);
  }

  /* setup pointers into work */
  r = work;
  rhat = work + n;
  p = work + 2*n;
  v = work + 3*n;
  w = work + 4*n;
  z = work + 5*n;
  bb = work + 6*n;

  c = -1.0;
  old_rho = 1.0;
  alpha = 1.0;
  old_omega = 1.0;

  printf("initial solution norm in bicgstab: %e\n", F77(dnrm2)(&n, x, &ONE));

  
  if (prec_obj) {
    SpMatrix_PRECON(prec_obj, n, b, bb);
  } else {
    F77(dcopy)(&n, b, &ONE, bb, &ONE);
  }

  /* compute initial residual */
  SpMatrix_MATVEC(mat_obj,n,x,n,w);
  if (prec_obj) {
    SpMatrix_PRECON(prec_obj, n, w, z);
  } else {
    F77(dcopy)(&n, w, &ONE, z, &ONE);
  }

  v_plus_cw(n,bb,z,c,r);
  F77(dcopy)(&n, r, &ONE, rhat, &ONE);

  res0 = F77(dnrm2)(&n, bb, &ONE); 
  printf("initial residual in bicgstab: %e\n", res0);

  *iter = 0;
  do {

    (*iter)++;
    rho = F77(ddot)(&n, r, &ONE, rhat, &ONE);
    beta = (rho/old_rho)*(alpha/old_omega);
#if 0
    printf("iter: %d rho = %e, beta = %e\n",*iter,rho,beta);
#endif

    /* compute new p */
    v_plus_cw(n,p,v,-old_omega,z);
    v_plus_cw(n,r,z,beta,p);

    /* compute new v, r, and alpha */
    SpMatrix_MATVEC(mat_obj,n,p,n,w);
    if (prec_obj) {
      SpMatrix_PRECON(prec_obj, n, w, v);
    } else {
      F77(dcopy)(&n, w, &ONE, v, &ONE);
    }

    alpha = rho/F77(ddot)(&n, rhat, &ONE, v, &ONE);
    v_plus_cw(n,r,v,-alpha,w);
    F77(dcopy)(&n, w, &ONE, r, &ONE);

    /* compute new omega */
    SpMatrix_MATVEC(mat_obj,n,r,n,w);
    if (prec_obj) {
      SpMatrix_PRECON(prec_obj, n, w, z);
    } else {
      F77(dcopy)(&n, w, &ONE, z, &ONE);
    }
    omega = F77(ddot)(&n, z, &ONE, r, &ONE)/F77(ddot)(&n, z, &ONE, z, &ONE);

    /* compute new x and new r */
    v_plus_cw(n,x,p,alpha,w);
    v_plus_cw(n,w,r,omega,x);
    v_plus_cw(n,r,z,-omega,w);
    F77(dcopy)(&n, w, &ONE, r, &ONE);
    old_rho = rho;
    old_omega = omega;

    /* compute exact residual -> w */
    SpMatrix_MATVEC(mat_obj,n,x,n,w);
    if (prec_obj) {
      SpMatrix_PRECON(prec_obj, n, w, z);
    } else {
      F77(dcopy)(&n, w, &ONE, z, &ONE);
    }
    v_plus_cw(n,bb,z,c,w);
    res = F77(dnrm2)(&n, w, &ONE); 

  }  while ((res/res0  > tol) && (*iter < maxit));

  *relres = res/res0;

  if (*relres>=tol) {
    *info = -1;
  } else {
    *info = 0;
  }

  if (clvl)
    itermsg(tol,maxit,*info,*iter,*relres);
  
  return(0);

}








/* BICGSTAB - stabilized BiConjugate Gradients algorithm
 */
int Itsolvers_bicgstab_kernel(int n, 
			      double *x, 
			      double *b,
			      double tol, 
			      int maxit,
			      int clvl,
			      int *iter, 
			      double *relres, 
			      int *info,
			      double *work,
			      PyObject *mat_obj,
			      PyObject *prec_obj)
{
  int ONE = 1;			/* to BLAS routines */
  int i;

  double n2b; /* 2 nrm of residual */
  double *r,*rhat,*p,*phat,*v,*s,*shat,*t; 
  double alpha,omega,rho_im1,rho_im2,beta;
  double res,res0;

  *info = -6;

  /* Check for all zero right hand side vector => all zero solution */
  n2b = F77(dnrm2)(&n, b, &ONE);/* Norm of rhs vector, b */
  if (n2b == 0.0) {		/* if rhs vector is all zeros */
    for (i = 0; i < n; i ++)	/* then  solution is all zeros */
      x[i] = 0.0;
    *info = 0;			/* a valid solution has been obtained */
    *relres = 0.0;		/* the relative residual is actually 0/0 */
    *iter = 0;			/* no iterations need be performed */
    if (clvl)
      itermsg(tol,maxit,*info,*iter,*relres);
    return(0);
  }

  /* setup pointers into work */
  r = work;
  rhat = work + n;
  p = work + 2*n;
  phat = work + 3*n;
  v = work + 4*n;
  s = work + 5*n;
  shat = work + 6*n;
  t = work + 7*n;

  omega   = 0.0;
  beta    = 0.0;
  alpha   = 0.0;
  rho_im2 = 0.0;

  /* compute residual */

  SpMatrix_MATVEC(mat_obj,n,x,n,r);
  for (i=0;i<n;i++) r[i] = b[i] + -r[i];
  res0 = F77(dnrm2)(&n, r, &ONE);

  F77(dcopy)(&n, r, &ONE, rhat, &ONE);

  *iter = 0;

  do {

    (*iter)++;
    
    rho_im1 = F77(ddot)(&n, rhat, &ONE, r, &ONE);
    if (rho_im1==0.0) return(-1);
    
    if (*iter==1) {
      F77(dcopy)(&n, r, &ONE, p, &ONE);
    } else {
      beta = (rho_im1/rho_im2)*(alpha/omega);
      for (i=0;i<n;i++)	p[i] = r[i] + beta*(p[i]-omega*v[i]);
    }

    if (prec_obj) {
      SpMatrix_PRECON(prec_obj, n, p, phat);
    } else {
      F77(dcopy)(&n, p, &ONE, phat, &ONE);
    }
    
    SpMatrix_MATVEC(mat_obj,n,phat,n,v);
    alpha = rho_im1/F77(ddot)(&n, rhat, &ONE, v, &ONE);
    v_plus_cw(n,r,v,-alpha,s);

    if (prec_obj) {
      SpMatrix_PRECON(prec_obj, n, s, shat);
    } else {
      F77(dcopy)(&n, s, &ONE, shat, &ONE);
    }
    
    SpMatrix_MATVEC(mat_obj,n,shat,n,t);
    omega = F77(ddot)(&n, t, &ONE, s, &ONE)/F77(ddot)(&n, t, &ONE, t, &ONE);
    for (i=0;i<n;i++) x[i] = x[i] + alpha*phat[i] + omega*shat[i];
    for (i=0;i<n;i++) r[i] = s[i] - omega*t[i];
    
    res = F77(dnrm2)(&n, r, &ONE); 

    if (omega==0.0) return(-1);    

    rho_im2 = rho_im1;

  }  while ((res/res0  > tol) && (*iter < maxit));


  *relres = res/res0;

  if (*relres>=tol) {
    *info = -1;
  } else {
    *info = 0;
  }

  if (clvl)
    itermsg(tol,maxit,*info,*iter,*relres);

  return(0);

}
















/* ITERMSG - Displays the final message for BICGSTAB method
 */
static void itermsg(double tol,
		    int maxit,
		    int info,
		    int iter,
		    double relres) {
  if (info != 0) {
    printf("BICGSTAB stopped at iteration %d without converging to the desired tolerance %0.2g", iter, tol);
  }
  
  switch(info) {
  case 0:
    if (iter == 0)
      printf("The initial guess has relative residual %0.2g which is within\nthe desired tolerance %0.2g so BICGSTAB returned it without iterating.",
	     relres, tol);
    else
      printf("BICGSTAB converged at iteration %d to a solution with relative residual %0.2g", iter, relres);
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
  
  if (info != 0)
    printf("\nThe iterate returned (number %d) has relative residual %0.2g",iter,relres);

  printf("\n");
}


void v_plus_cw(int n,
	       double *v,
	       double *w,
	       double c,
	       double *z)
{
  int i;
  for (i=0; i<n; i++)
    z[i] = v[i] + c*w[i];
}
