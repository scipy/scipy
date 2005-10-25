#include <math.h>
#include "Python.h"
#include "pysparse/blas.h"
#define SPMATRIX_UNIQUE_SYMBOL itsolvers_spmatrix
#include "pysparse/spmatrix.h"
#include "pysparse/cgs.h"

#define SpMatrix_PRECON(prec_obj, n, p, tmp) \
        {if (SpMatrix_Precon((prec_obj),(n),(p),(tmp))) return -1;}
#define SpMatrix_MATVEC(mat_obj, n, x, m, y) \
        {if (SpMatrix_Matvec((mat_obj), (n), (x), (m), (y))) return -1;}

int Itsolvers_cgs_kernel(int n, 
			 double *b, 
			 double *x, 
			 int maxit, 
			 double tol,
			 double *work,
			 int *iter,
			 double *res,
			 PyObject *mat_obj,
			 PyObject *prec_obj) {
  
  double *r0, *r, *p, *q, *u, *v, *tmp, *tmp2, alpha, beta, rho, rho_new, tol_sq, bnrm_sq, ddummy;

  double DMONE = -1.0, DONE = 1.0;
  int ONE = 1;

  /* Place vectors in WorkSpace */
  r0 = work;
  r = work + n;
  p = work + 2*n;
  q = work + 3*n;
  u = work + 4*n;
  v = work + 5*n;
  tmp = work + 6*n;
  tmp2 = work + 7*n;

  /* For ease */
  tol_sq = tol*tol;

  /* Initializing Values : r0 = b - A*x */
  *iter = 0;
  SpMatrix_MATVEC(mat_obj, n, x, n, tmp);
  F77(dcopy)(&n, b, &ONE, r0, &ONE);
  F77(daxpy)(&n, &DMONE, tmp, &ONE, r0, &ONE);

  /* p := u := r := r0 */
  F77(dcopy)(&n, r0, &ONE, r, &ONE);
  F77(dcopy)(&n, r0, &ONE, u, &ONE);
  F77(dcopy)(&n, r0, &ONE, p, &ONE);

  /* rho = r0'*r0 , bnrm_sq = b'*b */
  rho = F77(ddot)(&n, r0, &ONE, r0, &ONE);
  bnrm_sq = F77(ddot)(&n, b, &ONE, b, &ONE);

  /* Starting-vector is already good enough */
  if (rho < bnrm_sq*tol_sq) {
    *res = sqrt(rho / bnrm_sq);
    return 0;
  }

  /* Iterate for at most maxit steps ...*/
  for (; *iter < maxit; (*iter) ++) {

    if (prec_obj != NULL) {
      SpMatrix_PRECON(prec_obj, n, p, tmp);
      SpMatrix_MATVEC(mat_obj, n, tmp, n, v);
    } else
      SpMatrix_MATVEC(mat_obj, n, p, n, v);

    alpha = rho/F77(ddot)(&n, v, &ONE, r0, &ONE);         /* alpha = rho/(v'*inv(K)*A*r0) */

    ddummy = -alpha;                                      /* ddummy = -alpha */
    F77(dcopy)(&n, u, &ONE, q, &ONE);
    F77(daxpy)(&n, &ddummy, v, &ONE, q, &ONE);            /* q = u - alpha*v */

    F77(dcopy)(&n, u, &ONE, tmp, &ONE);
    F77(daxpy)(&n, &DONE, q, &ONE, tmp, &ONE);            /* tmp = u + q */
    if (prec_obj != NULL)	                          /* x = x + alpha*tmp2 = x + alpha*inv(K)*(u + q) */
      SpMatrix_PRECON(prec_obj, n, tmp, tmp2)
    else
      F77(dcopy)(&n, tmp, &ONE, tmp2, &ONE);
    F77(daxpy)(&n, &alpha, tmp2, &ONE, x, &ONE);        

    SpMatrix_MATVEC(mat_obj, n, tmp2, n, tmp);
    F77(daxpy)(&n, &ddummy, tmp, &ONE, r, &ONE);          /* r = r - alpha*A*tmp2 = r - alpha*A*inv(K)*(u + q) */

    *res = F77(ddot)(&n, r, &ONE, r, &ONE);
    if (*res < bnrm_sq*tol_sq) {                          /* Test covergence ... */
      *res = sqrt(*res / bnrm_sq);
      return 0;
    }

    rho_new = F77(ddot)(&n, r, &ONE, r0, &ONE);           /* rho_new = r_(i+1)'*r0 */
    beta = rho_new/rho;                                   /* beta = rho_new/rho */
    rho = rho_new;                                        /* rho = rho_new */

    F77(dcopy)(&n, r, &ONE, u, &ONE);
    F77(daxpy)(&n, &beta, q, &ONE, u, &ONE);              /* u = r + beta*q */

    F77(dcopy)(&n, q, &ONE, tmp, &ONE);
    F77(daxpy)(&n, &beta, p, &ONE, tmp, &ONE);
    F77(dcopy)(&n, u, &ONE, p, &ONE);
    F77(daxpy)(&n, &beta, tmp, &ONE, p, &ONE);            /* p = u + beta*(q + beta*p) */
  }

  /* CGS did not converge */
  *res = sqrt(*res / bnrm_sq) ;
  return -1;
}
