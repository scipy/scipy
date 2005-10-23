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
*  Routines for Orthogonalization and B-Orthogonalization
*
*  $Id: orthopack.c,v 1.1.1.1 2003/11/09 12:35:47 geus Exp $
*
*
**************************************************************************/


/******************************************************************************
 *                                                                            *
 * ICGS -- Iterated Classical Gram-Schmidt Orthogonalization                  *
 *                                                                            *
 * Orthogonalizes u against the columns of Q.                                 *
 *                                                                            *
 ******************************************************************************/

static void
icgs(double *u, double *unrm, int n, int m, double *A, 
     double *um) {

  const int maxcgsit = 5;
  const double alpha = 0.5;

  double unrm_old;
  int i, isorth = 0;

  *unrm = F77(dnrm2)(&n, u, &ONE);

  if (m == 0)
    return;

  for (i = 0; !isorth && i < maxcgsit; i ++) {
    
    F77(dgemv)("t", &n, &m,  &DONE, A, &n, u, &ONE, &DZER, um, &ONE, 1);
    F77(dgemv)("n", &n, &m, &DMONE, A, &n, um, &ONE, &DONE, u, &ONE, 1);
    
    unrm_old = (*unrm);
    *unrm = F77(dnrm2)(&n, u, &ONE);
    
    isorth=((*unrm) > alpha*unrm_old);
  }
  if (i >= maxcgsit) {
    printf("warning: loss of orthogonality. ");
    printf("icgs() not converged after %d steps.\n", maxcgsit);
  }
}

/******************************************************************************
 *                                                                            *
 * ICGSM -- Iterated Classical M-orthogonal Gram-Schmidt                      *
 *                                                                            *
 * M-orthogonalizes u against the columns of Q.                               *
 *                                                                            *
 ******************************************************************************/

static void 
icgsm(double *u, double *unrm, int n, int m, double *Q, 
      PyObject *mmat,
      double *um, double *temp) {

  const int maxcgsit = 5;
  const double alpha = 0.5;

  double unrm_old;
  int ret, i, isorth = 0;

  ret = SpMatrix_Matvec(mmat, n, u, n, um);
  assert(ret == 0);
  *unrm = sqrt(F77(ddot)(&n, u, &ONE, um, &ONE));

  if (m == 0)
    return;
	  
  for (i = 0; !isorth && i < maxcgsit; i ++) {

    F77(dgemv)("t", &n, &m,  &DONE, Q, &n, um, &ONE, &DZER, temp, &ONE, 1);
    F77(dgemv)("n", &n, &m, &DMONE, Q, &n, temp, &ONE, &DONE, u, &ONE, 1);

    ret = SpMatrix_Matvec(mmat, n, u, n, um);
    assert(ret == 0);

    unrm_old = (*unrm);
    *unrm = sqrt(F77(ddot)(&n, u, &ONE, um, &ONE));

    isorth=((*unrm) > alpha*unrm_old);
  }
  if (i >= maxcgsit) {
    printf("warning: loss of orthogonality. ");
    printf("icgsm() not converged after %d steps.\n", maxcgsit);
  }
}



/******************************************************************************
 *                                                                            *
 * MGS -- Modified Gram-Schmidt                                               *
 *                                                                            *
 * Orthogonlaizes v with respect to span{A[:,1:m]}                            *
 *                                                                            *
 ******************************************************************************/

static void 
mgs(double *u, int n, int m, double *Q) {

  int i;
  double s;
  
  for (i = 0; i < m; i ++) {
    s = - F77(ddot)(&n, Q+i*n, &ONE, u, &ONE);
    F77(daxpy)(&n, &s, Q+i*n, &ONE, u, &ONE);
  }
}

/******************************************************************************
 *                                                                            *
 * MGSM -- Modified M-orthogonal GramSchmidt                                  *
 *                                                                            *
 * M-Orthogonalizes u against the columns of Q.                               *
 *                                                                            *
 ******************************************************************************/

static void 
mgsm(double *u, int n, int m, double *Q, double *QM) {

  int i;
  double s;

  for (i = 0; i < m; i ++) {
    s = - F77(ddot)(&n, QM+i*n, &ONE, u, &ONE);
    F77(daxpy)(&n, &s, Q+i*n, &ONE, u, &ONE);
  }
}
