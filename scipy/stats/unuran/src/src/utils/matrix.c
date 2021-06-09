/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: matrix.c                                                          *
 *                                                                           *
 *   Routines for computations with square matrices.                         *
 *   Matrices are stored as double array (i.e.: double *matrix).             *
 *   The rows of the matrix have to be stored consecutively in this array,   *
 *   i.e. the entry with index [i,j] can be entered via (i*dim+j), where     *
 *   dim is the dimension of the dim x dim - matrix.                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Routines are mainly adapted from the GSL (GNU Scientifiy Library)       *
 *   Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough  *
 *                                                                           *
 *   adapted by Wolfgang Hoermann and Josef Leydold                          *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/
				                                                                                    
/*--------------------------------------------------------------------------*/

#include <unur_source.h>
#include "matrix_source.h"

/*---------------------------------------------------------------------------*/

static int _unur_matrix_swap_rows (int dim, double *A, int i, int j);
/* Swap rows i and j in square matrix A of dimension dim.                    */

static int _unur_matrix_permutation_swap (int dim, int *p, int i, int j);
/* Swap entries i and j of the integer array p.                              */

static int _unur_matrix_LU_decomp (int dim, double *A, int *P, int *signum);
/* Factorise a general dim x dim matrix A into P A = L U.                    */

static int _unur_matrix_backsubstitution_dtrsv(int dim, double *LU, double *X);
/* Backsubstitution used for inversion alg. _unur_matrix_LU_invert().        */

static int _unur_matrix_forwardsubstitution_dtrsv(int dim, double *LU, double *X);
/* Forwardsubstitution used for inversion alg. _unur_matrix_LU_invert().     */

static int _unur_matrix_LU_invert (int dim, double *LU, int *p, double *inverse);
/* Compute inverse of matrix with given LU decomposition.                    */

/*---------------------------------------------------------------------------*/

int
_unur_matrix_transform_diagonal (int dim, const double *M, const double *D, double *res)
     /*----------------------------------------------------------------------*/
     /* Computes the transformation M^t . D . M of diagonal matrix D         */
     /* and stores it in matrix 'res'.                                       */
     /*                                                                      */
     /* input:								     */
     /*    dim  ... number of columns and rows of m and diag                 */
     /*    M    ... square matrix                                            */
     /*    D    ... diagonal entries of a diagonal matrix                    */
     /*                                                                      */
     /* output:								     */
     /*   res   ... square matrix to store result                            */
     /*                                                                      */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   error code      otherwise                                          */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))

  int i,j,k;
  double sum;

  /* check arguments */
  CHECK_NULL(M,UNUR_ERR_NULL);
  CHECK_NULL(D,UNUR_ERR_NULL);
  CHECK_NULL(res,UNUR_ERR_NULL);

  for (i=0; i<dim; i++)
    for (j=0; j<dim; j++) {
      for (sum=0., k=0; k<dim; k++)
	sum += D[k] * M[idx(k,i)] * M[idx(k,j)];
      res[idx(i,j)] = sum;
    }
  
  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_matrix_transform_diagonal() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_swap_rows (int dim, double *A, int i, int j)
     /*----------------------------------------------------------------------*/
     /* Swap rows i and j in square matrix A of dimension dim.               */
     /* The resulting matrix is again stored in A.                           */ 
     /*                                                                      */
     /* input:								     */
     /*   dim ... number of columns and rows of				     */
     /*   A   ... dim x dim -matrix					     */
     /*   i,j ... row numbers that should be swapped			     */
     /*                                                                      */
     /* output:								     */
     /*   the given matrix A contains the new matrix with swapped rows       */
     /*                                                                      */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   error code      otherwise                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(A,UNUR_ERR_NULL);

  if (i != j)
  {
    double *row1 = A + i * dim;
    double *row2 = A + j * dim;
      
    int k;
      
    for (k = 0; k < dim; k++)
    {
       double tmp = row1[k] ;
       row1[k] = row2[k] ;
       row2[k] = tmp ;
    }
  }

  return UNUR_SUCCESS;
} /* end of _unur_matrix_swap_rows() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_permutation_swap (int dim ATTRIBUTE__UNUSED, int *p, int i, int j)
     /*----------------------------------------------------------------------*/
     /* Swap entries i and j of the integer array p.                         */
     /*                                                                      */
     /* input:							             */
     /*   dim ... number of elements				             */
     /*   p   ... permutation vector of length dim		             */
     /*           (contains all integers between 0 and dim-1) 	             */
     /*   i,j ... indexes that should be swapped		             */
     /*                                                                      */
     /* output:							             */
     /*   permutation p with swapped elements  			             */
     /*----------------------------------------------------------------------*/
{
  if (i != j)
    {
      int tmp = p[i];
      p[i] = p[j];
      p[j] = tmp;
    }
  
  return UNUR_SUCCESS;
} /* end of _unur_matrix_permutation_swap() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_LU_decomp (int dim, double *A, int *p, int *signum)
     /*----------------------------------------------------------------------*/
     /* Factorise a general dim x dim matrix A into,			     */
     /*									     */
     /*   P A = L U							     */
     /*									     */
     /* where P is a permutation matrix, L is unit lower triangular and U    */
     /* is upper triangular.						     */
     /*									     */
     /* L is stored in the strict lower triangular part of the input	     */
     /* matrix. The diagonal elements of L are unity and are not stored.     */
     /*									     */
     /* U is stored in the diagonal and upper triangular part of the	     */
     /* input matrix.  							     */
     /* 								     */
     /* P is stored in the permutation p. Column j of P is column k of the   */
     /* identity matrix, where k = permutation->data[j]			     */
     /*									     */
     /* signum gives the sign of the permutation, (-1)^n, where n is the     */
     /* number of interchanges in the permutation. 			     */
     /*									     */
     /* See Golub & Van Loan, Matrix Computations, Algorithm 3.4.1 (Gauss    */
     /* Elimination with Partial Pivoting).				     */
     /*									     */
     /* input:								     */
     /*   dim    ... number of columns and rows of			     */
     /*   A      ... dim x dim matrix					     */
     /*   p      ... pointer to array where permutation should be stored     */
     /*   signum ... pointer where sign of permutation should be stored      */
     /*									     */
     /* output:								     */
     /*   A      ... L (strict lower triangular part) and U (upper part)     */
     /*   p      ... permutation vector of length dim			     */
     /*              (contains all integers between 0 and dim-1) 	     */
     /*   signum ... sign of the permutation (1 or -1)			     */
     /*                                                                      */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   error code      otherwise                                          */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))
  int i, j, k;

  /* check arguments */
  CHECK_NULL(A,UNUR_ERR_NULL);
  CHECK_NULL(p,UNUR_ERR_NULL);
  CHECK_NULL(signum,UNUR_ERR_NULL);

  /* initialize permutation vector */
  *signum = 1;
  for(i=0;i<dim;i++) p[i]=i;

  for (j = 0; j < dim - 1; j++){
    /* Find maximum in the j-th column */

    double ajj, max = fabs (A[idx(j,j)]);
    int i_pivot = j;
    for (i = j + 1; i < dim; i++){
      double aij = fabs (A[idx(i,j)]);
      if (aij > max){
	max = aij;
	i_pivot = i;
      }
    }
    if (i_pivot != j){
      _unur_matrix_swap_rows (dim, A, j, i_pivot);
      _unur_matrix_permutation_swap (dim, p, j, i_pivot);
      *signum = -(*signum);
    }

    ajj = A[idx(j,j)];

    if ( !_unur_iszero(ajj) ){
      for (i = j + 1; i < dim; i++){
	double aij = A[idx(i,j)] / ajj;
	A[idx(i,j)] = aij;
	for (k = j + 1; k < dim; k++){
	  double aik = A[idx(i,k)];
	  double ajk = A[idx(j,k)];
	  A[idx(i,k)]= aik - aij * ajk;
	}
      }
    }
  }
      
  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_matrix_LU_decomp() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_backsubstitution_dtrsv(int dim, double *LU, double *X)
     /*----------------------------------------------------------------------*/
     /* Backsubstitution used for inversion alg. _unur_matrix_LU_invert()    */
     /*                                                                      */
     /* input:								     */
     /*   dim ... number of columns and rows of				     */
     /*   LU  ... dim x dim -matrix that contains LU decomposition of matrix */
     /*   X   ... vector 						     */
     /*                                                                      */
     /* output:								     */
     /*   X 								     */
     /*                                                                      */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   error code      otherwise                                          */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))

  int ix,jx,i,j;

  /* check arguments */
  CHECK_NULL(LU,UNUR_ERR_NULL);
  CHECK_NULL(X,UNUR_ERR_NULL);

  /* backsubstitution */
  ix = (dim - 1);
   
  X[ix] = X[ix] / LU[idx(ix,ix)];
 
  ix--;
  for (i = dim - 1; i > 0 && i--;) {
    double tmp = X[ix];
    jx = ix + 1;
    for (j = i + 1; j < dim; j++) {
      tmp -= LU[idx(i,j)] * X[jx];
      jx ++;
    }
      
    X[ix] = tmp / LU[idx(i,i)];
    ix --;
  }

  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_matrix_backsubstitution_dtrsv() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_forwardsubstitution_dtrsv(int dim, double *LU, double *X)
     /*----------------------------------------------------------------------*/
     /* Forwardsubstitution used for inversion alg. _unur_matrix_LU_invert() */
     /*                                                                      */
     /* input:								     */
     /*   dim ... number of columns and rows of				     */
     /*   LU  ... dim x dim -matrix that contains LU decomposition of matrix */
     /*   X   ... vector 						     */
     /*                                                                      */
     /* output:								     */
     /*   X 								     */
     /*                                                                      */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   error code      otherwise                                          */
     /*----------------------------------------------------------------------*/
{ 
#define idx(a,b) ((a)*dim+(b))

  int ix,jx,i,j;

  /* check arguments */
  CHECK_NULL(LU,UNUR_ERR_NULL);
  CHECK_NULL(X,UNUR_ERR_NULL);

  /* forward substitution */
  ix = 0;
  ix++;
  for (i = 1; i < dim; i++) {
    double tmp = X[ix];
    jx = 0;
    for (j = 0; j < i; j++) {
      tmp -= LU[idx(i,j)] * X[jx];
      jx += 1;
    }
    X[ix] = tmp;
    ix ++;
  }
  
  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_matrix_forwardsubstitution_dtrsv() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_LU_invert (int dim, double *LU, int *p, double *inverse)
     /*----------------------------------------------------------------------*/
     /* Compute inverse of matrix with given LU decomposition.               */
     /*                                                                      */
     /* input:								     */
     /*   dim ... number of columns and rows of				     */
     /*   LU  ... dim x dim -matrix, LU decomposition			     */
     /*   p   ... permutation vector of length dim			     */
     /*           (contains all integers between 0 and dim-1)     	     */
     /*   inverse ... pointer to array where inverse should be stored        */
     /*                                                                      */
     /* output								     */
     /*   inverse ... the inverse matrix 				     */
     /*                                                                      */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   error code      otherwise                                          */
     /*----------------------------------------------------------------------*/
{ 
#define idx(a,b) ((a)*dim+(b))

  double *vector;
  int i,j;

  /* check arguments */
  CHECK_NULL(LU,UNUR_ERR_NULL);
  CHECK_NULL(p,UNUR_ERR_NULL);
  CHECK_NULL(inverse,UNUR_ERR_NULL);

  /* allocate working array */
  vector = _unur_xmalloc(dim*sizeof(double));

  for (i = 0; i < dim; i++){
    for(j=0;j<dim;j++) vector[j] = 0.;
    vector[i] = 1.;
    /* Solve for c using forward-substitution, L c = P b */
    _unur_matrix_forwardsubstitution_dtrsv (dim, LU, vector);
    /* Perform back-substitution, U x = c */
    _unur_matrix_backsubstitution_dtrsv (dim, LU, vector);

    for(j=0;j<dim;j++){
      inverse[idx(j,p[i])] = vector[j]; /* set column vector of inverse matrix */
    }
  }

  /* free working arrays */
  free(vector);

  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_matrix_LU_invert() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_invert_matrix(int dim, const double *A, double *Ainv, double *det)
     /*----------------------------------------------------------------------*/
     /* Calculates the inverse matrix (by means of LU decomposition).        */
     /* As a side effect det(A) is comuted.                                  */
     /*									     */
     /* input:                                                               */
     /*   dim    ... dimension of the square matrix A                        */
     /*   A      ... dim x dim -matrix                                       */
     /*   Ainv   ... pointer to array where inverse matrix should be stored  */
     /*   det    ... pointer where det(A) should be stored                   */
     /*                                                                      */
     /* output:                                                              */
     /*   Ainv   ... inverse matrix of A	                             */
     /*   det    ... determinant of A                                        */
     /*									     */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   other error code, otherwise                                        */
     /*----------------------------------------------------------------------*/
{ 
#define idx(a,b) ((a)*dim+(b))

  int *p, s, i;
  double *LU;             /* array for storing LU decomposition of matrix A */
  
  /* check arguments */
  CHECK_NULL(A,UNUR_ERR_NULL);
  CHECK_NULL(Ainv,UNUR_ERR_NULL);
  CHECK_NULL(det,UNUR_ERR_NULL);
  if (dim<1) {
    _unur_error("matrix",UNUR_ERR_GENERIC,"dimension < 1");
    return UNUR_ERR_GENERIC;
  }
  
  /* allocate working space */
  p = _unur_xmalloc(dim*sizeof(int));
  LU = _unur_xmalloc(dim*dim*sizeof(double));
  
  /* compute LU decomposition */
  memcpy(LU, A, dim*dim*sizeof(double));
  _unur_matrix_LU_decomp(dim, LU, p, &s);
  
  /* compute determinant */
  *det = s;
  for(i=0;i<dim;i++)
    *det *= LU[idx(i,i)];
      
  /* compute inverse by means of LU factors */
  _unur_matrix_LU_invert(dim, LU, p, Ainv);   
  
  /* free working space */
  free(LU);
  free(p);
  
  return UNUR_SUCCESS;
  
#undef idx
} /* end of _unur_matrix_invert_matrix() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_multiplication(int dim, const double *A, const double *B, double *AB)
     /*----------------------------------------------------------------------*/
     /* Calculates the matrix multiplication of two matrices A and B         */
     /*									     */
     /* input:                                                               */
     /*   dim    ... dimension of the square matrix A                        */
     /*   A      ... dim x dim -matrix                                       */
     /*   B      ... dim x dim -matrix                                       */
     /*                                                                      */
     /* output:                                                              */
     /*   AB     ... A*B                                                     */
     /*									     */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   other error code, otherwise                                        */
     /*----------------------------------------------------------------------*/
{ 
#define idx(a,b) ((a)*dim+(b))

  int i, j, k;
  
  /* check arguments */
  CHECK_NULL(A,UNUR_ERR_NULL);
  CHECK_NULL(B,UNUR_ERR_NULL);
  CHECK_NULL(AB,UNUR_ERR_NULL);
  
  if (dim<1) {
    _unur_error("matrix",UNUR_ERR_GENERIC,"dimension < 1");
    return UNUR_ERR_GENERIC;
  }
  
  /* compute product */
  for(i=0;i<dim;i++) 
  for(j=0;j<dim;j++) {
    AB[idx(i,j)]=0.;
    for(k=0;k<dim;k++) {
      AB[idx(i,j)] += A[idx(i,k)]*B[idx(k,j)];
    }
  }
  
  return UNUR_SUCCESS;
  
#undef idx
} /* end of _unur_matrix_multiplication() */

/*---------------------------------------------------------------------------*/

double
_unur_matrix_determinant ( int dim, const double *A )
     /*----------------------------------------------------------------------*/
     /* Calculates the determinant of the matrix A                           */
     /* (by means of LU decomposition).                                      */
     /* input:                                                               */
     /*   dim    ... dimension of the square matrix A                        */
     /*   A      ... dim x dim -matrix                                       */
     /*									     */
     /* return:                                                              */
     /*   determinant of A                                                   */
     /*									     */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))
  
  int *p, s, i;
  double *LU;     /* array for storing LU decomposition of matrix A */
  double det;
  
  /* check arguments */
  CHECK_NULL(A,  INFINITY);
  
  /* one-dimensional case */
  if (dim==1) return A[0];
  
  /* allocate working space */
  p = _unur_xmalloc(dim*sizeof(int));
  LU = _unur_xmalloc(dim*dim*sizeof(double));
  
  /* compute LU decomposition */
  memcpy(LU, A, dim*dim*sizeof(double));
  _unur_matrix_LU_decomp(dim, LU, p, &s);
  
  /* compute determinant */
  det = s;
  for(i=0;i<dim;i++)
    det *= LU[idx(i,i)];
  
  /* free working space */
  free(LU);
  free(p);
  
  return det;
  
#undef idx
} /* end of _unur_matrix_determinant() */

/*---------------------------------------------------------------------------*/

double 
_unur_matrix_qf (int dim, double *x, double *A)
     /*----------------------------------------------------------------------*/
     /* Compute quadratic form x'Ax                                          */
     /*									     */
     /* input:                                                               */
     /*   dim ... number of columns and rows of                              */
     /*   x   ... vector                                                     */
     /*   A   ... dim x dim matrix                                           */
     /*									     */
     /* return:								     */
     /*   returns the result of x'Ax                                         */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))
  
  int i,j;
  double sum,outersum;
  
  /* check arguments */
  CHECK_NULL(x,INFINITY);
  CHECK_NULL(A,INFINITY);
  if (dim<1) {
    _unur_error("matrix",UNUR_ERR_GENERIC,"dimension < 1");
    return INFINITY;
  }
  
  outersum=0.;
  for(i=0;i<dim;i++){
    sum=0.;
    for(j=0;j<dim;j++)
      sum+=A[idx(i,j)]*x[j];
    outersum+=sum*x[i];
  }
  
  return outersum;
  
#undef idx
} /* end of _unur_matrix_qf() */

/*---------------------------------------------------------------------------*/

int
_unur_matrix_cholesky_decomposition (int dim, const double *S, double *L )
     /*----------------------------------------------------------------------*/
     /* The Colesky factor L of a variance-covariance matrix S is computed:  */
     /*    S = LL'                                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   dim ... dimension of covariance matrix S                           */
     /*   S   ... variance-covariance matrix                                 */
     /*   L   ... pointer to array where cholesky factor should be stored    */
     /*                                                                      */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   UNUR_FAILURE when matrix is not positive definite                  */
     /*   other error code, otherwise                                        */
     /*----------------------------------------------------------------------*/
{ 
#define idx(a,b) ((a)*dim+(b))
  
  int i,j,k;
  double sum1,sum2;
  
  /* check arguments */
  CHECK_NULL(S,UNUR_ERR_NULL);
  if (dim<1) {
    _unur_error("matrix",UNUR_ERR_GENERIC,"dimension < 1");
    return UNUR_ERR_GENERIC;
  }
  
  L[idx(0,0)] = sqrt( S[idx(0,0)] );
  
  for(j=1; j<dim; j++) {
    L[idx(j,0)] = S[idx(j,0)] / L[idx(0,0)];

    sum1 = L[idx(j,0)] * L[idx(j,0)];
    for(k=1; k<j; k++) {
      sum2 = 0.;
      for(i=0; i<k; i++)
	sum2 += L[idx(j,i)] * L[idx(k,i)];
      
      L[idx(j,k)] = (S[idx(j,k)] - sum2) / L[idx(k,k)];
      sum1 += L[idx(j,k)] * L[idx(j,k)];
    }

    if (!(S[idx(j,j)] > sum1)) {
      /* covariance matrix not positive definite */
      return UNUR_FAILURE;
    }

    L[idx(j,j)] = sqrt( S[idx(j,j)] - sum1 );
  }

  /* although not necessary upper triangular of L - matrix is set to 0 */
  for(j=0; j<dim; j++)
    for(k=j+1; k<dim; k++)
      L[idx(j,k)]=0.;

  /* o.k. */
  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_matrix_cholesky_decomposition() */

/*--------------------------------------------------------------------------*/

void
_unur_matrix_print_vector ( int dim, const double *vec, const char *info,
			    FILE *LOG, const char *genid, const char *indent )
     /*----------------------------------------------------------------------*/
     /* Print elements of vector in a single row enclosed by parenthesis     */
     /* into LOG file.                                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   dim    ... dimension                                               */
     /*   vec    ... vector with <dim> entries                               */
     /*   info   ... additional info-string to be printed as first line      */
     /*   LOG    ... output stream                                           */
     /*   genid  ... id string                                               */
     /*   indent ... left margin of printed vector                           */
     /*----------------------------------------------------------------------*/
{
  int i;

  if (vec) {
    fprintf(LOG,"%s: %s\n", genid, info );
    fprintf(LOG,"%s: %s( %g", genid, indent, vec[0]);
    for (i=1; i<dim; i++) 
      fprintf(LOG,", %g", vec[i]);
    fprintf(LOG," )\n");
  }
  else {
    fprintf(LOG,"%s: %s [unknown]\n", genid, info );
  }

  fprintf(LOG,"%s:\n",genid);
} /* end of _unur_matrix_print_vector() */

/*--------------------------------------------------------------------------*/

void
_unur_matrix_print_matrix ( int dim, const double *mat, const char *info,
			   FILE *LOG, const char *genid, const char *indent )
     /*----------------------------------------------------------------------*/
     /* Print elements of the given <dim>x<dim> square matrix into LOG file. */
     /* The matrix is stored row-wise in <mat>.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   dim    ... dimension                                               */
     /*   mat    ... square matrix with <dim> rows and columns               */
     /*   info   ... additional info-string to be printed as first line      */
     /*   LOG    ... output stream                                           */
     /*   genid  ... id string                                               */
     /*   indent ... left margin of printed vector                           */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))
  int i,j;

  if (mat) {
    fprintf(LOG,"%s: %s\n", genid, info); 
    for (i=0; i<dim; i++) {
      fprintf(LOG,"%s: %s(% e", genid, indent, mat[idx(i,0)]);
      for (j=1; j<dim; j++)
	fprintf(LOG,",% e",mat[idx(i,j)]);
      fprintf(LOG," )\n");
    }
  }
  else {
    fprintf(LOG,"%s: %s [unknown]\n", genid, info );
  }

  fprintf(LOG,"%s:\n",genid);

#undef idx
} /* end of _unur_matrix_print_matrix() */

/*--------------------------------------------------------------------------*/
