/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file cgsequ.c
 * \brief Computes row and column scalings
 *
 * <pre>
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 * Modified from LAPACK routine CGEEQU
 * </pre>
 */
/*
 * File name:	cgsequ.c
 * History:     Modified from LAPACK routine CGEEQU
 */
#include <math.h>
#include "slu_cdefs.h"



/*! \brief
 *
 * <pre>
 * Purpose   
 *   =======   
 *
 *   CGSEQU computes row and column scalings intended to equilibrate an   
 *   M-by-N sparse matrix A and reduce its condition number. R returns the row
 *   scale factors and C the column scale factors, chosen to try to make   
 *   the largest element in each row and column of the matrix B with   
 *   elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.   
 *
 *   R(i) and C(j) are restricted to be between SMLNUM = smallest safe   
 *   number and BIGNUM = largest safe number.  Use of these scaling   
 *   factors is not guaranteed to reduce the condition number of A but   
 *   works well in practice.   
 *
 *   See supermatrix.h for the definition of 'SuperMatrix' structure.
 *
 *   Arguments   
 *   =========   
 *
 *   A       (input) SuperMatrix*
 *           The matrix of dimension (A->nrow, A->ncol) whose equilibration
 *           factors are to be computed. The type of A can be:
 *           Stype = SLU_NC; Dtype = SLU_C; Mtype = SLU_GE.
 *	    
 *   R       (output) float*, size A->nrow
 *           If INFO = 0 or INFO > M, R contains the row scale factors   
 *           for A.
 *	    
 *   C       (output) float*, size A->ncol
 *           If INFO = 0,  C contains the column scale factors for A.
 *	    
 *   ROWCND  (output) float*
 *           If INFO = 0 or INFO > M, ROWCND contains the ratio of the   
 *           smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and   
 *           AMAX is neither too large nor too small, it is not worth   
 *           scaling by R.
 *	    
 *   COLCND  (output) float*
 *           If INFO = 0, COLCND contains the ratio of the smallest   
 *           C(i) to the largest C(i).  If COLCND >= 0.1, it is not   
 *           worth scaling by C.
 *	    
 *   AMAX    (output) float*
 *           Absolute value of largest matrix element.  If AMAX is very   
 *           close to overflow or very close to underflow, the matrix   
 *           should be scaled.
 *	    
 *   INFO    (output) int*
 *           = 0:  successful exit   
 *           < 0:  if INFO = -i, the i-th argument had an illegal value   
 *           > 0:  if INFO = i,  and i is   
 *                 <= A->nrow:  the i-th row of A is exactly zero   
 *                 >  A->ncol:  the (i-M)-th column of A is exactly zero   
 *
 *   ===================================================================== 
 * </pre>
 */
void
cgsequ(SuperMatrix *A, float *r, float *c, float *rowcnd,
	float *colcnd, float *amax, int *info)
{


    /* Local variables */
    NCformat *Astore;
    singlecomplex   *Aval;
    int_t i, j;
    int   irow;
    float rcmin, rcmax;
    float bignum, smlnum;
    extern float smach(char *);
    
    /* Test the input parameters. */
    *info = 0;
    if ( A->nrow < 0 || A->ncol < 0 ||
	 A->Stype != SLU_NC || A->Dtype != SLU_C || A->Mtype != SLU_GE )
	*info = -1;
    if (*info != 0) {
	int ii = -(*info);
	input_error("cgsequ", &ii);
	return;
    }

    /* Quick return if possible */
    if ( A->nrow == 0 || A->ncol == 0 ) {
	*rowcnd = 1.;
	*colcnd = 1.;
	*amax = 0.;
	return;
    }

    Astore = A->Store;
    Aval = Astore->nzval;
    
    /* Get machine constants. */
    smlnum = smach("S");  /* slamch_("S"); */
    bignum = 1. / smlnum;

    /* Compute row scale factors. */
    for (i = 0; i < A->nrow; ++i) r[i] = 0.;

    /* Find the maximum element in each row. */
    for (j = 0; j < A->ncol; ++j)
	for (i = Astore->colptr[j]; i < Astore->colptr[j+1]; ++i) {
	    irow = Astore->rowind[i];
	    r[irow] = SUPERLU_MAX( r[irow], c_abs1(&Aval[i]) );
	}

    /* Find the maximum and minimum scale factors. */
    rcmin = bignum;
    rcmax = 0.;
    for (i = 0; i < A->nrow; ++i) {
	rcmax = SUPERLU_MAX(rcmax, r[i]);
	rcmin = SUPERLU_MIN(rcmin, r[i]);
    }
    *amax = rcmax;

    if (rcmin == 0.) {
	/* Find the first zero scale factor and return an error code. */
	for (i = 0; i < A->nrow; ++i)
	    if (r[i] == 0.) {
		*info = i + 1;
		return;
	    }
    } else {
	/* Invert the scale factors. */
	for (i = 0; i < A->nrow; ++i)
	    r[i] = 1. / SUPERLU_MIN( SUPERLU_MAX( r[i], smlnum ), bignum );
	/* Compute ROWCND = min(R(I)) / max(R(I)) */
	*rowcnd = SUPERLU_MAX( rcmin, smlnum ) / SUPERLU_MIN( rcmax, bignum );
    }

    /* Compute column scale factors */
    for (j = 0; j < A->ncol; ++j) c[j] = 0.;

    /* Find the maximum element in each column, assuming the row
       scalings computed above. */
    for (j = 0; j < A->ncol; ++j)
	for (i = Astore->colptr[j]; i < Astore->colptr[j+1]; ++i) {
	    irow = Astore->rowind[i];
	    c[j] = SUPERLU_MAX( c[j], c_abs1(&Aval[i]) * r[irow] );
	}

    /* Find the maximum and minimum scale factors. */
    rcmin = bignum;
    rcmax = 0.;
    for (j = 0; j < A->ncol; ++j) {
	rcmax = SUPERLU_MAX(rcmax, c[j]);
	rcmin = SUPERLU_MIN(rcmin, c[j]);
    }

    if (rcmin == 0.) {
	/* Find the first zero scale factor and return an error code. */
	for (j = 0; j < A->ncol; ++j)
	    if ( c[j] == 0. ) {
		*info = A->nrow + j + 1;
		return;
	    }
    } else {
	/* Invert the scale factors. */
	for (j = 0; j < A->ncol; ++j)
	    c[j] = 1. / SUPERLU_MIN( SUPERLU_MAX( c[j], smlnum ), bignum);
	/* Compute COLCND = min(C(J)) / max(C(J)) */
	*colcnd = SUPERLU_MAX( rcmin, smlnum ) / SUPERLU_MIN( rcmax, bignum );
    }

    return;

} /* cgsequ */


