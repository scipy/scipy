//------------------------------------------------------------------------------
// AMD/Source/amd_preprocess: sort, remove duplicates, transpose a matrix
//------------------------------------------------------------------------------

// AMD, Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and
// Iain S. Duff.  All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

/* Sorts, removes duplicate entries, and transposes from the nonzero pattern of
 * a column-form matrix A, to obtain the matrix R.  The input matrix can have
 * duplicate entries and/or unsorted columns (AMD_valid (n,Ap,Ai) must not be
 * AMD_INVALID).
 *
 * This input condition is NOT checked.  This routine is not user-callable.
 */

#include "amd_internal.h"

/* AMD_preprocess does not check its input for errors or allocate workspace.
 * On input, the condition (AMD_valid (n,n,Ap,Ai) != AMD_INVALID) must hold.
 */

void Highs_amd_preprocess
(
    amd_int n,		/* input matrix: A is n-by-n */
    const amd_int Ap [ ],	/* size n+1 */
    const amd_int Ai [ ],	/* size nz = Ap [n] */

    /* output matrix R: */
    amd_int Rp [ ],		/* size n+1 */
    amd_int Ri [ ],		/* size nz (or less, if duplicates present) */

    amd_int W [ ],		/* workspace of size n */
    amd_int Flag [ ]	/* workspace of size n */
)
{

    /* --------------------------------------------------------------------- */
    /* local variables */
    /* --------------------------------------------------------------------- */

    amd_int i, j, p, p2 ;

    

    /* --------------------------------------------------------------------- */
    /* count the entries in each row of A (excluding duplicates) */
    /* --------------------------------------------------------------------- */

    for (i = 0 ; i < n ; i++)
    {
	W [i] = 0 ;		/* # of nonzeros in row i (excl duplicates) */
	Flag [i] = EMPTY ;	/* Flag [i] = j if i appears in column j */
    }
    for (j = 0 ; j < n ; j++)
    {
	p2 = Ap [j+1] ;
	for (p = Ap [j] ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    if (Flag [i] != j)
	    {
		/* row index i has not yet appeared in column j */
		W [i]++ ;	    /* one more entry in row i */
		Flag [i] = j ;	    /* flag row index i as appearing in col j*/
	    }
	}
    }

    /* --------------------------------------------------------------------- */
    /* compute the row pointers for R */
    /* --------------------------------------------------------------------- */

    Rp [0] = 0 ;
    for (i = 0 ; i < n ; i++)
    {
	Rp [i+1] = Rp [i] + W [i] ;
    }
    for (i = 0 ; i < n ; i++)
    {
	W [i] = Rp [i] ;
	Flag [i] = EMPTY ;
    }

    /* --------------------------------------------------------------------- */
    /* construct the row form matrix R */
    /* --------------------------------------------------------------------- */

    /* R = row form of pattern of A */
    for (j = 0 ; j < n ; j++)
    {
	p2 = Ap [j+1] ;
	for (p = Ap [j] ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    if (Flag [i] != j)
	    {
		/* row index i has not yet appeared in column j */
		Ri [W [i]++] = j ;  /* put col j in row i */
		Flag [i] = j ;	    /* flag row index i as appearing in col j*/
	    }
	}
    }

}
