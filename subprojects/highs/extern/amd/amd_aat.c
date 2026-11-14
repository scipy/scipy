//------------------------------------------------------------------------------
// AMD/Source/amd_aat: compute symmetry of A and nnz in each column of A+A'
//------------------------------------------------------------------------------

// AMD, Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and
// Iain S. Duff.  All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

/* AMD_aat:  compute the symmetry of the pattern of A, and count the number of
 * nonzeros each column of A+A' (excluding the diagonal).  Assumes the input
 * matrix has no errors, with sorted columns and no duplicates
 * (AMD_valid (n, n, Ap, Ai) must be AMD_OK, but this condition is not
 * checked).
 */

#include "amd_internal.h"

size_t Highs_amd_aat	/* returns nz in A+A' */
(
    amd_int n,
    const amd_int Ap [ ],
    const amd_int Ai [ ],
    amd_int Len [ ],	/* Len [j]: length of column j of A+A', excl diagonal*/
    amd_int Tp [ ],		/* workspace of size n */
    double Info [ ]
)
{
    amd_int p1, p2, p, i, j, pj, pj2, k, nzdiag, nzboth, nz ;
    double sym ;
    size_t nzaat ;


    if (Info != (double *) NULL)
    {
	/* clear the Info array, if it exists */
	for (i = 0 ; i < AMD_INFO ; i++)
	{
	    Info [i] = EMPTY ;
	}
	Info [AMD_STATUS] = AMD_OK ;
    }

    for (k = 0 ; k < n ; k++)
    {
	Len [k] = 0 ;
    }

    nzdiag = 0 ;
    nzboth = 0 ;
    nz = Ap [n] ;

    for (k = 0 ; k < n ; k++)
    {
	p1 = Ap [k] ;
	p2 = Ap [k+1] ;
	

	/* construct A+A' */
	for (p = p1 ; p < p2 ; )
	{
	    /* scan the upper triangular part of A */
	    j = Ai [p] ;
	    if (j < k)
	    {
		/* entry A (j,k) is in the strictly upper triangular part,
		 * add both A (j,k) and A (k,j) to the matrix A+A' */
		Len [j]++ ;
		Len [k]++ ;
		
		p++ ;
	    }
	    else if (j == k)
	    {
		/* skip the diagonal */
		p++ ;
		nzdiag++ ;
		break ;
	    }
	    else /* j > k */
	    {
		/* first entry below the diagonal */
		break ;
	    }
	    /* scan lower triangular part of A, in column j until reaching
	     * row k.  Start where last scan left off. */
	    
	    
	    pj2 = Ap [j+1] ;
	    for (pj = Tp [j] ; pj < pj2 ; )
	    {
		i = Ai [pj] ;
		if (i < k)
		{
		    /* A (i,j) is only in the lower part, not in upper.
		     * add both A (i,j) and A (j,i) to the matrix A+A' */
		    Len [i]++ ;
		    Len [j]++ ;
		    

		    pj++ ;
		}
		else if (i == k)
		{
		    /* entry A (k,j) in lower part and A (j,k) in upper */
		    pj++ ;
		    nzboth++ ;
		    break ;
		}
		else /* i > k */
		{
		    /* consider this entry later, when k advances to i */
		    break ;
		}
	    }
	    Tp [j] = pj ;
	}
	/* Tp [k] points to the entry just below the diagonal in column k */
	Tp [k] = p ;
    }

    /* clean up, for remaining mismatched entries */
    for (j = 0 ; j < n ; j++)
    {
	for (pj = Tp [j] ; pj < Ap [j+1] ; pj++)
	{
	    i = Ai [pj] ;
	    /* A (i,j) is only in the lower part, not in upper.
	     * add both A (i,j) and A (j,i) to the matrix A+A' */
	    Len [i]++ ;
	    Len [j]++ ;
	    

	}
    }

    /* --------------------------------------------------------------------- */
    /* compute the symmetry of the nonzero pattern of A */
    /* --------------------------------------------------------------------- */

    /* Given a matrix A, the symmetry of A is:
     *	B = tril (spones (A), -1) + triu (spones (A), 1) ;
     *  sym = nnz (B & B') / nnz (B) ;
     *  or 1 if nnz (B) is zero.
     */

    if (nz == nzdiag)
    {
	sym = 1 ;
    }
    else
    {
	sym = (2 * (double) nzboth) / ((double) (nz - nzdiag)) ;
    }

    nzaat = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	nzaat += Len [k] ;
    }

    

    


    if (Info != (double *) NULL)
    {
	Info [AMD_STATUS] = AMD_OK ;
	Info [AMD_N] = n ;
	Info [AMD_NZ] = nz ;
	Info [AMD_SYMMETRY] = sym ;	    /* symmetry of pattern of A */
	Info [AMD_NZDIAG] = nzdiag ;	    /* nonzeros on diagonal of A */
	Info [AMD_NZ_A_PLUS_AT] = nzaat ;   /* nonzeros in A+A' */
    }

    return (nzaat) ;
}
