/* ========================================================================== */
/* === AMD_aat ============================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* AMD Version 1.0 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A. Davis,   */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README for License.          */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/amd                           */
/* -------------------------------------------------------------------------- */

/* AMD_aat:  compute the symmetry of the pattern of A, and count the number of
 * nonzeros each column of A+A' (excluding the diagonal).  Assume the input
 * matrix has no errors.
 */

#include "amd_internal.h"

GLOBAL Int AMD_aat	/* returns nz in A+A' */
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int Len [ ],	/* Len [j]: length of column j of A+A', excl. diagonal*/
    Int Tp [ ],		/* workspace of size n */
    double Info [ ]
) 
{
    Int p1, p2, p, i, j, pj, pj2, k, nzdiag, nzboth, nz, nzaat ;
    double sym ;

#ifndef NDEBUG
    AMD_debug_init ("AMD AAT") ;
    for (k = 0 ; k < n ; k++) Tp [k] = EMPTY ;
    ASSERT (AMD_valid (n, n, Ap, Ai)) ;
#endif

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
	AMD_DEBUG2 (("\nAAT Column: "ID" p1: "ID" p2: "ID"\n", k, p1, p2)) ;

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
		AMD_DEBUG3 (("    upper ("ID","ID") ("ID","ID")\n", j,k, k,j)) ;
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
	    ASSERT (Tp [j] != EMPTY) ;
	    ASSERT (Ap [j] <= Tp [j] && Tp [j] <= Ap [j+1]) ;
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
		    AMD_DEBUG3 (("    lower ("ID","ID") ("ID","ID")\n",
			i,j, j,i)) ;
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
	    AMD_DEBUG3 (("    lower cleanup ("ID","ID") ("ID","ID")\n",
		i,j, j,i)) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* compute the symmetry of the nonzero pattern of A */
    /* ---------------------------------------------------------------------- */

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
	sym = ((double) (2 * nzboth)) / ((double) (nz - nzdiag)) ;
    }

    nzaat = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	nzaat += Len [k] ;
    }
    AMD_DEBUG1 (("AMD nz in A+A', excluding diagonal (nzaat) = "ID"\n", nzaat));
    AMD_DEBUG1 (("   nzboth: "ID" nz: "ID" nzdiag: "ID" symmetry: %g\n",
		nzboth, nz, nzdiag, sym)) ;

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
