/* ========================================================================== */
/* === AMD_order ============================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* AMD Version 1.0 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A. Davis,   */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README for License.          */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/amd                           */
/* -------------------------------------------------------------------------- */

/* User-callable AMD minimum degree ordering routine.  See amd.h for
 * documentation.
 */

#include "amd_internal.h"

GLOBAL Int AMD_order
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int P [ ],
    double Control [ ],
    double Info [ ]
)
{
    Int slen, *Len, *S, nz, nzaat, i, *Pinv, info ;

#ifndef NDEBUG
    AMD_debug_init ("amd") ;
#endif

    /* clear the Info array, if it exists */
    info = Info != (double *) NULL ;
    if (info)
    {
	for (i = 0 ; i < AMD_INFO ; i++)
	{
	    Info [i] = EMPTY ;
	}
	Info [AMD_N] = n ;
	Info [AMD_STATUS] = AMD_OK ;
    }

    /* make sure inputs exist and n is >= 0 */
    if (Ai == (Int *) NULL || Ap == (Int *) NULL || P == (Int *) NULL || n < 0)
    {
	if (info) Info [AMD_STATUS] = AMD_INVALID ;
	return (AMD_INVALID) ;	    /* arguments are invalid */
    }

    if (n == 0)
    {
	return (AMD_OK) ;	    /* n is 0 so there's nothing to do */
    }

    nz = Ap [n] ;
    if (info)
    {
	Info [AMD_NZ] = nz ;
    }
    if (nz < 0)
    {
	if (info) Info [AMD_STATUS] = AMD_INVALID ;
	return (AMD_INVALID) ;
    }

    /* Avoid integer overflow in memory size calculations.  The space required
     * by AMD is at most 2.4nz + 8n for S, and n for Len.
     * Note nz - n <= nzaat <= 2*nz, below. */
    if ((2.4 * (double) nz + 8 * (double) n) > (double) Int_MAX / sizeof (Int))
    {
	/* :: int overflow :: */
	if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	return (AMD_OUT_OF_MEMORY) ;
    }

    if (!AMD_valid (n, n, Ap, Ai))
    {
	if (info) Info [AMD_STATUS] = AMD_INVALID ;
	return (AMD_INVALID) ;	    /* matrix is invalid */
    }

    /* ---------------------------------------------------------------------- */
    /* determine the symmetry and count off-diagonal nonzeros in A+A' */
    /* ---------------------------------------------------------------------- */

    /* allocate size-n integer workspace */
    Len = (Int *) ALLOCATE (n * sizeof (Int)) ;
    if (!Len)
    {
	/* :: out of memory :: */
	if (info) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	return (AMD_OUT_OF_MEMORY) ;
    }
    nzaat = AMD_aat (n, Ap, Ai, Len, P, Info) ;
    AMD_DEBUG1 (("nzaat: "ID"\n", nzaat)) ;
    ASSERT (nz-n <= nzaat && nzaat <= 2*nz) ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace for matrix, elbow room, and 7 size-n vectors */
    /* ---------------------------------------------------------------------- */

    slen = (nzaat + nzaat/5 + n) + 7*n ;
    if (info)
    {
	/* memory usage (Len and S), in bytes. */
	Info [AMD_MEMORY] = ((double) slen + n) * sizeof (Int) ;
    }
    S = (Int *) ALLOCATE (slen * sizeof (Int)) ;
    AMD_DEBUG1 ((" S "ID" Len "ID" n "ID" nzaat "ID" slen "ID"\n",
	(Int) S, (Int) Len, n, nzaat, slen)) ;
    if (S == (Int *) NULL)
    {
	/* :: out of memory :: */
	FREE (Len) ;
	if (Info != (double *) NULL) Info [AMD_STATUS] = AMD_OUT_OF_MEMORY ;
	return (AMD_OUT_OF_MEMORY) ;
    }

    /* allocate space from S for Pinv */
    Pinv = S + slen - n ;
    slen -= n ;

    /* ---------------------------------------------------------------------- */
    /* order the matrix */
    /* ---------------------------------------------------------------------- */

    AMD_1 (n, Ap, Ai, P, Pinv, Len, slen, S, Control, Info) ;

    /* ---------------------------------------------------------------------- */
    /* free the workspace */
    /* ---------------------------------------------------------------------- */

    FREE (Len) ;
    FREE (S) ;
    return (AMD_OK) ;	    /* successful ordering */
}
