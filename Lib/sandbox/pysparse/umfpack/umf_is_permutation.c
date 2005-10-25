/* ========================================================================== */
/* === UMF_is_permutation =================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Return TRUE if P is a r-permutation vector, FALSE otherwise */
/* P [0..r-1] must be an r-permutation of 0..n-1 */

#include "umf_internal.h"

GLOBAL Int UMF_is_permutation
(
    const Int P [ ],	/* permutation of size r */
    Int W [ ],		/* workspace of size n */
    Int n,
    Int r
)
{
    Int i, k ;

    if (!P)
    {
	/* if P is (Int *) NULL, this is the identity permutation */
	return (TRUE) ;
    }

    ASSERT (W != (Int *) NULL) ;

    for (i = 0 ; i < n ; i++)
    {
	W [i] = FALSE ;
    }
    for (k = 0 ; k < r ; k++)
    {
	i = P [k] ;
	DEBUG5 (("k "ID" i "ID"\n", k, i)) ;
	if (i < 0 || i >= n)
	{
	    DEBUG0 (("i out of range "ID" "ID"\n", i, n)) ;
	    return (FALSE) ;
	}
	if (W [i])
	{
	    DEBUG0 (("i duplicate "ID"\n", i)) ;
	    return (FALSE) ;
	}
	W [i] = TRUE ;
    }
    return (TRUE) ;
}
