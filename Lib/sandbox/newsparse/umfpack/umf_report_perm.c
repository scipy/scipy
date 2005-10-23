/* ========================================================================== */
/* === UMF_report_perm ====================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

#include "umf_internal.h"

#define PRINTF4U(params) { if (user || prl >= 4) PRINTF (params) ; }

GLOBAL Int UMF_report_perm
(
    Int n,
    const Int P [ ],
    Int W [ ],		/* workspace of size n */
    Int prl,
    Int user
)
{
    Int i, k, valid, prl1 ;

    ASSERT (prl >= 3) ;

    PRINTF4U (("permutation vector, n = "ID". ", n)) ;

    if (n <= 0)
    {
	PRINTF (("ERROR: length of permutation is <= 0\n\n")) ;
	return (UMFPACK_ERROR_n_nonpositive) ;
    }

    if (!P)
    {
	/* if P is (Int *) NULL, this is the identity permutation */
	PRINTF (("(not present)\n\n")) ;
	return (UMFPACK_OK) ;
    }

    if (!W)
    {
	PRINTF (("ERROR: out of memory\n\n")) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

    PRINTF4 (("\n")) ;

    for (i = 0 ; i < n ; i++)
    {
	W [i] = TRUE ;
    }

    prl1 = prl ;
    for (k = 0 ; k < n ; k++)
    {
	i = P [k] ;
	PRINTF4 (("    "ID" : "ID" ", INDEX (k), INDEX (i))) ;
	valid = (i >= 0 && i < n) ;
	if (valid)
	{
	    valid = W [i] ;
	    W [i] = FALSE ;
	}
	if (!valid)
	{
	    /* out of range or duplicate entry */
	    PRINTF (("ERROR: invalid\n\n")) ;
	    return (UMFPACK_ERROR_invalid_permutation) ;
	}
	PRINTF4 (("\n")) ;
	if (prl == 4 && k == 9 && n > 10)
	{
	    PRINTF (("    ...\n")) ;
	    prl-- ;
	}
    }
    prl = prl1 ;

    PRINTF4 (("    permutation vector ")) ;
    PRINTF4U (("OK\n\n")) ;
    return (UMFPACK_OK) ;
}
