/* ========================================================================== */
/* === UMF_fsize ============================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Determine the largest frontal matrix size for each subtree.   Called by
 * UMF_colamd and UMF_analyze.  Only required to sort the children of each
 * node prior to AMD_postorder. */

#include "umf_internal.h"

GLOBAL void UMF_fsize
(
    Int nn,
    Int Fsize [ ],
    Int Fnrows [ ],
    Int Fncols [ ],
    Int Parent [ ],
    Int Npiv [ ]
)
{
    Int j, parent, frsize, r, c ;

    for (j = 0 ; j < nn ; j++)
    {
	Fsize [j] = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* find max front size for tree rooted at node j, for each front j */
    /* ---------------------------------------------------------------------- */

    DEBUG1 (("\n\n========================================FRONTS:\n")) ;
    for (j = 0 ; j < nn ; j++)
    {
	if (Npiv [j] > 0)
	{
	    /* this is a frontal matrix */
	    parent = Parent [j] ;
	    r = Fnrows [j] ;
	    c = Fncols [j] ;
	    frsize = r * c ;
	    /* avoid integer overflow */
	    if (INT_OVERFLOW (((double) r) * ((double) c)))
	    {
		/* :: frsize int overflow :: */
		frsize = Int_MAX ;
	    }
	    DEBUG1 ((""ID" : npiv "ID" size "ID" parent "ID" ",
		j, Npiv [j], frsize, parent)) ;
	    Fsize [j] = MAX (Fsize [j], frsize) ;
	    DEBUG1 (("Fsize [j = "ID"] = "ID"\n", j, Fsize [j])) ;
	    if (parent != EMPTY)
	    {
		/* find the maximum frontsize of self and children */
		ASSERT (Npiv [parent] > 0) ;
		ASSERT (parent > j) ;
		Fsize [parent] = MAX (Fsize [parent], Fsize [j]) ;
		DEBUG1 (("Fsize [parent = "ID"] = "ID"\n",
		    parent, Fsize [parent]));
	    }
	}
    }
}
