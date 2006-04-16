/* ========================================================================== */
/* === UMFPACK_report_triplet =============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Prints a matrix in triplet form.  See
    umfpack_report_triplet.h for details.
*/

#include "umf_internal.h"

GLOBAL Int UMFPACK_report_triplet
(
    Int n_row,
    Int n_col,
    Int nz,
    const Int Ti [ ],
    const Int Tj [ ],
    const double Tx [ ],
#ifdef COMPLEX
    const double Tz [ ],
#endif
    const double Control [UMFPACK_CONTROL]
)
{
    Int prl, prl1, k, i, j, do_values ;
    Entry t ;

    prl = GET_CONTROL (UMFPACK_PRL, UMFPACK_DEFAULT_PRL) ;

    if (prl <= 2)
    {
	return (UMFPACK_OK) ;
    }

    PRINTF (("triplet-form matrix, n_row = "ID", n_col = "ID" nz = "ID". ",
	n_row, n_col, nz)) ;

    if (!Ti || !Tj)
    {
	PRINTF (("ERROR: indices not present\n\n")) ;
	return (UMFPACK_ERROR_argument_missing) ;
    }

    if (n_row <= 0 || n_col <= 0)
    {
	PRINTF (("ERROR: n_row or n_col is <= 0\n\n")) ;
	return (UMFPACK_ERROR_n_nonpositive) ;
    }

    if (nz < 0)
    {
	PRINTF (("ERROR: nz is < 0\n\n")) ;
	return (UMFPACK_ERROR_invalid_matrix) ;
    }

    PRINTF4 (("\n")) ;

#ifdef COMPLEX
    do_values = Tx && Tz ;
#else
    do_values = Tx != (double *) NULL ;
#endif

    prl1 = prl ;
    for (k = 0 ; k < nz ; k++)
    {
	i = Ti [k] ;
	j = Tj [k] ;
	PRINTF4 (("    "ID" : "ID" "ID" ", INDEX (k), INDEX (i), INDEX (j))) ;
	if (do_values && prl >= 4)
	{
	    ASSIGN (t, Tx [k], Tz [k]) ;
	    PRINT_ENTRY (t) ;
	}
	PRINTF4 (("\n")) ;
	if (i < 0 || i >= n_row || j < 0 || j >= n_col)
	{
	    /* invalid triplet */
	    PRINTF (("ERROR: invalid triplet\n\n")) ;
	    return (UMFPACK_ERROR_invalid_matrix) ;
	}
	if (prl == 4 && k == 9 && nz > 10)
	{
	    PRINTF (("    ...\n")) ;
	    prl-- ;
	}
    }
    prl = prl1 ;

    PRINTF4 (("    triplet-form matrix ")) ;
    PRINTF (("OK\n\n")) ;
    return (UMFPACK_OK) ;
}
