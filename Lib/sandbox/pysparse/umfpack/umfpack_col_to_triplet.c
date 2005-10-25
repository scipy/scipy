/* ========================================================================== */
/* === UMFPACK_col_to_triplet =============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User callable.  Converts a column-oriented input matrix to triplet form by
    constructing the column indices Tj from the column pointers Ap.  The matrix
    may be singular.  See umfpack_col_to_triplet.h for details.

*/

#include "umf_internal.h"

GLOBAL Int UMFPACK_col_to_triplet
(
    Int n_col,
    const Int Ap [ ],
    Int Tj [ ]
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int nz, j, p, p1, p2, length ;

    /* ---------------------------------------------------------------------- */
    /* construct the column indices */
    /* ---------------------------------------------------------------------- */

    if (!Ap || !Tj)
    {
	return (UMFPACK_ERROR_argument_missing) ;
    }
    if (n_col <= 0)
    {
	return (UMFPACK_ERROR_n_nonpositive) ;
    }
    if (Ap [0] != 0)
    {
	return (UMFPACK_ERROR_invalid_matrix) ;
    }
    nz = Ap [n_col] ;
    if (nz < 0)
    {
	return (UMFPACK_ERROR_invalid_matrix) ;
    }

    for (j = 0 ; j < n_col ; j++)
    {
	p1 = Ap [j] ;
	p2 = Ap [j+1] ;
	length = p2 - p1 ;
	if (length < 0 || p2 > nz)
	{
	    return (UMFPACK_ERROR_invalid_matrix) ;
	}
	for (p = p1 ; p < p2 ; p++)
	{
	    Tj [p] = j ;
	}
    }

    return (UMFPACK_OK) ;
}
