/* ========================================================================== */
/* === UMFPACK_get_lunz ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Determines the number of nonzeros in L and U, and the size
    of L and U.
*/

#include "umf_internal.h"
#include "umf_valid_numeric.h"

GLOBAL Int UMFPACK_get_lunz
(
    Int *lnz,
    Int *unz,
    Int *n_row,
    Int *n_col,
    Int *nz_udiag,
    void *NumericHandle
)
{
    NumericType *Numeric ;

    Numeric = (NumericType *) NumericHandle ;

    if (!UMF_valid_numeric (Numeric))
    {
	return (UMFPACK_ERROR_invalid_Numeric_object) ;
    }
    if (!lnz || !unz || !n_row || !n_col || !nz_udiag)
    {
	return (UMFPACK_ERROR_argument_missing) ;
    }

    *n_row = Numeric->n_row ;
    *n_col = Numeric->n_col ;

    /* number of nz's in L below diagonal, plus the unit diagonal of L */
    *lnz = Numeric->lnz + MIN (Numeric->n_row, Numeric->n_col) ;

    /* number of nz's in U above diagonal, plus nz's on diagaonal of U */
    *unz = Numeric->unz + Numeric->nnzpiv ;

    /* number of nz's on the diagonal */
    *nz_udiag = Numeric->nnzpiv ;

    return (UMFPACK_OK) ;
}
