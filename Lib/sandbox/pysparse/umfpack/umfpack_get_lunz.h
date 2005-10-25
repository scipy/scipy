/* ========================================================================== */
/* === umfpack_get_lunz ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_get_lunz
(
    int *lnz,
    int *unz,
    int *n_row,
    int *n_col,
    int *nz_udiag,
    void *Numeric
) ;

long umfpack_dl_get_lunz
(
    long *lnz,
    long *unz,
    long *n_row,
    long *n_col,
    long *nz_udiag,
    void *Numeric
) ;

int umfpack_zi_get_lunz
(
    int *lnz,
    int *unz,
    int *n_row,
    int *n_col,
    int *nz_udiag,
    void *Numeric
) ;

long umfpack_zl_get_lunz
(
    long *lnz,
    long *unz,
    long *n_row,
    long *n_col,
    long *nz_udiag,
    void *Numeric
) ;

/*
double int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    int status, lnz, unz, n_row, n_col ;
    status = umfpack_di_get_lunz (&lnz, &unz, &n_row, &n_col, Numeric) ;

double long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    long status, lnz, unz, n_row, n_col ;
    status = umfpack_dl_get_lunz (&lnz, &unz, &n_row, &n_col, Numeric) ;

complex int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    int status, lnz, unz, n_row, n_col ;
    status = umfpack_zi_get_lunz (&lnz, &unz, &n_row, &n_col, Numeric) ;

complex long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    long status, lnz, unz, n_row, n_col ;
    status = umfpack_zl_get_lunz (&lnz, &unz, &n_row, &n_col, Numeric) ;

Purpose:

    Determines the size and number of nonzeros in the LU factors held by the
    Numeric object.  These are also the sizes of the output arrays required
    by umfpack_*_get_numeric.

    The matrix L is n_row -by- min(n_row,n_col), with lnz nonzeros, including
    the entries on the unit diagonal of L.

    The matrix U is min(n_row,n_col) -by- n_col, with unz nonzeros, including
    nonzeros on the diagonal of U.

Returns:

    UMFPACK_OK if successful.
    UMFPACK_ERROR_invalid_Numeric_object if Numeric is not a valid object.
    UMFPACK_ERROR_argument_missing if any other argument is (Int *) NULL.

Arguments:

    Int *lnz ;		Output argument.

	The number of nonzeros in L, including the diagonal (which is all
	one's).  This value is the required size of the Lj and Lx arrays as
	computed by umfpack_*_get_numeric.  The value of lnz is identical to
	Info [UMFPACK_LNZ], if that value was returned by umfpack_*_numeric.

    Int *unz ;		Output argument.

	The number of nonzeros in U, including the diagonal.  This value is the
	required size of the Ui and Ux arrays as computed by
	umfpack_*_get_numeric.  The value of unz is identical to
	Info [UMFPACK_UNZ], if that value was returned by umfpack_*_numeric.

    Int *n_row ;	Output argument.
    Int *n_col ;	Output argument.

	The order of the L and U matrices.  L is n_row -by- min(n_row,n_col)
	and U is min(n_row,n_col) -by- n_col.

    Int *nz_udiag ;	Output argument.

	The number of numerically nonzero values on the diagonal of U.  The
	matrix is singular if nz_diag < min(n_row,n_col).  A divide-by-zero
	will occur if nz_diag < n_row == n_col when solving a sparse system
	involving the matrix U in umfpack_*_*solve.  The value of nz_udiag is
	identical to Info [UMFPACK_UDIAG_NZ] if that value was returned by
	umfpack_*_numeric.

    void *Numeric ;	Input argument, not modified.

	Numeric must point to a valid Numeric object, computed by
	umfpack_*_numeric.
*/
