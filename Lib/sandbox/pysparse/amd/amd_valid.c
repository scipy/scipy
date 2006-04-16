/* ========================================================================== */
/* === AMD_valid ============================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* AMD Version 1.0 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A. Davis,   */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README for License.          */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/amd                           */
/* -------------------------------------------------------------------------- */

/* Check if a column-form matrix is valid or not.  The matrix A is
 * n_row-by-n_col.  The row indices of entries in column j are in
 * Ai [Ap [j] ... Ap [j+1]-1].  Required conditions are:
 *
 *	n_row >= 0
 *	n_col >= 0
 *	nz = Ap [n_col] >= 0	    number of entries in the matrix
 *	Ap [0] == 0
 *	Ap [j] <= Ap [j+1] for all j in the range 0 to n_col.
 *	row indices in Ai [Ap [j] ... Ap [j+1]-1] must be sorted in ascending
 *	    order, must be in the range 0 to n_row-1, and no duplicate entries
 *	    can exist.
 *
 * Not user-callable.
 */

#include "amd_internal.h"

GLOBAL Int AMD_valid
(
    /* inputs, not modified on output: */
    Int n_row,		/* A is n_row-by-n_col */
    Int n_col,
    const Int Ap [ ],	/* column pointers of A, of size n_col+1 */
    const Int Ai [ ]	/* row indices of A, of size nz = Ap [n_col] */
)
{
    Int nz, j, p1, p2, ilast, i, p ;
    if (n_row < 0 || n_col < 0)
    {
	AMD_DEBUG0 (("n must be >= 0: "ID" "ID"\n", n_row, n_col)) ;
	return (FALSE) ;
    }
    nz = Ap [n_col] ;
    if (Ap [0] != 0 || nz < 0)
    {
	/* column pointers must start at Ap [0] = 0, and Ap [n] must be >= 0 */
	AMD_DEBUG0 (("column 0 pointer bad or nz < 0\n")) ;
	return (FALSE) ;
    }
    for (j = 0 ; j < n_col ; j++)
    {
	p1 = Ap [j] ;
	p2 = Ap [j+1] ;
	AMD_DEBUG2 (("\nColumn: "ID" p1: "ID" p2: "ID"\n", j, p1, p2)) ;
	if (p1 > p2)
	{
	    /* column pointers must be ascending */
	    AMD_DEBUG0 (("column "ID" pointer bad\n", j)) ;
	    return (FALSE) ;
	}
	ilast = EMPTY ;
	for (p = p1 ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    AMD_DEBUG3 (("row: "ID"\n", i)) ;
	    if (i <= ilast || i >= n_row)
	    {
		/* row index out of range, or unsorted */
		AMD_DEBUG0 (("index out of range, col "ID" row "ID"\n", j, i)) ;
		return (FALSE) ;
	    }
	    ilast = i ;
	}
    }
    return (TRUE) ;
}
