//------------------------------------------------------------------------------
// AMD/Source/amd_valid: check if a matrix is valid for AMD
//------------------------------------------------------------------------------

// AMD, Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and
// Iain S. Duff.  All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

/* Check if a column-form matrix is valid or not.  The matrix A is
 * n_row-by-n_col.  The row indices of entries in column j are in
 * Ai [Ap [j] ... Ap [j+1]-1].  Required conditions are:
 *
 *	n_row >= 0
 *	n_col >= 0
 *	nz = Ap [n_col] >= 0	    number of entries in the matrix
 *	Ap [0] == 0
 *	Ap [j] <= Ap [j+1] for all j in the range 0 to n_col.
 *      Ai [0 ... nz-1] must be in the range 0 to n_row-1.
 *
 * If any of the above conditions hold, AMD_INVALID is returned.  If the
 * following condition holds, AMD_OK_BUT_JUMBLED is returned (a warning,
 * not an error):
 *
 *	row indices in Ai [Ap [j] ... Ap [j+1]-1] are not sorted in ascending
 *	    order, and/or duplicate entries exist.
 *
 * Otherwise, AMD_OK is returned.
 *
 * In v1.2 and earlier, this function returned TRUE if the matrix was valid
 * (now returns AMD_OK), or FALSE otherwise (now returns AMD_INVALID or
 * AMD_OK_BUT_JUMBLED).
 */

#include "amd_internal.h"

int Highs_amd_valid
(
    /* inputs, not modified on output: */
    amd_int n_row,		/* A is n_row-by-n_col */
    amd_int n_col,
    const amd_int Ap [ ],	/* column pointers of A, of size n_col+1 */
    const amd_int Ai [ ]	/* row indices of A, of size nz = Ap [n_col] */
)
{
    amd_int nz, j, p1, p2, ilast, i, p ;
    int result = AMD_OK ;

    if (n_row < 0 || n_col < 0 || Ap == NULL || Ai == NULL)
    {
	return (AMD_INVALID) ;
    }
    nz = Ap [n_col] ;
    if (Ap [0] != 0 || nz < 0)
    {
	/* column pointers must start at Ap [0] = 0, and Ap [n] must be >= 0 */
	
	return (AMD_INVALID) ;
    }
    for (j = 0 ; j < n_col ; j++)
    {
	p1 = Ap [j] ;
	p2 = Ap [j+1] ;
	
	if (p1 > p2)
	{
	    /* column pointers must be ascending */
	    
	    return (AMD_INVALID) ;
	}
	ilast = EMPTY ;
	for (p = p1 ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    
	    if (i < 0 || i >= n_row)
	    {
		/* row index out of range */
		
		return (AMD_INVALID) ;
	    }
	    if (i <= ilast)
	    {
		/* row index unsorted, or duplicate entry present */
		
		result = AMD_OK_BUT_JUMBLED ;
	    }
	    ilast = i ;
	}
    }
    return (result) ;
}
