/* ========================================================================== */
/* === umfpack_triplet_to_col =============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_triplet_to_col
(
    int n_row,
    int n_col,
    int nz,
    const int Ti [ ],
    const int Tj [ ],
    const double Tx [ ],
    int Ap [ ],
    int Ai [ ],
    double Ax [ ],
    int Map [ ]
) ;

long umfpack_dl_triplet_to_col
(
    long n_row,
    long n_col,
    long nz,
    const long Ti [ ],
    const long Tj [ ],
    const double Tx [ ],
    long Ap [ ],
    long Ai [ ],
    double Ax [ ],
    long Map [ ]
) ;

int umfpack_zi_triplet_to_col
(
    int n_row,
    int n_col,
    int nz,
    const int Ti [ ],
    const int Tj [ ],
    const double Tx [ ], const double Tz [ ],
    int Ap [ ],
    int Ai [ ],
    double Ax [ ], double Az [ ],
    int Map [ ]
) ;

long umfpack_zl_triplet_to_col
(
    long n_row,
    long n_col,
    long nz,
    const long Ti [ ],
    const long Tj [ ],
    const double Tx [ ], const double Tz [ ],
    long Ap [ ],
    long Ai [ ],
    double Ax [ ], double Az [ ],
    long Map [ ]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    int n_row, n_col, nz, *Ti, *Tj, *Ap, *Ai, status, *Map ;
    double *Tx, *Ax ;
    status = umfpack_di_triplet_to_col (n_row, n_col, nz, Ti, Tj, Tx,
	Ap, Ai, Ax, Map) ;

double long Syntax:

    #include "umfpack.h"
    long n_row, n_col, nz, *Ti, *Tj, *Ap, *Ai, status, *Map ;
    double *Tx, *Ax ;
    status = umfpack_dl_triplet_to_col (n_row, n_col, nz, Ti, Tj, Tx,
	Ap, Ai, Ax, Map) ;

complex int Syntax:

    #include "umfpack.h"
    int n_row, n_col, nz, *Ti, *Tj, *Ap, *Ai, status, *Map ;
    double *Tx, *Tz, *Ax, *Az ;
    status = umfpack_zi_triplet_to_col (n_row, n_col, nz, Ti, Tj, Tx, Tz,
	Ap, Ai, Ax, Az, Map) ;

long Syntax:

    #include "umfpack.h"
    long n_row, n_col, nz, *Ti, *Tj, *Ap, *Ai, status, *Map ;
    double *Tx, *Tz, *Ax, *Az ;
    status = umfpack_zl_triplet_to_col (n_row, n_col, nz, Ti, Tj, Tx, Tz,
	Ap, Ai, Ax, Az, Map) ;

Purpose:

    Converts a sparse matrix from "triplet" form to compressed-column form.
    Analogous to A = spconvert (Ti, Tj, Tx + Tx*1i) in MATLAB, except that
    zero entries present in the triplet form are present in A.

    The triplet form of a matrix is a very simple data structure for basic
    sparse matrix operations.  For example, suppose you wish to factorize a
    matrix A coming from a finite element method, in which A is a sum of
    dense submatrices, A = E1 + E2 + E3 + ... .  The entries in each element
    matrix Ei can be concatenated together in the three triplet arrays, and
    any overlap between the elements will be correctly summed by
    umfpack_*_triplet_to_col.

    Transposing a matrix in triplet form is simple; just interchange the
    use of Ti and Tj.  You can construct the complex conjugate transpose by
    negating Tz, for the complex versions.

    Permuting a matrix in triplet form is also simple.  If you want the matrix
    PAQ, or A (P,Q) in MATLAB notation, where P [k] = i means that row i of
    A is the kth row of PAQ and Q [k] = j means that column j of A is the kth
    column of PAQ, then do the following.  First, create inverse permutations
    Pinv and Qinv such that Pinv [i] = k if P [k] = i and Qinv [j] = k if
    Q [k] = j.  Next, for the mth triplet (Ti [m], Tj [m], Tx [m], Tz [m]),
    replace Ti [m] with Pinv [Ti [m]] and replace Tj [m] with Qinv [Tj [m]].

    If you have a column-form matrix with duplicate entries or unsorted
    columns, you can sort it and sum up the duplicates by first converting it
    to triplet form with umfpack_*_col_to_triplet, and then converting it back
    with umfpack_*_triplet_to_col.

    Constructing a submatrix is also easy.  Just scan the triplets and remove
    those entries outside the desired subset of 0...n_row-1 and 0...n_col-1,
    and renumber the indices according to their position in the subset.

    You can do all these operations on a column-form matrix by first
    converting it to triplet form with umfpack_*_col_to_triplet, doing the
    operation on the triplet form, and then converting it back with
    umfpack_*_triplet_to_col.

    The only operation not supported easily in the triplet form is the
    multiplication of two sparse matrices (UMFPACK does not provide this
    operation).

    You can print the input triplet form with umfpack_*_report_triplet, and
    the output matrix with umfpack_*_report_matrix.

    The matrix may be singular (nz can be zero, and empty rows and/or columns
    may exist).  It may also be rectangular and/or complex.

Returns:

    UMFPACK_OK if successful.
    UMFPACK_ERROR_argument_missing if Ap, Ai, Ti, and/or Tj are missing.
    UMFPACK_ERROR_n_nonpositive if n_row <= 0 or n_col <= 0.
    UMFPACK_ERROR_invalid_matrix if nz < 0, or if for any k, Ti [k] and/or
	Tj [k] are not in the range 0 to n_row-1 or 0 to n_col-1, respectively.
    UMFPACK_ERROR_out_of_memory if unable to allocate sufficient workspace.

Arguments:

    Int n_row ;		Input argument, not modified.
    Int n_col ;		Input argument, not modified.

	A is an n_row-by-n_col matrix.  Restriction: n_row > 0 and n_col > 0.
	All row and column indices in the triplet form must be in the range
	0 to n_row-1 and 0 to n_col-1, respectively.

    Int nz ;		Input argument, not modified.

	The number of entries in the triplet form of the matrix.  Restriction:
	nz >= 0.

    Int Ti [nz] ;	Input argument, not modified.
    Int Tj [nz] ;	Input argument, not modified.
    double Tx [nz] ;	Input argument, not modified.
    double Tz [nz] ;	Input argument, not modified, for complex versions.

	Ti, Tj, Tx, and Tz hold the "triplet" form of a sparse matrix.  The kth
	nonzero entry is in row i = Ti [k], column j = Tj [k], and the real part
	of a_ij is Tx [k].  The imaginary part of a_ij is Tz [k], for complex
	versions.  The row and column indices i and j must be in the range 0 to
	n_row-1 and 0 to n_col-1, respectively.  Duplicate entries may be
	present; they are summed in the output matrix.  This is not an error
	condition.  The "triplets" may be in any order.  Tx, Tz, Ax, and Az
	are optional.  For the real version, Ax is computed only if both Ax
	and Tx are present (not (double *) NULL).  For the complex version, Ax
	and Az are computed only if Tx, Tz, Ax, and Az are all present.  These
	are not error conditions; the routine can create just the pattern of
	the output matrix from the pattern of the triplets.

	Future complex version:  if Tx is present and Tz is NULL, then both real
	and imaginary parts will be contained in Tx[0..2*nz-1], with Tx[2*k]
	and Tx[2*k+1] being the real and imaginary part of the kth entry.

    Int Ap [n_col+1] ;	Output argument.

	Ap is an integer array of size n_col+1 on input.  On output, Ap holds
	the "pointers" for the column form of the sparse matrix A.  Column j of
	the matrix A is held in Ai [(Ap [j]) ... (Ap [j+1]-1)].  The first
	entry, Ap [0], is zero, and Ap [j] <= Ap [j+1] holds for all j in the
	range 0 to n_col-1.  The value nz2 = Ap [n_col] is thus the total
	number of entries in the pattern of the matrix A.  Equivalently, the
	number of duplicate triplets is nz - Ap [n_col].

    Int Ai [nz] ;	Output argument.

	Ai is an integer array of size nz on input.  Note that only the first
	Ap [n_col] entries are used.

	The nonzero pattern (row indices) for column j is stored in
	Ai [(Ap [j]) ... (Ap [j+1]-1)].  The row indices in a given column j
	are in ascending order, and no duplicate row indices are present.
	Row indices are in the range 0 to n_col-1 (the matrix is 0-based).

    double Ax [nz] ;	Output argument.
    double Az [nz] ;	Output argument for complex versions.

	Ax and Az (for the complex versions) are double arrays of size nz on
	input.  Note that only the first Ap [n_col] entries are used
	in both arrays.

	Ax is optional; if Tx and/or Ax are not present (a (double *) NULL
	pointer), then Ax is not computed.  Az is also optional; if Tz and/or
	Az are not present, then Az is not computed.  If present, Ax holds the
	numerical values of the the real part of the sparse matrix A and Az
	holds the imaginary parts.  The nonzero pattern (row indices) for
	column j is stored in Ai [(Ap [j]) ... (Ap [j+1]-1)], and the
	corresponding numerical values are stored in
	Ax [(Ap [j]) ... (Ap [j+1]-1)].  The imaginary parts are stored in
	Az [(Ap [j]) ... (Ap [j+1]-1)], for the complex versions.

	Future complex version:  if Ax is present and Az is NULL, then both real
	and imaginary parts will be returned in Ax[0..2*nz2-1], with Ax[2*k]
	and Ax[2*k+1] being the real and imaginary part of the kth entry.

    int Map [nz] ;	Optional output argument.

	If Map is present (a non-NULL pointer to an Int array of size nz), then
	on output it holds the position of the triplets in the column-form
	matrix.  That is, suppose p = Map [k], and the k-th triplet is i=Ti[k],
	j=Tj[k], and aij=Tx[k].  Then i=Ai[p], and aij will have been summed
	into Ax[p] (or simply aij=Ax[p] if there were no duplicate entries also
	in row i and column j).  Also, Ap[j] <= p < Ap[j+1].  The Map array is
	not computed if it is (Int *) NULL.  The Map array is useful for
	converting a subsequent triplet form matrix with the same pattern as the
	first one, without calling this routine.  If Ti and Tj do not change,
	then Ap, and Ai can be reused from the prior call to
	umfpack_*_triplet_to_col.  You only need to recompute Ax (and Az for the
	complex version).  This code excerpt properly sums up all duplicate
	values (for the real version):

	    for (p = 0 ; p < Ap [n_col] ; p++) Ax [p] = 0 ;
	    for (k = 0 ; k < nz ; k++) Ax [Map [k]] += Tx [k] ;

	This feature is useful (along with the reuse of the Symbolic object) if
	you need to factorize a sequence of triplet matrices with identical
	nonzero pattern (the order of the triplets in the Ti,Tj,Tx arrays must
	also remain unchanged).  It is faster than calling this routine for
	each matrix, and requires no workspace.
*/
