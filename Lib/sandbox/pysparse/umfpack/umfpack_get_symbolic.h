/* ========================================================================== */
/* === umfpack_get_symbolic ================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_get_symbolic
(
    int *n_row,
    int *n_col,
    int *n1,
    int *nz,
    int *nfr,
    int *nchains,
    int P [ ],
    int Q [ ],
    int Front_npivcol [ ],
    int Front_parent [ ],
    int Front_1strow [ ],
    int Front_leftmostdesc [ ],
    int Chain_start [ ],
    int Chain_maxrows [ ],
    int Chain_maxcols [ ],
    void *Symbolic
) ;

long umfpack_dl_get_symbolic
(
    long *n_row,
    long *n_col,
    long *n1,
    long *nz,
    long *nfr,
    long *nchains,
    long P [ ],
    long Q [ ],
    long Front_npivcol [ ],
    long Front_parent [ ],
    long Front_1strow [ ],
    long Front_leftmostdesc [ ],
    long Chain_start [ ],
    long Chain_maxrows [ ],
    long Chain_maxcols [ ],
    void *Symbolic
) ;

int umfpack_zi_get_symbolic
(
    int *n_row,
    int *n_col,
    int *n1,
    int *nz,
    int *nfr,
    int *nchains,
    int P [ ],
    int Q [ ],
    int Front_npivcol [ ],
    int Front_parent [ ],
    int Front_1strow [ ],
    int Front_leftmostdesc [ ],
    int Chain_start [ ],
    int Chain_maxrows [ ],
    int Chain_maxcols [ ],
    void *Symbolic
) ;

long umfpack_zl_get_symbolic
(
    long *n_row,
    long *n_col,
    long *n1,
    long *nz,
    long *nfr,
    long *nchains,
    long P [ ],
    long Q [ ],
    long Front_npivcol [ ],
    long Front_parent [ ],
    long Front_1strow [ ],
    long Front_leftmostdesc [ ],
    long Chain_start [ ],
    long Chain_maxrows [ ],
    long Chain_maxcols [ ],
    void *Symbolic
) ;

/*

double int Syntax:

    #include "umfpack.h"
    int status, n_row, n_col, nz, nfr, nchains, *P, *Q,
	*Front_npivcol, *Front_parent, *Front_1strow, *Front_leftmostdesc,
	*Chain_start, *Chain_maxrows, *Chain_maxcols ;
    void *Symbolic ;
    status = umfpack_di_get_symbolic (&n_row, &n_col, &nz, &nfr, &nchains,
	P, Q, Front_npivcol, Front_parent, Front_1strow,
	Front_leftmostdesc, Chain_start, Chain_maxrows, Chain_maxcols,
	Symbolic) ;

double long Syntax:

    #include "umfpack.h"
    long status, n_row, n_col, nz, nfr, nchains, *P, *Q,
	*Front_npivcol, *Front_parent, *Front_1strow, *Front_leftmostdesc,
	*Chain_start, *Chain_maxrows, *Chain_maxcols ;
    void *Symbolic ;
    status = umfpack_dl_get_symbolic (&n_row, &n_col, &nz, &nfr, &nchains,
	P, Q, Front_npivcol, Front_parent, Front_1strow,
	Front_leftmostdesc, Chain_start, Chain_maxrows, Chain_maxcols,
	Symbolic) ;

complex int Syntax:

    #include "umfpack.h"
    int status, n_row, n_col, nz, nfr, nchains, *P, *Q,
	*Front_npivcol, *Front_parent, *Front_1strow, *Front_leftmostdesc,
	*Chain_start, *Chain_maxrows, *Chain_maxcols ;
    void *Symbolic ;
    status = umfpack_zi_get_symbolic (&n_row, &n_col, &nz, &nfr, &nchains,
	P, Q, Front_npivcol, Front_parent, Front_1strow,
	Front_leftmostdesc, Chain_start, Chain_maxrows, Chain_maxcols,
	Symbolic) ;

complex long Syntax:

    #include "umfpack.h"
    long status, n_row, n_col, nz, nfr, nchains, *P, *Q,
	*Front_npivcol, *Front_parent, *Front_1strow, *Front_leftmostdesc,
	*Chain_start, *Chain_maxrows, *Chain_maxcols ;
    void *Symbolic ;
    status = umfpack_zl_get_symbolic (&n_row, &n_col, &nz, &nfr, &nchains,
	P, Q, Front_npivcol, Front_parent, Front_1strow,
	Front_leftmostdesc, Chain_start, Chain_maxrows, Chain_maxcols,
	Symbolic) ;

Purpose:

    Copies the contents of the Symbolic object into simple integer arrays
    accessible to the user.  This routine is not needed to factorize and/or
    solve a sparse linear system using UMFPACK.  Note that the output arrays
    P, Q, Front_npivcol, Front_parent, Front_1strow, Front_leftmostdesc,
    Chain_start, Chain_maxrows, and Chain_maxcols are not allocated by
    umfpack_*_get_symbolic; they must exist on input.

    All output arguments are optional.  If any of them are NULL
    on input, then that part of the symbolic analysis is not copied.  You can
    use this routine to extract just the parts of the symbolic analysis that
    you want.  For example, to retrieve just the column permutation Q, use:

    #define noI (int *) NULL
    status = umfpack_di_get_symbolic (noI, noI, noI, noI, noI, noI, noI,
	    Q, noI, noI, noI, noI, noI, noI, noI, Symbolic) ;

    The only required argument the last one, the pointer to the Symbolic object.

    The Symbolic object is small.  Its size for an n-by-n square matrix varies
    from 4*n to 13*n, depending on the matrix.  The object holds the initial
    column permutation, the supernodal column elimination tree, and information
    about each frontal matrix.  You can print it with umfpack_*_report_symbolic.

Returns:

    Returns UMFPACK_OK if successful, UMFPACK_ERROR_invalid_Symbolic_object
    if Symbolic is an invalid object.

Arguments:

    Int *n_row ;	Output argument.
    Int *n_col ;	Output argument.

	The dimensions of the matrix A analyzed by the call to
	umfpack_*_symbolic that generated the Symbolic object.

    Int *n1 ;		Output argument.

	The number of pivots with zero Markowitz cost (they have just one entry
	in the pivot row, or the pivot column, or both).  These appear first in
	the output permutations P and Q.

	NOTE: this argument is new to version 4.1.

    Int *nz ;		Output argument.

	The number of nonzeros in A.

    Int *nfr ;	Output argument.

	The number of frontal matrices that will be used by umfpack_*_numeric
	to factorize the matrix A.  It is in the range 0 to n_col.

    Int *nchains ;	Output argument.

	The frontal matrices are related to one another by the supernodal
	column elimination tree.  Each node in this tree is one frontal matrix.
	The tree is partitioned into a set of disjoint paths, and a frontal
	matrix chain is one path in this tree.  Each chain is factorized using
	a unifrontal technique, with a single working array that holds each
	frontal matrix in the chain, one at a time.  nchains is in the range
	0 to nfr.

    Int P [n_row] ;	Output argument.

	The initial row permutation.  If P [k] = i, then this means that
	row i is the kth row in the pre-ordered matrix.  In general, this P is
	not the same as the final row permutation computed by umfpack_*_numeric.

	For the unsymmetric strategy, P defines the row-merge order.  Let j be
	the column index of the leftmost nonzero entry in row i of A*Q.  Then
	P defines a sort of the rows according to this value.  A row can appear
	earlier in this ordering if it is aggressively absorbed before it can
	become a pivot row.  If P [k] = i, row i typically will not be the kth
	pivot row.

	For the symmetric strategy, P = Q.  For the 2-by-2 strategy, P is the
	row permutation that places large entries on the diagonal of P*A*Q.
	If no pivoting occurs during numerical factorization, P [k] = i also
	defines the final permutation of umfpack_*_numeric, for either the
	symmetric or 2-by-2 strategies.

    Int Q [n_col] ;	Output argument.

	The initial column permutation.  If Q [k] = j, then this means that
	column j is the kth pivot column in the pre-ordered matrix.  Q is
	not necessarily the same as the final column permutation Q, computed by
	umfpack_*_numeric.  The numeric factorization may reorder the pivot
	columns within each frontal matrix to reduce fill-in.  If the matrix is
	structurally singular, and if the symmetric or 2-by-2 strategies or
	used (or if Control [UMFPACK_FIXQ] > 0), then this Q will be the same
	as the final column permutation computed in umfpack_*_numeric.

    Int Front_npivcol [n_col+1] ;	Output argument.

	This array should be of size at least n_col+1, in order to guarantee
	that it will be large enough to hold the output.  Only the first nfr+1
	entries are used, however.

	The kth frontal matrix holds Front_npivcol [k] pivot columns.  Thus, the
	first frontal matrix, front 0, is used to factorize the first
	Front_npivcol [0] columns; these correspond to the original columns
	Q [0] through Q [Front_npivcol [0]-1].  The next frontal matrix
	is used to factorize the next Front_npivcol [1] columns, which are thus
	the original columns Q [Front_npivcol [0]] through
	Q [Front_npivcol [0] + Front_npivcol [1] - 1], and so on.  Columns
	with no entries at all are put in a placeholder "front",
	Front_npivcol [nfr].  The sum of Front_npivcol [0..nfr] is equal to
	n_col.

	Any modifications that umfpack_*_numeric makes to the initial column
	permutation are constrained to within each frontal matrix.  Thus, for
	the first frontal matrix, Q [0] through Q [Front_npivcol [0]-1] is some
	permutation of the columns Q [0] through
	Q [Front_npivcol [0]-1].  For second frontal matrix,
	Q [Front_npivcol [0]] through Q [Front_npivcol [0] + Front_npivcol[1]-1]
	is some permutation of the same portion of Q, and so on.  All pivot
	columns are numerically factorized within the frontal matrix originally
	determined by the symbolic factorization; there is no delayed pivoting
	across frontal matrices.

    Int Front_parent [n_col+1] ;	Output argument.

	This array should be of size at least n_col+1, in order to guarantee
	that it will be large enough to hold the output.  Only the first nfr+1
	entries are used, however.

	Front_parent [0..nfr] holds the supernodal column elimination tree
	(including the placeholder front nfr, which may be empty).  Each node in
	the tree corresponds to a single frontal matrix.  The parent of node f
	is Front_parent [f].

    Int Front_1strow [n_col+1] ;	Output argument.

	This array should be of size at least n_col+1, in order to guarantee
	that it will be large enough to hold the output.  Only the first nfr+1
	entries are used, however.

	Front_1strow [k] is the row index of the first row in A (P,Q)
	whose leftmost entry is in a pivot column for the kth front.  This is
	necessary only to properly factorize singular matrices.  It is new to
	Version 4.0.  Rows in the range Front_1strow [k] to
	Front_1strow [k+1]-1 first become pivot row candidates at the kth front.
	Any rows not eliminated in the kth front may be selected as pivot rows
	in the parent of k (Front_parent [k]) and so on up the tree.

    Int Front_leftmostdesc [n_col+1] ;	Output argument.

	This array should be of size at least n_col+1, in order to guarantee
	that it will be large enough to hold the output.  Only the first nfr+1
	entries are used, however.

	Front_leftmostdesc [k] is the leftmost descendant of front k, or k
	if the front has no children in the tree.  Since the rows and columns
	(P and Q) have been post-ordered via a depth-first-search of
	the tree, rows in the range Front_1strow [Front_leftmostdesc [k]] to
	Front_1strow [k+1]-1 form the entire set of candidate pivot rows for
	the kth front (some of these will typically have already been selected
	by fronts in the range Front_leftmostdesc [k] to front k-1, before
	the factorization reaches front k).

    Chain_start [n_col+1] ;	Output argument.

	This array should be of size at least n_col+1, in order to guarantee
	that it will be large enough to hold the output.  Only the first
	nchains+1 entries are used, however.

	The kth frontal matrix chain consists of frontal matrices Chain_start[k]
	through Chain_start [k+1]-1.  Thus, Chain_start [0] is always 0, and
	Chain_start [nchains] is the total number of frontal matrices, nfr.  For
	two adjacent fronts f and f+1 within a single chain, f+1 is always the
	parent of f (that is, Front_parent [f] = f+1).

    Int Chain_maxrows [n_col+1] ;	Output argument.
    Int Chain_maxcols [n_col+1] ;	Output argument.

	These arrays should be of size at least n_col+1, in order to guarantee
	that they will be large enough to hold the output.  Only the first
	nchains entries are used, however.

	The kth frontal matrix chain requires a single working array of
	dimension Chain_maxrows [k] by Chain_maxcols [k], for the unifrontal
	technique that factorizes the frontal matrix chain.  Since the symbolic
	factorization only provides an upper bound on the size of each frontal
	matrix, not all of the working array is necessarily used during the
	numerical factorization.

	Note that the upper bound on the number of rows and columns of each
	frontal matrix is computed by umfpack_*_symbolic, but all that is
	required by umfpack_*_numeric is the maximum of these two sets of
	values for each frontal matrix chain.  Thus, the size of each
	individual frontal matrix is not preserved in the Symbolic object.

    void *Symbolic ;			Input argument, not modified.

	The Symbolic object, which holds the symbolic factorization computed by
	umfpack_*_symbolic.  The Symbolic object is not modified by
	umfpack_*_get_symbolic.
*/
