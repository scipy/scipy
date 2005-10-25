/* ========================================================================== */
/* === umfpack_symbolic ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_symbolic
(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

long umfpack_dl_symbolic
(
    long n_row,
    long n_col,
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

int umfpack_zi_symbolic
(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ], const double Az [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

long umfpack_zl_symbolic
(
    long n_row,
    long n_col,
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ], const double Az [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    int n_row, n_col, *Ap, *Ai, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO], *Ax ;
    status = umfpack_di_symbolic (n_row, n_col, Ap, Ai, Ax,
	&Symbolic, Control, Info) ;

double long Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    long n_row, n_col, *Ap, *Ai, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO], *Ax ;
    status = umfpack_dl_symbolic (n_row, n_col, Ap, Ai, Ax,
	&Symbolic, Control, Info) ;

complex int Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    int n_row, n_col, *Ap, *Ai, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO], *Ax, *Az ;
    status = umfpack_zi_symbolic (n_row, n_col, Ap, Ai, Ax, Az,
	&Symbolic, Control, Info) ;

complex long Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    long n_row, n_col, *Ap, *Ai, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO], *Ax, *Az ;
    status = umfpack_zl_symbolic (n_row, n_col, Ap, Ai, Ax, Az,
	&Symbolic, Control, Info) ;

Purpose:

    Given nonzero pattern of a sparse matrix A in column-oriented form,
    umfpack_*_symbolic performs a column pre-ordering to reduce fill-in
    (using COLAMD or AMD) and a symbolic factorization.  This is required
    before the matrix can be numerically factorized with umfpack_*_numeric.
    If you wish to bypass the COLAMD or AMD pre-ordering and provide your own
    ordering, use umfpack_*_qsymbolic instead.

    Since umfpack_*_symbolic and umfpack_*_qsymbolic are very similar, options
    for both routines are discussed below.

    For the following discussion, let S be the submatrix of A obtained after
    eliminating all pivots of zero Markowitz cost.  S has dimension
    (n_row-n1-nempty_row) -by- (n_col-n1-nempty_col), where
    n1 = Info [UMFPACK_COL_SINGLETONS] + Info [UMFPACK_ROW_SINGLETONS],
    nempty_row = Info [UMFPACK_NEMPTY_ROW] and
    nempty_col = Info [UMFPACK_NEMPTY_COL].

Returns:

    The status code is returned.  See Info [UMFPACK_STATUS], below.

Arguments:

    Int n_row ;		Input argument, not modified.
    Int n_col ;		Input argument, not modified.

	A is an n_row-by-n_col matrix.  Restriction: n_row > 0 and n_col > 0.

    Int Ap [n_col+1] ;	Input argument, not modified.

	Ap is an integer array of size n_col+1.  On input, it holds the
	"pointers" for the column form of the sparse matrix A.  Column j of
	the matrix A is held in Ai [(Ap [j]) ... (Ap [j+1]-1)].  The first
	entry, Ap [0], must be zero, and Ap [j] <= Ap [j+1] must hold for all
	j in the range 0 to n_col-1.  The value nz = Ap [n_col] is thus the
	total number of entries in the pattern of the matrix A.  nz must be
	greater than or equal to zero.

    Int Ai [nz] ;	Input argument, not modified, of size nz = Ap [n_col].

	The nonzero pattern (row indices) for column j is stored in
	Ai [(Ap [j]) ... (Ap [j+1]-1)].  The row indices in a given column j
	must be in ascending order, and no duplicate row indices may be present.
	Row indices must be in the range 0 to n_row-1 (the matrix is 0-based).
	See umfpack_*_triplet_to_col for how to sort the columns of a matrix
	and sum up the duplicate entries.  See umfpack_*_report_matrix for how
	to print the matrix A.

    double Ax [nz] ;	Optional input argument, not modified.

	The numerical values of the sparse matrix A.  The nonzero pattern (row
	indices) for column j is stored in Ai [(Ap [j]) ... (Ap [j+1]-1)], and
	the corresponding numerical values are stored in
	Ax [(Ap [j]) ... (Ap [j+1]-1)].  Used only by the 2-by-2 strategy to
	determine whether entries are "large" or "small".  You do not have to
	pass the same numerical values to umfpack_*_numeric.  If Ax is not
	present (a (double *) NULL pointer), then any entry in A is assumed to
	be "large".

    double Az [nz] ;	Optional input argument, not modified, for complex
			versions.

	For the complex versions, this holds the imaginary part of A.  The
	imaginary part of column j is held in Az [(Ap [j]) ... (Ap [j+1]-1)].

	Future complex version:  if Ax is present and Az is NULL, then both real
	and imaginary parts will be contained in Ax[0..2*nz-1], with Ax[2*k]
	and Ax[2*k+1] being the real and imaginary part of the kth entry.

	Used by the 2-by-2 strategy only.  See the description of Ax, above.

    void **Symbolic ;	Output argument.

	**Symbolic is the address of a (void *) pointer variable in the user's
	calling routine (see Syntax, above).  On input, the contents of this
	variable are not defined.  On output, this variable holds a (void *)
	pointer to the Symbolic object (if successful), or (void *) NULL if
	a failure occurred.

    double Control [UMFPACK_CONTROL] ;	Input argument, not modified.

	If a (double *) NULL pointer is passed, then the default control
	settings are used (the defaults are suitable for all matrices,
	ranging from those with highly unsymmetric nonzero pattern, to
	symmetric matrices).  Otherwise, the settings are determined from the
	Control array.  See umfpack_*_defaults on how to fill the Control
	array with the default settings.  If Control contains NaN's, the
	defaults are used.  The following Control parameters are used:

	Control [UMFPACK_STRATEGY]:  This is the most important control
	    parameter.  It determines what kind of ordering and pivoting
	    strategy that UMFPACK should use.  It is new to Version 4.1
	    There are 4 options:

	    UMFPACK_STRATEGY_AUTO:  This is the default.  The input matrix is
		analyzed to determine how symmetric the nonzero pattern is, and
		how many entries there are on the diagonal.  It then selects one
		of the following strategies.  Refer to the User Guide for a
		description of how the strategy is automatically selected.

	    UMFPACK_STRATEGY_UNSYMMETRIC:  Use the unsymmetric strategy.  COLAMD
		is used to order the columns of A, followed by a postorder of
		the column elimination tree.  No attempt is made to perform
		diagonal pivoting.  The column ordering is refined during
		factorization.  This strategy was the only one provided with
		UMFPACK V4.0.

		In the numerical factorization, the
		Control [UMFPACK_SYM_PIVOT_TOLERANCE] parameter is ignored.  A
		pivot is selected if its magnitude is >=
		Control [UMFPACK_PIVOT_TOLERANCE] (default 0.1) times the
		largest entry in its column.

	    UMFPACK_STRATEGY_SYMMETRIC:  Use the symmetric strategy (new to
		Version 4.1).  In this method, the approximate minimum degree
		ordering (AMD) is applied to A+A', followed by a postorder of
		the elimination tree of A+A'.  UMFPACK attempts to perform
		diagonal pivoting during numerical factorization.  No refinement
		of the column pre-ordering is performed during factorization.

		In the numerical factorization, a nonzero entry on the diagonal
		is selected as the pivot if its magnitude is >= Control
		[UMFPACK_SYM_PIVOT_TOLERANCE] (default 0.001) times the largest
		entry in its column.  If this is not acceptable, then an
		off-diagonal pivot is selected with magnitude >= Control
		[UMFPACK_PIVOT_TOLERANCE] (default 0.1) times the largest entry
		in its column.

	    UMFPACK_STRATEGY_2BY2:  a row permutation P2 is found that places
		large entries on the diagonal.  The matrix P2*A is then
		factorized using the symmetric strategy, described above.
		Refer to the User Guide for more information.

	Control [UMFPACK_DENSE_COL]:
	    If COLAMD is used, columns with more than
	    max (16, Control [UMFPACK_DENSE_COL] * 16 * sqrt (n_row)) entries
	    are placed placed last in the column pre-ordering.  Default: 0.2.

	Control [UMFPACK_DENSE_ROW]:
	    Rows with more than max (16, Control [UMFPACK_DENSE_ROW] * 16 *
	    sqrt (n_col)) entries are treated differently in the COLAMD
	    pre-ordering, and in the internal data structures during the
	    subsequent numeric factorization.  Default: 0.2.

	Control [UMFPACK_AMD_DENSE]:  rows/columns in A+A' with more than
	    max (16, Control [UMFPACK_AMD_DENSE] * sqrt (n)) entries
	    (where n = n_row = n_col) are ignored in the AMD pre-ordering.
	    Default: 10.

	Control [UMFPACK_BLOCK_SIZE]:  the block size to use for Level-3 BLAS
	    in the subsequent numerical factorization (umfpack_*_numeric).
	    A value less than 1 is treated as 1.  Default: 32.  Modifying this
	    parameter affects when updates are applied to the working frontal
	    matrix, and can indirectly affect fill-in and operation count.
	    As long as the block size is large enough (8 or so), this parameter
	    has a modest effect on performance. 

	Control [UMFPACK_2BY2_TOLERANCE]:  a diagonal entry S (k,k) is
	    considered "small" if it is < tol * max (abs (S (:,k))), where S a
	    submatrix of the scaled input matrix, with pivots of zero Markowitz
	    cost removed.

	Control [UMFPACK_SCALE]:  This parameter is new to V4.1.  See
	    umfpack_numeric.h for a description.  Only affects the 2-by-2
	    strategy.  Default: UMFPACK_SCALE_SUM.

	Control [UMFPACK_FIXQ]:  If > 0, then the pre-ordering Q is not modified
	    during numeric factorization.  If < 0, then Q may be modified.  If
	    zero, then this is controlled automatically (the unsymmetric
	    strategy modifies Q, the others do not).  Default: 0.

	Control [UMFPACK_AGGRESSIVE]:  If nonzero, aggressive absorption is used
	    in COLAMD and AMD.  Default: 1.

    double Info [UMFPACK_INFO] ;	Output argument, not defined on input.

	Contains statistics about the symbolic analysis.  If a (double *) NULL
	pointer is passed, then no statistics are returned in Info (this is not
	an error condition).  The entire Info array is cleared (all entries set
	to -1) and then the following statistics are computed:

	Info [UMFPACK_STATUS]: status code.  This is also the return value,
	    whether or not Info is present.

	    UMFPACK_OK

		Each column of the input matrix contained row indices
		in increasing order, with no duplicates.  Only in this case
		does umfpack_*_symbolic compute a valid symbolic factorization.
		For the other cases below, no Symbolic object is created
		(*Symbolic is (void *) NULL).

	    UMFPACK_ERROR_n_nonpositive

		n is less than or equal to zero.

	    UMFPACK_ERROR_invalid_matrix

		Number of entries in the matrix is negative, Ap [0] is nonzero,
		a column has a negative number of entries, a row index is out of
		bounds, or the columns of input matrix were jumbled (unsorted
		columns or duplicate entries).

	    UMFPACK_ERROR_out_of_memory

		Insufficient memory to perform the symbolic analysis.  If the
		analysis requires more than 2GB of memory and you are using
		the 32-bit ("int") version of UMFPACK, then you are guaranteed
		to run out of memory.  Try using the 64-bit version of UMFPACK.

	    UMFPACK_ERROR_argument_missing

		One or more required arguments is missing.

	    UMFPACK_ERROR_internal_error

		Something very serious went wrong.  This is a bug.
		Please contact the author (davis@cise.ufl.edu).

	    Note that the UMFPACK_ERROR_problem_too_large error code is no
	    longer returned (it was in Version 4.0).

	Info [UMFPACK_NROW]:  the value of the input argument n_row.

	Info [UMFPACK_NCOL]:  the value of the input argument n_col.

	Info [UMFPACK_NZ]:  the number of entries in the input matrix
	    (Ap [n_col]).

	Info [UMFPACK_SIZE_OF_UNIT]:  the number of bytes in a Unit,
	    for memory usage statistics below.

	Info [UMFPACK_SIZE_OF_INT]:  the number of bytes in an int.

	Info [UMFPACK_SIZE_OF_LONG]:  the number of bytes in a long.

	Info [UMFPACK_SIZE_OF_POINTER]:  the number of bytes in a void *
	    pointer.

	Info [UMFPACK_SIZE_OF_ENTRY]:  the number of bytes in a numerical entry.

	Info [UMFPACK_NDENSE_ROW]:  number of "dense" rows in A.  These rows are
	    ignored when the column pre-ordering is computed in COLAMD.  They
	    are also treated differently during numeric factorization.  If > 0,
	    then the matrix had to be re-analyzed by UMF_analyze, which does
	    not ignore these rows.

	Info [UMFPACK_NEMPTY_ROW]:  number of "empty" rows in A, as determined
	    These are rows that either have no entries, or whose entries are
	    all in pivot columns of zero-Markowitz-cost pivots.

	Info [UMFPACK_NDENSE_COL]:  number of "dense" columns in A.  COLAMD
	    orders these columns are ordered last in the factorization, but
	    before "empty" columns.

	Info [UMFPACK_NEMPTY_COL]:  number of "empty" columns in A.  These are
	    columns that either have no entries, or whose entries are all in
	    pivot rows of zero-Markowitz-cost pivots.  These columns are
	    ordered last in the factorization, to the right of "dense" columns.

	Info [UMFPACK_SYMBOLIC_DEFRAG]:  number of garbage collections
	    performed during ordering and symbolic pre-analysis.

	Info [UMFPACK_SYMBOLIC_PEAK_MEMORY]:  the amount of memory (in Units)
	    required for umfpack_*_symbolic to complete.  This count includes
	    the size of the Symbolic object itself, which is also reported in
	    Info [UMFPACK_SYMBOLIC_SIZE].

	Info [UMFPACK_SYMBOLIC_SIZE]: the final size of the Symbolic object (in
	    Units).  This is fairly small, roughly 2*n to 13*n integers,
	    depending on the matrix.

	Info [UMFPACK_VARIABLE_INIT_ESTIMATE]: the Numeric object contains two
	    parts.  The first is fixed in size (O (n_row+n_col)).  The
	    second part holds the sparse LU factors and the contribution blocks
	    from factorized frontal matrices.  This part changes in size during
	    factorization.  Info [UMFPACK_VARIABLE_INIT_ESTIMATE] is the exact
	    size (in Units) required for this second variable-sized part in
	    order for the numerical factorization to start.

	Info [UMFPACK_VARIABLE_PEAK_ESTIMATE]: the estimated peak size (in
	    Units) of the variable-sized part of the Numeric object.  This is
	    usually an upper bound, but that is not guaranteed. 

	Info [UMFPACK_VARIABLE_FINAL_ESTIMATE]: the estimated final size (in
	    Units) of the variable-sized part of the Numeric object.  This is
	    usually an upper bound, but that is not guaranteed.  It holds just
	    the sparse LU factors.

	Info [UMFPACK_NUMERIC_SIZE_ESTIMATE]:  an estimate of the final size (in
	    Units) of the entire Numeric object (both fixed-size and variable-
	    sized parts), which holds the LU factorization (including the L, U,
	    P and Q matrices).

	Info [UMFPACK_PEAK_MEMORY_ESTIMATE]:  an estimate of the total amount of
	    memory (in Units) required by umfpack_*_symbolic and
	    umfpack_*_numeric to perform both the symbolic and numeric
	    factorization.  This is the larger of the amount of memory needed
	    in umfpack_*_numeric itself, and the amount of memory needed in
	    umfpack_*_symbolic (Info [UMFPACK_SYMBOLIC_PEAK_MEMORY]).  The
	    count includes the size of both the Symbolic and Numeric objects
	    themselves.  It can be a very loose upper bound, particularly when
	    the symmetric or 2-by-2 strategies are used.

	Info [UMFPACK_FLOPS_ESTIMATE]:  an estimate of the total floating-point
	    operations required to factorize the matrix.  This is a "true"
	    theoretical estimate of the number of flops that would be performed
	    by a flop-parsimonious sparse LU algorithm.  It assumes that no
	    extra flops are performed except for what is strictly required to
	    compute the LU factorization.  It ignores, for example, the flops
            performed by umfpack_di_numeric to add contribution blocks of
	    frontal matrices together.  If L and U are the upper bound on the
	    pattern of the factors, then this flop count estimate can be
	    represented in MATLAB (for real matrices, not complex) as:

		Lnz = full (sum (spones (L))) - 1 ;	% nz in each col of L
		Unz = full (sum (spones (U')))' - 1 ;	% nz in each row of U
		flops = 2*Lnz*Unz + sum (Lnz) ;

	    The actual "true flop" count found by umfpack_*_numeric will be
	    less than this estimate.

	    For the real version, only (+ - * /) are counted.  For the complex
	    version, the following counts are used:

		operation	flops
	    	c = 1/b		6
		c = a*b		6
		c -= a*b	8

	Info [UMFPACK_LNZ_ESTIMATE]:  an estimate of the number of nonzeros in
	    L, including the diagonal.  Since L is unit-diagonal, the diagonal
	    of L is not stored.  This estimate is a strict upper bound on the
	    actual nonzeros in L to be computed by umfpack_*_numeric.

	Info [UMFPACK_UNZ_ESTIMATE]:  an estimate of the number of nonzeros in
	    U, including the diagonal.  This estimate is a strict upper bound on
	    the actual nonzeros in U to be computed by umfpack_*_numeric.

	Info [UMFPACK_MAX_FRONT_SIZE_ESTIMATE]: estimate of the size of the
	    largest frontal matrix (# of entries), for arbitrary partial
	    pivoting during numerical factorization.

	Info [UMFPACK_SYMBOLIC_TIME]:  The CPU time taken, in seconds.

	------------------------------------------------------------------------
	The rest of the statistics are new to Version 4.1:
	------------------------------------------------------------------------

	Info [UMFPACK_SYMBOLIC_WALLTIME]:  The wallclock time taken, in seconds.

	Info [UMFPACK_STRATEGY_USED]: The ordering strategy used:
	    UMFPACK_STRATEGY_SYMMETRIC, UMFPACK_STRATEGY_UNSYMMETRIC, or
	    UMFPACK_STRATEGY_2BY2.

	Info [UMFPACK_ORDERING_USED]:  The ordering method used:
	    UMFPACK_ORDERING_COLAMD or UMFPACK_ORDERING_AMD.  It can be
	    UMFPACK_ORDERING_GIVEN for umfpack_*_qsymbolic.

	Info [UMFPACK_QFIXED]: 1 if the column pre-ordering will be refined
	    during numerical factorization, 0 if not.

	Info [UMFPACK_DIAG_PREFERED]: 1 if diagonal pivoting will be attempted,
	    0 if not.

	Info [UMFPACK_COL_SINGLETONS]:  the matrix A is analyzed by first
	    eliminating all pivots with zero Markowitz cost.  This count is the
	    number of these pivots with exactly one nonzero in their pivot
	    column.

	Info [UMFPACK_ROW_SINGLETONS]:  the number of zero-Markowitz-cost
	    pivots with exactly one nonzero in their pivot row.

	Info [UMFPACK_PATTERN_SYMMETRY]: the symmetry of the pattern of S.

	Info [UMFPACK_NZ_A_PLUS_AT]: the number of off-diagonal entries in S+S'.

	Info [UMFPACK_NZDIAG]:  the number of entries on the diagonal of S.

	Info [UMFPACK_N2]:  if S is square, and nempty_row = nempty_col, this
	    is equal to n_row - n1 - nempty_row.

	Info [UMFPACK_S_SYMMETRIC]: 1 if S is square and its diagonal has been
	    preserved, 0 otherwise.


	Info [UMFPACK_MAX_FRONT_NROWS_ESTIMATE]: estimate of the max number of
	    rows in any frontal matrix, for arbitrary partial pivoting.

	Info [UMFPACK_MAX_FRONT_NCOLS_ESTIMATE]: estimate of the max number of
	    columns in any frontal matrix, for arbitrary partial pivoting.

	------------------------------------------------------------------------
	The next four statistics are computed only if AMD is used:
	------------------------------------------------------------------------

	Info [UMFPACK_SYMMETRIC_LUNZ]: The number of nonzeros in L and U,
	    assuming no pivoting during numerical factorization, and assuming a
	    zero-free diagonal of U.  Excludes the entries on the diagonal of
	    L.  If the matrix has a purely symmetric nonzero pattern, this is
	    often a lower bound on the nonzeros in the actual L and U computed
	    in the numerical factorization, for matrices that fit the criteria
	    for the "symmetric" strategy.

	Info [UMFPACK_SYMMETRIC_FLOPS]: The floating-point operation count in
	    the numerical factorization phase, assuming no pivoting.  If the
	    pattern of the matrix is symmetric, this is normally a lower bound
	    on the floating-point operation count in the actual numerical
	    factorization, for matrices that fit the criteria for the symmetric
	    or 2-by-2 strategies

	Info [UMFPACK_SYMMETRIC_NDENSE]: The number of "dense" rows/columns of
	    S+S' that were ignored during the AMD ordering.  These are placed
	    last in the output order.  If > 0, then the
	    Info [UMFPACK_SYMMETRIC_*] statistics, above are rough upper bounds.

	Info [UMFPACK_SYMMETRIC_DMAX]: The maximum number of nonzeros in any
	    column of L, if no pivoting is performed during numerical
	    factorization.  Excludes the part of the LU factorization for
	    pivots with zero Markowitz cost.

	------------------------------------------------------------------------
	The following statistics are computed only if the 2-by-2 strategy is
	used or attempted:
	------------------------------------------------------------------------

	Info [UMFPACK_2BY2_NWEAK]: the number of small diagonal entries in S.

	Info [UMFPACK_2BY2_UNMATCHED]: the number of small diagonal entries
	    in P2*S.

	Info [UMFPACK_2BY2_PATTERN_SYMMETRY]: the symmetry of P2*S.

	Info [UMFPACK_2BY2_NZ_PA_PLUS_AT]:  the number of off-diagonal entries
	    in (P2*S)+(P2*S)'.

	Info [UMFPACK_2BY2_NZDIAG]:  the number of nonzero entries on the
	    diagonal of P2*S.


	At the start of umfpack_*_symbolic, all of Info is set of -1, and then
	after that only the above listed Info [...] entries are accessed.
	Future versions might modify different parts of Info.
*/
