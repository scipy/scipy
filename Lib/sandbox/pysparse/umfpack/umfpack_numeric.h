/* ========================================================================== */
/* === umfpack_numeric ====================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_numeric
(
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    void *Symbolic,
    void **Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

long umfpack_dl_numeric
(
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ],
    void *Symbolic,
    void **Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

int umfpack_zi_numeric
(
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ], const double Az [ ],
    void *Symbolic,
    void **Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

long umfpack_zl_numeric
(
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ], const double Az [ ],
    void *Symbolic,
    void **Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    void *Symbolic, *Numeric ;
    int *Ap, *Ai, status ;
    double *Ax, Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);

double long Syntax:

    #include "umfpack.h"
    void *Symbolic, *Numeric ;
    long *Ap, *Ai, status ;
    double *Ax, Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    status = umfpack_dl_numeric (Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);

complex int Syntax:

    #include "umfpack.h"
    void *Symbolic, *Numeric ;
    int *Ap, *Ai, status ;
    double *Ax, *Az, Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    status = umfpack_zi_numeric (Ap, Ai, Ax, Az, Symbolic, &Numeric,
	Control, Info) ;

complex long Syntax:

    #include "umfpack.h"
    void *Symbolic, *Numeric ;
    long *Ap, *Ai, status ;
    double *Ax, *Az, Control [UMFPACK_CONTROL], Info [UMFPACK_INFO] ;
    status = umfpack_zl_numeric (Ap, Ai, Ax, Symbolic, &Numeric,
	Control, Info) ;

Purpose:

    Given a sparse matrix A in column-oriented form, and a symbolic analysis
    computed by umfpack_*_*symbolic, the umfpack_*_numeric routine performs the
    numerical factorization, PAQ=LU, PRAQ=LU, or P(R\A)Q=LU, where P and Q are
    permutation matrices (represented as permutation vectors), R is the row
    scaling, L is unit-lower triangular, and U is upper triangular.  This is
    required before the system Ax=b (or other related linear systems) can be
    solved.  umfpack_*_numeric can be called multiple times for each call to
    umfpack_*_*symbolic, to factorize a sequence of matrices with identical
    nonzero pattern.  Simply compute the Symbolic object once, with
    umfpack_*_*symbolic, and reuse it for subsequent matrices.  This routine
    safely detects if the pattern changes, and sets an appropriate error code.

Returns:

    The status code is returned.  See Info [UMFPACK_STATUS], below.

Arguments:

    Int Ap [n_col+1] ;	Input argument, not modified.

	This must be identical to the Ap array passed to umfpack_*_*symbolic.
	The value of n_col is what was passed to umfpack_*_*symbolic (this is
	held in the Symbolic object).

    Int Ai [nz] ;	Input argument, not modified, of size nz = Ap [n_col].

	This must be identical to the Ai array passed to umfpack_*_*symbolic.

    double Ax [nz] ;	Input argument, not modified, of size nz = Ap [n_col].

	The numerical values of the sparse matrix A.  The nonzero pattern (row
	indices) for column j is stored in Ai [(Ap [j]) ... (Ap [j+1]-1)], and
	the corresponding numerical values are stored in
	Ax [(Ap [j]) ... (Ap [j+1]-1)].

    double Az [nz] ;	Input argument, not modified, for complex versions.

	For the complex versions, this holds the imaginary part of A.  The
	imaginary part of column j is held in Az [(Ap [j]) ... (Ap [j+1]-1)].

	Future complex version:  if Ax is present and Az is NULL, then both real
	and imaginary parts will be contained in Ax[0..2*nz-1], with Ax[2*k]
	and Ax[2*k+1] being the real and imaginary part of the kth entry.

    void *Symbolic ;	Input argument, not modified.

	The Symbolic object, which holds the symbolic factorization computed by
	umfpack_*_*symbolic.  The Symbolic object is not modified by
	umfpack_*_numeric.

    void **Numeric ;	Output argument.

	**Numeric is the address of a (void *) pointer variable in the user's
	calling routine (see Syntax, above).  On input, the contents of this
	variable are not defined.  On output, this variable holds a (void *)
	pointer to the Numeric object (if successful), or (void *) NULL if
	a failure occurred.

    double Control [UMFPACK_CONTROL] ;   Input argument, not modified.

	If a (double *) NULL pointer is passed, then the default control
	settings are used.  Otherwise, the settings are determined from the
	Control array.  See umfpack_*_defaults on how to fill the Control
	array with the default settings.  If Control contains NaN's, the
	defaults are used.  The following Control parameters are used:

	Control [UMFPACK_PIVOT_TOLERANCE]:  relative pivot tolerance for
	    threshold partial pivoting with row interchanges.  In any given
	    column, an entry is numerically acceptable if its absolute value is
	    greater than or equal to Control [UMFPACK_PIVOT_TOLERANCE] times
	    the largest absolute value in the column.  A value of 1.0 gives true
	    partial pivoting.  If less than or equal to zero, then any nonzero
	    entry is numerically acceptable as a pivot (this is changed from
	    Version 4.0).  Default: 0.1.

	    Smaller values tend to lead to sparser LU factors, but the solution
	    to the linear system can become inaccurate.  Larger values can lead
	    to a more accurate solution (but not always), and usually an
	    increase in the total work.

	    For complex matrices, a cheap approximate of the absolute value
	    is used for the threshold partial pivoting test (|a_real| + |a_imag|
	    instead of the more expensive-to-compute exact absolute value
	    sqrt (a_real^2 + a_imag^2)).

	Control [UMFPACK_SYM_PIVOT_TOLERANCE]:  This parameter is new to V4.1.
	    If diagonal pivoting is attempted (the symmetric or symmetric-2by2
	    strategies are used) then this parameter is used to control when the
	    diagonal entry is selected in a given pivot column.  The absolute
	    value of the entry must be >= Control [UMFPACK_SYM_PIVOT_TOLERANCE]
	    times the largest absolute value in the column.  A value of zero
	    will ensure that no off-diagonal pivoting is performed, except that
	    zero diagonal entries are not selected if there are any off-diagonal
	    nonzero entries.

	    If an off-diagonal pivot is selected, an attempt is made to restore
	    symmetry later on.  Suppose A (i,j) is selected, where i != j.
	    If column i has not yet been selected as a pivot column, then
	    the entry A (j,i) is redefined as a "diagonal" entry, except that
	    the tighter tolerance (Control [UMFPACK_PIVOT_TOLERANCE]) is
	    applied.  This strategy has an effect similar to 2-by-2 pivoting
	    for symmetric indefinite matrices.  If a 2-by-2 block pivot with
	    nonzero structure

		       i j
		    i: 0 x
		    j: x 0

	    is selected in a symmetric indefinite factorization method, the
	    2-by-2 block is inverted and a rank-2 update is applied.  In
	    UMFPACK, this 2-by-2 block would be reordered as

		       j i
		    i: x 0
		    j: 0 x

	    In both cases, the symmetry of the Schur complement is preserved.

	Control [UMFPACK_SCALE]:  This parameter is new to V4.1.  Version 4.0
	    did not scale the matrix.  Note that the user's input matrix is
	    never modified, only an internal copy is scaled.

	    There are three valid settings for this parameter.  If any other
	    value is provided, the default is used.

	    UMFPACK_SCALE_NONE:  no scaling is performed.

	    UMFPACK_SCALE_SUM:  each row of the input matrix A is divided by
		the sum of the absolute values of the entries in that row.
		The scaled matrix has an infinity norm of 1.

	    UMFPACK_SCALE_MAX:  each row of the input matrix A is divided by
		the maximum the absolute values of the entries in that row.
		In the scaled matrix the largest entry in each row has
		a magnitude exactly equal to 1.

	    Note that for complex matrices, a cheap approximate absolute value
	    is used, |a_real| + |a_imag|, instead of the exact absolute value
	    sqrt ((a_real)^2 + (a_imag)^2).

	    Scaling is very important for the "symmetric" strategy when
	    diagonal pivoting is attempted.  It also improves the performance
	    of the "unsymmetric" strategy.

	    Default: UMFPACK_SCALE_SUM.

	Control [UMFPACK_ALLOC_INIT]:  This parameter has changed in V4.1.

	    When umfpack_*_numeric starts, it allocates memory for the Numeric
	    object.  Part of this is of fixed size (approximately n double's +
	    12*n integers).  The remainder is of variable size, which grows to
	    hold the LU factors and the frontal matrices created during
	    factorization.  A estimate of the upper bound is computed by
	    umfpack_*_*symbolic, and returned by umfpack_*_*symbolic in
	    Info [UMFPACK_VARIABLE_PEAK_ESTIMATE] (in Units).

	    If Control [UMFPACK_ALLOC_INIT] is >= 0, umfpack_*_numeric initially
	    allocates space for the variable-sized part equal to this estimate
	    times Control [UMFPACK_ALLOC_INIT].  Typically, for matrices for
	    which the "unsymmetric" strategy applies, umfpack_*_numeric needs
	    only about half the estimated memory space, so a setting of 0.5 or
	    0.6 often provides enough memory for umfpack_*_numeric to factorize
	    the matrix with no subsequent increases in the size of this block.

	    If the matrix is ordered via AMD, then this non-negative parameter
	    is ignored.  The initial allocation ratio computed automatically,
	    as 1.2 * (nz + Info [UMFPACK_SYMMETRIC_LUNZ]) /
	    (Info [UMFPACK_LNZ_ESTIMATE] + Info [UMFPACK_UNZ_ESTIMATE] -
	    min (n_row, n_col)).

	    If Control [UMFPACK_ALLOC_INIT] is negative, then umfpack_*_numeric
	    allocates a space with initial size (in Units) equal to
	    (-Control [UMFPACK_ALLOC_INIT]).

	    Regardless of the value of this parameter, a space equal to or
	    greater than the the bare minimum amount of memory needed to start
	    the factorization is always initially allocated.  The bare initial
	    memory required is returned by umfpack_*_*symbolic in
	    Info [UMFPACK_VARIABLE_INIT_ESTIMATE] (an exact value, not an
	    estimate).

	    If the variable-size part of the Numeric object is found to be too
	    small sometime after numerical factorization has started, the memory
	    is increased in size by a factor of 1.2.   If this fails, the
	    request is reduced by a factor of 0.95 until it succeeds, or until
	    it determines that no increase in size is possible.  Garbage
	    collection then occurs.

	    The strategy of attempting to "malloc" a working space, and
	    re-trying with a smaller space, may not work under MATLAB, since
	    mxMalloc aborts the mexFunction if it fails.  The built-in umfpack
	    routine (version 4.0) in MATLAB 6.5 uses utMalloc instead, which
	    avoids this problem.  As a mexFunction, utMalloc is used unless
	    -DNUTIL is defined at compile time.  The utMalloc routine, and
	    utFree and utRealloc, are not documented.  If the mexFunction
	    doesn't work, then compile it with -DNUTIL instead.

	    If you are using the umfpack mexFunction, decrease the magnitude of
	    Control [UMFPACK_ALLOC_INIT] if you run out of memory in MATLAB.

	    Default initial allocation size: 0.7.  Thus, with the default
	    control settings and the "unsymmetric" strategy, the upper-bound is
	    reached after two reallocations (0.7 * 1.2 * 1.2 = 1.008).

	    Changing this parameter has little effect on fill-in or operation
	    count.  It has a small impact on run-time (the extra time required
	    to do the garbage collection and memory reallocation).

	Control [UMFPACK_FRONT_ALLOC_INIT]:  This parameter is new to V4.1.

	    When UMFPACK starts the factorization of each "chain" of frontal
	    matrices, it allocates a working array to hold the frontal matrices
	    as they are factorized.  The symbolic factorization computes the
	    size of the largest possible frontal matrix that could occur during
	    the factorization of each chain.

	    If Control [UMFPACK_FRONT_ALLOC_INIT] is >= 0, the following
	    strategy is used.  If the AMD ordering was used, this non-negative
	    parameter is ignored.  A front of size (d+2)*(d+2) is allocated,
	    where d = Info [UMFPACK_SYMMETRIC_DMAX].  Otherwise, a front of
	    size Control [UMFPACK_FRONT_ALLOC_INIT] times the largest front
	    possible for this chain is allocated.

	    If Control [UMFPACK_FRONT_ALLOC_INIT] is negative, then a front of
	    size (-Control [UMFPACK_FRONT_ALLOC_INIT]) is allocated (where the
	    size is in terms of the number of numerical entries).  This is done
	    regardless of the ordering method or ordering strategy used.

	    Default: 0.5.

    double Info [UMFPACK_INFO] ;	Output argument.

	Contains statistics about the numeric factorization.  If a
	(double *) NULL pointer is passed, then no statistics are returned in
	Info (this is not an error condition).  The following statistics are
	computed in umfpack_*_numeric:

	Info [UMFPACK_STATUS]: status code.  This is also the return value,
	    whether or not Info is present.

	    UMFPACK_OK

		Numeric factorization was successful.  umfpack_*_numeric
		computed a valid numeric factorization.

	    UMFPACK_WARNING_singular_matrix

		Numeric factorization was successful, but the matrix is
		singular.  umfpack_*_numeric computed a valid numeric
		factorization, but you will get a divide by zero in
		umfpack_*_*solve.  For the other cases below, no Numeric object
		is created (*Numeric is (void *) NULL).

	    UMFPACK_ERROR_out_of_memory

		Insufficient memory to complete the numeric factorization.

	    UMFPACK_ERROR_argument_missing

		One or more required arguments are missing.

	    UMFPACK_ERROR_invalid_Symbolic_object

		Symbolic object provided as input is invalid.

	    UMFPACK_ERROR_different_pattern

		The pattern (Ap and/or Ai) has changed since the call to
		umfpack_*_*symbolic which produced the Symbolic object.

	Info [UMFPACK_NROW]:  the value of n_row stored in the Symbolic object.

	Info [UMFPACK_NCOL]:  the value of n_col stored in the Symbolic object.

	Info [UMFPACK_NZ]:  the number of entries in the input matrix.
	    This value is obtained from the Symbolic object.

	Info [UMFPACK_SIZE_OF_UNIT]:  the number of bytes in a Unit, for memory
	    usage statistics below.

	Info [UMFPACK_VARIABLE_INIT]: the initial size (in Units) of the
	    variable-sized part of the Numeric object.  If this differs from
	    Info [UMFPACK_VARIABLE_INIT_ESTIMATE], then the pattern (Ap and/or
	    Ai) has changed since the last call to umfpack_*_*symbolic, which is
	    an error condition.

	Info [UMFPACK_VARIABLE_PEAK]: the peak size (in Units) of the
	    variable-sized part of the Numeric object.  This size is the amount
	    of space actually used inside the block of memory, not the space
	    allocated via UMF_malloc.  You can reduce UMFPACK's memory
	    requirements by setting Control [UMFPACK_ALLOC_INIT] to the ratio
	    Info [UMFPACK_VARIABLE_PEAK] / Info[UMFPACK_VARIABLE_PEAK_ESTIMATE].
	    This will ensure that no memory reallocations occur (you may want to
	    add 0.001 to make sure that integer roundoff does not lead to a
	    memory size that is 1 Unit too small; otherwise, garbage collection
	    and reallocation will occur).

	Info [UMFPACK_VARIABLE_FINAL]: the final size (in Units) of the
	    variable-sized part of the Numeric object.  It holds just the
	    sparse LU factors.

	Info [UMFPACK_NUMERIC_SIZE]:  the actual final size (in Units) of the
	    entire Numeric object, including the final size of the variable
	    part of the object.  Info [UMFPACK_NUMERIC_SIZE_ESTIMATE],
	    an estimate, was computed by umfpack_*_*symbolic.  The estimate is
	    normally an upper bound on the actual final size, but this is not
	    guaranteed.

	Info [UMFPACK_PEAK_MEMORY]:  the actual peak memory usage (in Units) of
	    both umfpack_*_*symbolic and umfpack_*_numeric.  An estimate,
	    Info [UMFPACK_PEAK_MEMORY_ESTIMATE], was computed by
	    umfpack_*_*symbolic.  The estimate is normally an upper bound on the
	    actual peak usage, but this is not guaranteed.  With testing on
	    hundreds of matrix arising in real applications, I have never
	    observed a matrix where this estimate or the Numeric size estimate
	    was less than the actual result, but this is theoretically possible.
	    Please send me one if you find such a matrix.

	Info [UMFPACK_FLOPS]:  the actual count of the (useful) floating-point
	    operations performed.  An estimate, Info [UMFPACK_FLOPS_ESTIMATE],
	    was computed by umfpack_*_*symbolic.  The estimate is guaranteed to
	    be an upper bound on this flop count.  The flop count excludes
	    "useless" flops on zero values, flops performed during the pivot
	    search (for tentative updates and assembly of candidate columns),
	    and flops performed to add frontal matrices together.

	    For the real version, only (+ - * /) are counted.  For the complex
	    version, the following counts are used:

		operation	flops
	    	c = 1/b		6
		c = a*b		6
		c -= a*b	8

	Info [UMFPACK_LNZ]: the actual nonzero entries in final factor L,
	    including the diagonal.  This excludes any zero entries in L,
	    although some of these are stored in the Numeric object.  The
	    Info [UMFPACK_LU_ENTRIES] statistic does account for all
	    explicitly stored zeros, however.  Info [UMFPACK_LNZ_ESTIMATE],
	    an estimate, was computed by umfpack_*_*symbolic.  The estimate is
	    guaranteed to be an upper bound on Info [UMFPACK_LNZ].

	Info [UMFPACK_UNZ]: the actual nonzero entries in final factor U,
	    including the diagonal.  This excludes any zero entries in U,
	    although some of these are stored in the Numeric object.  The
	    Info [UMFPACK_LU_ENTRIES] statistic does account for all
	    explicitly stored zeros, however.  Info [UMFPACK_UNZ_ESTIMATE],
	    an estimate, was computed by umfpack_*_*symbolic.  The estimate is
	    guaranteed to be an upper bound on Info [UMFPACK_UNZ].

	Info [UMFPACK_NUMERIC_DEFRAG]:  The number of garbage collections
	    performed during umfpack_*_numeric, to compact the contents of the
	    variable-sized workspace used by umfpack_*_numeric.  No estimate was
	    computed by umfpack_*_*symbolic.  In the current version of UMFPACK,
	    garbage collection is performed and then the memory is reallocated,
	    so this statistic is the same as Info [UMFPACK_NUMERIC_REALLOC],
	    below.  It may differ in future releases.

	Info [UMFPACK_NUMERIC_REALLOC]:  The number of times that the Numeric
	    object was increased in size from its initial size.  A rough upper
	    bound on the peak size of the Numeric object was computed by
	    umfpack_*_*symbolic, so reallocations should be rare.  However, if
	    umfpack_*_numeric is unable to allocate that much storage, it
	    reduces its request until either the allocation succeeds, or until
	    it gets too small to do anything with.  If the memory that it
	    finally got was small, but usable, then the reallocation count
	    could be high.  No estimate of this count was computed by
	    umfpack_*_*symbolic.

	Info [UMFPACK_NUMERIC_COSTLY_REALLOC]:  The number of times that the
	    system realloc library routine (or mxRealloc for the mexFunction)
	    had to move the workspace.  Realloc can sometimes increase the size
	    of a block of memory without moving it, which is much faster.  This
	    statistic will always be <= Info [UMFPACK_NUMERIC_REALLOC].  If your
	    memory space is fragmented, then the number of "costly" realloc's
	    will be equal to Info [UMFPACK_NUMERIC_REALLOC].

	Info [UMFPACK_COMPRESSED_PATTERN]:  The number of integers used to
	    represent the pattern of L and U.

	Info [UMFPACK_LU_ENTRIES]:  The total number of numerical values that
	    are stored for the LU factors.  Some of the values may be explicitly
	    zero in order to save space (allowing for a smaller compressed
	    pattern).

	Info [UMFPACK_NUMERIC_TIME]:  The CPU time taken, in seconds.

	Info [UMFPACK_RCOND]:  A rough estimate of the condition number, equal
	    to min (abs (diag (U))) / max (abs (diag (U))), or zero if the
	    diagonal of U is all zero.

	Info [UMFPACK_UDIAG_NZ]:  The number of numerically nonzero values on
	    the diagonal of U.

	Info [UMFPACK_UMIN]:  the smallest absolute value on the diagonal of U.

	Info [UMFPACK_UMAX]:  the smallest absolute value on the diagonal of U.

	Info [UMFPACK_MAX_FRONT_SIZE]: the size of the
	    largest frontal matrix (number of entries).

	------------------------------------------------------------------------
	The following statistics were added to Version 4.1:
	------------------------------------------------------------------------

	Info [UMFPACK_NUMERIC_WALLTIME]:  The wallclock time taken, in seconds.

	Info [UMFPACK_MAX_FRONT_NROWS]: the max number of
	    rows in any frontal matrix.

	Info [UMFPACK_MAX_FRONT_NCOLS]: the max number of
	    columns in any frontal matrix.

	Info [UMFPACK_WAS_SCALED]:  the scaling used, either UMFPACK_SCALE_NONE,
	    UMFPACK_SCALE_SUM, or UMFPACK_SCALE_MAX.

	Info [UMFPACK_RSMIN]: if scaling is performed, the smallest scale factor
	    for any row (either the smallest sum of absolute entries, or the
	    smallest maximum of absolute entries).

	Info [UMFPACK_RSMAX]: if scaling is performed, the largest scale factor
	    for any row (either the largest sum of absolute entries, or the
	    largest maximum of absolute entries).

	Info [UMFPACK_ALLOC_INIT_USED]:  the initial allocation parameter used.

	Info [UMFPACK_FORCED_UPDATES]:  the number of BLAS-3 updates to the
	    frontal matrices that were required because the frontal matrix
	    grew larger than its current working array.

	Info [UMFPACK_NOFF_DIAG]: number of off-diagonal pivots selected, if the
	    symmetric or 2-by-2 strategies are used.

	Only the above listed Info [...] entries are accessed.  The remaining
	entries of Info are not accessed or modified by umfpack_*_numeric.
	Future versions might modify different parts of Info.
*/
