/* ========================================================================== */
/* === UMF_colamd =========================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
UMF_colamd:  an approximate minimum degree column ordering algorithm,
    used as a preordering for UMFPACK.

NOTE: if this routine is used outside of UMFPACK, for a sparse Cholesky
factorization of (AQ)'*(AQ) or a QR factorization of A, then one line should
be removed (the "&& pivot_row_thickness > 0" expression).  See the comment
regarding the Cholesky factorization, below.

Purpose:

    Colamd computes a permutation Q such that the Cholesky factorization of
    (AQ)'(AQ) has less fill-in and requires fewer floating point operations
    than A'A.  This also provides a good ordering for sparse partial
    pivoting methods, P(AQ) = LU, where Q is computed prior to numerical
    factorization, and P is computed during numerical factorization via
    conventional partial pivoting with row interchanges.  Colamd is the
    column ordering method used in SuperLU, part of the ScaLAPACK library.
    It is also available as built-in function in MATLAB Version 6,
    available from MathWorks, Inc. (http://www.mathworks.com).  This
    routine can be used in place of colmmd in MATLAB.  By default, the \
    and / operators in MATLAB perform a column ordering (using colmmd
    or colamd) prior to LU factorization using sparse partial pivoting,
    in the built-in MATLAB lu(A) routine.

    This code is derived from Colamd Version 2.0.

Authors:

    The authors of the COLAMD code itself are Stefan I. Larimore and Timothy A.
    Davis, University of Florida.  The algorithm was developed in collaboration
    with John Gilbert, Xerox PARC, and Esmond Ng, Oak Ridge National Laboratory.
    The AMD metric on which this is based is by Patrick Amestoy, T. Davis,
    and Iain Duff.

Date:

    UMFPACK Version: see above.
    COLAMD Version 2.0 was released on January 31, 2000.

Acknowledgements:

    This work was supported by the National Science Foundation, under
    grants DMS-9504974 and DMS-9803599.

UMFPACK:  Copyright (c) 2003 by Timothy A. Davis.  All Rights Reserved.

See the UMFPACK README file for the License for your use of this code.

Availability:

    Both UMFPACK and the original unmodified colamd/symamd library are
    available at http://www.cise.ufl.edu/research/sparse.

Changes for inclusion in UMFPACK:

    * symamd, symamd_report, and colamd_report removed

    * additional terms added to RowInfo, ColInfo, and stats

    * Frontal matrix information computed for UMFPACK

    * routines renamed

    * column elimination tree post-ordering incorporated.  In the original
	version 2.0, this was performed in colamd.m.

For more information, see:

	Amestoy, P. R. and Davis, T. A. and Duff, I. S.,
	An approximate minimum degree ordering algorithm,
	SIAM J. Matrix Analysis and Applic, vol 17, no 4., pp 886-905, 1996.

	Davis, T. A. and Gilbert, J. R. and Larimore, S. I. and Ng, E. G.,
	A column approximate minimum degree ordering algorithm,
	Univ. of Florida, CISE Dept., TR-00-005, Gainesville, FL
	Oct. 2000.  Submitted to ACM Trans. Math. Softw.

*/

/* ========================================================================== */
/* === Description of user-callable routines ================================ */
/* ========================================================================== */

/*
    ----------------------------------------------------------------------------
    colamd_recommended: removed for UMFPACK
    ----------------------------------------------------------------------------

    ----------------------------------------------------------------------------
    colamd_set_defaults:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    colamd_set_defaults (double knobs [COLAMD_KNOBS]) ;

	Purpose:

	    Sets the default parameters.  The use of this routine is optional.

	Arguments:

	    double knobs [COLAMD_KNOBS] ;	Output only.

		Let c = knobs [COLAMD_DENSE_COL], r = knobs [COLAMD_DENSE_ROW].
		Colamd: rows with more than max (16, r*16*sqrt(n_col))
		entries are removed prior to ordering.  Columns with more than
		max (16, c*16*sqrt(n_row)) entries are removed prior to
		ordering, and placed last in the output column ordering.

		Symamd: removed for UMFPACK.

		COLAMD_DENSE_ROW and COLAMD_DENSE_COL are defined as 0 and 1,
		respectively, in colamd.h.  Default values of these two knobs
		are both 0.5.  Currently, only knobs [0] and knobs [1] are
		used, but future versions may use more knobs.  If so, they will
		be properly set to their defaults by the future version of
		colamd_set_defaults, so that the code that calls colamd will
		not need to change, assuming that you either use
		colamd_set_defaults, or pass a (double *) NULL pointer as the
		knobs array to colamd or symamd.

		knobs [COLAMD_AGGRESSIVE]: if nonzero, then perform aggressive
		absorption.  Otherwise, do not.  This version does aggressive
		absorption by default.  COLAMD v2.1 (in MATLAB) always
		does aggressive absorption (it doesn't have an option to turn
		it off).

    ----------------------------------------------------------------------------
    colamd:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    Int UMF_colamd (Int n_row, Int n_col, Int Alen, Int *A, Int *p,
		double knobs [COLAMD_KNOBS], Int stats [COLAMD_STATS]) ;

	Purpose:

	    Computes a column ordering (Q) of A such that P(AQ)=LU or
	    (AQ)'AQ=LL' have less fill-in and require fewer floating point
	    operations than factorizing the unpermuted matrix A or A'A,
	    respectively.

	Returns:

	    TRUE (1) if successful, FALSE (0) otherwise.

	Arguments:

	    Int n_row ;		Input argument.

		Number of rows in the matrix A.
		Restriction:  n_row >= 0.
		Colamd returns FALSE if n_row is negative.

	    Int n_col ;		Input argument.

		Number of columns in the matrix A.
		Restriction:  n_col >= 0.
		Colamd returns FALSE if n_col is negative.

	    Int Alen ;		Input argument.

		Restriction (see note):
		Alen >= 2*nnz + 8*(n_col+1) + 6*(n_row+1) + n_col
		Colamd returns FALSE if these conditions are not met.

		Note:  this restriction makes an modest assumption regarding
		the size of the two typedef's structures in colamd.h.
		We do, however, guarantee that

			Alen >= UMF_COLAMD_RECOMMENDED (nnz, n_row, n_col)

		will be sufficient.

	    Int A [Alen] ;	Input and output argument.

		A is an integer array of size Alen.  Alen must be at least as
		large as the bare minimum value given above, but this is very
		low, and can result in excessive run time.  For best
		performance, we recommend that Alen be greater than or equal to
		UMF_COLAMD_RECOMMENDED (nnz, n_row, n_col), which adds
		nnz/5 to the bare minimum value given above.

		On input, the row indices of the entries in column c of the
		matrix are held in A [(p [c]) ... (p [c+1]-1)].  The row indices
		in a given column c need not be in ascending order, and
		duplicate row indices may be be present.  However, colamd will
		work a little faster if both of these conditions are met
		(Colamd puts the matrix into this format, if it finds that the
		the conditions are not met).

		The matrix is 0-based.  That is, rows are in the range 0 to
		n_row-1, and columns are in the range 0 to n_col-1.  Colamd
		returns FALSE if any row index is out of range.

		A holds the inverse permutation on output.

	    Int p [n_col+1] ;	Both input and output argument.

		p is an integer array of size n_col+1.  On input, it holds the
		"pointers" for the column form of the matrix A.  Column c of
		the matrix A is held in A [(p [c]) ... (p [c+1]-1)].  The first
		entry, p [0], must be zero, and p [c] <= p [c+1] must hold
		for all c in the range 0 to n_col-1.  The value p [n_col] is
		thus the total number of entries in the pattern of the matrix A.
		Colamd returns FALSE if these conditions are not met.

		On output, if colamd returns TRUE, the array p holds the column
		permutation (Q, for P(AQ)=LU or (AQ)'(AQ)=LL'), where p [0] is
		the first column index in the new ordering, and p [n_col-1] is
		the last.  That is, p [k] = j means that column j of A is the
		kth pivot column, in AQ, where k is in the range 0 to n_col-1
		(p [0] = j means that column j of A is the first column in AQ).

		If colamd returns FALSE, then no permutation is returned, and
		p is undefined on output.

	    double knobs [COLAMD_KNOBS] ;	Input argument.

		See colamd_set_defaults for a description.
		The behavior is undefined if knobs contains NaN's.
		(UMFPACK does not call umf_colamd with NaN-valued knobs).

	    Int stats [COLAMD_STATS] ;		Output argument.

		Statistics on the ordering, and error status.
		See colamd.h for related definitions.
		Colamd returns FALSE if stats is not present.

		stats [0]:  number of dense or empty rows ignored.

		stats [1]:  number of dense or empty columns ignored (and
				ordered last in the output permutation p)
				Note that a row can become "empty" if it
				contains only "dense" and/or "empty" columns,
				and similarly a column can become "empty" if it
				only contains "dense" and/or "empty" rows.

		stats [2]:  number of garbage collections performed.
				This can be excessively high if Alen is close
				to the minimum required value.

		stats [3]:  status code.  < 0 is an error code.
			    > 1 is a warning or notice.

			0	OK.  Each column of the input matrix contained
				row indices in increasing order, with no
				duplicates.

			-11	Columns of input matrix jumbled
				(unsorted columns or duplicate entries).

					stats [4]: the bad column index
					stats [5]: the bad row index

			-1	A is a null pointer

			-2	p is a null pointer

			-3	n_row is negative

					stats [4]: n_row

			-4	n_col is negative

					stats [4]: n_col

			-5	number of nonzeros in matrix is negative

					stats [4]: number of nonzeros, p [n_col]

			-6	p [0] is nonzero

					stats [4]: p [0]

			-7	A is too small

					stats [4]: required size
					stats [5]: actual size (Alen)

			-8	a column has a zero or negative number of
				entries (changed for UMFPACK)

					stats [4]: column with <= 0 entries
					stats [5]: number of entries in col

			-9	a row index is out of bounds

					stats [4]: column with bad row index
					stats [5]: bad row index
					stats [6]: n_row, # of rows of matrx

			-10	unused

			-999	(unused; see symamd.c)

		Future versions may return more statistics in the stats array.

	Example:

	    See http://www.cise.ufl.edu/~davis/colamd/example.c
	    for a complete example.

	    To order the columns of a 5-by-4 matrix with 11 nonzero entries in
	    the following nonzero pattern

		x 0 x 0
		x 0 x x
		0 x x 0
		0 0 x x
		x x 0 0

	    with default knobs and no output statistics, do the following:

		#include "colamd.h"
		#define ALEN UMF_COLAMD_RECOMMENDED (11, 5, 4)
		Int A [ALEN] = {1, 2, 5, 3, 5, 1, 2, 3, 4, 2, 4} ;
		Int p [ ] = {0, 3, 5, 9, 11} ;
		Int stats [COLAMD_STATS] ;
		UMF_colamd (5, 4, ALEN, A, p, (double *) NULL, stats) ;

	    The permutation is returned in the array p, and A is destroyed.


    ----------------------------------------------------------------------------
    symamd:  does not appear in this version for UMFPACK
    ----------------------------------------------------------------------------

    ----------------------------------------------------------------------------
    colamd_report: does not appear in this version for UMFPACK
    ----------------------------------------------------------------------------

    ----------------------------------------------------------------------------
    symamd_report: does not appear in this version for UMFPACK
    ----------------------------------------------------------------------------

*/

/* ========================================================================== */
/* === Scaffolding code definitions  ======================================== */
/* ========================================================================== */

/* UMFPACK debugging control moved to amd_internal.h */

/*
   Our "scaffolding code" philosophy:  In our opinion, well-written library
   code should keep its "debugging" code, and just normally have it turned off
   by the compiler so as not to interfere with performance.  This serves
   several purposes:

   (1) assertions act as comments to the reader, telling you what the code
	expects at that point.  All assertions will always be true (unless
	there really is a bug, of course).

   (2) leaving in the scaffolding code assists anyone who would like to modify
	the code, or understand the algorithm (by reading the debugging output,
	one can get a glimpse into what the code is doing).

   (3) (gasp!) for actually finding bugs.  This code has been heavily tested
	and "should" be fully functional and bug-free ... but you never know...

    To enable debugging, comment out the "#define NDEBUG" above.  For a MATLAB
    mexFunction, you will also need to modify mexopts.sh to remove the -DNDEBUG
    definition.  The code will become outrageously slow when debugging is
    enabled.  To control the level of debugging output, set an environment
    variable D to 0 (little), 1 (some), 2, 3, or 4 (lots).  When debugging,
    you should see the following message on the standard output:

	colamd: debug version, D = 1 (THIS WILL BE SLOW!)

    or a similar message for symamd.  If you don't, then debugging has not
    been enabled.

*/

/* ========================================================================== */
/* === Include files ======================================================== */
/* ========================================================================== */

/* ------------------ */
/* modified for UMFPACK: */
#include "umf_internal.h"
#include "umf_colamd.h"
#include "umf_apply_order.h"
#include "umf_fsize.h"
/* ------------------ */

/* ========================================================================== */
/* === Definitions ========================================================== */
/* ========================================================================== */

/* ------------------ */
/* UMFPACK: duplicate definitions moved to umf_internal.h */
/* ------------------ */

/* Row and column status */
#define ALIVE	(0)
#define DEAD	(-1)

/* Column status */
#define DEAD_PRINCIPAL		(-1)
#define DEAD_NON_PRINCIPAL	(-2)

/* Macros for row and column status update and checking. */
#define ROW_IS_DEAD(r)			ROW_IS_MARKED_DEAD (Row[r].shared2.mark)
#define ROW_IS_MARKED_DEAD(row_mark)	(row_mark < ALIVE)
#define ROW_IS_ALIVE(r)			(Row [r].shared2.mark >= ALIVE)
#define COL_IS_DEAD(c)			(Col [c].start < ALIVE)
#define COL_IS_ALIVE(c)			(Col [c].start >= ALIVE)
#define COL_IS_DEAD_PRINCIPAL(c)	(Col [c].start == DEAD_PRINCIPAL)
#define KILL_ROW(r)			{ Row [r].shared2.mark = DEAD ; }
#define KILL_PRINCIPAL_COL(c)		{ Col [c].start = DEAD_PRINCIPAL ; }
#define KILL_NON_PRINCIPAL_COL(c)	{ Col [c].start = DEAD_NON_PRINCIPAL ; }

/* ------------------ */
/* UMFPACK: Colamd reporting mechanism moved to umf_internal.h */
/* ------------------ */

/* ========================================================================== */
/* === Prototypes of PRIVATE routines ======================================= */
/* ========================================================================== */

PRIVATE Int init_rows_cols
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int p []
    /* Int stats [COLAMD_STATS] */
) ;

PRIVATE void init_scoring
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int head [],
    double knobs [COLAMD_KNOBS],
    Int *p_n_row2,
    Int *p_n_col2,
    Int *p_max_deg
    /* ------------------ */
    /* added for UMFPACK */
    , Int *p_ndense_row		/* number of dense rows */
    , Int *p_nempty_row		/* number of original empty rows */
    , Int *p_nnewlyempty_row	/* number of newly empty rows */
    , Int *p_ndense_col		/* number of dense cols (excl "empty" cols) */
    , Int *p_nempty_col		/* number of original empty cols */
    , Int *p_nnewlyempty_col	/* number of newly empty cols */
) ;

PRIVATE Int find_ordering
(
    Int n_row,
    Int n_col,
    Int Alen,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int head [],
    Int n_col2,
    Int max_deg,
    Int pfree
    /* ------------------ */
    /* added for UMFPACK: */
    , Int Front_npivcol [ ]
    , Int Front_nrows [ ]
    , Int Front_ncols [ ]
    , Int Front_parent [ ]
    , Int Front_cols [ ]
    , Int *p_nfr
    , Int aggressive
    , Int InFront [ ]
    /* ------------------ */
) ;

/* ------------------ */
/* order_children deleted for UMFPACK: */
/* ------------------ */

PRIVATE void detect_super_cols
(

#ifndef NDEBUG
    Int n_col,
    Colamd_Row Row [],
#endif /* NDEBUG */

    Colamd_Col Col [],
    Int A [],
    Int head [],
    Int row_start,
    Int row_length
) ;

PRIVATE Int garbage_collection
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int *pfree
) ;

PRIVATE Int clear_mark
(
    Int n_row,
    Colamd_Row Row []
) ;

/* ------------------ */
/* print_report deleted for UMFPACK */
/* ------------------ */

/* ========================================================================== */
/* === Debugging prototypes and definitions ================================= */
/* ========================================================================== */

#ifndef NDEBUG

/* ------------------ */
/* debugging macros moved for UMFPACK */
/* ------------------ */

PRIVATE void debug_deg_lists
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int head [],
    Int min_score,
    Int should,
    Int max_deg
) ;

PRIVATE void debug_mark
(
    Int n_row,
    Colamd_Row Row [],
    Int tag_mark,
    Int max_mark
) ;

PRIVATE void debug_matrix
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A []
) ;

PRIVATE void debug_structures
(
    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int n_col2
) ;

/* ------------------ */
/* dump_super added for UMFPACK: */
PRIVATE void dump_super
(
    Int super_c,
    Colamd_Col Col [],
    Int n_col
) ;
/* ------------------ */

#endif /* NDEBUG */

/* ========================================================================== */



/* ========================================================================== */
/* === USER-CALLABLE ROUTINES: ============================================== */
/* ========================================================================== */


/* ========================================================================== */
/* === colamd_set_defaults ================================================== */
/* ========================================================================== */

/*
    The colamd_set_defaults routine sets the default values of the user-
    controllable parameters for colamd:

	knobs [0]	rows with knobs[0]*n_col entries or more are removed
			prior to ordering in colamd.  Rows and columns with
			knobs[0]*n_col entries or more are removed prior to
			ordering in symamd and placed last in the output
			ordering.

	knobs [1]	columns with knobs[1]*n_row entries or more are removed
			prior to ordering in colamd, and placed last in the
			column permutation.  Symamd ignores this knob.

	knobs [2]	if nonzero, then perform aggressive absorption.

	knobs [3..19]	unused, but future versions might use this
*/

GLOBAL void UMF_colamd_set_defaults
(
    /* === Parameters ======================================================= */

    double knobs [COLAMD_KNOBS]		/* knob array */
)
{
    /* === Local variables ================================================== */

    Int i ;

#if 0
    if (!knobs)
    {
	return ;			/* UMFPACK always passes knobs array */
    }
#endif
    for (i = 0 ; i < COLAMD_KNOBS ; i++)
    {
	knobs [i] = 0 ;
    }
    knobs [COLAMD_DENSE_ROW] = 0.2 ;	/* default changed for UMFPACK */
    knobs [COLAMD_DENSE_COL] = 0.2 ;	/* default changed for UMFPACK */
    knobs [COLAMD_AGGRESSIVE] = TRUE ;	/* default is to do aggressive
					 * absorption */
}


/* ========================================================================== */
/* === symamd removed for UMFPACK =========================================== */
/* ========================================================================== */



/* ========================================================================== */
/* === colamd =============================================================== */
/* ========================================================================== */

/*
    The colamd routine computes a column ordering Q of a sparse matrix
    A such that the LU factorization P(AQ) = LU remains sparse, where P is
    selected via partial pivoting.   The routine can also be viewed as
    providing a permutation Q such that the Cholesky factorization
    (AQ)'(AQ) = LL' remains sparse.
*/

/* For UMFPACK: colamd always returns TRUE */

GLOBAL Int UMF_colamd		/* returns TRUE if successful, FALSE otherwise*/
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows in A */
    Int n_col,			/* number of columns in A */
    Int Alen,			/* length of A */
    Int A [],			/* row indices of A */
    Int p [],			/* pointers to columns in A */
    double knobs [COLAMD_KNOBS],/* parameters (uses defaults if NULL) */
    Int stats [COLAMD_STATS]	/* output statistics and error codes */

    /* ------------------ */
    /* added for UMFPACK: each Front_ array is of size n_col+1 */
    , Int Front_npivcol [ ]	/* # pivot cols in each front */
    , Int Front_nrows [ ]	/* # of rows in each front (incl. pivot rows) */
    , Int Front_ncols [ ]	/* # of cols in each front (incl. pivot cols) */
    , Int Front_parent [ ]	/* parent of each front */
    , Int Front_cols [ ]	/* link list of pivot columns for each front */
    , Int *p_nfr		/* total number of frontal matrices */
    , Int InFront [ ]		/* InFront [row] = f if the original row was
				 * absorbed into front f.  EMPTY if the row was
				 * empty, dense, or not absorbed.  This array
				 * has size n_row+1 */
    /* ------------------ */
)
{
    /* === Local variables ================================================== */

    Int row ;			/* row index */
    Int i ;			/* loop index */
    Int nnz ;			/* nonzeros in A */
    Int Row_size ;		/* size of Row [], in integers */
    Int Col_size ;		/* size of Col [], in integers */
#if 0
    Int need ;			/* minimum required length of A */
#endif
    Colamd_Row *Row ;		/* pointer into A of Row [0..n_row] array */
    Colamd_Col *Col ;		/* pointer into A of Col [0..n_col] array */
    Int n_col2 ;		/* number of non-dense, non-empty columns */
    Int n_row2 ;		/* number of non-dense, non-empty rows */
    Int ngarbage ;		/* number of garbage collections performed */
    Int max_deg ;		/* maximum row degree */
    Int aggressive ;		/* TRUE if doing aggressive absorption */
#if 0
    double default_knobs [COLAMD_KNOBS] ;	/* default knobs array */
#endif

    /* ------------------ */
    /* debugging initializations moved for UMFPACK */
    /* ------------------ */

    /* ------------------ */
    /* added for UMFPACK: */
    Int ndense_row, nempty_row, parent, ndense_col,
	nempty_col, k, col, nfr, *Front_child, *Front_sibling, *Front_stack,
	*Front_order, *Front_size ;
    Int nnewlyempty_col, nnewlyempty_row ;
    /* ------------------ */

    /* === Check the input arguments ======================================== */

#if 0
    if (!stats)
    {
	DEBUG0 (("colamd: stats not present\n")) ;
	return (FALSE) ;	/* UMFPACK:  always passes stats [ ] */
    }
#endif

    ASSERT (stats != (Int *) NULL) ;

    for (i = 0 ; i < COLAMD_STATS ; i++)
    {
	stats [i] = 0 ;
    }
    stats [COLAMD_STATUS] = COLAMD_OK ;
    stats [COLAMD_INFO1] = -1 ;
    stats [COLAMD_INFO2] = -1 ;

#if 0
    if (!A)		/* A is not present */
    {
	/* UMFPACK:  always passes A [ ] */
	DEBUG0 (("colamd: A not present\n")) ;
	stats [COLAMD_STATUS] = COLAMD_ERROR_A_not_present ;
	return (FALSE) ;
    }

    if (!p)		/* p is not present */
    {
	/* UMFPACK:  always passes p [ ] */
	DEBUG0 (("colamd: p not present\n")) ;
	stats [COLAMD_STATUS] = COLAMD_ERROR_p_not_present ;
	return (FALSE) ;
    }

    if (n_row < 0)	/* n_row must be >= 0 */
    {
	/* UMFPACK:  does not call UMF_colamd if n <= 0 */
	DEBUG0 (("colamd: nrow negative "ID"\n", n_row)) ;
	stats [COLAMD_STATUS] = COLAMD_ERROR_nrow_negative ;
	stats [COLAMD_INFO1] = n_row ;
	return (FALSE) ;
    }

    if (n_col < 0)	/* n_col must be >= 0 */
    {
	/* UMFPACK:  does not call UMF_colamd if n <= 0 */
	DEBUG0 (("colamd: ncol negative "ID"\n", n_col)) ;
	stats [COLAMD_STATUS] = COLAMD_ERROR_ncol_negative ;
	stats [COLAMD_INFO1] = n_col ;
	return (FALSE) ;
    }
#endif

    ASSERT (A != (Int *) NULL) ;
    ASSERT (p != (Int *) NULL) ;
    ASSERT (n_row >= 0) ;
    ASSERT (n_col >= 0) ;

    nnz = p [n_col] ;

#if 0
    if (nnz < 0)	/* nnz must be >= 0 */
    {
	/* UMFPACK:  does not call UMF_colamd if nnz < 0 */
	DEBUG0 (("colamd: number of entries negative "ID"\n", nnz)) ;
	stats [COLAMD_STATUS] = COLAMD_ERROR_nnz_negative ;
	stats [COLAMD_INFO1] = nnz ;
	return (FALSE) ;
    }

    if (p [0] != 0)	/* p [0] must be exactly zero */
    {
	DEBUG0 (("colamd: p[0] not zero "ID"\n", p [0])) ;
	stats [COLAMD_STATUS] = COLAMD_ERROR_p0_nonzero	;
	stats [COLAMD_INFO1] = p [0] ;
	return (FALSE) ;
    }
#endif

    ASSERT (nnz >= 0) ;
    ASSERT (p [0] == 0) ;

    /* === If no knobs, set default knobs =================================== */

#if 0
    if (!knobs)
    {
	/* UMFPACK:  always passes the knobs */
	UMF_colamd_set_defaults (default_knobs) ;
	knobs = default_knobs ;
    }
#endif

    ASSERT (knobs != (double *) NULL) ;

    /* --------------------- */
    /* added for UMFPACK v4.1: */
    aggressive = (knobs [COLAMD_AGGRESSIVE] != 0) ;
    /* --------------------- */

    /* === Allocate the Row and Col arrays from array A ===================== */

    Col_size = UMF_COLAMD_C (n_col) ;
    Row_size = UMF_COLAMD_R (n_row) ;

#if 0
    need = MAX (2*nnz, 4*n_col) + n_col + Col_size + Row_size ;
    if (need > Alen)
    {
	/* UMFPACK: always passes enough space */
	/* not enough space in array A to perform the ordering */
	DEBUG0 (("colamd: Need Alen >= "ID", given only Alen = "ID"\n",
	    need, Alen)) ;
	stats [COLAMD_STATUS] = COLAMD_ERROR_A_too_small ;
	stats [COLAMD_INFO1] = need ;
	stats [COLAMD_INFO2] = Alen ;
	return (FALSE) ;
    }
#endif

    Alen -= Col_size + Row_size ;
    Col = (Colamd_Col *) &A [Alen] ;
    Row = (Colamd_Row *) &A [Alen + Col_size] ;

    /* Size of A is now Alen >= MAX (2*nnz, 4*n_col) + n_col.  The ordering
     * requires Alen >= 2*nnz + n_col, and the postorder requires
     * Alen >= 5*n_col. */

    /* === Construct the row and column data structures ===================== */

    i = init_rows_cols (n_row, n_col, Row, Col, A, p) ;

#if 0
    if (!i)
    {
	/* input matrix is invalid */
	DEBUG0 (("colamd: Matrix invalid\n")) ;
	return (FALSE) ;
    }
#endif

    ASSERT (i) ;

    /* === UMFPACK: Initialize front info =================================== */

    for (col = 0 ; col < n_col ; col++)
    {
	Front_npivcol [col] = 0 ;
	Front_nrows [col] = 0 ;
	Front_ncols [col] = 0 ;
	Front_parent [col] = EMPTY ;
	Front_cols [col] = EMPTY ;
    }

    /* === Initialize scores, kill dense rows/columns ======================= */

    init_scoring (n_row, n_col, Row, Col, A, p, knobs,
	&n_row2, &n_col2, &max_deg
	/* ------------------ */
	/* added for UMFPACK: */
	, &ndense_row, &nempty_row, &nnewlyempty_row
	, &ndense_col, &nempty_col, &nnewlyempty_col
	/* ------------------ */
	) ;
    ASSERT (n_row2 == n_row - nempty_row - nnewlyempty_row - ndense_row) ;
    ASSERT (n_col2 == n_col - nempty_col - nnewlyempty_col - ndense_col) ;

    /* === Order the supercolumns =========================================== */

    ngarbage = find_ordering (n_row, n_col, Alen, Row, Col, A, p,
	n_col2, max_deg, 2*nnz
	/* ------------------ */
	/* added for UMFPACK: */
	, Front_npivcol, Front_nrows, Front_ncols, Front_parent, Front_cols
	, &nfr, aggressive, InFront
	/* ------------------ */
	) ;

    /* ------------------ */
    /* changed for UMFPACK: */

    /* A is no longer needed, so use A [0..5*nfr-1] as workspace [ [ */
    /* This step requires Alen >= 5*n_col */
    Front_child   = A ;
    Front_sibling = Front_child + nfr ;
    Front_stack   = Front_sibling + nfr ;
    Front_order   = Front_stack + nfr ;
    Front_size    = Front_order + nfr ;

    UMF_fsize (nfr, Front_size, Front_nrows, Front_ncols,
	    Front_parent, Front_npivcol) ;

    AMD_postorder (nfr, Front_parent, Front_npivcol, Front_size,
	Front_order, Front_child, Front_sibling, Front_stack) ;

    /* Front_size, Front_stack, Front_child, Front_sibling no longer needed ] */

    /* use A [0..nfr-1] as workspace */
    UMF_apply_order (Front_npivcol, Front_order, A, nfr, nfr) ;
    UMF_apply_order (Front_nrows,   Front_order, A, nfr, nfr) ;
    UMF_apply_order (Front_ncols,   Front_order, A, nfr, nfr) ;
    UMF_apply_order (Front_parent,  Front_order, A, nfr, nfr) ;
    UMF_apply_order (Front_cols,    Front_order, A, nfr, nfr) ;

    /* fix the parent to refer to the new numbering */
    for (i = 0 ; i < nfr ; i++)
    {
	parent = Front_parent [i] ;
	if (parent != EMPTY)
	{
	    Front_parent [i] = Front_order [parent] ;
	}
    }

    /* fix InFront to refer to the new numbering */
    for (row = 0 ; row < n_row ; row++)
    {
	i = InFront [row] ;
	ASSERT (i >= EMPTY && i < nfr) ;
	if (i != EMPTY)
	{
	    InFront [row] = Front_order [i] ;
	}
    }

    /* Front_order longer needed ] */

    /* === Order the columns in the fronts ================================== */

    /* use A [0..n_col-1] as inverse permutation */
    for (i = 0 ; i < n_col ; i++)
    {
	A [i] = EMPTY ;
    }
    k = 0 ;
    for (i = 0 ; i < nfr ; i++)
    {
	ASSERT (Front_npivcol [i] > 0) ;
	for (col = Front_cols [i] ; col != EMPTY ; col = Col [col].nextcol)
	{
	    ASSERT (col >= 0 && col < n_col) ;
	    DEBUG1 (("Colamd output ordering: k "ID" col "ID"\n", k, col)) ;
	    p [k] = col ;
	    ASSERT (A [col] == EMPTY) ;
	    A [col] = k ;
	    k++ ;
	}
    }

    /* === Order the "dense" and null columns =============================== */

    ASSERT (k == n_col2) ;
    if (n_col2 < n_col)
    {
	for (col = 0 ; col < n_col ; col++)
	{
	    if (A [col] == EMPTY)
	    {
		k = Col [col].shared2.order ;
		ASSERT (k >= n_col2 && k < n_col) ;
		DEBUG1 (("Colamd output ordering: k "ID" col "ID
		    " (dense or null col)\n", k, col)) ;
		p [k] = col ;
		A [col] = k ;
	    }
	}
    }

    /* ------------------ */

    /* === Return statistics in stats ======================================= */

    /* ------------------ */
    /* modified for UMFPACK */
    stats [COLAMD_DENSE_ROW] = ndense_row ;
    stats [COLAMD_EMPTY_ROW] = nempty_row ;
    stats [COLAMD_NEWLY_EMPTY_ROW] = nnewlyempty_row ;
    stats [COLAMD_DENSE_COL] = ndense_col ;
    stats [COLAMD_EMPTY_COL] = nempty_col ;
    stats [COLAMD_NEWLY_EMPTY_COL] = nnewlyempty_col ;
    ASSERT (ndense_col + nempty_col + nnewlyempty_col == n_col - n_col2) ;
    /* ------------------ */
    stats [COLAMD_DEFRAG_COUNT] = ngarbage ;
    *p_nfr = nfr ;
    DEBUG1 (("colamd: done.\n")) ;
    return (TRUE) ;
}




/* ========================================================================== */
/* === colamd_report removed for UMFPACK ==================================== */
/* ========================================================================== */

/* ========================================================================== */
/* === symamd_report removed for UMFPACK ==================================== */
/* ========================================================================== */



/* ========================================================================== */
/* === NON-USER-CALLABLE ROUTINES: ========================================== */
/* ========================================================================== */

/* There are no user-callable routines beyond this point in the file */


/* ========================================================================== */
/* === init_rows_cols ======================================================= */
/* ========================================================================== */

/*
    Takes the column form of the matrix in A and creates the row form of the
    matrix.  Also, row and column attributes are stored in the Col and Row
    structs.  If the columns are un-sorted or contain duplicate row indices,
    this routine will also sort and remove duplicate row indices from the
    column form of the matrix.  Returns FALSE if the matrix is invalid,
    TRUE otherwise.  Not user-callable.
*/

/* For UMFPACK, this always returns TRUE */

PRIVATE Int init_rows_cols	/* returns TRUE if OK, or FALSE otherwise */
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows of A */
    Int n_col,			/* number of columns of A */
    Colamd_Row Row [],		/* of size n_row+1 */
    Colamd_Col Col [],		/* of size n_col+1 */
    Int A [],			/* row indices of A, of size Alen */
    Int p []			/* pointers to columns in A, of size n_col+1 */
/*
    Int stats [COLAMD_STATS]	colamd statistics, removed for UMFPACK
*/
)
{
    /* === Local variables ================================================== */

    Int col ;			/* a column index */
    Int row ;			/* a row index */
    Int *cp ;			/* a column pointer */
    Int *cp_end ;		/* a pointer to the end of a column */

    /* === Initialize columns, and check column pointers ==================== */

    for (col = 0 ; col < n_col ; col++)
    {
	Col [col].start = p [col] ;
	Col [col].length = p [col+1] - p [col] ;

#if 0
	if (Col [col].length < 0)
	{
	    /* column pointers must be non-decreasing */
	    stats [COLAMD_STATUS] = COLAMD_ERROR_col_length_negative ;
	    stats [COLAMD_INFO1] = col ;
	    stats [COLAMD_INFO2] = Col [col].length ;
	    DEBUG0 (("colamd: col "ID" length "ID" <= 0\n",
		col, Col [col].length));
	    return (FALSE) ;
	}
#endif

	ASSERT (Col [col].length >= 0) ;

	/* added for UMFPACK v4.1 */
	ASSERT (Col [col].length > 0) ;

	Col [col].shared1.thickness = 1 ;
	Col [col].shared2.score = 0 ;
	Col [col].shared3.prev = EMPTY ;
	Col [col].shared4.degree_next = EMPTY ;

	/* ------------------ */
	/* added for UMFPACK: */
	Col [col].nextcol = EMPTY ;
	Col [col].lastcol = col ;
	/* ------------------ */
    }

    /* p [0..n_col] no longer needed, used as "head" in subsequent routines */

    /* === Scan columns, compute row degrees, and check row indices ========= */

    /* ------------------ */
    /* stats [COLAMD_INFO3] = 0 ; */
    /* number of duplicate or unsorted row indices - not computed in UMFPACK */
    /* ------------------ */

    for (row = 0 ; row < n_row ; row++)
    {
	Row [row].length = 0 ;
	/* ------------------ */
	/* removed for UMFPACK */
	/* Row [row].shared2.mark = -1 ; */
	/* ------------------ */
	/* ------------------ */
	/* added for UMFPACK: */
	Row [row].thickness = 1 ;
	Row [row].front = EMPTY ;
	/* ------------------ */
    }

    for (col = 0 ; col < n_col ; col++)
    {
#ifndef NDEBUG
	Int last_row = -1 ;
#endif

	cp = &A [p [col]] ;
	cp_end = &A [p [col+1]] ;

	while (cp < cp_end)
	{
	    row = *cp++ ;

#if 0
	    /* make sure row indices within range */
	    if (row < 0 || row >= n_row)
	    {
		stats [COLAMD_STATUS] = COLAMD_ERROR_row_index_out_of_bounds ;
		stats [COLAMD_INFO1] = col ;
		stats [COLAMD_INFO2] = row ;
		/* ------------------ */
		/* not needed in UMFPACK: */
		/* stats [COLAMD_INFO3] = n_row ; */
		/* ------------------ */
		DEBUG0 (("colamd: row "ID" col "ID" out of bounds\n", row,col));
		return (FALSE) ;
	    }
#endif

	    ASSERT (row >= 0 && row < n_row) ;

#if 0
	    /* ------------------ */
	    /* changed for UMFPACK */
	    if (row <= last_row)
	    {
		/* row index are unsorted or repeated (or both), thus col */
		/* is jumbled.  This is an error condition for UMFPACK */
		stats [COLAMD_STATUS] = COLAMD_ERROR_jumbled_matrix ;
		stats [COLAMD_INFO1] = col ;
		stats [COLAMD_INFO2] = row ;
		DEBUG1 (("colamd: row "ID" col "ID" unsorted/duplicate\n",
		    row, col)) ;
		return (FALSE) ;
	    }
	    /* ------------------ */
#endif

	    ASSERT (row > last_row) ;

	    /* ------------------ */
	    /* changed for UMFPACK - jumbled columns not tolerated */
	    Row [row].length++ ;
	    /* ------------------ */

#ifndef NDEBUG
	    last_row = row ;
#endif
	}
    }

    /* === Compute row pointers ============================================= */

    /* row form of the matrix starts directly after the column */
    /* form of matrix in A */
    Row [0].start = p [n_col] ;
    Row [0].shared1.p = Row [0].start ;
    /* ------------------ */
    /* removed for UMFPACK */
    /* Row [0].shared2.mark = -1 ; */
    /* ------------------ */
    for (row = 1 ; row < n_row ; row++)
    {
	Row [row].start = Row [row-1].start + Row [row-1].length ;
	Row [row].shared1.p = Row [row].start ;
	/* ------------------ */
	/* removed for UMFPACK */
	/* Row [row].shared2.mark = -1 ; */
	/* ------------------ */
    }

    /* === Create row form ================================================== */

    /* ------------------ */
    /* jumbled matrix case removed for UMFPACK */
    /* ------------------ */

	for (col = 0 ; col < n_col ; col++)
	{
	    cp = &A [p [col]] ;
	    cp_end = &A [p [col+1]] ;
	    while (cp < cp_end)
	    {
		A [(Row [*cp++].shared1.p)++] = col ;
	    }
	}

    /* === Clear the row marks and set row degrees ========================== */

    for (row = 0 ; row < n_row ; row++)
    {
	Row [row].shared2.mark = 0 ;
	Row [row].shared1.degree = Row [row].length ;
    }

    /* ------------------ */
    /* recreate columns for jumbled matrix case removed for UMFPACK */
    /* ------------------ */

    return (TRUE) ;
}


/* ========================================================================== */
/* === init_scoring ========================================================= */
/* ========================================================================== */

/*
    Kills dense or empty columns and rows, calculates an initial score for
    each column, and places all columns in the degree lists.  Not user-callable.
*/

PRIVATE void init_scoring
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows of A */
    Int n_col,			/* number of columns of A */
    Colamd_Row Row [],		/* of size n_row+1 */
    Colamd_Col Col [],		/* of size n_col+1 */
    Int A [],			/* column form and row form of A */
    Int head [],		/* of size n_col+1 */
    double knobs [COLAMD_KNOBS],/* parameters */
    Int *p_n_row2,		/* number of non-dense, non-empty rows */
    Int *p_n_col2,		/* number of non-dense, non-empty columns */
    Int *p_max_deg		/* maximum row degree */
    /* ------------------ */
    /* added for UMFPACK */
    , Int *p_ndense_row		/* number of dense rows */
    , Int *p_nempty_row		/* number of original empty rows */
    , Int *p_nnewlyempty_row	/* number of newly empty rows */
    , Int *p_ndense_col		/* number of dense cols (excl "empty" cols) */
    , Int *p_nempty_col		/* number of original empty cols */
    , Int *p_nnewlyempty_col	/* number of newly empty cols */
    /* ------------------ */
)
{
    /* === Local variables ================================================== */

    Int c ;			/* a column index */
    Int r, row ;		/* a row index */
    Int *cp ;			/* a column pointer */
    Int deg ;			/* degree of a row or column */
    Int *cp_end ;		/* a pointer to the end of a column */
    Int *new_cp ;		/* new column pointer */
    Int col_length ;		/* length of pruned column */
    Int score ;			/* current column score */
    Int n_col2 ;		/* number of non-dense, non-empty columns */
    Int n_row2 ;		/* number of non-dense, non-empty rows */
    Int dense_row_count ;	/* remove rows with more entries than this */
    Int dense_col_count ;	/* remove cols with more entries than this */
    Int min_score ;		/* smallest column score */
    Int max_deg ;		/* maximum row degree */
    Int next_col ;		/* Used to add to degree list.*/

    /* ------------------ */
    /* added for UMFPACK */
    Int ndense_row ;		/* number of dense rows */
    Int nempty_row ;		/* number of empty rows */
    Int nnewlyempty_row ;	/* number of newly empty rows */
    Int ndense_col ;		/* number of dense cols (excl "empty" cols) */
    Int nempty_col ;		/* number of original empty cols */
    Int nnewlyempty_col ;	/* number of newly empty cols */
    Int ne ;
    /* ------------------ */

#ifndef NDEBUG
    Int debug_count ;		/* debug only. */
#endif /* NDEBUG */

    /* === Extract knobs ==================================================== */

    /* --------------------- */
    /* old dense row/column knobs:
    dense_row_count = MAX (0, MIN (knobs [COLAMD_DENSE_ROW] * n_col, n_col)) ;
    dense_col_count = MAX (0, MIN (knobs [COLAMD_DENSE_COL] * n_row, n_row)) ;
    */
    /* new, for UMFPACK: */
    /* Note: if knobs contains a NaN, this is undefined: */
    dense_row_count =
	UMFPACK_DENSE_DEGREE_THRESHOLD (knobs [COLAMD_DENSE_ROW], n_col) ;
    dense_col_count =
	UMFPACK_DENSE_DEGREE_THRESHOLD (knobs [COLAMD_DENSE_COL], n_row) ;
    /* Make sure dense_*_count is between 0 and n: */
    dense_row_count = MAX (0, MIN (dense_row_count, n_col)) ;
    dense_col_count = MAX (0, MIN (dense_col_count, n_row)) ;
    /* --------------------- */

    DEBUG1 (("colamd: densecount: "ID" "ID"\n",
	dense_row_count, dense_col_count)) ;
    max_deg = 0 ;
    n_col2 = n_col ;
    n_row2 = n_row ;

    /* --------------------- */
    /* added for UMFPACK */
    ndense_col = 0 ;
    nempty_col = 0 ;
    nnewlyempty_col = 0 ;
    ndense_row = 0 ;
    nempty_row = 0 ;
    nnewlyempty_row = 0 ;
    /* --------------------- */

    /* === Kill empty columns =============================================== */

    /* removed for UMFPACK v4.1.  prune_singletons has already removed empty
     * columns and empty rows */

#if 0
    /* Put the empty columns at the end in their natural order, so that LU */
    /* factorization can proceed as far as possible. */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	deg = Col [c].length ;
	if (deg == 0)
	{
	    /* this is a empty column, kill and order it last */
	    Col [c].shared2.order = --n_col2 ;
	    KILL_PRINCIPAL_COL (c) ;
	    /* --------------------- */
	    /* added for UMFPACK */
	    nempty_col++ ;
	    /* --------------------- */
	}
    }
    DEBUG1 (("colamd: null columns killed: "ID"\n", n_col - n_col2)) ;
#endif

#ifndef NDEBUG
    for (c = 0 ; c < n_col ; c++)
    {
	ASSERT (Col [c].length > 0) ;
    }
#endif

    /* === Count null rows ================================================== */

#if 0
    for (r = 0 ; r < n_row ; r++)
    {
	deg = Row [r].shared1.degree ;
	if (deg == 0)
	{
	    /* this is an original empty row */
	    nempty_row++ ;
	}
    }
#endif

#ifndef NDEBUG
    for (r = 0 ; r < n_row ; r++)
    {
	ASSERT (Row [r].shared1.degree > 0) ;
	ASSERT (Row [r].length > 0) ;
    }
#endif

    /* === Kill dense columns =============================================== */

    /* Put the dense columns at the end, in their natural order */
    for (c = n_col-1 ; c >= 0 ; c--)
    {

	/* ----------------------------------------------------------------- */
#if 0
	/* removed for UMFPACK v4.1: no empty columns */
	/* skip any dead columns */
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
#endif
	ASSERT (COL_IS_ALIVE (c)) ;
	ASSERT (Col [c].length > 0) ;
	/* ----------------------------------------------------------------- */

	deg = Col [c].length ;
	if (deg > dense_col_count)
	{
	    /* this is a dense column, kill and order it last */
	    Col [c].shared2.order = --n_col2 ;
	    /* --------------------- */
	    /* added for UMFPACK */
	    ndense_col++ ;
	    /* --------------------- */
	    /* decrement the row degrees */
	    cp = &A [Col [c].start] ;
	    cp_end = cp + Col [c].length ;
	    while (cp < cp_end)
	    {
		Row [*cp++].shared1.degree-- ;
	    }
	    KILL_PRINCIPAL_COL (c) ;
	}
    }
    DEBUG1 (("colamd: Dense and null columns killed: "ID"\n", n_col - n_col2)) ;

    /* === Kill dense and empty rows ======================================== */

    /* Note that there can now be empty rows, since dense columns have
     * been deleted.  These are "newly" empty rows. */

    ne = 0 ;
    for (r = 0 ; r < n_row ; r++)
    {
	deg = Row [r].shared1.degree ;
	ASSERT (deg >= 0 && deg <= n_col) ;
	/* --------------------- */
	/* added for UMFPACK */
	if (deg > dense_row_count)
	{
	    /* There is at least one dense row.  Continue ordering, but */
	    /* symbolic factorization will be redone after UMF_colamd is done.*/
	    ndense_row++ ;
	}
	if (deg == 0)
	{
	    /* this is a newly empty row, or original empty row */
	    ne++ ;
	}
	/* --------------------- */
	if (deg > dense_row_count || deg == 0)
	{
	    /* kill a dense or empty row */
	    KILL_ROW (r) ;
	    /* --------------------- */
	    /* added for UMFPACK */
	    Row [r].thickness = 0 ;
	    /* --------------------- */
	    --n_row2 ;
	}
	else
	{
	    /* keep track of max degree of remaining rows */
	    max_deg = MAX (max_deg, deg) ;
	}
    }
    nnewlyempty_row = ne - nempty_row ;
    DEBUG1 (("colamd: Dense rows killed: "ID"\n", ndense_row)) ;
    DEBUG1 (("colamd: Dense and null rows killed: "ID"\n", n_row - n_row2)) ;

    /* === Compute initial column scores ==================================== */

    /* At this point the row degrees are accurate.  They reflect the number */
    /* of "live" (non-dense) columns in each row.  No empty rows exist. */
    /* Some "live" columns may contain only dead rows, however.  These are */
    /* pruned in the code below. */

    /* now find the initial matlab score for each column */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	/* skip dead column */
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
	score = 0 ;
	cp = &A [Col [c].start] ;
	new_cp = cp ;
	cp_end = cp + Col [c].length ;
	while (cp < cp_end)
	{
	    /* get a row */
	    row = *cp++ ;
	    /* skip if dead */
	    if (ROW_IS_DEAD (row))
	    {
		continue ;
	    }
	    /* compact the column */
	    *new_cp++ = row ;
	    /* add row's external degree */
	    score += Row [row].shared1.degree - 1 ;
	    /* guard against integer overflow */
	    score = MIN (score, n_col) ;
	}
	/* determine pruned column length */
	col_length = (Int) (new_cp - &A [Col [c].start]) ;
	if (col_length == 0)
	{
	    /* a newly-made null column (all rows in this col are "dense" */
	    /* and have already been killed) */
	    DEBUG2 (("Newly null killed: "ID"\n", c)) ;
	    Col [c].shared2.order = --n_col2 ;
	    KILL_PRINCIPAL_COL (c) ;
	    /* --------------------- */
	    /* added for UMFPACK */
	    nnewlyempty_col++ ;
	    /* --------------------- */
	}
	else
	{
	    /* set column length and set score */
	    ASSERT (score >= 0) ;
	    ASSERT (score <= n_col) ;
	    Col [c].length = col_length ;
	    Col [c].shared2.score = score ;
	}
    }
    DEBUG1 (("colamd: Dense, null, and newly-null columns killed: "ID"\n",
	n_col-n_col2)) ;

    /* At this point, all empty rows and columns are dead.  All live columns */
    /* are "clean" (containing no dead rows) and simplicial (no supercolumns */
    /* yet).  Rows may contain dead columns, but all live rows contain at */
    /* least one live column. */

#ifndef NDEBUG
    debug_structures (n_row, n_col, Row, Col, A, n_col2) ;
#endif /* NDEBUG */

    /* === Initialize degree lists ========================================== */

#ifndef NDEBUG
    debug_count = 0 ;
#endif /* NDEBUG */

    /* clear the hash buckets */
    for (c = 0 ; c <= n_col ; c++)
    {
	head [c] = EMPTY ;
    }
    min_score = n_col ;
    /* place in reverse order, so low column indices are at the front */
    /* of the lists.  This is to encourage natural tie-breaking */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	/* only add principal columns to degree lists */
	if (COL_IS_ALIVE (c))
	{
	    DEBUG4 (("place "ID" score "ID" minscore "ID" ncol "ID"\n",
		c, Col [c].shared2.score, min_score, n_col)) ;

	    /* === Add columns score to DList =============================== */

	    score = Col [c].shared2.score ;

	    ASSERT (min_score >= 0) ;
	    ASSERT (min_score <= n_col) ;
	    ASSERT (score >= 0) ;
	    ASSERT (score <= n_col) ;
	    ASSERT (head [score] >= EMPTY) ;

	    /* now add this column to dList at proper score location */
	    next_col = head [score] ;
	    Col [c].shared3.prev = EMPTY ;
	    Col [c].shared4.degree_next = next_col ;

	    /* if there already was a column with the same score, set its */
	    /* previous pointer to this new column */
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = c ;
	    }
	    head [score] = c ;

	    /* see if this score is less than current min */
	    min_score = MIN (min_score, score) ;

#ifndef NDEBUG
	    debug_count++ ;
#endif /* NDEBUG */

	}
    }

#ifndef NDEBUG
    DEBUG1 (("colamd: Live cols "ID" out of "ID", non-princ: "ID"\n",
	debug_count, n_col, n_col-debug_count)) ;
    ASSERT (debug_count == n_col2) ;
    debug_deg_lists (n_row, n_col, Row, Col, head, min_score, n_col2, max_deg) ;
#endif /* NDEBUG */

    /* === Return number of remaining columns, and max row degree =========== */

    *p_n_col2 = n_col2 ;
    *p_n_row2 = n_row2 ;
    *p_max_deg = max_deg ;

    /* --------------------- */
    /* added for UMFPACK */
    *p_ndense_row = ndense_row ;
    *p_nempty_row = nempty_row ;	/* original empty rows */
    *p_nnewlyempty_row = nnewlyempty_row ;
    *p_ndense_col = ndense_col ;
    *p_nempty_col = nempty_col ;	/* original empty cols */
    *p_nnewlyempty_col = nnewlyempty_col ;
    /* --------------------- */
}


/* ========================================================================== */
/* === find_ordering ======================================================== */
/* ========================================================================== */

/*
    Order the principal columns of the supercolumn form of the matrix
    (no supercolumns on input).  Uses a minimum approximate column minimum
    degree ordering method.  Not user-callable.
*/

PRIVATE Int find_ordering	/* return the number of garbage collections */
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows of A */
    Int n_col,			/* number of columns of A */
    Int Alen,			/* size of A, 2*nnz + n_col or larger */
    Colamd_Row Row [],		/* of size n_row+1 */
    Colamd_Col Col [],		/* of size n_col+1 */
    Int A [],			/* column form and row form of A */
    Int head [],		/* of size n_col+1 */
    Int n_col2,			/* Remaining columns to order */
    Int max_deg,		/* Maximum row degree */
    Int pfree			/* index of first free slot (2*nnz on entry) */
    /* ------------------ */
    /* added for UMFPACK: */
    , Int Front_npivcol [ ]
    , Int Front_nrows [ ]
    , Int Front_ncols [ ]
    , Int Front_parent [ ]
    , Int Front_cols [ ]
    , Int *p_nfr		/* number of fronts */
    , Int aggressive
    , Int InFront [ ]
    /* ------------------ */
)
{
    /* === Local variables ================================================== */

    Int k ;			/* current pivot ordering step */
    Int pivot_col ;		/* current pivot column */
    Int *cp ;			/* a column pointer */
    Int *rp ;			/* a row pointer */
    Int pivot_row ;		/* current pivot row */
    Int *new_cp ;		/* modified column pointer */
    Int *new_rp ;		/* modified row pointer */
    Int pivot_row_start ;	/* pointer to start of pivot row */
    Int pivot_row_degree ;	/* number of columns in pivot row */
    Int pivot_row_length ;	/* number of supercolumns in pivot row */
    Int pivot_col_score ;	/* score of pivot column */
    Int needed_memory ;		/* free space needed for pivot row */
    Int *cp_end ;		/* pointer to the end of a column */
    Int *rp_end ;		/* pointer to the end of a row */
    Int row ;			/* a row index */
    Int col ;			/* a column index */
    Int max_score ;		/* maximum possible score */
    Int cur_score ;		/* score of current column */
    unsigned Int hash ;		/* hash value for supernode detection */
    Int head_column ;		/* head of hash bucket */
    Int first_col ;		/* first column in hash bucket */
    Int tag_mark ;		/* marker value for mark array */
    Int row_mark ;		/* Row [row].shared2.mark */
    Int set_difference ;	/* set difference size of row with pivot row */
    Int min_score ;		/* smallest column score */
    Int col_thickness ;		/* "thickness" (no. of columns in a supercol) */
    Int max_mark ;		/* maximum value of tag_mark */
    Int pivot_col_thickness ;	/* number of columns represented by pivot col */
    Int prev_col ;		/* Used by Dlist operations. */
    Int next_col ;		/* Used by Dlist operations. */
    Int ngarbage ;		/* number of garbage collections performed */

#ifndef NDEBUG
    Int debug_d ;		/* debug loop counter */
    Int debug_step = 0 ;	/* debug loop counter */
#endif /* NDEBUG */

    /* ------------------ */
    /* added for UMFPACK: */
    Int pivot_row_thickness ;	/* number of rows represented by pivot row */
    Int nfr = 0 ;		/* number of fronts */
    Int child ;
    /* ------------------ */

    /* === Initialization and clear mark ==================================== */

    max_mark = MAX_MARK (n_col) ;	/* defined in umfpack.h */
    tag_mark = clear_mark (n_row, Row) ;
    min_score = 0 ;
    ngarbage = 0 ;
    DEBUG1 (("colamd: Ordering, n_col2="ID"\n", n_col2)) ;

    for (row = 0 ; row < n_row ; row++)
    {
	InFront [row] = EMPTY ;
    }

    /* === Order the columns ================================================ */

    for (k = 0 ; k < n_col2 ; /* 'k' is incremented below */)
    {

#ifndef NDEBUG
	if (debug_step % 100 == 0)
	{
	    DEBUG2 (("\n...  Step k: "ID" out of n_col2: "ID"\n", k, n_col2)) ;
	}
	else
	{
	    DEBUG3 (("\n-----Step k: "ID" out of n_col2: "ID"\n", k, n_col2)) ;
	}
	debug_step++ ;
	debug_deg_lists (n_row, n_col, Row, Col, head,
		min_score, n_col2-k, max_deg) ;
	debug_matrix (n_row, n_col, Row, Col, A) ;
#endif /* NDEBUG */

	/* === Select pivot column, and order it ============================ */

	/* make sure degree list isn't empty */
	ASSERT (min_score >= 0) ;
	ASSERT (min_score <= n_col) ;
	ASSERT (head [min_score] >= EMPTY) ;

#ifndef NDEBUG
	for (debug_d = 0 ; debug_d < min_score ; debug_d++)
	{
	    ASSERT (head [debug_d] == EMPTY) ;
	}
#endif /* NDEBUG */

	/* get pivot column from head of minimum degree list */
	while (head [min_score] == EMPTY && min_score < n_col)
	{
	    min_score++ ;
	}
	pivot_col = head [min_score] ;
	ASSERT (pivot_col >= 0 && pivot_col <= n_col) ;
	next_col = Col [pivot_col].shared4.degree_next ;
	head [min_score] = next_col ;
	if (next_col != EMPTY)
	{
	    Col [next_col].shared3.prev = EMPTY ;
	}

	ASSERT (COL_IS_ALIVE (pivot_col)) ;
	DEBUG3 (("Pivot col: "ID"\n", pivot_col)) ;

	/* remember score for defrag check */
	pivot_col_score = Col [pivot_col].shared2.score ;

	/* the pivot column is the kth column in the pivot order */
	Col [pivot_col].shared2.order = k ;

	/* increment order count by column thickness */
	pivot_col_thickness = Col [pivot_col].shared1.thickness ;
	/* ------------------ */
	/* changed for UMFPACK: */
	k += pivot_col_thickness ;
	/* ------------------ */
	ASSERT (pivot_col_thickness > 0) ;

	/* === Garbage_collection, if necessary ============================= */

	needed_memory = MIN (pivot_col_score, n_col - k) ;
	if (pfree + needed_memory >= Alen)
	{
	    pfree = garbage_collection (n_row, n_col, Row, Col, A, &A [pfree]) ;
	    ngarbage++ ;
	    /* after garbage collection we will have enough */
	    ASSERT (pfree + needed_memory < Alen) ;
	    /* garbage collection has wiped out the Row[].shared2.mark array */
	    tag_mark = clear_mark (n_row, Row) ;

#ifndef NDEBUG
	    debug_matrix (n_row, n_col, Row, Col, A) ;
#endif /* NDEBUG */
	}

	/* === Compute pivot row pattern ==================================== */

	/* get starting location for this new merged row */
	pivot_row_start = pfree ;

	/* initialize new row counts to zero */
	pivot_row_degree = 0 ;

	/* ------------------ */
	/* added for UMFPACK: */
	pivot_row_thickness = 0 ;
	/* ------------------ */

	/* [ tag pivot column as having been visited so it isn't included */
	/* in merged pivot row */
	Col [pivot_col].shared1.thickness = -pivot_col_thickness ;

	/* pivot row is the union of all rows in the pivot column pattern */
	cp = &A [Col [pivot_col].start] ;
	cp_end = cp + Col [pivot_col].length ;
	while (cp < cp_end)
	{
	    /* get a row */
	    row = *cp++ ;
	    DEBUG4 (("Pivot col pattern %d "ID"\n", ROW_IS_ALIVE(row), row)) ;
	    /* skip if row is dead */
	    if (ROW_IS_DEAD (row))
	    {
		continue ;
	    }

	    /* ------------------ */
	    /* added for UMFPACK: */
	    /* sum the thicknesses of all the rows */
	    /* ASSERT (Row [row].thickness > 0) ; */
	    pivot_row_thickness += Row [row].thickness ;
	    /* ------------------ */

	    rp = &A [Row [row].start] ;
	    rp_end = rp + Row [row].length ;
	    while (rp < rp_end)
	    {
		/* get a column */
		col = *rp++ ;
		/* add the column, if alive and untagged */
		col_thickness = Col [col].shared1.thickness ;
		if (col_thickness > 0 && COL_IS_ALIVE (col))
		{
		    /* tag column in pivot row */
		    Col [col].shared1.thickness = -col_thickness ;
		    ASSERT (pfree < Alen) ;
		    /* place column in pivot row */
		    A [pfree++] = col ;
		    pivot_row_degree += col_thickness ;
		    /* ------------------ */
		    /* added for UMFPACK: */
		    DEBUG4 (("\t\t\tNew live column in pivot row: "ID"\n",col));
		    /* ------------------ */
		}
		/* ------------------ */
		/* added for UMFPACK */
#ifndef NDEBUG
		if (col_thickness < 0 && COL_IS_ALIVE (col))
		{
		    DEBUG4 (("\t\t\tOld live column in pivot row: "ID"\n",col));
		}
#endif
		/* ------------------ */
	    }
	}

	/* ------------------ */
	/* added for UMFPACK: */
	/* pivot_row_thickness is the number of rows in frontal matrix */
	/* both pivotal rows and nonpivotal rows */
	/* ------------------ */

	/* clear tag on pivot column */
	Col [pivot_col].shared1.thickness = pivot_col_thickness ;	/* ] */
	max_deg = MAX (max_deg, pivot_row_degree) ;

#ifndef NDEBUG
	DEBUG3 (("check2\n")) ;
	debug_mark (n_row, Row, tag_mark, max_mark) ;
#endif /* NDEBUG */

	/* === Kill all rows used to construct pivot row ==================== */

	/* also kill pivot row, temporarily */
	cp = &A [Col [pivot_col].start] ;
	cp_end = cp + Col [pivot_col].length ;
	while (cp < cp_end)
	{
	    /* may be killing an already dead row */
	    row = *cp++ ;

	    DEBUG2 (("Kill row in pivot col: "ID" alive? %d, front "ID"\n",
		row, ROW_IS_ALIVE (row), Row [row].front)) ;

	    /* added for UMFPACK: */
	    if (ROW_IS_ALIVE (row))
	    {
		if (Row [row].front != EMPTY)
		{
		    /* This row represents a frontal matrix. */
		    /* Row [row].front is a child of current front */
		    child = Row [row].front ;
		    Front_parent [child] = nfr ;
		    DEBUG1 (("Front "ID" => front "ID", normal\n", child, nfr));
		}
		else
		{
		    /* This is an original row.  Keep track of which front
		     * is its parent in the row-merge tree. */
		    InFront [row] = nfr ;
		    DEBUG1 (("Row "ID" => front "ID", normal\n", row, nfr)) ;
		}
	    }

	    KILL_ROW (row) ;

	    /* ------------------ */
	    /* added for UMFPACK: */
	    Row [row].thickness = 0 ;
	    /* ------------------ */
	}

	/* === Select a row index to use as the new pivot row =============== */

	pivot_row_length = pfree - pivot_row_start ;
	if (pivot_row_length > 0)
	{
	    /* pick the "pivot" row arbitrarily (first row in col) */
	    pivot_row = A [Col [pivot_col].start] ;
	    DEBUG3 (("Pivotal row is "ID"\n", pivot_row)) ;
	}
	else
	{
	    /* there is no pivot row, since it is of zero length */
	    pivot_row = EMPTY ;
	    ASSERT (pivot_row_length == 0) ;
	}
	ASSERT (Col [pivot_col].length > 0 || pivot_row_length == 0) ;

	/* === Approximate degree computation =============================== */

	/* Here begins the computation of the approximate degree.  The column */
	/* score is the sum of the pivot row "length", plus the size of the */
	/* set differences of each row in the column minus the pattern of the */
	/* pivot row itself.  The column ("thickness") itself is also */
	/* excluded from the column score (we thus use an approximate */
	/* external degree). */

	/* The time taken by the following code (compute set differences, and */
	/* add them up) is proportional to the size of the data structure */
	/* being scanned - that is, the sum of the sizes of each column in */
	/* the pivot row.  Thus, the amortized time to compute a column score */
	/* is proportional to the size of that column (where size, in this */
	/* context, is the column "length", or the number of row indices */
	/* in that column).  The number of row indices in a column is */
	/* monotonically non-decreasing, from the length of the original */
	/* column on input to colamd. */

	/* === Compute set differences ====================================== */

	DEBUG3 (("** Computing set differences phase. **\n")) ;

	/* pivot row is currently dead - it will be revived later. */

	DEBUG3 (("Pivot row: \n")) ;
	/* for each column in pivot row */
	rp = &A [pivot_row_start] ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    col = *rp++ ;
	    ASSERT (COL_IS_ALIVE (col) && col != pivot_col) ;
	    DEBUG3 (("    Col: "ID"\n", col)) ;

	    /* clear tags used to construct pivot row pattern */
	    col_thickness = -Col [col].shared1.thickness ;
	    ASSERT (col_thickness > 0) ;
	    Col [col].shared1.thickness = col_thickness ;

	    /* === Remove column from degree list =========================== */

	    cur_score = Col [col].shared2.score ;
	    prev_col = Col [col].shared3.prev ;
	    next_col = Col [col].shared4.degree_next ;
	    ASSERT (cur_score >= 0) ;
	    ASSERT (cur_score <= n_col) ;
	    ASSERT (cur_score >= EMPTY) ;
	    if (prev_col == EMPTY)
	    {
		head [cur_score] = next_col ;
	    }
	    else
	    {
		Col [prev_col].shared4.degree_next = next_col ;
	    }
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = prev_col ;
	    }

	    /* === Scan the column ========================================== */

	    cp = &A [Col [col].start] ;
	    cp_end = cp + Col [col].length ;
	    while (cp < cp_end)
	    {
		/* get a row */
		row = *cp++ ;
		row_mark = Row [row].shared2.mark ;
		/* skip if dead */
		if (ROW_IS_MARKED_DEAD (row_mark))
		{
		    continue ;
		}
		ASSERT (row != pivot_row) ;
		set_difference = row_mark - tag_mark ;
		/* check if the row has been seen yet */
		if (set_difference < 0)
		{
		    ASSERT (Row [row].shared1.degree <= max_deg) ;
		    set_difference = Row [row].shared1.degree ;
		}
		/* subtract column thickness from this row's set difference */
		set_difference -= col_thickness ;
		ASSERT (set_difference >= 0) ;
		ASSERT (ROW_IS_ALIVE (row)) ;

		/* absorb this row if the set difference becomes zero */
		if (set_difference == 0 && aggressive)
		{
		    /* v4.1: do aggressive absorption */
		    DEBUG3 (("aggressive absorption. Row: "ID"\n", row)) ;

		    if (Row [row].front != EMPTY)
		    {
			/* Row [row].front is a child of current front. */
			child = Row [row].front ;
			Front_parent [child] = nfr ;
			DEBUG1 (("Front "ID" => front "ID", aggressive\n",
				    child, nfr)) ;
		    }
		    else
		    {
			/* this is an original row.  Keep track of which front
			 * assembles it, for the row-merge tree */
			InFront [row] = nfr ;
			DEBUG1 (("Row "ID" => front "ID", aggressive\n",
				    row, nfr)) ;
		    }

		    KILL_ROW (row) ;

		    /* sum the thicknesses of all the rows */
		    /* ASSERT (Row [row].thickness > 0) ; */
		    pivot_row_thickness += Row [row].thickness ;
		    Row [row].thickness = 0 ;

		}
		else
		{
		    /* save the new mark */
		    Row [row].shared2.mark = set_difference + tag_mark ;
		}
	    }
	}

#ifndef NDEBUG
	debug_deg_lists (n_row, n_col, Row, Col, head,
		min_score, n_col2-k-pivot_row_degree, max_deg) ;
#endif /* NDEBUG */

	/* === Add up set differences for each column ======================= */

	DEBUG3 (("** Adding set differences phase. **\n")) ;

	/* for each column in pivot row */
	rp = &A [pivot_row_start] ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    /* get a column */
	    col = *rp++ ;
	    ASSERT (COL_IS_ALIVE (col) && col != pivot_col) ;
	    hash = 0 ;
	    cur_score = 0 ;
	    cp = &A [Col [col].start] ;
	    /* compact the column */
	    new_cp = cp ;
	    cp_end = cp + Col [col].length ;

	    DEBUG4 (("Adding set diffs for Col: "ID".\n", col)) ;

	    while (cp < cp_end)
	    {
		/* get a row */
		row = *cp++ ;
		ASSERT(row >= 0 && row < n_row) ;
		row_mark = Row [row].shared2.mark ;
		/* skip if dead */
		if (ROW_IS_MARKED_DEAD (row_mark))
		{
		    /* ------------------ */
		    /* changed for UMFPACK: */
		    DEBUG4 ((" Row "ID", dead\n", row)) ;
		    /* ------------------ */
		    continue ;
		}
		/* ------------------ */
		/* changed for UMFPACK: */
		/* ASSERT (row_mark > tag_mark) ; */
		DEBUG4 ((" Row "ID", set diff "ID"\n", row, row_mark-tag_mark));
		ASSERT (row_mark >= tag_mark) ;
		/* ------------------ */
		/* compact the column */
		*new_cp++ = row ;
		/* compute hash function */
		hash += row ;
		/* add set difference */
		cur_score += row_mark - tag_mark ;
		/* integer overflow... */
		cur_score = MIN (cur_score, n_col) ;
	    }

	    /* recompute the column's length */
	    Col [col].length = (Int) (new_cp - &A [Col [col].start]) ;

	    /* === Further mass elimination ================================= */

	    if (Col [col].length == 0)
	    {
		DEBUG4 (("further mass elimination. Col: "ID"\n", col)) ;
		/* nothing left but the pivot row in this column */
		KILL_PRINCIPAL_COL (col) ;
		pivot_row_degree -= Col [col].shared1.thickness ;
		ASSERT (pivot_row_degree >= 0) ;
		/* order it */
		Col [col].shared2.order = k ;
		/* increment order count by column thickness */
		k += Col [col].shared1.thickness ;

		/* ------------------ */
		/* added for UMFPACK: */
		pivot_col_thickness += Col [col].shared1.thickness ;

		/* add to column list of front ... */
#ifndef NDEBUG
		DEBUG1 (("Mass")) ;
		dump_super (col, Col, n_col) ;
#endif
		Col [Col [col].lastcol].nextcol = Front_cols [nfr] ;
		Front_cols [nfr] = col ;
		/* ------------------ */

	    }
	    else
	    {
		/* === Prepare for supercolumn detection ==================== */

		DEBUG4 (("Preparing supercol detection for Col: "ID".\n", col));

		/* save score so far */
		Col [col].shared2.score = cur_score ;

		/* add column to hash table, for supercolumn detection */
		/* NOTE: hash is an unsigned Int to avoid a problem in ANSI C.
		 * The sign of the expression a % b is not defined when a and/or
		 * b are negative.  Since hash is unsigned and n_col >= 0,
		 * this problem is avoided. */
		hash %= n_col + 1 ;

		DEBUG4 ((" Hash = "ID", n_col = "ID".\n", (Int) hash, n_col)) ;
		ASSERT (((Int) hash) <= n_col) ;

		head_column = head [hash] ;
		if (head_column > EMPTY)
		{
		    /* degree list "hash" is non-empty, use prev (shared3) of */
		    /* first column in degree list as head of hash bucket */
		    first_col = Col [head_column].shared3.headhash ;
		    Col [head_column].shared3.headhash = col ;
		}
		else
		{
		    /* degree list "hash" is empty, use head as hash bucket */
		    first_col = - (head_column + 2) ;
		    head [hash] = - (col + 2) ;
		}
		Col [col].shared4.hash_next = first_col ;

		/* save hash function in Col [col].shared3.hash */
		Col [col].shared3.hash = (Int) hash ;
		ASSERT (COL_IS_ALIVE (col)) ;
	    }
	}

	/* The approximate external column degree is now computed.  */

	/* === Supercolumn detection ======================================== */

	DEBUG3 (("** Supercolumn detection phase. **\n")) ;

	detect_super_cols (

#ifndef NDEBUG
		n_col, Row,
#endif /* NDEBUG */

		Col, A, head, pivot_row_start, pivot_row_length) ;

	/* === Kill the pivotal column ====================================== */

	KILL_PRINCIPAL_COL (pivot_col) ;

	/* ------------------ */
	/* added for UMFPACK: */
	/* add columns to column list of front */
#ifndef NDEBUG
	DEBUG1 (("Pivot")) ;
	dump_super (pivot_col, Col, n_col) ;
#endif
	Col [Col [pivot_col].lastcol].nextcol = Front_cols [nfr] ;
	Front_cols [nfr] = pivot_col ;
	/* ------------------ */

	/* === Clear mark =================================================== */

	tag_mark += (max_deg + 1) ;
	if (tag_mark >= max_mark)
	{
	    DEBUG2 (("clearing tag_mark\n")) ;
	    tag_mark = clear_mark (n_row, Row) ;
	}

#ifndef NDEBUG
	DEBUG3 (("check3\n")) ;
	debug_mark (n_row, Row, tag_mark, max_mark) ;
#endif /* NDEBUG */

	/* === Finalize the new pivot row, and column scores ================ */

	DEBUG3 (("** Finalize scores phase. **\n")) ;
	DEBUG3 (("pivot_row_degree "ID"\n", pivot_row_degree)) ;

	/* for each column in pivot row */
	rp = &A [pivot_row_start] ;
	/* compact the pivot row */
	new_rp = rp ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    col = *rp++ ;
	    DEBUG3 (("Col "ID" \n", col)) ;
	    /* skip dead columns */
	    if (COL_IS_DEAD (col))
	    {
		DEBUG3 (("dead\n")) ;
		continue ;
	    }
	    *new_rp++ = col ;
	    /* add new pivot row to column */
	    A [Col [col].start + (Col [col].length++)] = pivot_row ;

	    /* retrieve score so far and add on pivot row's degree. */
	    /* (we wait until here for this in case the pivot */
	    /* row's degree was reduced due to mass elimination). */
	    cur_score = Col [col].shared2.score + pivot_row_degree ;
	    DEBUG3 ((" cur_score "ID" ", cur_score)) ;

	    /* calculate the max possible score as the number of */
	    /* external columns minus the 'k' value minus the */
	    /* columns thickness */
	    max_score = n_col - k - Col [col].shared1.thickness ;
	    DEBUG3 ((" max_score "ID" ", max_score)) ;

	    /* make the score the external degree of the union-of-rows */
	    cur_score -= Col [col].shared1.thickness ;
	    DEBUG3 ((" cur_score "ID" ", cur_score)) ;

	    /* make sure score is less or equal than the max score */
	    cur_score = MIN (cur_score, max_score) ;
	    ASSERT (cur_score >= 0) ;

	    /* store updated score */
	    Col [col].shared2.score = cur_score ;
	    DEBUG3 ((" "ID"\n", cur_score)) ;

	    /* === Place column back in degree list ========================= */

	    ASSERT (min_score >= 0) ;
	    ASSERT (min_score <= n_col) ;
	    ASSERT (cur_score >= 0) ;
	    ASSERT (cur_score <= n_col) ;
	    ASSERT (head [cur_score] >= EMPTY) ;
	    next_col = head [cur_score] ;
	    Col [col].shared4.degree_next = next_col ;
	    Col [col].shared3.prev = EMPTY ;
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = col ;
	    }
	    head [cur_score] = col ;

	    /* see if this score is less than current min */
	    min_score = MIN (min_score, cur_score) ;

	}

#ifndef NDEBUG
	debug_deg_lists (n_row, n_col, Row, Col, head,
		min_score, n_col2-k, max_deg) ;
#endif /* NDEBUG */

	/* ------------------ */
	/* added for UMFPACK: */
	/* frontal matrix can have more pivot cols than pivot rows for */
	/* singular matrices. */

	/* number of candidate pivot columns */
	Front_npivcol [nfr] = pivot_col_thickness ;

	/* all rows (not just size of contrib. block) */
	Front_nrows [nfr] = pivot_row_thickness ;

	/* all cols */
	Front_ncols [nfr] = pivot_col_thickness + pivot_row_degree ;

	Front_parent [nfr] = EMPTY ;

	pivot_row_thickness -= pivot_col_thickness ;
	DEBUG1 (("Front "ID" Pivot_row_thickness after pivot cols elim: "ID"\n",
	    nfr, pivot_row_thickness)) ;
	pivot_row_thickness = MAX (0, pivot_row_thickness) ;
	/* ------------------ */

	/* === Resurrect the new pivot row ================================== */

	if (pivot_row_degree > 0
	/* ------------------ */
	/* added for UMFPACK.  Note that this part of the expression should be
	 * removed if this routine is used outside of UMFPACK, for a Cholesky
	 * factorization of (AQ)'(AQ) */
	&& pivot_row_thickness > 0
	/* ------------------ */
	)
	{
	    /* update pivot row length to reflect any cols that were killed */
	    /* during super-col detection and mass elimination */
	    Row [pivot_row].start  = pivot_row_start ;
	    Row [pivot_row].length = (Int) (new_rp - &A[pivot_row_start]) ;
	    ASSERT (Row [pivot_row].length > 0) ;
	    Row [pivot_row].shared1.degree = pivot_row_degree ;
	    Row [pivot_row].shared2.mark = 0 ;
	    /* ------------------ */
	    /* added for UMFPACK: */
	    Row [pivot_row].thickness = pivot_row_thickness ;
	    Row [pivot_row].front = nfr ;
	    /* ------------------ */
	    /* pivot row is no longer dead */
	}

	/* ------------------ */
	/* added for UMFPACK: */

#ifndef NDEBUG
	DEBUG1 (("Front "ID" : "ID" "ID" "ID" ", nfr,
		Front_npivcol [nfr], Front_nrows [nfr], Front_ncols [nfr])) ;
	DEBUG1 ((" cols:[ ")) ;
	debug_d = 0 ;
	for (col = Front_cols [nfr] ; col != EMPTY ; col = Col [col].nextcol)
	{
	    DEBUG1 ((" "ID, col)) ;
	    ASSERT (col >= 0 && col < n_col) ;
	    ASSERT (COL_IS_DEAD (col)) ;
	    debug_d++ ;
	    ASSERT (debug_d <= pivot_col_thickness) ;
	}
	ASSERT (debug_d == pivot_col_thickness) ;
	DEBUG1 ((" ]\n ")) ;
#endif
	nfr++ ; /* one more front */
	/* ------------------ */

    }

    /* === All principal columns have now been ordered ====================== */

    /* ------------------ */
    /* added for UMFPACK: */
    *p_nfr = nfr ;
    /* ------------------ */

    return (ngarbage) ;
}


/* ========================================================================== */
/* === order_children deleted for UMFPACK =================================== */
/* ========================================================================== */

/* ========================================================================== */
/* === detect_super_cols ==================================================== */
/* ========================================================================== */

/*
    Detects supercolumns by finding matches between columns in the hash buckets.
    Check amongst columns in the set A [row_start ... row_start + row_length-1].
    The columns under consideration are currently *not* in the degree lists,
    and have already been placed in the hash buckets.

    The hash bucket for columns whose hash function is equal to h is stored
    as follows:

	if head [h] is >= 0, then head [h] contains a degree list, so:

		head [h] is the first column in degree bucket h.
		Col [head [h]].headhash gives the first column in hash bucket h.

	otherwise, the degree list is empty, and:

		-(head [h] + 2) is the first column in hash bucket h.

    For a column c in a hash bucket, Col [c].shared3.prev is NOT a "previous
    column" pointer.  Col [c].shared3.hash is used instead as the hash number
    for that column.  The value of Col [c].shared4.hash_next is the next column
    in the same hash bucket.

    Assuming no, or "few" hash collisions, the time taken by this routine is
    linear in the sum of the sizes (lengths) of each column whose score has
    just been computed in the approximate degree computation.
    Not user-callable.
*/

PRIVATE void detect_super_cols
(
    /* === Parameters ======================================================= */

#ifndef NDEBUG
    /* these two parameters are only needed when debugging is enabled: */
    Int n_col,			/* number of columns of A */
    Colamd_Row Row [],		/* of size n_row+1 */
#endif /* NDEBUG */

    Colamd_Col Col [],		/* of size n_col+1 */
    Int A [],			/* row indices of A */
    Int head [],		/* head of degree lists and hash buckets */
    Int row_start,		/* pointer to set of columns to check */
    Int row_length		/* number of columns to check */
)
{
    /* === Local variables ================================================== */

    Int hash ;			/* hash value for a column */
    Int *rp ;			/* pointer to a row */
    Int c ;			/* a column index */
    Int super_c ;		/* column index of the column to absorb into */
    Int *cp1 ;			/* column pointer for column super_c */
    Int *cp2 ;			/* column pointer for column c */
    Int length ;		/* length of column super_c */
    Int prev_c ;		/* column preceding c in hash bucket */
    Int i ;			/* loop counter */
    Int *rp_end ;		/* pointer to the end of the row */
    Int col ;			/* a column index in the row to check */
    Int head_column ;		/* first column in hash bucket or degree list */
    Int first_col ;		/* first column in hash bucket */

    /* === Consider each column in the row ================================== */

    rp = &A [row_start] ;
    rp_end = rp + row_length ;
    while (rp < rp_end)
    {
	col = *rp++ ;
	if (COL_IS_DEAD (col))
	{
	    continue ;
	}

	/* get hash number for this column */
	hash = Col [col].shared3.hash ;
	ASSERT (hash <= n_col) ;

	/* === Get the first column in this hash bucket ===================== */

	head_column = head [hash] ;
	if (head_column > EMPTY)
	{
	    first_col = Col [head_column].shared3.headhash ;
	}
	else
	{
	    first_col = - (head_column + 2) ;
	}

	/* === Consider each column in the hash bucket ====================== */

	for (super_c = first_col ; super_c != EMPTY ;
	    super_c = Col [super_c].shared4.hash_next)
	{
	    ASSERT (COL_IS_ALIVE (super_c)) ;
	    ASSERT (Col [super_c].shared3.hash == hash) ;
	    length = Col [super_c].length ;

	    /* prev_c is the column preceding column c in the hash bucket */
	    prev_c = super_c ;

	    /* === Compare super_c with all columns after it ================ */

	    for (c = Col [super_c].shared4.hash_next ;
		c != EMPTY ; c = Col [c].shared4.hash_next)
	    {
		ASSERT (c != super_c) ;
		ASSERT (COL_IS_ALIVE (c)) ;
		ASSERT (Col [c].shared3.hash == hash) ;

		/* not identical if lengths or scores are different */
		if (Col [c].length != length ||
		    Col [c].shared2.score != Col [super_c].shared2.score)
		{
		    prev_c = c ;
		    continue ;
		}

		/* compare the two columns */
		cp1 = &A [Col [super_c].start] ;
		cp2 = &A [Col [c].start] ;

		for (i = 0 ; i < length ; i++)
		{
		    /* the columns are "clean" (no dead rows) */
		    ASSERT (ROW_IS_ALIVE (*cp1))  ;
		    ASSERT (ROW_IS_ALIVE (*cp2))  ;
		    /* row indices will same order for both supercols, */
		    /* no gather scatter nessasary */
		    if (*cp1++ != *cp2++)
		    {
			break ;
		    }
		}

		/* the two columns are different if the for-loop "broke" */
		if (i != length)
		{
		    prev_c = c ;
		    continue ;
		}

		/* === Got it!  two columns are identical =================== */

		ASSERT (Col [c].shared2.score == Col [super_c].shared2.score) ;

		Col [super_c].shared1.thickness += Col [c].shared1.thickness ;
		Col [c].shared1.parent = super_c ;
		KILL_NON_PRINCIPAL_COL (c) ;

		Col [c].shared2.order = EMPTY ;
		/* remove c from hash bucket */
		Col [prev_c].shared4.hash_next = Col [c].shared4.hash_next ;

		/* ------------------ */
		/* added for UMFPACK: */
		/* add c to end of list of super_c */
		ASSERT (Col [super_c].lastcol >= 0) ;
		ASSERT (Col [super_c].lastcol < n_col) ;
		Col [Col [super_c].lastcol].nextcol = c ;
		Col [super_c].lastcol = Col [c].lastcol ;
#ifndef NDEBUG
		/* dump the supercolumn */
		DEBUG1 (("Super")) ;
		dump_super (super_c, Col, n_col) ;
#endif
		/* ------------------ */

	    }
	}

	/* === Empty this hash bucket ======================================= */

	if (head_column > EMPTY)
	{
	    /* corresponding degree list "hash" is not empty */
	    Col [head_column].shared3.headhash = EMPTY ;
	}
	else
	{
	    /* corresponding degree list "hash" is empty */
	    head [hash] = EMPTY ;
	}
    }
}


/* ========================================================================== */
/* === garbage_collection =================================================== */
/* ========================================================================== */

/*
    Defragments and compacts columns and rows in the workspace A.  Used when
    all avaliable memory has been used while performing row merging.  Returns
    the index of the first free position in A, after garbage collection.  The
    time taken by this routine is linear is the size of the array A, which is
    itself linear in the number of nonzeros in the input matrix.
    Not user-callable.
*/

PRIVATE Int garbage_collection  /* returns the new value of pfree */
(
    /* === Parameters ======================================================= */

    Int n_row,			/* number of rows */
    Int n_col,			/* number of columns */
    Colamd_Row Row [],		/* row info */
    Colamd_Col Col [],		/* column info */
    Int A [],			/* A [0 ... Alen-1] holds the matrix */
    Int *pfree			/* &A [0] ... pfree is in use */
)
{
    /* === Local variables ================================================== */

    Int *psrc ;			/* source pointer */
    Int *pdest ;		/* destination pointer */
    Int j ;			/* counter */
    Int r ;			/* a row index */
    Int c ;			/* a column index */
    Int length ;		/* length of a row or column */

#ifndef NDEBUG
    Int debug_rows ;
    DEBUG2 (("Defrag..\n")) ;
    for (psrc = &A[0] ; psrc < pfree ; psrc++) ASSERT (*psrc >= 0) ;
    debug_rows = 0 ;
#endif /* NDEBUG */

    /* === Defragment the columns =========================================== */

    pdest = &A[0] ;
    for (c = 0 ; c < n_col ; c++)
    {
	if (COL_IS_ALIVE (c))
	{
	    psrc = &A [Col [c].start] ;

	    /* move and compact the column */
	    ASSERT (pdest <= psrc) ;
	    Col [c].start = (Int) (pdest - &A [0]) ;
	    length = Col [c].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		r = *psrc++ ;
		if (ROW_IS_ALIVE (r))
		{
		    *pdest++ = r ;
		}
	    }
	    Col [c].length = (Int) (pdest - &A [Col [c].start]) ;
	}
    }

    /* === Prepare to defragment the rows =================================== */

    for (r = 0 ; r < n_row ; r++)
    {
	if (ROW_IS_ALIVE (r))
	{
	    if (Row [r].length == 0)
	    {
		/* :: defrag row kill :: */
		/* This row is of zero length.  cannot compact it, so kill it.
		 * NOTE: in the current version, there are no zero-length live
		 * rows when garbage_collection is called.  So this code will
		 * never trigger.  However, if the code is modified, or if
		 * garbage_collection is called at a different place, then rows
		 * can be of zero length.  So this test is kept, just in case.
		 */
		DEBUGm4 (("Defrag row kill\n")) ;
		KILL_ROW (r) ;
	    }
	    else
	    {
		/* save first column index in Row [r].shared2.first_column */
		psrc = &A [Row [r].start] ;
		Row [r].shared2.first_column = *psrc ;
		ASSERT (ROW_IS_ALIVE (r)) ;
		/* flag the start of the row with the one's complement of row */
		*psrc = ONES_COMPLEMENT (r) ;
#ifndef NDEBUG
		debug_rows++ ;
#endif /* NDEBUG */
	    }
	}
    }

    /* === Defragment the rows ============================================== */

    psrc = pdest ;
    while (psrc < pfree)
    {
	/* find a negative number ... the start of a row */
	if (*psrc++ < 0)
	{
	    psrc-- ;
	    /* get the row index */
	    r = ONES_COMPLEMENT (*psrc) ;
	    ASSERT (r >= 0 && r < n_row) ;
	    /* restore first column index */
	    *psrc = Row [r].shared2.first_column ;
	    ASSERT (ROW_IS_ALIVE (r)) ;

	    /* move and compact the row */
	    ASSERT (pdest <= psrc) ;
	    Row [r].start = (Int) (pdest - &A [0]) ;
	    length = Row [r].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		c = *psrc++ ;
		if (COL_IS_ALIVE (c))
		{
		    *pdest++ = c ;
		}
	    }
	    Row [r].length = (Int) (pdest - &A [Row [r].start]) ;

#ifndef NDEBUG
	    debug_rows-- ;
#endif /* NDEBUG */

	}
    }
    /* ensure we found all the rows */
    ASSERT (debug_rows == 0) ;

    /* === Return the new value of pfree ==================================== */

    return ((Int) (pdest - &A [0])) ;
}


/* ========================================================================== */
/* === clear_mark =========================================================== */
/* ========================================================================== */

/*
    Clears the Row [].shared2.mark array, and returns the new tag_mark.
    Return value is the new tag_mark.  Not user-callable.
*/

PRIVATE Int clear_mark	/* return the new value for tag_mark */
(
    /* === Parameters ======================================================= */

    Int n_row,		/* number of rows in A */
    Colamd_Row Row []	/* Row [0 ... n-1].shared2.mark is set to zero */
)
{
    /* === Local variables ================================================== */

    Int r ;

    for (r = 0 ; r < n_row ; r++)
    {
	if (ROW_IS_ALIVE (r))
	{
	    Row [r].shared2.mark = 0 ;
	}
    }

    /* ------------------ */
    return (1) ;
    /* ------------------ */

}


/* ========================================================================== */
/* === print_report removed for UMFPACK ===================================== */
/* ========================================================================== */



/* ========================================================================== */
/* === colamd debugging routines ============================================ */
/* ========================================================================== */

/* When debugging is disabled, the remainder of this file is ignored. */

#ifndef NDEBUG


/* ========================================================================== */
/* === debug_structures ===================================================== */
/* ========================================================================== */

/*
    At this point, all empty rows and columns are dead.  All live columns
    are "clean" (containing no dead rows) and simplicial (no supercolumns
    yet).  Rows may contain dead columns, but all live rows contain at
    least one live column.
*/

PRIVATE void debug_structures
(
    /* === Parameters ======================================================= */

    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A [],
    Int n_col2
)
{
    /* === Local variables ================================================== */

    Int i ;
    Int c ;
    Int *cp ;
    Int *cp_end ;
    Int len ;
    Int score ;
    Int r ;
    Int *rp ;
    Int *rp_end ;
    Int deg ;

    /* === Check A, Row, and Col ============================================ */

    for (c = 0 ; c < n_col ; c++)
    {
	if (COL_IS_ALIVE (c))
	{
	    len = Col [c].length ;
	    score = Col [c].shared2.score ;
	    DEBUG4 (("initial live col "ID" "ID" "ID"\n", c, len, score)) ;
	    ASSERT (len > 0) ;
	    ASSERT (score >= 0) ;
	    ASSERT (Col [c].shared1.thickness == 1) ;
	    cp = &A [Col [c].start] ;
	    cp_end = cp + len ;
	    while (cp < cp_end)
	    {
		r = *cp++ ;
		ASSERT (ROW_IS_ALIVE (r)) ;
	    }
	}
	else
	{
	    i = Col [c].shared2.order ;
	    ASSERT (i >= n_col2 && i < n_col) ;
	}
    }

    for (r = 0 ; r < n_row ; r++)
    {
	if (ROW_IS_ALIVE (r))
	{
	    i = 0 ;
	    len = Row [r].length ;
	    deg = Row [r].shared1.degree ;
	    ASSERT (len > 0) ;
	    ASSERT (deg > 0) ;
	    rp = &A [Row [r].start] ;
	    rp_end = rp + len ;
	    while (rp < rp_end)
	    {
		c = *rp++ ;
		if (COL_IS_ALIVE (c))
		{
		    i++ ;
		}
	    }
	    ASSERT (i > 0) ;
	}
    }
}


/* ========================================================================== */
/* === debug_deg_lists ====================================================== */
/* ========================================================================== */

/*
    Prints the contents of the degree lists.  Counts the number of columns
    in the degree list and compares it to the total it should have.  Also
    checks the row degrees.
*/

PRIVATE void debug_deg_lists
(
    /* === Parameters ======================================================= */

    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int head [],
    Int min_score,
    Int should,
    Int max_deg
)
{
    /* === Local variables ================================================== */

    Int deg ;
    Int col ;
    Int have ;
    Int row ;

    /* === Check the degree lists =========================================== */

    if (n_col > 10000 && UMF_debug <= 0)
    {
	return ;
    }
    have = 0 ;
    DEBUG4 (("Degree lists: "ID"\n", min_score)) ;
    for (deg = 0 ; deg <= n_col ; deg++)
    {
	col = head [deg] ;
	if (col == EMPTY)
	{
	    continue ;
	}
	DEBUG4 ((ID":", deg)) ;
	while (col != EMPTY)
	{
	    DEBUG4 ((" "ID, col)) ;
	    have += Col [col].shared1.thickness ;
	    ASSERT (COL_IS_ALIVE (col)) ;
	    col = Col [col].shared4.degree_next ;
	}
	DEBUG4 (("\n")) ;
    }
    DEBUG4 (("should "ID" have "ID"\n", should, have)) ;
    ASSERT (should == have) ;

    /* === Check the row degrees ============================================ */

    if (n_row > 10000 && UMF_debug <= 0)
    {
	return ;
    }
    for (row = 0 ; row < n_row ; row++)
    {
	if (ROW_IS_ALIVE (row))
	{
	    ASSERT (Row [row].shared1.degree <= max_deg) ;
	}
    }
}


/* ========================================================================== */
/* === debug_mark =========================================================== */
/* ========================================================================== */

/*
    Ensures that the tag_mark is less that the maximum and also ensures that
    each entry in the mark array is less than the tag mark.
*/

PRIVATE void debug_mark
(
    /* === Parameters ======================================================= */

    Int n_row,
    Colamd_Row Row [],
    Int tag_mark,
    Int max_mark
)
{
    /* === Local variables ================================================== */

    Int r ;

    /* === Check the Row marks ============================================== */

    ASSERT (tag_mark > 0 && tag_mark <= max_mark) ;
    if (n_row > 10000 && UMF_debug <= 0)
    {
	return ;
    }
    for (r = 0 ; r < n_row ; r++)
    {
	ASSERT (Row [r].shared2.mark < tag_mark) ;
    }
}


/* ========================================================================== */
/* === debug_matrix ========================================================= */
/* ========================================================================== */

/*
    Prints out the contents of the columns and the rows.
*/

PRIVATE void debug_matrix
(
    /* === Parameters ======================================================= */

    Int n_row,
    Int n_col,
    Colamd_Row Row [],
    Colamd_Col Col [],
    Int A []
)
{
    /* === Local variables ================================================== */

    Int r ;
    Int c ;
    Int *rp ;
    Int *rp_end ;
    Int *cp ;
    Int *cp_end ;

    /* === Dump the rows and columns of the matrix ========================== */

    if (UMF_debug < 3)
    {
	return ;
    }
    DEBUG3 (("DUMP MATRIX:\n")) ;
    for (r = 0 ; r < n_row ; r++)
    {
	DEBUG3 (("Row "ID" alive? %d\n", r, ROW_IS_ALIVE (r))) ;
	if (ROW_IS_DEAD (r))
	{
	    continue ;
	}

	/* ------------------ */
	/* changed for UMFPACK: */
	DEBUG3 (("start "ID" length "ID" degree "ID" thickness "ID"\n",
		Row [r].start, Row [r].length, Row [r].shared1.degree,
		Row [r].thickness)) ;
	/* ------------------ */

	rp = &A [Row [r].start] ;
	rp_end = rp + Row [r].length ;
	while (rp < rp_end)
	{
	    c = *rp++ ;
	    DEBUG4 (("	%d col "ID"\n", COL_IS_ALIVE (c), c)) ;
	}
    }

    for (c = 0 ; c < n_col ; c++)
    {
	DEBUG3 (("Col "ID" alive? %d\n", c, COL_IS_ALIVE (c))) ;
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
	/* ------------------ */
	/* changed for UMFPACK: */
	DEBUG3 (("start "ID" length "ID" shared1[thickness,parent] "ID
		" shared2 [order,score] "ID"\n", Col [c].start, Col [c].length,
		Col [c].shared1.thickness, Col [c].shared2.score));
	/* ------------------ */
	cp = &A [Col [c].start] ;
	cp_end = cp + Col [c].length ;
	while (cp < cp_end)
	{
	    r = *cp++ ;
	    DEBUG4 (("	%d row "ID"\n", ROW_IS_ALIVE (r), r)) ;
	}

	/* ------------------ */
	/* added for UMFPACK: */
	DEBUG1 (("Col")) ;
	dump_super (c, Col, n_col) ;
	/* ------------------ */

    }
}

/* ------------------ */
/* dump_super added for UMFPACK: */
PRIVATE void dump_super
(
    Int super_c,
    Colamd_Col Col [],
    Int n_col
)
{
    Int col, ncols ;

    DEBUG1 ((" =[ ")) ;
    ncols = 0 ;
    for (col = super_c ; col != EMPTY ; col = Col [col].nextcol)
    {
	DEBUG1 ((" "ID, col)) ;
	ASSERT (col >= 0 && col < n_col) ;
	if (col != super_c)
	{
	    ASSERT (COL_IS_DEAD (col)) ;
	}
	if (Col [col].nextcol == EMPTY)
	{
	    ASSERT (col == Col [super_c].lastcol) ;
	}
	ncols++ ;
	ASSERT (ncols <= Col [super_c].shared1.thickness) ;
    }
    ASSERT (ncols == Col [super_c].shared1.thickness) ;
    DEBUG1 (("]\n")) ;
}
/* ------------------ */


#endif /* NDEBUG */
