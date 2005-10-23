/* ========================================================================== */
/* === UMFPACK_qsymbolic ==================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Performs a symbolic factorization.
    See umfpack_qsymbolic.h and umfpack_symbolic.h for details.

    Dynamic memory usage:  about (3.4nz + 8n + n) integers and n double's as
    workspace (via UMF_malloc, for a square matrix).  All of it is free'd via
    UMF_free if an error occurs.  If successful, the Symbolic object contains
    12 to 14 objects allocated by UMF_malloc, with a total size of no more
    than about 13*n integers.
*/

#include "umf_internal.h"
#include "umf_symbolic_usage.h"
#include "umf_colamd.h"
#include "umf_set_stats.h"
#include "umf_analyze.h"
#include "umf_transpose.h"
#include "umf_is_permutation.h"
#include "umf_malloc.h"
#include "umf_free.h"
#include "umf_2by2.h"
#include "umf_singletons.h"

typedef struct	/* SWType */
{
    Int *Front_npivcol ;    /* size n_col + 1 */
    Int *Front_nrows ;	    /* size n_col */
    Int *Front_ncols ;	    /* size n_col */
    Int *Front_parent ;	    /* size n_col */
    Int *Front_cols ;	    /* size n_col */
    Int *InFront ;	    /* size n_row */
    Int *Ci ;		    /* size Clen */
    Int *Cperm1 ;	    /* size n_col */
    Int *Rperm1 ;	    /* size n_row */
    Int *InvRperm1 ;	    /* size n_row */
    Int *Si ;		    /* size nz */
    Int *Sp ;		    /* size n_col + 1 */
    double *Rs ;	    /* size n_row */
    Int *Rperm_2by2 ;	    /* size n_row */

} SWType ;

PRIVATE void free_work
(
    SWType *SW
) ;

PRIVATE void error
(
    SymbolicType **Symbolic,
    SWType *SW
) ;

/* worst-case usage for SW object */
#define SYM_WORK_USAGE(n_col,n_row,Clen) \
    (DUNITS (Int, Clen) + \
     DUNITS (Int, nz) + \
     4 * DUNITS (Int, n_row) + \
     4 * DUNITS (Int, n_col) + \
     2 * DUNITS (Int, n_col + 1) + \
     DUNITS (double, n_row))

/* required size of Ci for code that calls UMF_transpose and UMF_analyze below*/
#define UMF_ANALYZE_CLEN(nz,n_row,n_col,nn) \
    ((n_col) + MAX ((nz),(n_col)) + 3*(nn)+1 + (n_col))

/* size of an element (in Units), including tuples */
#define ELEMENT_SIZE(r,c) \
    (DGET_ELEMENT_SIZE (r, c) + 1 + (r + c) * UNITS (Tuple, 1))

#ifndef NDEBUG
PRIVATE Int init_count ;
#endif

/* ========================================================================== */
/* === do_amd =============================================================== */
/* ========================================================================== */

PRIVATE void do_amd
(
    Int n,
    const Int Ap [ ],		/* size n+1 */
    const Int Ai [ ],		/* size nz = Ap [n] */
    Int Q [ ],			/* output permutation, j = Q [k] */
    Int Qinv [ ],		/* output inverse permutation, Qinv [j] = k */
    Int Sdeg [ ],		/* degree of A+A', from AMD_aat */
    Int Clen,			/* size of Ci */
    Int Ci [ ],			/* size Ci workspace */
    double amd_Control [ ],	/* AMD control parameters */
    double amd_Info [ ],	/* AMD info */
    SymbolicType *Symbolic,	/* Symbolic object */
    double Info [ ]		/* UMFPACK info */
)
{

    if (n == 0)
    {
	Symbolic->amd_dmax = 0 ;
	Symbolic->amd_lunz = 0 ;
	Info [UMFPACK_SYMMETRIC_LUNZ] = 0 ;
	Info [UMFPACK_SYMMETRIC_FLOPS] = 0 ;
	Info [UMFPACK_SYMMETRIC_DMAX] = 0 ;
	Info [UMFPACK_SYMMETRIC_NDENSE] = 0 ;
    }
    else
    {
	AMD_1 (n, Ap, Ai, Q, Qinv, Sdeg, Clen, Ci, amd_Control, amd_Info) ;

	/* return estimates computed from AMD on PA+PA' */
	Symbolic->amd_dmax = amd_Info [AMD_DMAX] ;
	Symbolic->amd_lunz = 2 * amd_Info [AMD_LNZ] + n ;
	Info [UMFPACK_SYMMETRIC_LUNZ] = Symbolic->amd_lunz ;
	Info [UMFPACK_SYMMETRIC_FLOPS] = DIV_FLOPS * amd_Info [AMD_NDIV] +
	    MULTSUB_FLOPS * amd_Info [AMD_NMULTSUBS_LU] ;
	Info [UMFPACK_SYMMETRIC_DMAX] = Symbolic->amd_dmax ;
	Info [UMFPACK_SYMMETRIC_NDENSE] = amd_Info [AMD_NDENSE] ;
	Info [UMFPACK_SYMBOLIC_DEFRAG] += amd_Info [AMD_NCMPA] ;
    }
}

/* ========================================================================== */
/* === prune_singletons ===================================================== */
/* ========================================================================== */

/* Create the submatrix after removing the n1 singletons.  The matrix has
 * row and column indices in the range 0 to n_row-n1 and 0 to n_col-n1,
 * respectively.  */

PRIVATE Int prune_singletons
(
    Int n1,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    Int Cperm1 [ ],
    Int InvRperm1 [ ],
    Int Si [ ],
    Int Sp [ ]
#ifndef NDEBUG
    , Int Rperm1 [ ]
    , Int n_row
#endif
)
{
    Int row, k, pp, p, oldcol, newcol, newrow, nzdiag, do_nzdiag ;

    nzdiag = 0 ;
    do_nzdiag = (Ax != (double *) NULL)
#ifdef COMPLEX
	&& (Az != (double *) NULL)
#endif
	;

#ifndef NDEBUG
    DEBUGm4 (("Prune : S = A (Cperm1 (n1+1:end), Rperm1 (n1+1:end))\n")) ;
    for (k = 0 ; k < n_row ; k++)
    {
	ASSERT (Rperm1 [k] >= 0 && Rperm1 [k] < n_row) ;
	ASSERT (InvRperm1 [Rperm1 [k]] == k) ;
    }
#endif

    /* create the submatrix after removing singletons */

    pp = 0 ;
    for (k = n1 ; k < n_col ; k++)
    {
	oldcol = Cperm1 [k] ;
	newcol = k - n1 ;
	DEBUG5 (("Prune singletons k "ID" oldcol "ID" newcol "ID": "ID"\n",
	    k, oldcol, newcol, pp)) ;
	Sp [newcol] = pp ;  /* load column pointers */
	for (p = Ap [oldcol] ; p < Ap [oldcol+1] ; p++)
	{
	    row = Ai [p] ;
	    DEBUG5 (("  "ID":  row "ID, pp, row)) ;
	    ASSERT (row >= 0 && row < n_row) ;
	    newrow = InvRperm1 [row] - n1 ;
	    ASSERT (newrow < n_row - n1) ;
	    if (newrow >= 0)
	    {
		DEBUG5 (("  newrow "ID, newrow)) ;
		Si [pp++] = newrow ;
		if (do_nzdiag)
		{
		    /* count the number of truly nonzero entries on the
		     * diagonal of S, excluding entries that are present,
		     * but numerically zero */
		    if (newrow == newcol)
		    {
			/* this is the diagonal entry */
			if (SCALAR_IS_NONZERO (Ax [p])
#ifdef COMPLEX
			 || SCALAR_IS_NONZERO (Az [p])
#endif
			)
			{
			    nzdiag++ ;
			}
		    }
		}
	    }
	    DEBUG5 (("\n")) ;
	}
    }
    Sp [n_col - n1] = pp ;

    ASSERT (AMD_valid (n_row - n1, n_col - n1, Sp, Si)) ;
    return (nzdiag) ;
}

/* ========================================================================== */
/* === combine_ordering ===================================================== */
/* ========================================================================== */

PRIVATE void combine_ordering
(
    Int n1,
    Int nempty_col,
    Int n_col,
    Int Cperm_init [ ],	    /* output permutation */
    Int Cperm1 [ ],	    /* singleton and empty column ordering */
    Int Qinv [ ]	    /* Qinv from AMD or COLAMD */
)
{
    Int k, oldcol, newcol, knew ;

    /* combine the singleton ordering with Qinv */
#ifndef NDEBUG
    for (k = 0 ; k < n_col ; k++)
    {
	Cperm_init [k] = EMPTY ;
    }
#endif
    for (k = 0 ; k < n1 ; k++)
    {
	DEBUG1 ((ID" Initial singleton: "ID"\n", k, Cperm1 [k])) ;
	Cperm_init [k] = Cperm1 [k] ;
    }
    for (k = n1 ; k < n_col - nempty_col ; k++)
    {
	/* this is a non-singleton column */
	oldcol = Cperm1 [k] ;	/* user's name for this column */
	newcol = k - n1 ;	/* Qinv's name for this column */
	knew = Qinv [newcol] ;	/* Qinv's ordering for this column */
	knew += n1 ;		/* shift order, after singletons */
	DEBUG1 ((" k "ID" oldcol "ID" newcol "ID" knew "ID"\n",
	    k, oldcol, newcol, knew)) ;
	ASSERT (knew >= 0 && knew < n_col - nempty_col) ;
	ASSERT (Cperm_init [knew] == EMPTY) ;
	Cperm_init [knew] = oldcol ;
    }
    for (k = n_col - nempty_col ; k < n_col ; k++)
    {
	Cperm_init [k] = Cperm1 [k] ;
    }
#ifndef NDEBUG
    {
	Int *W = (Int *) malloc ((n_col + 1) * sizeof (Int)) ;
	ASSERT (UMF_is_permutation (Cperm_init, W, n_col, n_col)) ;
	free (W) ;
    }
#endif

}

/* ========================================================================== */
/* === UMFPACK_qsymbolic ==================================================== */
/* ========================================================================== */

GLOBAL Int UMFPACK_qsymbolic
(
    Int n_row,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    const Int Quser [ ],
    void **SymbolicHandle,
    const double Control [UMFPACK_CONTROL],
    double User_Info [UMFPACK_INFO]
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int i, nz, j, newj, status, f1, f2, maxnrows, maxncols, nfr, col,
	nchains, maxrows, maxcols, p, nb, nn, *Chain_start, *Chain_maxrows,
	*Chain_maxcols, *Front_npivcol, *Ci, Clen, colamd_stats [COLAMD_STATS],
	fpiv, n_inner, child, parent, *Link, row, *Front_parent,
	analyze_compactions, k, chain, is_sym, *Si, *Sp, n2, do_UMF_analyze,
	fpivcol, fallrows, fallcols, *InFront, *F1, snz, *Front_1strow, f1rows,
	kk, *Cperm_init, *Rperm_init, newrow, *InvRperm1, *Front_leftmostdesc,
	Clen_analyze, strategy, Clen_amd, fixQ, prefer_diagonal, nzdiag, nzaat,
	*Wq, *Sdeg, *Fr_npivcol, nempty, *Fr_nrows, *Fr_ncols, *Fr_parent,
	*Fr_cols, nempty_row, nempty_col, user_auto_strategy, fail, max_rdeg,
	head_usage, tail_usage, lnz, unz, esize, *Esize, rdeg, *Cdeg, *Rdeg,
	*Cperm1, *Rperm1, n1, oldcol, newcol, n1c, n1r, *Rperm_2by2, oldrow,
	dense_row_threshold, tlen, aggressive ;
    double knobs [COLAMD_KNOBS], flops, f, r, c, *Info, force_fixQ,
	Info2 [UMFPACK_INFO], drow, dcol, dtail_usage, dlf, duf, dmax_usage,
	dhead_usage, dlnz, dunz, dmaxfrsize, dClen, dClen_analyze, sym,
	amd_Info [AMD_INFO], dClen_amd, dr, dc, cr, cc, cp,
	amd_Control [AMD_CONTROL], stats [2], tol, scale ;
    SymbolicType *Symbolic ;
    SWType SWspace, *SW ;

#ifndef NDEBUG
    UMF_dump_start ( ) ;
    init_count = UMF_malloc_count ;
    PRINTF ((
"**** Debugging enabled (UMFPACK will be exceedingly slow!) *****************\n"
	)) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get the amount of time used by the process so far */
    /* ---------------------------------------------------------------------- */

    umfpack_tic (stats) ;

    /* ---------------------------------------------------------------------- */
    /* get control settings and check input parameters */
    /* ---------------------------------------------------------------------- */

    drow = GET_CONTROL (UMFPACK_DENSE_ROW, UMFPACK_DEFAULT_DENSE_ROW) ;
    dcol = GET_CONTROL (UMFPACK_DENSE_COL, UMFPACK_DEFAULT_DENSE_COL) ;
    nb = GET_CONTROL (UMFPACK_BLOCK_SIZE, UMFPACK_DEFAULT_BLOCK_SIZE) ;
    strategy = GET_CONTROL (UMFPACK_STRATEGY, UMFPACK_DEFAULT_STRATEGY) ;
    tol = GET_CONTROL (UMFPACK_2BY2_TOLERANCE, UMFPACK_DEFAULT_2BY2_TOLERANCE) ;
    scale = GET_CONTROL (UMFPACK_SCALE, UMFPACK_DEFAULT_SCALE) ;
    force_fixQ = GET_CONTROL (UMFPACK_FIXQ, UMFPACK_DEFAULT_FIXQ) ;
    AMD_defaults (amd_Control) ;
    amd_Control [AMD_DENSE] =
	GET_CONTROL (UMFPACK_AMD_DENSE, UMFPACK_DEFAULT_AMD_DENSE) ;
    aggressive =
	(GET_CONTROL (UMFPACK_AGGRESSIVE, UMFPACK_DEFAULT_AGGRESSIVE) != 0) ;
    amd_Control [AMD_AGGRESSIVE] = aggressive ;

    nb = MAX (2, nb) ;
    nb = MIN (nb, MAXNB) ;
    ASSERT (nb >= 0) ;
    if (nb % 2 == 1) nb++ ;	/* make sure nb is even */
    DEBUG0 (("UMFPACK_qsymbolic: nb = "ID" aggressive = "ID"\n", nb,
	aggressive)) ;

    tol = MAX (0.0, MIN (tol,  1.0)) ;
    if (scale != UMFPACK_SCALE_NONE && scale != UMFPACK_SCALE_MAX)
    {
	scale = UMFPACK_DEFAULT_SCALE ;
    }

    if (User_Info != (double *) NULL)
    {
	/* return Info in user's array */
	Info = User_Info ;
    }
    else
    {
	/* no Info array passed - use local one instead */
	Info = Info2 ;
    }
    /* clear all of Info */
    for (i = 0 ; i < UMFPACK_INFO ; i++)
    {
	Info [i] = EMPTY ;
    }

    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;

    Info [UMFPACK_STATUS] = UMFPACK_OK ;
    Info [UMFPACK_NROW] = n_row ;
    Info [UMFPACK_NCOL] = n_col ;
    Info [UMFPACK_SIZE_OF_UNIT] = (double) (sizeof (Unit)) ;
    Info [UMFPACK_SIZE_OF_INT] = (double) (sizeof (int)) ;
    Info [UMFPACK_SIZE_OF_LONG] = (double) (sizeof (long)) ;
    Info [UMFPACK_SIZE_OF_POINTER] = (double) (sizeof (void *)) ;
    Info [UMFPACK_SIZE_OF_ENTRY] = (double) (sizeof (Entry)) ;
    Info [UMFPACK_SYMBOLIC_DEFRAG] = 0 ;

    if (!Ai || !Ap || !SymbolicHandle)
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_argument_missing ;
	return (UMFPACK_ERROR_argument_missing) ;
    }

    *SymbolicHandle = (void *) NULL ;

    if (n_row <= 0 || n_col <= 0)	/* n_row, n_col must be > 0 */
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_n_nonpositive ;
	return (UMFPACK_ERROR_n_nonpositive) ;
    }

    nz = Ap [n_col] ;
    DEBUG0 (("n_row "ID" n_col "ID" nz "ID"\n", n_row, n_col, nz)) ;
    Info [UMFPACK_NZ] = nz ;
    if (nz < 0)
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_invalid_matrix ;
	return (UMFPACK_ERROR_invalid_matrix) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get the requested strategy */
    /* ---------------------------------------------------------------------- */

    if (n_row != n_col)
    {
	/* if the matrix is rectangular, the only available strategy is
	 *  unsymmetric */
	strategy = UMFPACK_STRATEGY_UNSYMMETRIC ;
	DEBUGm3 (("Rectangular: forcing unsymmetric strategy\n")) ;
    }

    if (strategy < UMFPACK_STRATEGY_AUTO
     || strategy > UMFPACK_STRATEGY_SYMMETRIC)
    {
	/* unrecognized strategy */
	strategy = UMFPACK_STRATEGY_AUTO ;
    }

    user_auto_strategy = (strategy == UMFPACK_STRATEGY_AUTO) ;

    /* ---------------------------------------------------------------------- */
    /* determine amount of memory required for UMFPACK_symbolic */
    /* ---------------------------------------------------------------------- */

    /* The size of Clen required for UMF_colamd is always larger than */
    /* UMF_analyze, but the max is included here in case that changes in */
    /* future versions. */

    /* This is about 2.2*nz + 9*n_col + 6*n_row, or nz/5 + 13*n_col + 6*n_row,
     * whichever is bigger.  For square matrices, it works out to
     * 2.2nz + 15n, or nz/5 + 19n, whichever is bigger (typically 2.2nz+15n). */
    dClen = UMF_COLAMD_RECOMMENDED ((double) nz, (double) n_row,
	(double) n_col) ;

    /* This is defined above, as max (nz,n_col) + 3*nn+1 + 2*n_col, where
     * nn = max (n_row,n_col).  It is always smaller than the space required
     * for colamd or amd. */
    dClen_analyze = UMF_ANALYZE_CLEN ((double) nz, (double) n_row,
	(double) n_col, (double) nn) ;
    dClen = MAX (dClen, dClen_analyze) ;

    /* The space for AMD can be larger than what's required for colamd: */
    dClen_amd = 2.4 * (double) nz + 8 * (double) n_inner ;
    /* additional space for the 2-by-2 strategy */
    dClen_amd += (double) MAX (nn, nz) ;
    dClen = MAX (dClen, dClen_amd) ;

    /* worst case total memory usage for UMFPACK_symbolic (revised below) */
    Info [UMFPACK_SYMBOLIC_PEAK_MEMORY] =
	SYM_WORK_USAGE (n_col, n_row, dClen) +
	UMF_symbolic_usage (n_row, n_col, n_col, n_col, n_col, TRUE) ;

    if (INT_OVERFLOW (dClen * sizeof (Int)))
    {
	/* :: int overflow, Clen too large :: */
	/* Problem is too large for array indexing (Ci [i]) with an Int i. */
	/* Cannot even analyze the problem to determine upper bounds on */
	/* memory usage. Need to use the long integer version, umfpack_*l_*. */
	DEBUGm4 (("out of memory: symbolic int overflow\n")) ;
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

    /* repeat the size calculations, in integers */
    Clen = UMF_COLAMD_RECOMMENDED (nz, n_row, n_col) ;
    Clen_analyze = UMF_ANALYZE_CLEN (nz, n_row, n_col, nn) ;
    Clen = MAX (Clen, Clen_analyze) ;
    Clen_amd = 2.4 * nz + 8 * n_inner ;
    Clen_amd += MAX (nn, nz) ;			/* for Ri, in UMF_2by2 */
    Clen = MAX (Clen, Clen_amd) ;

    /* ---------------------------------------------------------------------- */
    /* allocate the first part of the Symbolic object (header and Cperm_init) */
    /* ---------------------------------------------------------------------- */

    /* (1) Five calls to UMF_malloc are made, for a total space of
     * 2 * (n_row + n_col) + 4 integers + sizeof (SymbolicType).
     * sizeof (SymbolicType) is a small constant.  This space is part of the
     * Symbolic object and is not freed unless an error occurs.  If A is square
     * then this is about 4*n integers.
     */

    Symbolic = (SymbolicType *) UMF_malloc (1, sizeof (SymbolicType)) ;

    if (!Symbolic)
    {
	/* If we fail here, Symbolic is NULL and thus it won't be */
	/* dereferenced by UMFPACK_free_symbolic, as called by error ( ). */
	DEBUGm4 (("out of memory: symbolic object\n")) ;
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	error (&Symbolic, (SWType *) NULL) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

    /* We now know that Symbolic has been allocated */
    Symbolic->valid = 0 ;
    Symbolic->Chain_start = (Int *) NULL ;
    Symbolic->Chain_maxrows = (Int *) NULL ;
    Symbolic->Chain_maxcols = (Int *) NULL ;
    Symbolic->Front_npivcol = (Int *) NULL ;
    Symbolic->Front_parent = (Int *) NULL ;
    Symbolic->Front_1strow = (Int *) NULL ;
    Symbolic->Front_leftmostdesc = (Int *) NULL ;
    Symbolic->Esize = (Int *) NULL ;
    Symbolic->esize = 0 ;

    Symbolic->Cperm_init   = (Int *) UMF_malloc (n_col+1, sizeof (Int)) ;
    Symbolic->Rperm_init   = (Int *) UMF_malloc (n_row+1, sizeof (Int)) ;
    Symbolic->Cdeg	   = (Int *) UMF_malloc (n_col+1, sizeof (Int)) ;
    Symbolic->Rdeg	   = (Int *) UMF_malloc (n_row+1, sizeof (Int)) ;
    Symbolic->Diagonal_map = (Int *) NULL ;

    Cperm_init = Symbolic->Cperm_init ;
    Rperm_init = Symbolic->Rperm_init ;
    Cdeg = Symbolic->Cdeg ;
    Rdeg = Symbolic->Rdeg ;

    if (!Cperm_init || !Rperm_init || !Cdeg || !Rdeg)
    {
	DEBUGm4 (("out of memory: symbolic perm\n")) ;
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	error (&Symbolic, (SWType *) NULL) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

    Symbolic->n_row = n_row ;
    Symbolic->n_col = n_col ;
    Symbolic->nz = nz ;
    Symbolic->nb = nb ;

    /* ---------------------------------------------------------------------- */
    /* check user's input permutation */
    /* ---------------------------------------------------------------------- */

    if (Quser != (Int *) NULL)
    {
	/* use Cperm_init as workspace to check input permutation */
	if (!UMF_is_permutation (Quser, Cperm_init, n_col, n_col))
	{
	    Info [UMFPACK_STATUS] = UMFPACK_ERROR_invalid_permutation ;
	    error (&Symbolic, (SWType *) NULL) ;
	    return (UMFPACK_ERROR_invalid_permutation) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* (2) Eleven calls to UMF_malloc are made, for workspace of size
     * Clen + nz + 7*n_col + 2*n_row + 2 integers.  Clen is the larger of
     *     MAX (2*nz, 4*n_col) + 8*n_col + 6*n_row + n_col + nz/5 and
     *     2.4*nz + 8 * MIN (n_row, n_col) + MAX (n_row, n_col, nz)
     * If A is square and non-singular, then Clen is
     *     MAX (MAX (2*nz, 4*n) + 7*n + nz/5,  3.4*nz) + 8*n
     * If A has at least 4*n nonzeros then Clen is
     *     MAX (2.2*nz + 7*n,  3.4*nz) + 8*n
     * If A has at least (7/1.2)*n nonzeros, (about 5.8*n), then Clen is
     *     3.4*nz + 8*n
     * This space will be free'd when this routine finishes.
     *
     * Total space thus far is about 3.4nz + 12n integers.
     * For the double precision, 32-bit integer version, the user's matrix
     * requires an equivalent space of 3*nz + n integers.  So this space is just
     * slightly larger than the user's input matrix (including the numerical
     * values themselves).
     */

    SW = &SWspace ;	/* used for UMFPACK_symbolic only */

    /* Note that SW->Front_* does not include the dummy placeholder front. */
    /* This space is accounted for by the SYM_WORK_USAGE macro. */

    /* this is free'd early */
    SW->Si	      = (Int *) UMF_malloc (nz, sizeof (Int)) ;
    SW->Sp	      = (Int *) UMF_malloc (n_col + 1, sizeof (Int)) ;
    SW->InvRperm1     = (Int *) UMF_malloc (n_row, sizeof (Int)) ;
    SW->Cperm1	      = (Int *) UMF_malloc (n_col, sizeof (Int)) ;

    /* this is free'd late */
    SW->Ci	      = (Int *) UMF_malloc (Clen, sizeof (Int)) ;
    SW->Front_npivcol = (Int *) UMF_malloc (n_col + 1, sizeof (Int)) ;
    SW->Front_nrows   = (Int *) UMF_malloc (n_col, sizeof (Int)) ;
    SW->Front_ncols   = (Int *) UMF_malloc (n_col, sizeof (Int)) ;
    SW->Front_parent  = (Int *) UMF_malloc (n_col, sizeof (Int)) ;
    SW->Front_cols    = (Int *) UMF_malloc (n_col, sizeof (Int)) ;
    SW->Rperm1	      = (Int *) UMF_malloc (n_row, sizeof (Int)) ;
    SW->InFront	      = (Int *) UMF_malloc (n_row, sizeof (Int)) ;

    /* this is allocated later, and free'd after Cperm1 but before Ci */
    SW->Rperm_2by2    = (Int *) NULL ;   /* will be nn Int's */

    /* this is allocated last, and free'd first */
    SW->Rs	      = (double *) NULL ;	/* will be n_row double's */

    Ci	       = SW->Ci ;
    Fr_npivcol = SW->Front_npivcol ;
    Fr_nrows   = SW->Front_nrows ;
    Fr_ncols   = SW->Front_ncols ;
    Fr_parent  = SW->Front_parent ;
    Fr_cols    = SW->Front_cols ;
    Cperm1     = SW->Cperm1 ;
    Rperm1     = SW->Rperm1 ;
    Si	       = SW->Si ;
    Sp	       = SW->Sp ;
    InvRperm1  = SW->InvRperm1 ;
    Rperm_2by2 = (Int *) NULL ;
    InFront    = SW->InFront ;

    if (!Ci || !Fr_npivcol || !Fr_nrows || !Fr_ncols || !Fr_parent || !Fr_cols
	|| !Cperm1 || !Rperm1 || !Si || !Sp || !InvRperm1 || !InFront)
    {
	DEBUGm4 (("out of memory: symbolic work\n")) ;
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	error (&Symbolic, SW) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

    DEBUG0 (("Symbolic UMF_malloc_count - init_count = "ID"\n",
	UMF_malloc_count - init_count)) ;
    ASSERT (UMF_malloc_count == init_count + 17) ;

    /* ---------------------------------------------------------------------- */
    /* find the row and column singletons */
    /* ---------------------------------------------------------------------- */

    /* [ use first nz + n_row + MAX (n_row, n_col) entries in Ci as workspace,
     * and use Rperm_init as workspace */
    ASSERT (Clen >= nz + n_row + MAX (n_row, n_col)) ;

    status = UMF_singletons (n_row, n_col, Ap, Ai, Quser, Cdeg, Cperm1, Rdeg,
	Rperm1, InvRperm1, &n1, &n1c, &n1r, &nempty_col, &nempty_row, &is_sym,
	&max_rdeg, /* workspace: */ Rperm_init, Ci, Ci + nz, Ci + nz + n_row) ;

    /* ] done using Rperm_init and Ci as workspace */

    /* InvRperm1 is now the inverse of Rperm1 */

    if (status != UMFPACK_OK)
    {
	DEBUGm4 (("matrix invalid: UMF_singletons\n")) ;
	Info [UMFPACK_STATUS] = status ;
	error (&Symbolic, SW) ;
	return (status) ;
    }
    Info [UMFPACK_NEMPTY_COL] = nempty_col ;
    Info [UMFPACK_NEMPTY_ROW] = nempty_row ;
    Info [UMFPACK_NDENSE_COL] = 0 ;	/* # dense rows/cols recomputed below */
    Info [UMFPACK_NDENSE_ROW] = 0 ;
    Info [UMFPACK_COL_SINGLETONS] = n1c ;
    Info [UMFPACK_ROW_SINGLETONS] = n1r ;
    Info [UMFPACK_S_SYMMETRIC] = is_sym ;

    nempty = MIN (nempty_col, nempty_row) ;
    Symbolic->nempty_row = nempty_row ;
    Symbolic->nempty_col = nempty_col ;

    /* UMF_singletons has verified that the user's input matrix is valid */
    ASSERT (AMD_valid (n_row, n_col, Ap, Ai)) ;

    Symbolic->n1 = n1 ;
    Symbolic->nempty = nempty ;
    ASSERT (n1 <= n_inner) ;
    n2 = nn - n1 - nempty ;

    dense_row_threshold =
	UMFPACK_DENSE_DEGREE_THRESHOLD (drow, n_col - n1 - nempty_col) ;
    Symbolic->dense_row_threshold = dense_row_threshold ;

    if (!is_sym)
    {
	/* either the pruned submatrix rectangular, or it is square and
	 * Rperm [n1 .. n-nempty-1] is not the same as Cperm [n1 .. n-nempty-1].
	 * For the auto strategy, switch to the unsymmetric strategy.
	 * Otherwise, if the strategy selected by the user is symmetric or
	 * 2-by-2, then the singletons will be discarded. */
	strategy = UMFPACK_STRATEGY_UNSYMMETRIC ;
	DEBUGm4 (("Strategy: Unsymmetric singletons\n")) ;
    }

    /* ---------------------------------------------------------------------- */
    /* determine symmetry, nzdiag, and degrees of S+S' */
    /* ---------------------------------------------------------------------- */

    /* S is the matrix obtained after removing singletons
     *   = A (Cperm1 [n1..n_col-nempty_col-1], Rperm1 [n1..n_row-nempty_row-1])
     */

    Wq = Rperm_init ;	    /* use Rperm_init as workspace for Wq [ */
    Sdeg = Cperm_init ;	    /* use Cperm_init as workspace for Sdeg [ */
    sym = EMPTY ;
    nzaat = EMPTY ;
    nzdiag = EMPTY ;
    for (i = 0 ; i < AMD_INFO ; i++)
    {
	amd_Info [i] = EMPTY ;
    }

    if (strategy != UMFPACK_STRATEGY_UNSYMMETRIC)
    {
	/* This also determines the degree of each node in S+S' (Sdeg), which
	 * is needed by the 2-by-2 strategy, the symmetry of S, and the number
	 * of nonzeros on the diagonal of S. */
	ASSERT (n_row == n_col) ;
	ASSERT (nempty_row == nempty_col) ;

	/* get the count of nonzeros on the diagonal of S, excluding explicitly
	 * zero entries.  nzdiag = amd_Info [AMD_NZDIAG] counts the zero entries
	 * in S. */

	nzdiag = prune_singletons (n1, nn, Ap, Ai, Ax,
#ifdef COMPLEX
	    Az,
#endif
	    Cperm1, InvRperm1, Si, Sp
#ifndef NDEBUG
	    , Rperm1, nn
#endif
	    ) ;
	nzaat = AMD_aat (n2, Sp, Si, Sdeg, Wq, amd_Info) ;
	sym = amd_Info [AMD_SYMMETRY] ;
	Info [UMFPACK_N2] = n2 ;
	/* nzdiag = amd_Info [AMD_NZDIAG] counts the zero entries of S too */

#ifndef NDEBUG
	for (k = 0 ; k < n2 ; k++)
	{
	    ASSERT (Sdeg [k] >= 0 && Sdeg [k] < n2) ;
	}
	ASSERT (Sp [n2] - n2 <= nzaat && nzaat <= 2 * Sp [n2]) ;
	DEBUG0 (("Explicit zeros: "ID" %g\n", nzdiag, amd_Info [AMD_NZDIAG])) ;
#endif
    }

    /* get statistics from amd_aat, if computed */
    Symbolic->sym = sym ;
    Symbolic->nzaat = nzaat ;
    Symbolic->nzdiag = nzdiag ;
    Symbolic->amd_dmax = EMPTY ;

    Info [UMFPACK_PATTERN_SYMMETRY] = sym ;
    Info [UMFPACK_NZ_A_PLUS_AT] = nzaat ;
    Info [UMFPACK_NZDIAG] = nzdiag ;

    /* ---------------------------------------------------------------------- */
    /* determine the initial strategy based on symmetry and nnz (diag (S)) */
    /* ---------------------------------------------------------------------- */

    if (strategy == UMFPACK_STRATEGY_AUTO)
    {
	if (sym < 0.10)
	{
	    /* highly unsymmetric: use the unsymmetric strategy */
	    strategy = UMFPACK_STRATEGY_UNSYMMETRIC ;
	    DEBUGm4 (("Strategy: select unsymmetric\n")) ;
	}
	else if (sym >= 0.7 && nzdiag == n2)
	{
	    /* mostly symmetric, zero-free diagonal: use symmetric strategy */
	    strategy = UMFPACK_STRATEGY_SYMMETRIC ;
	    DEBUGm4 (("Strategy: select symmetric\n")) ;
	}
	else
	{
	    /* Evaluate the symmetric 2-by-2 strategy, and select it, or
	     * the unsymmetric strategy if the 2-by-2 strategy doesn't look
	     * promising. */
	    strategy = UMFPACK_STRATEGY_2BY2 ;
	    DEBUGm4 (("Strategy: try 2-by-2\n")) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* try the 2-by-2 strategy */
    /* ---------------------------------------------------------------------- */

    /* (3) If the 2-by-2 strategy is attempted, additional workspace of size
     * nn integers and nn double's is allocated, where nn = n_row = n_col.
     * The real workspace is immediately free'd.  The integer workspace of
     * size nn remains until the end of umfpack_qsymbolic. */

    /* If the resulting matrix S (Rperm_2by2, :) is too unsymmetric, then the
     * unsymmetric strategy will be used instead. */

    if (strategy == UMFPACK_STRATEGY_2BY2)
    {
	Int *Rp, *Ri, *Blen, *W, nz_papat, nzd2, nweak, unmatched,
	    Clen3 ;
	double sym2 ;

	/* ------------------------------------------------------------------ */
	/* get workspace for UMF_2by2 */
	/* ------------------------------------------------------------------ */

	ASSERT (n_row == n_col && nn == n_row) ;

#ifndef NDEBUG
	for (k = 0 ; k < n2 ; k++)
	{
	    ASSERT (Sdeg [k] >= 0 && Sdeg [k] < n2) ;
	}
#endif

	/* allocate Rperm_2by2 */
	SW->Rperm_2by2 = (Int *) UMF_malloc (nn, sizeof (Int)) ;
	Rperm_2by2 = SW->Rperm_2by2 ;
	if (Rperm_2by2 == (Int *) NULL)
	{
	    DEBUGm4 (("out of memory: Rperm_2by2\n")) ;
	    Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	    error (&Symbolic, SW) ;
	    return (UMFPACK_ERROR_out_of_memory) ;
	}

	/* allocate Ri from the tail end of Ci [ */
	Clen3 = Clen - (MAX (nn, nz) + 1) ;
	Ri = Ci + Clen3 ;
	ASSERT (Clen3 >= nz) ;	/* space required for UMF_2by2 */

	/* use Fr_* as workspace for Rp, Blen, and W [ */
	Rp = Fr_npivcol ;
	Blen = Fr_ncols ;
	W = Fr_cols ;

	if (scale != UMFPACK_SCALE_NONE)
	{
	    SW->Rs = (double *) UMF_malloc (nn, sizeof (double)) ;
	    if (SW->Rs == (double *) NULL)
	    {
		DEBUGm4 (("out of memory: scale factors for 2-by-2\n")) ;
		Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
		error (&Symbolic, SW) ;
		return (UMFPACK_ERROR_out_of_memory) ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* find the 2-by-2 row permutation */
	/* ------------------------------------------------------------------ */

	/* find a row permutation Rperm_2by2 such that S (Rperm_2by2, :)
	 * has a healthy diagonal */

	UMF_2by2 (nn, Ap, Ai, Ax,
#ifdef COMPLEX
		Az,
#endif
		tol, scale, Cperm1,
#ifndef NDEBUG
		Rperm1,
#endif
		InvRperm1, n1, nempty, Sdeg, Rperm_2by2, &nweak, &unmatched,
		Ri, Rp, SW->Rs, Blen, W, Ci, Wq) ;
	DEBUGm3 (("2by2: nweak "ID" unmatched "ID"\n", nweak, unmatched)) ;
	Info [UMFPACK_2BY2_NWEAK] = nweak ;
	Info [UMFPACK_2BY2_UNMATCHED] = unmatched ;

	SW->Rs = (double *) UMF_free ((void *) SW->Rs) ;

	/* R = S (Rperm_2by2,:)' */
	(void) UMF_transpose (n2, n2, Sp, Si, (double *) NULL, Rperm_2by2,
	    (Int *) NULL, 0, Rp, Ri, (double *) NULL, W, FALSE
#ifdef COMPLEX
	    , (double *) NULL, (double *) NULL, FALSE
#endif
	    ) ;
	ASSERT (AMD_valid (n2, n2, Rp, Ri)) ;

	/* contents of Si and Sp no longer needed, but the space is
	 * still needed */

	/* ------------------------------------------------------------------ */
	/* find symmetry of S (Rperm_2by2, :)', and prepare to order with AMD */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < AMD_INFO ; i++)
	{
	    amd_Info [i] = EMPTY ;
	}
	nz_papat = AMD_aat (n2, Rp, Ri, Sdeg, Wq, amd_Info) ;
	sym2 = amd_Info [AMD_SYMMETRY] ;
	nzd2 = amd_Info [AMD_NZDIAG] ;

	Info [UMFPACK_2BY2_PATTERN_SYMMETRY] = sym2 ;
	Info [UMFPACK_2BY2_NZ_PA_PLUS_PAT] = nz_papat ;
	Info [UMFPACK_2BY2_NZDIAG] = nzd2 ;

	DEBUG0 (("2by2: sym2 %g nzd2 "ID" n2 "ID"\n", sym2, nzd2, n2)) ;

	/* ------------------------------------------------------------------ */
	/* evaluate the 2-by-2 results */
	/* ------------------------------------------------------------------ */

	if (user_auto_strategy)
	{
	    if ((sym2 > 1.1 * sym) && (nzd2 > 0.9 * n2))
	    {
		/* 2-by-2 made it much more symmetric */
		DEBUGm4 (("eval Strategy 2by2: much more symmetric:  2by2\n")) ;
		strategy = UMFPACK_STRATEGY_2BY2 ;
	    }
	    else if (sym2 < 0.7 * sym)
	    {
		/* 2-by-2 made it much more unsymmetric */
		DEBUGm4 (("eval Strategy 2by2: much more UNsymmetric:unsym\n"));
		strategy = UMFPACK_STRATEGY_UNSYMMETRIC ;
	    }
	    else if (sym2 < 0.25)
	    {
		DEBUGm4 (("eval Strategy 2by2: is UNsymmetric: unsym\n"));
		strategy = UMFPACK_STRATEGY_UNSYMMETRIC ;
	    }
	    else if (sym2 >= 0.51)
	    {
		DEBUGm4 (("eval Strategy 2by2: sym2 >= 0.51: 2by2\n")) ;
		strategy = UMFPACK_STRATEGY_2BY2 ;
	    }
	    else if (sym2 >= 0.999 * sym)
	    {
		/* 2-by-2 improved symmetry, or made it only slightly worse */
		DEBUGm4 (("eval Strategy 2by2: sym2 >= 0.999 sym: 2by2\n")) ;
		strategy = UMFPACK_STRATEGY_2BY2 ;
	    }
	    else
	    {
		/* can't decide what to do, so pick the unsymmetric strategy */
		DEBUGm4 (("eval Strategy 2by2: punt: unsym\n"));
		strategy = UMFPACK_STRATEGY_UNSYMMETRIC ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* if the 2-by-2 strategy is selected: */
	/* ------------------------------------------------------------------ */

	if (strategy == UMFPACK_STRATEGY_2BY2)
	{
	    if (Quser == (Int *) NULL)
	    {
		/* 2-by-2 strategy is successful */
		/* compute amd (S) */
		Int *Qinv = Fr_npivcol ;
		ASSERT (Clen3 >= (nz_papat + nz_papat/5 + nn) + 7*nn) ;
		do_amd (n2, Rp, Ri, Wq, Qinv, Sdeg, Clen3, Ci,
		    amd_Control, amd_Info, Symbolic, Info) ;
		/* combine the singleton ordering and the AMD ordering */
		combine_ordering (n1, nempty, nn, Cperm_init, Cperm1, Qinv) ;
	    }
	    /* fix Rperm_2by2 to reflect A, not S */
	    for (k = 0 ; k < n1 ; k++)
	    {
		oldcol = Cperm1 [k] ;
		i = k ;
		oldrow = Rperm1 [k] ;
		W [oldcol] = oldrow ;
	    }
	    for (k = n1 ; k < nn - nempty ; k++)
	    {
		oldcol = Cperm1 [k] ;
		i = Rperm_2by2 [k - n1] + n1 ;
		oldrow = Rperm1 [i] ;
		W [oldcol] = oldrow ;
	    }
	    for (k = nn - nempty ; k < nn ; k++)
	    {
		oldcol = Cperm1 [k] ;
		i = k ;
		oldrow = Rperm1 [k] ;
		W [oldcol] = oldrow ;
	    }
	    for (k = 0 ; k < nn ; k++)
	    {
		Rperm_2by2 [k] = W [k] ;
	    }

	    /* Now, the "diagonal" entry in oldcol (where oldcol is the user's
	     * name for a column, is the entry in row oldrow (where oldrow is
	     * the user's name for a row, and oldrow = Rperm_2by2 [oldcol] */
	}

	/* Fr_* no longer needed for Rp, Blen, W ] */
    }

    /* ---------------------------------------------------------------------- */
    /* finalize the strategy, including fixQ and prefer_diagonal */
    /* ---------------------------------------------------------------------- */

    if (strategy == UMFPACK_STRATEGY_SYMMETRIC)
    {
	/* use given Quser or AMD (A+A'), fix Q during factorization,
	 * prefer diagonal */
	DEBUG0 (("\nStrategy: symmetric\n")) ;
	ASSERT (n_row == n_col) ;
	Symbolic->ordering = UMFPACK_ORDERING_AMD ;
	fixQ = TRUE ;
	prefer_diagonal = TRUE ;
    }
    else if (strategy == UMFPACK_STRATEGY_2BY2)
    {
	/* use Q = given Quser or Q = AMD (PA+PA'), fix Q during factorization,
	 * prefer diagonal, and factorize PAQ, where P is found by UMF_2by2. */
	DEBUG0 (("\nStrategy: symmetric 2-by-2\n")) ;
	ASSERT (n_row == n_col) ;
	Symbolic->ordering = UMFPACK_ORDERING_AMD ;
	fixQ = TRUE ;
	prefer_diagonal = TRUE ;
    }
    else
    {
	/* use given Quser or COLAMD (A), refine Q during factorization,
	 * no diagonal preference */
	ASSERT (strategy == UMFPACK_STRATEGY_UNSYMMETRIC) ;
	DEBUG0 (("\nStrategy: unsymmetric\n")) ;
	Symbolic->ordering = UMFPACK_ORDERING_COLAMD ;
	fixQ = FALSE ;
	prefer_diagonal = FALSE ;
    }

    if (Quser != (Int *) NULL)
    {
	Symbolic->ordering = UMFPACK_ORDERING_GIVEN ;
    }

    if (force_fixQ > 0)
    {
	fixQ = TRUE ;
	DEBUG0 (("Force fixQ true\n")) ;
    }
    else if (force_fixQ < 0)
    {
	fixQ = FALSE ;
	DEBUG0 (("Force fixQ false\n")) ;
    }

    DEBUG0 (("Strategy: ordering:   "ID"\n", Symbolic->ordering)) ;
    DEBUG0 (("Strategy: fixQ:       "ID"\n", fixQ)) ;
    DEBUG0 (("Strategy: prefer diag "ID"\n", prefer_diagonal)) ;

    /* get statistics from amd_aat, if computed */
    Symbolic->strategy = strategy ;
    Symbolic->fixQ = fixQ ;
    Symbolic->prefer_diagonal = prefer_diagonal ;

    Info [UMFPACK_STRATEGY_USED] = strategy ;
    Info [UMFPACK_ORDERING_USED] = Symbolic->ordering ;
    Info [UMFPACK_QFIXED] = fixQ ;
    Info [UMFPACK_DIAG_PREFERRED] = prefer_diagonal ;

    /* ---------------------------------------------------------------------- */
    /* get the AMD ordering for the symmetric strategy */
    /* ---------------------------------------------------------------------- */

    if (strategy == UMFPACK_STRATEGY_SYMMETRIC && Quser == (Int *) NULL)
    {
	/* symmetric strategy for a matrix with mostly symmetric pattern */
	Int *Qinv = Fr_npivcol ;
	ASSERT (n_row == n_col && nn == n_row) ;
	ASSERT (Clen >= (nzaat + nzaat/5 + nn) + 7*nn) ;
	do_amd (n2, Sp, Si, Wq, Qinv, Sdeg, Clen, Ci,
		amd_Control, amd_Info, Symbolic, Info) ;
	/* combine the singleton ordering and the AMD ordering */
	combine_ordering (n1, nempty, nn, Cperm_init, Cperm1, Qinv) ;
    }
    /* Sdeg no longer needed ] */
    /* done using Rperm_init as workspace for Wq ] */

    /* Contents of Si and Sp no longer needed, but the space is still needed */

    /* ---------------------------------------------------------------------- */
    /* use the user's input column ordering (already in Cperm1) */
    /* ---------------------------------------------------------------------- */

    if (Quser != (Int *) NULL)
    {
	for (k = 0 ; k < n_col ; k++)
	{
	    Cperm_init [k] = Cperm1 [k] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* use COLAMD to order the matrix */
    /* ---------------------------------------------------------------------- */

    if (strategy == UMFPACK_STRATEGY_UNSYMMETRIC && Quser == (Int *) NULL)
    {

	/* ------------------------------------------------------------------ */
	/* copy the matrix into colamd workspace (colamd destroys its input) */
	/* ------------------------------------------------------------------ */

	/* C = A (Cperm1 (n1+1:end), Rperm1 (n1+1:end)), where Ci is used as
	 * the row indices and Cperm_init (on input) is used as the column
	 * pointers. */

	(void) prune_singletons (n1, n_col, Ap, Ai,
	    (double *) NULL,
#ifdef COMPLEX
	    (double *) NULL,
#endif
	    Cperm1, InvRperm1, Ci, Cperm_init
#ifndef NDEBUG
	    , Rperm1, n_row
#endif
	    ) ;

	/* ------------------------------------------------------------------ */
	/* set UMF_colamd defaults */
	/* ------------------------------------------------------------------ */

	UMF_colamd_set_defaults (knobs) ;
	knobs [COLAMD_DENSE_ROW] = drow ;
	knobs [COLAMD_DENSE_COL] = dcol ;
	knobs [COLAMD_AGGRESSIVE] = aggressive ;

	/* ------------------------------------------------------------------ */
	/* check input matrix and find the initial column pre-ordering */
	/* ------------------------------------------------------------------ */

	/* NOTE: umf_colamd is not given any original empty rows or columns.
	 * Those have already been removed via prune_singletons, above.  The
	 * umf_colamd routine has been modified to assume that all rows and
	 * columns have at least one entry in them.  It will break if it is
	 * given empty rows or columns (an assertion is triggered when running
	 * in debug mode. */

	(void) UMF_colamd (
		n_row - n1 - nempty_row,
		n_col - n1 - nempty_col,
		Clen, Ci, Cperm_init, knobs, colamd_stats,
		Fr_npivcol, Fr_nrows, Fr_ncols, Fr_parent, Fr_cols, &nfr,
		InFront) ;
	ASSERT (colamd_stats [COLAMD_EMPTY_ROW] == 0) ;
	ASSERT (colamd_stats [COLAMD_EMPTY_COL] == 0) ;

	/* # of dense rows will be recomputed below */
	Info [UMFPACK_NDENSE_ROW]  = colamd_stats [COLAMD_DENSE_ROW] ;
	Info [UMFPACK_NDENSE_COL]  = colamd_stats [COLAMD_DENSE_COL] ;
	Info [UMFPACK_SYMBOLIC_DEFRAG] = colamd_stats [COLAMD_DEFRAG_COUNT] ;

	/* re-analyze if any "dense" rows or cols ignored by UMF_colamd */
	do_UMF_analyze =
	    colamd_stats [COLAMD_DENSE_ROW] > 0 ||
	    colamd_stats [COLAMD_DENSE_COL] > 0 ;

	/* Combine the singleton and colamd ordering into Cperm_init */
	/* Note that colamd returns its inverse permutation in Ci */
	combine_ordering (n1, nempty_col, n_col, Cperm_init, Cperm1, Ci) ;

	/* contents of Ci no longer needed */

#ifndef NDEBUG
	for (col = 0 ; col < n_col ; col++)
	{
	    DEBUG1 (("Cperm_init ["ID"] = "ID"\n", col, Cperm_init[col]));
	}
	/* make sure colamd returned a valid permutation */
	ASSERT (Cperm_init != (Int *) NULL) ;
	ASSERT (UMF_is_permutation (Cperm_init, Ci, n_col, n_col)) ;
#endif

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* do not call colamd - use input Quser or AMD instead */
	/* ------------------------------------------------------------------ */

	/* The ordering (Quser or Qamd) is already in Cperm_init */
	do_UMF_analyze = TRUE ;

    }

    Cperm_init [n_col] = EMPTY ;	/* unused in Cperm_init */

    /* ---------------------------------------------------------------------- */
    /* AMD ordering, if it exists, has been copied into Cperm_init */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    DEBUG3 (("Cperm_init column permutation:\n")) ;
    ASSERT (UMF_is_permutation (Cperm_init, Ci, n_col, n_col)) ;
    for (k = 0 ; k < n_col ; k++)
    {
	DEBUG3 ((ID"\n", Cperm_init [k])) ;
    }
    /* ensure that empty columns have been placed last in A (:,Cperm_init) */
    for (newj = 0 ; newj < n_col ; newj++)
    {
	/* empty columns will be last in A (:, Cperm_init (1:n_col)) */
	j = Cperm_init [newj] ;
	ASSERT (IMPLIES (newj >= n_col-nempty_col, Cdeg [j] == 0)) ;
	ASSERT (IMPLIES (newj <  n_col-nempty_col, Cdeg [j] > 0)) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* symbolic factorization (unless colamd has already done it) */
    /* ---------------------------------------------------------------------- */

    if (do_UMF_analyze)
    {

	Int *W, *Bp, *Bi, *Cperm2, ok, *P, Clen2, bsize, Clen0 ;

	/* ------------------------------------------------------------------ */
	/* construct column pre-ordered, pruned submatrix */
	/* ------------------------------------------------------------------ */

	/* S = column form submatrix after removing singletons and applying
	 * initial column ordering (includes singleton ordering) */
	(void) prune_singletons (n1, n_col, Ap, Ai,
	    (double *) NULL,
#ifdef COMPLEX
	    (double *) NULL,
#endif
	    Cperm_init, InvRperm1, Si, Sp
#ifndef NDEBUG
	    , Rperm1, n_row
#endif
	    ) ;

	/* ------------------------------------------------------------------ */
	/* Ci [0 .. Clen-1] holds the following work arrays:

		first Clen0 entries	empty space, where Clen0 =
					Clen - (nn+1 + 2*nn + n_col)
					and Clen0 >= nz + n_col
		next nn+1 entries	Bp [0..nn]
		next nn entries		Link [0..nn-1]
		next nn entries		W [0..nn-1]
		last n_col entries	Cperm2 [0..n_col-1]

	    We have Clen >= n_col + MAX (nz,n_col) + 3*nn+1 + n_col,
	    So  Clen0 >= 2*n_col as required for AMD_postorder
	    and Clen0 >= n_col + nz as required
	*/

	Clen0 = Clen - (nn+1 + 2*nn + n_col) ;
	Bp = Ci + Clen0 ;
	Link = Bp + (nn+1) ;
	W = Link + nn ;
	Cperm2 = W + nn ;
	ASSERT (Cperm2 + n_col == Ci + Clen) ;
	ASSERT (Clen0 >= nz + n_col) ;
	ASSERT (Clen0 >= 2*n_col) ;

	/* ------------------------------------------------------------------ */
	/* P = order that rows will be used in UMF_analyze */
	/* ------------------------------------------------------------------ */

	/* use W to mark rows, and use Link for row permutation P [ [ */
	for (row = 0 ; row < n_row - n1 ; row++)
	{
	    W [row] = FALSE ;
	}
	P = Link ;

	k = 0 ;

	for (col = 0 ; col < n_col - n1 ; col++)
	{
	    /* empty columns are last in S */
	    for (p = Sp [col] ; p < Sp [col+1] ; p++)
	    {
		row = Si [p] ;
		if (!W [row])
		{
		    /* this row has just been seen for the first time */
		    W [row] = TRUE ;
		    P [k++] = row ;
		}
	    }
	}

	/* If the matrix has truly empty rows, then P will not be */
	/* complete, and visa versa.  The matrix is structurally singular. */
	nempty_row = n_row - n1 - k ;
	if (k < n_row - n1)
	{
	    /* complete P by putting empty rows last in their natural order, */
	    /* rather than declaring an error (the matrix is singular) */
	    for (row = 0 ; row < n_row - n1 ; row++)
	    {
		if (!W [row])
		{
		    /* W [row] = TRUE ;  (not required) */
		    P [k++] = row ;
		}
	    }
	}

	/* contents of W no longer needed ] */

#ifndef NDEBUG
	DEBUG3 (("Induced row permutation:\n")) ;
	ASSERT (k == n_row - n1) ;
	ASSERT (UMF_is_permutation (P, W, n_row - n1, n_row - n1)) ;
	for (k = 0 ; k < n_row - n1 ; k++)
	{
	    DEBUG3 ((ID"\n", P [k])) ;
	}
#endif

	/* ------------------------------------------------------------------ */
	/* B = row-form of the pattern of S (excluding empty columns) */
	/* ------------------------------------------------------------------ */

	/* Ci [0 .. Clen-1] holds the following work arrays:

		first Clen2 entries	empty space, must be at least >= n_col
		next max (nz,1)		Bi [0..max (nz,1)-1]
		next nn+1 entries	Bp [0..nn]
		next nn entries		Link [0..nn-1]
		next nn entries		W [0..nn-1]
		last n_col entries	Cperm2 [0..n_col-1]

		This memory usage is accounted for by the UMF_ANALYZE_CLEN
		macro.
	*/

	Clen2 = Clen0 ;
	snz = Sp [n_col - n1] ;
	bsize = MAX (snz, 1) ;
	Clen2 -= bsize ;
	Bi = Ci + Clen2 ;
	ASSERT (Clen2 >= n_col) ;

	(void) UMF_transpose (n_row - n1, n_col - n1 - nempty_col,
	    Sp, Si, (double *) NULL,
	    P, (Int *) NULL, 0, Bp, Bi, (double *) NULL, W, FALSE
#ifdef COMPLEX
	    , (double *) NULL, (double *) NULL, FALSE
#endif
	    ) ;

	/* contents of Si and Sp no longer needed */

	/* contents of P (same as Link) and W not needed */
	/* still need Link and W as work arrays, though ] */

	ASSERT (Bp [0] == 0) ;
	ASSERT (Bp [n_row - n1] == snz) ;

	/* increment Bp to point into Ci, not Bi */
	for (i = 0 ; i <= n_row - n1 ; i++)
	{
	    Bp [i] += Clen2 ;
	}
	ASSERT (Bp [0] == Clen0 - bsize) ;
	ASSERT (Bp [n_row - n1] <= Clen0) ;

	/* Ci [0 .. Clen-1] holds the following work arrays:

		first Clen0 entries	Ci [0 .. Clen0-1], where the col indices
					of B are at the tail end of this part,
					and Bp [0] = Clen2 >= n_col.  Note
					that Clen0 = Clen2 + max (snz,1).
		next nn+1 entries	Bp [0..nn]
		next nn entries		Link [0..nn-1]
		next nn entries		W [0..nn-1]
		last n_col entries	Cperm2 [0..n_col-1]
	*/

	/* ------------------------------------------------------------------ */
	/* analyze */
	/* ------------------------------------------------------------------ */

	/* only analyze the non-empty, non-singleton part of the matrix */
	ok = UMF_analyze (
		n_row - n1 - nempty_row,
		n_col - n1 - nempty_col,
		Ci, Bp, Cperm2, fixQ, W, Link,
		Fr_ncols, Fr_nrows, Fr_npivcol,
		Fr_parent, &nfr, &analyze_compactions) ;
	if (!ok)
	{
	    /* :: internal error in umf_analyze :: */
	    Info [UMFPACK_STATUS] = UMFPACK_ERROR_internal_error ;
	    error (&Symbolic, SW) ;
	    return (UMFPACK_ERROR_internal_error) ;
	}
	Info [UMFPACK_SYMBOLIC_DEFRAG] += analyze_compactions ;

	/* ------------------------------------------------------------------ */
	/* combine the input permutation and UMF_analyze's permutation */
	/* ------------------------------------------------------------------ */

	if (!fixQ)
	{
	    /* Cperm2 is the column etree post-ordering */
	    ASSERT (UMF_is_permutation (Cperm2, W,
	    n_col-n1-nempty_col, n_col-n1-nempty_col)) ;

	    /* Note that the empty columns remain at the end of Cperm_init */
	    for (k = 0 ; k < n_col - n1 - nempty_col ; k++)
	    {
		W [k] = Cperm_init [n1 + Cperm2 [k]] ;
	    }

	    for (k = 0 ; k < n_col - n1 - nempty_col ; k++)
	    {
		Cperm_init [n1 + k] = W [k] ;
	    }
	}

	ASSERT (UMF_is_permutation (Cperm_init, W, n_col, n_col)) ;

    }

    /* ---------------------------------------------------------------------- */
    /* free some of the workspace */
    /* ---------------------------------------------------------------------- */

    /* (4) The real workspace, Rs, of size n_row doubles has already been
     * free'd.  An additional workspace of size nz + n_col+1 + n_col integers
     * is now free'd as well. */

    SW->Si = (Int *) UMF_free ((void *) SW->Si) ;
    SW->Sp = (Int *) UMF_free ((void *) SW->Sp) ;
    SW->Cperm1 = (Int *) UMF_free ((void *) SW->Cperm1) ;
    ASSERT (SW->Rs == (double *) NULL) ;

    /* ---------------------------------------------------------------------- */
    /* determine the size of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* determine the size of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    nchains = 0 ;
    for (i = 0 ; i < nfr ; i++)
    {
	if (Fr_parent [i] != i+1)
	{
	    nchains++ ;
	}
    }

    Symbolic->nchains = nchains ;
    Symbolic->nfr = nfr ;
    Symbolic->esize
	= (max_rdeg > dense_row_threshold) ? (n_col - n1 - nempty_col) : 0 ;

    /* true size of Symbolic object */
    Info [UMFPACK_SYMBOLIC_SIZE] = UMF_symbolic_usage (n_row, n_col, nchains,
	    nfr, Symbolic->esize, prefer_diagonal) ;

    /* actual peak memory usage for UMFPACK_symbolic (actual nfr, nchains) */
    Info [UMFPACK_SYMBOLIC_PEAK_MEMORY] =
	SYM_WORK_USAGE (n_col, n_row, Clen) + Info [UMFPACK_SYMBOLIC_SIZE] ;
    Symbolic->peak_sym_usage = Info [UMFPACK_SYMBOLIC_PEAK_MEMORY] ;

    DEBUG0 (("Number of fronts: "ID"\n", nfr)) ;

    /* ---------------------------------------------------------------------- */
    /* allocate the second part of the Symbolic object (Front_*, Chain_*) */
    /* ---------------------------------------------------------------------- */

    /* (5) UMF_malloc is called 7 or 8 times, for a total space of
     * (4*(nfr+1) + 3*(nchains+1) + esize) integers, where nfr is the total
     * number of frontal matrices and nchains is the total number of frontal
     * matrix chains, and where nchains <= nfr <= n_col.  esize is zero if there
     * are no dense rows, or n_col-n1-nempty_col otherwise (n1 is the number of
     * singletons and nempty_col is the number of empty columns).  This space is
     * part of the Symbolic object and is not free'd unless an error occurs.
     * This is between 7 and about 8n integers when A is square.
     */

    /* Note that Symbolic->Front_* does include the dummy placeholder front */
    Symbolic->Front_npivcol = (Int *) UMF_malloc (nfr+1, sizeof (Int)) ;
    Symbolic->Front_parent = (Int *) UMF_malloc (nfr+1, sizeof (Int)) ;
    Symbolic->Front_1strow = (Int *) UMF_malloc (nfr+1, sizeof (Int)) ;
    Symbolic->Front_leftmostdesc = (Int *) UMF_malloc (nfr+1, sizeof (Int)) ;
    Symbolic->Chain_start = (Int *) UMF_malloc (nchains+1, sizeof (Int)) ;
    Symbolic->Chain_maxrows = (Int *) UMF_malloc (nchains+1, sizeof (Int)) ;
    Symbolic->Chain_maxcols = (Int *) UMF_malloc (nchains+1, sizeof (Int)) ;

    fail = (!Symbolic->Front_npivcol || !Symbolic->Front_parent ||
	!Symbolic->Front_1strow || !Symbolic->Front_leftmostdesc ||
	!Symbolic->Chain_start || !Symbolic->Chain_maxrows ||
	!Symbolic->Chain_maxcols) ;

    if (Symbolic->esize > 0)
    {
	Symbolic->Esize = (Int *) UMF_malloc (Symbolic->esize, sizeof (Int)) ;
	fail = fail || !Symbolic->Esize ;
    }

    if (fail)
    {
	DEBUGm4 (("out of memory: rest of symbolic object\n")) ;
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	error (&Symbolic, SW) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }
    DEBUG0 (("Symbolic UMF_malloc_count - init_count = "ID"\n",
	UMF_malloc_count - init_count)) ;
    ASSERT (UMF_malloc_count == init_count + 21
	+ (SW->Rperm_2by2 != (Int *) NULL)
	+ (Symbolic->Esize != (Int *) NULL)) ;

    Front_npivcol = Symbolic->Front_npivcol ;
    Front_parent = Symbolic->Front_parent ;
    Front_1strow = Symbolic->Front_1strow ;
    Front_leftmostdesc = Symbolic->Front_leftmostdesc ;

    Chain_start = Symbolic->Chain_start ;
    Chain_maxrows = Symbolic->Chain_maxrows ;
    Chain_maxcols = Symbolic->Chain_maxcols ;

    Esize = Symbolic->Esize ;

    /* ---------------------------------------------------------------------- */
    /* assign rows to fronts */
    /* ---------------------------------------------------------------------- */

    /* find InFront, unless colamd has already computed it */
    if (do_UMF_analyze)
    {

	DEBUGm4 ((">>>>>>>>>Computing Front_1strow from scratch\n")) ;
	/* empty rows go to dummy front nfr */
	for (row = 0 ; row < n_row ; row++)
	{
	    InFront [row] = nfr ;	
	}
	/* assign the singleton pivot rows to the "empty" front */
	for (k = 0 ; k < n1 ; k++)
	{
	    row = Rperm1 [k] ;
	    InFront [row] = EMPTY ;
	}
	DEBUG1 (("Front (EMPTY), singleton nrows "ID" ncols "ID"\n", k, k)) ;
	newj = n1 ;
	for (i = 0 ; i < nfr ; i++)
	{
	    fpivcol = Fr_npivcol [i] ;
	    f1rows = 0 ;
	    /* for all pivot columns in front i */
	    for (kk = 0 ; kk < fpivcol ; kk++, newj++)
	    {
		j = Cperm_init [newj] ;
		ASSERT (IMPLIES (newj >= n_col-nempty_col,
				Ap [j+1] - Ap [j] == 0));
		for (p = Ap [j] ; p < Ap [j+1] ; p++)
		{
		    row = Ai [p] ;
		    if (InFront [row] == nfr)
		    {
			/* this row belongs to front i */
			DEBUG1 (("    Row "ID" in Front "ID"\n", row, i)) ;
			InFront [row] = i ;
			f1rows++ ;
		    }
		}
	    }
	    Front_1strow [i] = f1rows ;
	    DEBUG1 (("    Front "ID" has 1strows: "ID" pivcols "ID"\n",
		i, f1rows, fpivcol)) ;
	}

    }
    else
    {

	/* COLAMD has already computed InFront, but it is not yet
	 * InFront [row] = front i, where row is an original row.  It is
	 * InFront [k-n1] = i for k in the range n1 to n_row-nempty_row,
	 * and where row = Rperm1 [k].  Need to permute InFront.  Also compute
	 * # of original rows assembled into each front.
	 * [ use Ci as workspace */
	DEBUGm4 ((">>>>>>>>>Computing Front_1strow from colamd's InFront\n")) ;
	for (i = 0 ; i <= nfr ; i++)
	{
	    Front_1strow [i] = 0 ;
	}
	/* assign the singleton pivot rows to "empty" front */
	for (k = 0 ; k < n1 ; k++)
	{
	    row = Rperm1 [k] ;
	    Ci [row] = EMPTY ;
	}
	/* assign the non-empty rows to the front that assembled them */
	for ( ; k < n_row - nempty_row ; k++)
	{
	    row = Rperm1 [k] ;
	    i = InFront [k - n1] ;
	    ASSERT (i >= EMPTY && i < nfr) ;
	    if (i != EMPTY)
	    {
		Front_1strow [i]++ ;
	    }
	    /* use Ci as permuted version of InFront */
	    Ci [row] = i ;
	}
	/* empty rows go to the "dummy" front */
	for ( ; k < n_row ; k++)
	{
	    row = Rperm1 [k] ;
	    Ci [row] = nfr ;
	}
	/* permute InFront so that InFront [row] = i if the original row is
	 * in front i */
	for (row = 0 ; row < n_row ; row++)
	{
	    InFront [row] = Ci [row] ;
	}
	/* ] no longer need Ci as workspace */
    }

#ifndef NDEBUG
    for (row = 0 ; row < n_row ; row++)
    {
	if (InFront [row] == nfr)
	{
	    DEBUG1 (("    Row "ID" in Dummy Front "ID"\n", row, nfr)) ;
	}
	else if (InFront [row] == EMPTY)
	{
	    DEBUG1 (("    singleton Row "ID"\n", row)) ;
	}
	else
	{
	    DEBUG1 (("    Row "ID" in Front "ID"\n", row, nfr)) ;
	}
    }
    for (i = 0 ; i <= nfr ; i++)
    {
	DEBUG1 (("Front "ID" has 1strows: "ID" pivcols "ID"\n",
		i, f1rows, fpivcol)) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* copy front information into Symbolic object */
    /* ---------------------------------------------------------------------- */

    k = n1 ;
    for (i = 0 ; i < nfr ; i++)
    {
	fpivcol = Fr_npivcol [i] ;
	DEBUG1 (("Front "ID" k "ID" npivcol "ID" nrows "ID" ncols "ID"\n",
	    i, k, fpivcol, Fr_nrows [i], Fr_ncols [i])) ;
	k += fpivcol ;
	/* copy Front info into Symbolic object from SW */
	Front_npivcol [i] = fpivcol ;
	Front_parent [i] = Fr_parent [i] ;
    }

    /* assign empty columns to dummy placehold front nfr */
    DEBUG1 (("Dummy Cols in Front "ID" : "ID"\n", nfr, n_col-k)) ;
    Front_npivcol [nfr] = n_col - k ;
    Front_parent [nfr] = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* find initial row permutation */
    /* ---------------------------------------------------------------------- */

    /* order the singleton pivot rows */
    for (k = 0 ; k < n1 ; k++)
    {
	Rperm_init [k] = Rperm1 [k] ;
    }

    /* determine the first row in each front (in the new row ordering) */
    for (i = 0 ; i < nfr ; i++)
    {
	f1rows = Front_1strow [i] ;
	DEBUG1 (("Front "ID" : npivcol "ID" parent "ID,
	    i, Front_npivcol [i], Front_parent [i])) ;
	DEBUG1 (("    1st rows in Front "ID" : "ID"\n", i, f1rows)) ;
	Front_1strow [i] = k ;
	k += f1rows ;
    }

    /* assign empty rows to dummy placehold front nfr */
    DEBUG1 (("Rows in Front "ID" (dummy): "ID"\n", nfr, n_row-k)) ;
    Front_1strow [nfr] = k ;
    DEBUG1 (("nfr "ID" 1strow[nfr] "ID" nrow "ID"\n", nfr, k, n_row)) ;

    /* Use Ci as temporary workspace for F1 */
    F1 = Ci ;				/* [ of size nfr+1 */
    ASSERT (Clen >= 2*n_row + nfr+1) ;

    for (i = 0 ; i <= nfr ; i++)
    {
	F1 [i] = Front_1strow [i] ;
    }

    for (row = 0 ; row < n_row ; row++)
    {
	i = InFront [row] ;
	if (i != EMPTY)
	{
	    newrow = F1 [i]++ ;
	    ASSERT (newrow >= n1) ;
	    Rperm_init [newrow] = row ;
	}
    }
    Rperm_init [n_row] = EMPTY ;	/* unused */

#ifndef NDEBUG
    for (k = 0 ; k < n_row ; k++)
    {
	DEBUG2 (("Rperm_init ["ID"] = "ID"\n", k, Rperm_init [k])) ;
    }
#endif

    /* ] done using F1 */

    /* ---------------------------------------------------------------------- */
    /* find the diagonal map */
    /* ---------------------------------------------------------------------- */

    /* Rperm_init [newrow] = row gives the row permutation that is implied
     * by the column permutation, where "row" is a row index of the original
     * matrix A.  It is not dependent on the Rperm_2by2 permutation, which
     * only redefines the "diagonal".   Both are used to construct the
     * Diagonal_map.  Diagonal_map only needs to be defined for
     * k = n1 to nn - nempty, but go ahead and define it for all of
     * k = 0 to nn */

    if (prefer_diagonal)
    {
	Int *Diagonal_map ;
	ASSERT (n_row == n_col && nn == n_row) ;
	ASSERT (nempty_row == nempty_col && nempty == nempty_row) ;

	/* allocate the Diagonal_map */
	Symbolic->Diagonal_map = (Int *) UMF_malloc (n_col+1, sizeof (Int)) ;
	Diagonal_map = Symbolic->Diagonal_map ;
	if (Diagonal_map == (Int *) NULL)
	{
	    /* :: out of memory (diagonal map) :: */
	    DEBUGm4 (("out of memory: Diagonal map\n")) ;
	    Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	    error (&Symbolic, SW) ;
	    return (UMFPACK_ERROR_out_of_memory) ;
	}

	/* use Ci as workspace to compute the inverse of Rperm_init [ */
	for (newrow = 0 ; newrow < nn ; newrow++)
	{
	    oldrow = Rperm_init [newrow] ;
	    ASSERT (oldrow >= 0 && oldrow < nn) ;
	    Ci [oldrow] = newrow ;
	}
	if (strategy == UMFPACK_STRATEGY_2BY2)
	{
	    ASSERT (Rperm_2by2 != (Int *) NULL) ;
	    for (newcol = 0 ; newcol < nn ; newcol++)
	    {
		oldcol = Cperm_init [newcol] ;
		/* 2-by-2 pivoting done in S */
		oldrow = Rperm_2by2 [oldcol] ;
		newrow = Ci [oldrow] ;
		Diagonal_map [newcol] = newrow ;
	    }
	}
	else
	{
	    for (newcol = 0 ; newcol < nn ; newcol++)
	    {
		oldcol = Cperm_init [newcol] ;
		/* no 2-by-2 pivoting in S */
		oldrow = oldcol ;
		newrow = Ci [oldrow] ;
		Diagonal_map [newcol] = newrow ;
	    }
	}

#ifndef NDEBUG
	DEBUG1 (("\nDiagonal map:\n")) ;
	for (newcol = 0 ; newcol < nn ; newcol++)
	{
	    oldcol = Cperm_init [newcol] ;
	    DEBUG3 (("oldcol "ID" newcol "ID":\n", oldcol, newcol)) ;
	    for (p = Ap [oldcol] ; p < Ap [oldcol+1] ; p++)
	    {
		Entry aij ;
		oldrow = Ai [p] ;
		newrow = Ci [oldrow] ;
		if (Ax != (double *) NULL
#ifdef COMPLEX
		 && Az != (double *) NULL
#endif
			)
		{
		    ASSIGN (aij, Ax [p], Az [p]) ;
		}
		if (oldrow == oldcol)
		{
		    DEBUG2 (("     old diagonal : oldcol "ID" oldrow "ID" ",
			    oldcol, oldrow)) ;
		    EDEBUG2 (aij) ;
		    DEBUG2 (("\n")) ;
		}
		if (newrow == Diagonal_map [newcol])
		{
		    DEBUG2 (("     MAP diagonal : newcol "ID" MAProw "ID" ",
			    newcol, Diagonal_map [newrow])) ;
		    EDEBUG2 (aij) ;
		    DEBUG2 (("\n")) ;
		}
	    }
	}
#endif
	/* done using Ci as workspace ] */

    }

    /* ---------------------------------------------------------------------- */
    /* find the leftmost descendant of each front */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i <= nfr ; i++)
    {
	Front_leftmostdesc [i] = EMPTY ;
    }

    for (i = 0 ; i < nfr ; i++)
    {
	/* start at i and walk up the tree */
	DEBUG2 (("Walk up front tree from "ID"\n", i)) ;
	j = i ;
	while (j != EMPTY && Front_leftmostdesc [j] == EMPTY)
	{
	    DEBUG3 (("  Leftmost desc of "ID" is "ID"\n", j, i)) ;
	    Front_leftmostdesc [j] = i ;
	    j = Front_parent [j] ;
	    DEBUG3 (("  go to j = "ID"\n", j)) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* find the frontal matrix chains and max frontal matrix sizes */
    /* ---------------------------------------------------------------------- */

    maxnrows = 1 ;		/* max # rows in any front */
    maxncols = 1 ;		/* max # cols in any front */
    dmaxfrsize = 1 ;		/* max frontal matrix size */

    /* start the first chain */
    nchains = 0 ;		/* number of chains */
    Chain_start [0] = 0 ;	/* front 0 starts a new chain */
    maxrows = 1 ;		/* max # rows for any front in current chain */
    maxcols = 1 ;		/* max # cols for any front in current chain */
    DEBUG1 (("Constructing chains:\n")) ;

    for (i = 0 ; i < nfr ; i++)
    {
	/* get frontal matrix info */
	fpivcol  = Front_npivcol [i] ;	    /* # candidate pivot columns */
	fallrows = Fr_nrows [i] ;	    /* all rows (not just Schur comp) */
	fallcols = Fr_ncols [i] ;	    /* all cols (not just Schur comp) */
	parent = Front_parent [i] ;	    /* parent in column etree */
	fpiv = MIN (fpivcol, fallrows) ;    /* # pivot rows and cols */
	maxrows = MAX (maxrows, fallrows) ;
	maxcols = MAX (maxcols, fallcols) ;

	DEBUG1 (("Front: "ID", pivcol "ID", "ID"-by-"ID" parent "ID
	    ", npiv "ID" Chain: maxrows "ID" maxcols "ID"\n", i, fpivcol,
	    fallrows, fallcols, parent, fpiv, maxrows, maxcols)) ;

	if (parent != i+1)
	{
	    /* this is the end of a chain */
	    double s ;
	    DEBUG1 (("\nEnd of chain "ID"\n", nchains)) ;

	    /* make sure maxrows is an odd number */
	    ASSERT (maxrows >= 0) ;
	    if (maxrows % 2 == 0) maxrows++ ;

	    DEBUG1 (("Chain maxrows "ID" maxcols "ID"\n", maxrows, maxcols)) ;

	    Chain_maxrows [nchains] = maxrows ;
	    Chain_maxcols [nchains] = maxcols ;

	    /* keep track of the maximum front size for all chains */

	    /* for Info only: */
	    s = (double) maxrows * (double) maxcols ;
	    dmaxfrsize = MAX (dmaxfrsize, s) ;

	    /* for the subsequent numerical factorization */
	    maxnrows = MAX (maxnrows, maxrows) ;
	    maxncols = MAX (maxncols, maxcols) ;

	    DEBUG1 (("Chain dmaxfrsize %g\n\n", dmaxfrsize)) ;

	    /* start the next chain */
	    nchains++ ;
	    Chain_start [nchains] = i+1 ;
	    maxrows = 1 ;
	    maxcols = 1 ;
	}
    }

    /* for Info only: */
    dmaxfrsize = ceil (dmaxfrsize) ;
    DEBUGm1 (("dmaxfrsize %30.20g Int_MAX "ID"\n", dmaxfrsize, Int_MAX)) ;
    ASSERT (Symbolic->nchains == nchains) ;

    /* For allocating objects in umfpack_numeric (does not include all possible
     * pivots, particularly pivots from prior fronts in the chain.  Need to add
     * nb to these to get the # of columns in the L block, for example.  This
     * is the largest row dimension and largest column dimension of any frontal
     * matrix.  maxnrows is always odd. */
    Symbolic->maxnrows = maxnrows ;
    Symbolic->maxncols = maxncols ;
    DEBUGm3 (("maxnrows "ID" maxncols "ID"\n", maxnrows, maxncols)) ;

    /* ---------------------------------------------------------------------- */
    /* find the initial element sizes */
    /* ---------------------------------------------------------------------- */

    if (max_rdeg > dense_row_threshold)
    {
	/* there are one or more dense rows in the input matrix */
	/* count the number of dense rows in each column */
	/* use Ci as workspace for inverse of Rperm_init [ */
	ASSERT (Esize != (Int *) NULL) ;
	for (newrow = 0 ; newrow < n_row ; newrow++)
	{
	    oldrow = Rperm_init [newrow] ;
	    ASSERT (oldrow >= 0 && oldrow < nn) ;
	    Ci [oldrow] = newrow ;
	}
	for (col = n1 ; col < n_col - nempty_col ; col++)
	{
	    oldcol = Cperm_init [col] ;
	    esize = Cdeg [oldcol] ;
	    ASSERT (esize > 0) ;
	    for (p = Ap [oldcol] ; p < Ap [oldcol+1] ; p++)
	    {
		oldrow = Ai [p] ;
		newrow = Ci [oldrow] ;
		if (newrow >= n1 && Rdeg [oldrow] > dense_row_threshold)
		{
		    esize-- ;
		}
	    }
	    ASSERT (esize >= 0) ;
	    Esize [col - n1] = esize ;
	}
	/* done using Ci as workspace ] */
    }

    /* If there are no dense rows, then Esize [col-n1] is identical to
     * Cdeg [col], once Cdeg is permuted below */

    /* ---------------------------------------------------------------------- */
    /* permute Cdeg and Rdeg according to initial column and row permutation */
    /* ---------------------------------------------------------------------- */

    /* use Ci as workspace [ */
    for (k = 0 ; k < n_col ; k++)
    {
	Ci [k] = Cdeg [Cperm_init [k]] ;
    }
    for (k = 0 ; k < n_col ; k++)
    {
	Cdeg [k] = Ci [k] ;
    }
    for (k = 0 ; k < n_row ; k++)
    {
	Ci [k] = Rdeg [Rperm_init [k]] ;
    }
    for (k = 0 ; k < n_row ; k++)
    {
	Rdeg [k] = Ci [k] ;
    }
    /* done using Ci as workspace ] */

    /* ---------------------------------------------------------------------- */
    /* simulate UMF_kernel_init */
    /* ---------------------------------------------------------------------- */

    /* count elements and tuples at tail, LU factors of singletons, and
     * head and tail markers */

    dlnz = n_inner ;	/* upper limit of nz in L (incl diag) */
    dunz = dlnz ;	/* upper limit of nz in U (incl diag) */

    /* head marker */
    head_usage  = 1 ;
    dhead_usage = 1 ;

    /* tail markers: */
    tail_usage  = 2 ;
    dtail_usage = 2 ;

    /* allocate the Rpi and Rpx workspace for UMF_kernel_init (incl. headers) */
    tail_usage  +=  UNITS (Int *, n_row+1) +  UNITS (Entry *, n_row+1) + 2 ;
    dtail_usage += DUNITS (Int *, n_row+1) + DUNITS (Entry *, n_row+1) + 2 ;
    DEBUG1 (("Symbolic usage after Rpi/Rpx allocation: head "ID" tail "ID"\n",
	head_usage, tail_usage)) ;

    /* LU factors for singletons, at the head of memory */
    for (k = 0 ; k < n1 ; k++)
    {
	lnz = Cdeg [k] - 1 ;
	unz = Rdeg [k] - 1 ;
	dlnz += lnz ;
	dunz += unz ;
	DEBUG1 (("singleton k "ID" pivrow "ID" pivcol "ID" lnz "ID" unz "ID"\n",
	    k, Rperm_init [k], Cperm_init [k], lnz, unz)) ;
	head_usage  +=  UNITS (Int, lnz) +  UNITS (Entry, lnz)
		    +   UNITS (Int, unz) +  UNITS (Entry, unz) ;
	dhead_usage += DUNITS (Int, lnz) + DUNITS (Entry, lnz)
		    +  DUNITS (Int, unz) + DUNITS (Entry, unz) ;
    }
    DEBUG1 (("Symbolic init head usage: "ID" for LU singletons\n",head_usage)) ;

    /* column elements: */
    for (k = n1 ; k < n_col - nempty_col; k++)
    {
	esize = Esize ? Esize [k-n1] : Cdeg [k] ;
	DEBUG2 (("   esize: "ID"\n", esize)) ;
	ASSERT (esize >= 0) ;
	if (esize > 0)
	{
	    tail_usage  +=  GET_ELEMENT_SIZE (esize, 1) + 1 ;
	    dtail_usage += DGET_ELEMENT_SIZE (esize, 1) + 1 ;
	}
    }

    /* dense row elements */
    if (Esize)
    {
	Int nrow_elements = 0 ;
	for (k = n1 ; k < n_row - nempty_row ; k++)
	{
	    rdeg = Rdeg [k] ;
	    if (rdeg > dense_row_threshold)
	    {
		tail_usage  += GET_ELEMENT_SIZE (1, rdeg) + 1 ;
		dtail_usage += GET_ELEMENT_SIZE (1, rdeg) + 1 ;
		nrow_elements++ ;
	    }
	}
	Info [UMFPACK_NDENSE_ROW] = nrow_elements ;
    }

    DEBUG1 (("Symbolic usage: "ID" = head "ID" + tail "ID" after els\n",
	head_usage + tail_usage, head_usage, tail_usage)) ;

    /* compute the tuple lengths */
    if (Esize)
    {
	/* row tuples */
	for (row = n1 ; row < n_row ; row++)
	{
	    rdeg = Rdeg [row] ;
	    tlen = (rdeg > dense_row_threshold) ? 1 : rdeg ;
	    tail_usage  += 1 +  UNITS (Tuple, TUPLES (tlen)) ;
	    dtail_usage += 1 + DUNITS (Tuple, TUPLES (tlen)) ;
	}
	/* column tuples */
	for (col = n1 ; col < n_col - nempty_col ; col++)
	{
	    /* tlen is 1 plus the number of dense rows in this column */
	    esize = Esize [col - n1] ;
	    tlen = (esize > 0) + (Cdeg [col] - esize) ;
	    tail_usage  += 1 +  UNITS (Tuple, TUPLES (tlen)) ;
	    dtail_usage += 1 + DUNITS (Tuple, TUPLES (tlen)) ;
	}
	for ( ; col < n_col ; col++)
	{
	    tail_usage  += 1 +  UNITS (Tuple, TUPLES (0)) ;
	    dtail_usage += 1 + DUNITS (Tuple, TUPLES (0)) ;
	}
    }
    else
    {
	/* row tuples */
	for (row = n1 ; row < n_row ; row++)
	{
	    tlen = Rdeg [row] ;
	    tail_usage  += 1 +  UNITS (Tuple, TUPLES (tlen)) ;
	    dtail_usage += 1 + DUNITS (Tuple, TUPLES (tlen)) ;
	}
	/* column tuples */
	for (col = n1 ; col < n_col ; col++)
	{
	    tail_usage  += 1 +  UNITS (Tuple, TUPLES (1)) ;
	    dtail_usage += 1 + DUNITS (Tuple, TUPLES (1)) ;
	}
    }

    Symbolic->num_mem_init_usage = head_usage + tail_usage ;
    DEBUG1 (("Symbolic usage: "ID" = head "ID" + tail "ID" final\n",
	Symbolic->num_mem_init_usage, head_usage, tail_usage)) ;

    ASSERT (UMF_is_permutation (Rperm_init, Ci, n_row, n_row)) ;

    /* initial head and tail usage in Numeric->Memory */
    dmax_usage = dhead_usage + dtail_usage ;
    dmax_usage = MAX (Symbolic->num_mem_init_usage, ceil (dmax_usage)) ;
    Info [UMFPACK_VARIABLE_INIT_ESTIMATE] = dmax_usage ;

    /* In case Symbolic->num_mem_init_usage overflows, keep as a double, too */
    Symbolic->dnum_mem_init_usage = dmax_usage ;

    /* free the Rpi and Rpx workspace */
    tail_usage  -=  UNITS (Int *, n_row+1) +  UNITS (Entry *, n_row+1) ;
    dtail_usage -= DUNITS (Int *, n_row+1) + DUNITS (Entry *, n_row+1) ;

    /* ---------------------------------------------------------------------- */
    /* simulate UMF_kernel, assuming unsymmetric pivoting */
    /* ---------------------------------------------------------------------- */

    /* Use Ci as temporary workspace for link lists [ */
    Link = Ci ;
    for (i = 0 ; i < nfr ; i++)
    {
	Link [i] = EMPTY ;
    }

    flops = 0 ;			/* flop count upper bound */

    for (chain = 0 ; chain < nchains ; chain++)
    {
	double fsize ;
	f1 = Chain_start [chain] ;
	f2 = Chain_start [chain+1] - 1 ;

	/* allocate frontal matrix working array (C, L, and U) */
	dr = Chain_maxrows [chain] ;
	dc = Chain_maxcols [chain] ;
	fsize =
	      nb*nb	    /* LU is nb-by-nb */
	    + dr*nb	    /* L is dr-by-nb */
	    + nb*dc	    /* U is nb-by-dc, stored by rows */
	    + dr*dc ;	    /* C is dr by dc */
	dtail_usage += DUNITS (Entry, fsize) ;
	dmax_usage = MAX (dmax_usage, dhead_usage + dtail_usage) ;

	for (i = f1 ; i <= f2 ; i++)
	{

	    /* get frontal matrix info */
	    fpivcol  = Front_npivcol [i] ; /* # candidate pivot columns */
	    fallrows = Fr_nrows [i] ;   /* all rows (not just Schur comp*/
	    fallcols = Fr_ncols [i] ;   /* all cols (not just Schur comp*/
	    parent = Front_parent [i] ; /* parent in column etree */
	    fpiv = MIN (fpivcol, fallrows) ;	/* # pivot rows and cols */
	    f = (double) fpiv ;
	    r = fallrows - fpiv ;		/* # rows in Schur comp. */
	    c = fallcols - fpiv ;		/* # cols in Schur comp. */

	    /* assemble all children of front i in column etree */
	    for (child = Link [i] ; child != EMPTY ; child = Link [child])
	    {
		ASSERT (child >= 0 && child < i) ;
		ASSERT (Front_parent [child] == i) ;
		/* free the child element and remove it from tuple lists */
		cp = MIN (Front_npivcol [child], Fr_nrows [child]) ;
		cr = Fr_nrows [child] - cp ;
		cc = Fr_ncols [child] - cp ;
		ASSERT (cp >= 0 && cr >= 0 && cc >= 0) ;
		dtail_usage -= ELEMENT_SIZE (cr, cc) ;

	    }

	    /* The flop count computed here is "canonical". */

	    /* factorize the frontal matrix */
	    flops += DIV_FLOPS * (f*r + (f-1)*f/2)  /* scale pivot columns */
		/* f outer products: */
		+ MULTSUB_FLOPS * (f*r*c + (r+c)*(f-1)*f/2 + (f-1)*f*(2*f-1)/6);

	    /* count nonzeros and memory usage in double precision */
	    dlf = (f*f-f)/2 + f*r ;		/* nz in L below diagonal */
	    duf = (f*f-f)/2 + f*c ;		/* nz in U above diagonal */
	    dlnz += dlf ;
	    dunz += duf ;

	    /* store f columns of L and f rows of U */
	    dhead_usage +=
		DUNITS (Entry, dlf + duf)   /* numerical values (excl diag) */
		+ DUNITS (Int, r + c + f) ; /* indices (compressed) */

	    if (parent != EMPTY)
	    {
		/* create new element and place in tuple lists */
		dtail_usage += ELEMENT_SIZE (r, c) ;

		/* place in link list of parent */
		Link [i] = Link [parent] ;
		Link [parent] = i ;
	    }

	    /* keep track of peak Numeric->Memory usage */
	    dmax_usage = MAX (dmax_usage, dhead_usage + dtail_usage) ;

	}

	/* free the current frontal matrix */
	dtail_usage -= DUNITS (Entry, fsize) ;
    }

    dhead_usage = ceil (dhead_usage) ;
    dmax_usage = ceil (dmax_usage) ;
    Symbolic->num_mem_size_est = dhead_usage ;
    Symbolic->num_mem_usage_est = dmax_usage ;
    Symbolic->lunz_bound = dlnz + dunz - n_inner ;

    /* ] done using Ci as workspace for Link array */

    /* ---------------------------------------------------------------------- */
    /* estimate total memory usage in UMFPACK_numeric */
    /* ---------------------------------------------------------------------- */

    UMF_set_stats (
	Info,
	Symbolic,
	dmax_usage,		/* estimated peak size of Numeric->Memory */
	dhead_usage,		/* estimated final size of Numeric->Memory */
	flops,			/* estimated "true flops" */
	dlnz,			/* estimated nz in L */
	dunz,			/* estimated nz in U */
	dmaxfrsize,		/* estimated largest front size */
	(double) n_col,		/* worst case Numeric->Upattern size */
	(double) n_inner,	/* max possible pivots to be found */
	(double) maxnrows,	/* estimated largest #rows in front */
	(double) maxncols,	/* estimated largest #cols in front */
	TRUE,			/* assume scaling is to be performed */
	prefer_diagonal,
	ESTIMATE) ;

    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    for (i = 0 ; i < nchains ; i++)
    {
	DEBUG2 (("Chain "ID" start "ID" end "ID" maxrows "ID" maxcols "ID"\n",
		i, Chain_start [i], Chain_start [i+1] - 1,
		Chain_maxrows [i], Chain_maxcols [i])) ;
	UMF_dump_chain (Chain_start [i], Fr_parent, Fr_npivcol, Fr_nrows,
	    Fr_ncols, nfr) ;
    }
    fpivcol = 0 ;
    for (i = 0 ; i < nfr ; i++)
    {
	fpivcol = MAX (fpivcol, Front_npivcol [i]) ;
    }
    DEBUG0 (("Max pivot cols in any front: "ID"\n", fpivcol)) ;
    DEBUG1 (("Largest front: maxnrows "ID" maxncols "ID" dmaxfrsize %g\n",
	maxnrows, maxncols, dmaxfrsize)) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* UMFPACK_symbolic was successful, return the object handle */
    /* ---------------------------------------------------------------------- */

    Symbolic->valid = SYMBOLIC_VALID ;
    *SymbolicHandle = (void *) Symbolic ;

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    /* (6) The last of the workspace is free'd.  The final Symbolic object
     * consists of 12 to 14 allocated objects.  Its final total size is lies
     * roughly between 4*n and 13*n for a square matrix, which is all that is
     * left of the memory allocated by this routine.  If an error occurs, the
     * entire Symbolic object is free'd when this routine returns (the error
     * return routine, below).
     */

    free_work (SW) ;

    DEBUG0 (("(3)Symbolic UMF_malloc_count - init_count = "ID"\n",
	UMF_malloc_count - init_count)) ;
    ASSERT (UMF_malloc_count == init_count + 12
	+ (Symbolic->Esize != (Int *) NULL)
	+ (Symbolic->Diagonal_map != (Int *) NULL)) ;

    /* ---------------------------------------------------------------------- */
    /* get the time used by UMFPACK_*symbolic */
    /* ---------------------------------------------------------------------- */

    umfpack_toc (stats) ;
    Info [UMFPACK_SYMBOLIC_WALLTIME] = stats [0] ;
    Info [UMFPACK_SYMBOLIC_TIME] = stats [1] ;

    return (UMFPACK_OK) ;
}


/* ========================================================================== */
/* === free_work ============================================================ */
/* ========================================================================== */

PRIVATE void free_work
(
    SWType *SW
)
{
    if (SW)
    {
	SW->Rperm_2by2 = (Int *) UMF_free ((void *) SW->Rperm_2by2) ;
	SW->InvRperm1 = (Int *) UMF_free ((void *) SW->InvRperm1) ;
	SW->Rs = (double *) UMF_free ((void *) SW->Rs) ;
	SW->Si = (Int *) UMF_free ((void *) SW->Si) ;
	SW->Sp = (Int *) UMF_free ((void *) SW->Sp) ;
	SW->Ci = (Int *) UMF_free ((void *) SW->Ci) ;
	SW->Front_npivcol = (Int *) UMF_free ((void *) SW->Front_npivcol);
	SW->Front_nrows = (Int *) UMF_free ((void *) SW->Front_nrows) ;
	SW->Front_ncols = (Int *) UMF_free ((void *) SW->Front_ncols) ;
	SW->Front_parent = (Int *) UMF_free ((void *) SW->Front_parent) ;
	SW->Front_cols = (Int *) UMF_free ((void *) SW->Front_cols) ;
	SW->Cperm1 = (Int *) UMF_free ((void *) SW->Cperm1) ;
	SW->Rperm1 = (Int *) UMF_free ((void *) SW->Rperm1) ;
	SW->InFront = (Int *) UMF_free ((void *) SW->InFront) ;
    }
}


/* ========================================================================== */
/* === error ================================================================ */
/* ========================================================================== */

/* Error return from UMFPACK_symbolic.  Free all allocated memory. */

PRIVATE void error
(
    SymbolicType **Symbolic,
    SWType *SW
)
{

    free_work (SW) ;
    UMFPACK_free_symbolic ((void **) Symbolic) ;
    ASSERT (UMF_malloc_count == init_count) ;
}
