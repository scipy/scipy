/* ========================================================================== */
/* === UMF_set_stats ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Sets statistics in Info array.  Calculates everything in double precision,
    rather than Int or size_t, so that usage estimates can be computed even if
    the problem is so large that it would cause integer overflow.

    This routine has many double relop's, but the NaN case is ignored.
*/

#include "umf_internal.h"
#include "umf_symbolic_usage.h"

GLOBAL void UMF_set_stats
(
    double Info [ ],
    SymbolicType *Symbolic,
    double max_usage,		/* peak size of Numeric->Memory, in Units */
    double num_mem_size,	/* final size of Numeric->Memory, in Units */
    double flops,		/* "true flops" */
    double lnz,			/* nz in L */
    double unz,			/* nz in U */
    double maxfrsize,		/* largest front size */
    double ulen,		/* size of Numeric->Upattern */
    double npiv,		/* number of pivots found */
    double maxnrows,		/* largest #rows in front */
    double maxncols,		/* largest #cols in front */
    Int scale,			/* true if scaling the rows of A */
    Int prefer_diagonal,	/* true if diagonal pivoting (only square A) */
    Int what			/* ESTIMATE or ACTUAL */
)
{

    double sym_size, work_usage, nn, n_row, n_col, n_inner, num_On_size1,
	num_On_size2, num_usage, sym_maxncols, sym_maxnrows, elen, n1 ;

    n_col = Symbolic->n_col ;
    n_row = Symbolic->n_row ;
    n1 = Symbolic->n1 ;
    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;
    sym_maxncols = MIN (Symbolic->maxncols + Symbolic->nb, n_col) ;
    sym_maxnrows = MIN (Symbolic->maxnrows + Symbolic->nb, n_row) ;
    elen = (n_col - n1) + (n_row - n1) + MIN (n_col - n1, n_row - n1) + 1 ;

    /* final Symbolic object size */
    sym_size = UMF_symbolic_usage (Symbolic->n_row, Symbolic->n_col,
	Symbolic->nchains, Symbolic->nfr, Symbolic->esize, prefer_diagonal) ;

    /* size of O(n) part of Numeric object during factorization, */
    /* except Numeric->Memory and Numeric->Upattern */
    num_On_size1 =
	DUNITS (NumericType, 1)		/* Numeric structure */
	+ DUNITS (Entry, n_inner+1)	/* D */
	+ 4 * DUNITS (Int, n_row+1)	/* Rperm, Lpos, Uilen, Uip */
	+ 4 * DUNITS (Int, n_col+1)	/* Cperm, Upos, Lilen, Lip */
	+ (scale ? DUNITS (Entry, n_row) : 0) ;   /* Rs, row scale factors */

    /* size of O(n) part of Numeric object after factorization, */
    /* except Numeric->Memory and Numeric->Upattern */
    num_On_size2 =
	DUNITS (NumericType, 1)		/* Numeric structure */
	+ DUNITS (Entry, n_inner+1)	/* D */
	+ DUNITS (Int, n_row+1)		/* Rperm */
	+ DUNITS (Int, n_col+1)		/* Cperm */
	+ 6 * DUNITS (Int, npiv+1)	/* Lpos, Uilen, Uip, Upos, Lilen, Lip */
	+ (scale ? DUNITS (Entry, n_row) : 0) ;	    /* Rs, row scale factors */

    DEBUG1 (("num O(n) size2: %g\n", num_On_size2)) ;

    /* peak size of Numeric->Memory, including LU factors, current frontal
     * matrix, elements, and tuple lists.  */
    Info [UMFPACK_VARIABLE_PEAK + what] = max_usage ;

    /* final size of Numeric->Memory (LU factors only) */
    Info [UMFPACK_VARIABLE_FINAL + what] = num_mem_size ;

    /* final size of Numeric object, including Numeric->Memory and ->Upattern */
    Info [UMFPACK_NUMERIC_SIZE + what] =
	num_On_size2
	+ num_mem_size		/* final Numeric->Memory size */
	+ DUNITS (Int, ulen+1) ;/* Numeric->Upattern (from Work->Upattern) */

    DEBUG1 (("num mem size: %g\n", num_mem_size)) ;
    DEBUG1 (("ulen units %g\n", DUNITS (Int, ulen))) ;
    DEBUG1 (("numeric size %g\n", Info [UMFPACK_NUMERIC_SIZE + what])) ;

    /* largest front size (working array size, or actual size used) */
    Info [UMFPACK_MAX_FRONT_SIZE + what] = maxfrsize ;
    Info [UMFPACK_MAX_FRONT_NROWS + what] = maxnrows ;
    Info [UMFPACK_MAX_FRONT_NCOLS + what] = maxncols ;
    DEBUGm4 (("maxnrows %g maxncols %g\n", maxnrows, maxncols)) ;
    DEBUGm4 (("maxfrsize %g\n", maxfrsize)) ;

    /* UMF_kernel usage, from work_alloc routine in umf_kernel.c */
    work_usage =
	/* Work-> arrays, except for current frontal matrix which is allocated
	 * inside Numeric->Memory. */
	2 * DUNITS (Entry, sym_maxnrows + 1)	/* Wx, Wy */
	+ 2 * DUNITS (Int, n_row+1)		/* Frpos, Lpattern */
	+ 2 * DUNITS (Int, n_col+1)		/* Fcpos, Upattern */
	+ DUNITS (Int, nn + 1)			/* Wp */
	+ DUNITS (Int, MAX (n_col, sym_maxnrows) + 1)	/* Wrp */
	+ 2 * DUNITS (Int, sym_maxnrows + 1)	/* Frows, Wm */
	+ 3 * DUNITS (Int, sym_maxncols + 1)	/* Fcols, Wio, Woi */
	+ DUNITS (Int, MAX (sym_maxnrows, sym_maxncols) + 1)	/* Woo */
	+ DUNITS (Int, elen)			/* E */
	+ DUNITS (Int, Symbolic->nfr + 1)	/* Front_new1strow */
	+ ((n_row == n_col) ? (2 * DUNITS (Int, nn)) : 0) ;  /* Diag map,imap */

    /* Peak memory for just UMFPACK_numeric. */
    num_usage =
	sym_size	/* size of Symbolic object */
	+ num_On_size1	/* O(n) part of Numeric object (excl. Upattern) */
	+ work_usage	/* Work-> arrays (including Upattern) */
	+ max_usage ;	/* peak size of Numeric->Memory */

    /* peak memory usage for both UMFPACK_*symbolic and UMFPACK_numeric. */
    Info [UMFPACK_PEAK_MEMORY + what] =
	MAX (Symbolic->peak_sym_usage, num_usage) ;

    Info [UMFPACK_FLOPS + what] = flops ;
    Info [UMFPACK_LNZ + what] = lnz ;
    Info [UMFPACK_UNZ + what] = unz ;
}
