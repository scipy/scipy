/* ========================================================================== */
/* === UMFPACK_numeric ====================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Factorizes A into its LU factors, given a symbolic
    pre-analysis computed by UMFPACK_symbolic.  See umfpack_numeric.h for a
    description.

    Dynamic memory allocation:  substantial.  See comments (1) through (7),
    below.  If an error occurs, all allocated space is free'd by UMF_free.
    If successful, the Numeric object contains 11 to 13 objects allocated by
    UMF_malloc that hold the LU factors of the input matrix.
*/

#include "umf_internal.h"
#include "umf_valid_symbolic.h"
#include "umf_set_stats.h"
#include "umf_kernel.h"
#include "umf_malloc.h"
#include "umf_free.h"
#include "umf_realloc.h"

#ifndef NDEBUG
PRIVATE Int init_count ;
#endif

PRIVATE Int work_alloc
(
    WorkType *Work,
    SymbolicType *Symbolic
) ;

PRIVATE void free_work
(
    WorkType *Work
) ;

PRIVATE Int numeric_alloc
(
    NumericType **NumericHandle,
    SymbolicType *Symbolic,
    double alloc_init,
    Int scale
) ;

PRIVATE void error
(
    NumericType **Numeric,
    WorkType *Work
) ;


/* ========================================================================== */
/* === UMFPACK_numeric ====================================================== */
/* ========================================================================== */

GLOBAL Int UMFPACK_numeric
(
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    void *SymbolicHandle,
    void **NumericHandle,
    const double Control [UMFPACK_CONTROL],
    double User_Info [UMFPACK_INFO]
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    NumericType *Numeric ;
    SymbolicType *Symbolic ;
    WorkType WorkSpace, *Work ;
    Int n_row, n_col, n_inner, newsize, i, status, *inew, npiv, ulen, scale ;
    Unit *mnew ;
    double Info2 [UMFPACK_INFO], *Info, alloc_init, relpt, relpt2,
	front_alloc_init, stats [2] ;

    /* ---------------------------------------------------------------------- */
    /* get the amount of time used by the process so far */
    /* ---------------------------------------------------------------------- */

    umfpack_tic (stats) ;

    /* ---------------------------------------------------------------------- */
    /* initialize and check inputs */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    UMF_dump_start ( ) ;
    init_count = UMF_malloc_count ;
    DEBUGm4 (("\nUMFPACK numeric: U transpose version\n")) ;
#endif

    /* If front_alloc_init negative then allocate that size of front in
     * UMF_start_front.  If alloc_init negative, then allocate that initial
     * size of Numeric->Memory. */

    relpt = GET_CONTROL (UMFPACK_PIVOT_TOLERANCE,
	UMFPACK_DEFAULT_PIVOT_TOLERANCE) ;
    relpt2 = GET_CONTROL (UMFPACK_SYM_PIVOT_TOLERANCE,
	UMFPACK_DEFAULT_SYM_PIVOT_TOLERANCE) ;
    alloc_init = GET_CONTROL (UMFPACK_ALLOC_INIT, UMFPACK_DEFAULT_ALLOC_INIT) ;
    front_alloc_init = GET_CONTROL (UMFPACK_FRONT_ALLOC_INIT,
	UMFPACK_DEFAULT_FRONT_ALLOC_INIT) ;
    scale = GET_CONTROL (UMFPACK_SCALE, UMFPACK_DEFAULT_SCALE) ;

    relpt   = MAX (0.0, MIN (relpt,  1.0)) ;
    relpt2  = MAX (0.0, MIN (relpt2, 1.0)) ;
    front_alloc_init = MIN (1.0, front_alloc_init) ;

    if (scale != UMFPACK_SCALE_NONE && scale != UMFPACK_SCALE_MAX)
    {
	scale = UMFPACK_DEFAULT_SCALE ;
    }

    if (User_Info != (double *) NULL)
    {
	/* return Info in user's array */
	Info = User_Info ;
	/* clear the parts of Info that are set by UMFPACK_numeric */
	for (i = UMFPACK_NUMERIC_SIZE ; i <= UMFPACK_MAX_FRONT_NCOLS ; i++)
	{
	    Info [i] = EMPTY ;
	}
	for (i = UMFPACK_NUMERIC_DEFRAG ; i < UMFPACK_IR_TAKEN ; i++)
	{
	    Info [i] = EMPTY ;
	}
    }
    else
    {
	/* no Info array passed - use local one instead */
	Info = Info2 ;
	for (i = 0 ; i < UMFPACK_INFO ; i++)
	{
	    Info [i] = EMPTY ;
	}
    }

    Symbolic = (SymbolicType *) SymbolicHandle ;
    Numeric = (NumericType *) NULL ;
    if (!UMF_valid_symbolic (Symbolic))
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_invalid_Symbolic_object ;
	return (UMFPACK_ERROR_invalid_Symbolic_object) ;
    }

    /* compute alloc_init automatically for AMD ordering */
    if (Symbolic->ordering == UMFPACK_ORDERING_AMD && alloc_init >= 0)
    {
	alloc_init = (Symbolic->nz + Symbolic->amd_lunz) / Symbolic->lunz_bound;
	alloc_init = MIN (1.0, alloc_init) ;
	alloc_init *= UMF_REALLOC_INCREASE ;
    }

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;
    n_inner = MIN (n_row, n_col) ;

    /* check for integer overflow in Numeric->Memory minimum size */
    if (INT_OVERFLOW (Symbolic->dnum_mem_init_usage * sizeof (Unit)))
    {
	/* :: int overflow, initial Numeric->Memory size :: */
	/* There's no hope to allocate a Numeric object big enough simply to
	 * hold the initial matrix, so return an out-of-memory condition */
	DEBUGm4 (("out of memory: numeric int overflow\n")) ;
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

    Info [UMFPACK_STATUS] = UMFPACK_OK ;
    Info [UMFPACK_NROW] = n_row ;
    Info [UMFPACK_NCOL] = n_col ;
    Info [UMFPACK_SIZE_OF_UNIT] = (double) (sizeof (Unit)) ;

    if (!Ap || !Ai || !Ax || !NumericHandle
#ifdef COMPLEX
	|| !Az
#endif
    )
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_argument_missing ;
	return (UMFPACK_ERROR_argument_missing) ;
    }

    Info [UMFPACK_NZ] = Ap [n_col] ;
    *NumericHandle = (void *) NULL ;

    /* ---------------------------------------------------------------------- */
    /* allocate the Work object */
    /* ---------------------------------------------------------------------- */

    /* (1) calls UMF_malloc 15 or 17 times, to obtain temporary workspace of
     * size c+1 Entry's and 2*(n_row+1) + 3*(n_col+1) + (n_col+n_inner+1) +
     * (nn+1) + * 3*(c+1) + 2*(r+1) + max(r,c) + (nfr+1) integers plus 2*nn
     * more integers if diagonal pivoting is to be done.  r is the maximum
     * number of rows in any frontal matrix, c is the maximum number of columns
     * in any frontal matrix, n_inner is min (n_row,n_col), nn is
     * max (n_row,n_col), and nfr is the number of frontal matrices.  For a
     * square matrix, this is c+1 Entry's and about 8n + 3c + 2r + max(r,c) +
     * nfr integers, plus 2n more for diagonal pivoting.
     */

    Work = &WorkSpace ;
    Work->n_row = n_row ;
    Work->n_col = n_col ;
    Work->nfr = Symbolic->nfr ;
    Work->nb = Symbolic->nb ;
    Work->n1 = Symbolic->n1 ;

    if (!work_alloc (Work, Symbolic))
    {
	DEBUGm4 (("out of memory: numeric work\n")) ;
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	error (&Numeric, Work) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }
    ASSERT (UMF_malloc_count == init_count + 16 + 2*Symbolic->prefer_diagonal) ;

    /* ---------------------------------------------------------------------- */
    /* allocate Numeric object */
    /* ---------------------------------------------------------------------- */

    /* (2) calls UMF_malloc 10 or 11 times, for a total space of
     * sizeof (NumericType) bytes, 4*(n_row+1) + 4*(n_row+1) integers, and
     * (n_inner+1) Entry's, plus n_row Entry's if row scaling is to be done.
     * sizeof (NumericType) is a small constant.  Next, it calls UMF_malloc
     * once, for the variable-sized part of the Numeric object
     * (Numeric->Memory).  The size of this object is the larger of
     * (Control [UMFPACK_ALLOC_INIT]) *  (the approximate upper bound computed
     * by UMFPACK_symbolic), and the minimum required to start the numerical
     * factorization.  * This request is reduced if it fails.
     */

    if (!numeric_alloc (&Numeric, Symbolic, alloc_init, scale))
    {
	DEBUGm4 (("out of memory: initial numeric\n")) ;
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	error (&Numeric, Work) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }
    DEBUG0 (("malloc: init_count "ID" UMF_malloc_count "ID"\n",
	init_count, UMF_malloc_count)) ;
    ASSERT (UMF_malloc_count == init_count
	+ (16 + 2*Symbolic->prefer_diagonal)
	+ (11 + (scale != UMFPACK_SCALE_NONE))) ;

    /* set control parameters */
    Numeric->relpt = relpt ;
    Numeric->relpt2 = relpt2 ;
    Numeric->alloc_init = alloc_init ;
    Numeric->front_alloc_init = front_alloc_init ;
    Numeric->scale = scale ;

    DEBUG0 (("umf relpt %g %g init %g %g inc %g red %g\n",
	relpt, relpt2, alloc_init, front_alloc_init,
	UMF_REALLOC_INCREASE, UMF_REALLOC_REDUCTION)) ;

    /* ---------------------------------------------------------------------- */
    /* scale and factorize */
    /* ---------------------------------------------------------------------- */

    /* (3) During numerical factorization (inside UMF_kernel), the variable-size
     * block of memory is increased in size via a call to UMF_realloc if it is
     * found to be too small.  During factorization, this block holds the
     * pattern and values of L and U at the top end, and the elements
     * (contibution blocks) and the current frontal matrix (Work->F*) at the
     * bottom end.  The peak size of the variable-sized object is estimated in
     * UMFPACK_*symbolic (Info [UMFPACK_VARIABLE_PEAK_ESTIMATE]), although this
     * upper bound can be very loose.  The size of the Symbolic object
     * (which is currently allocated) is in Info [UMFPACK_SYMBOLIC_SIZE], and
     * is between 2*n and 13*n integers.
     */

    DEBUG0 (("Calling umf_kernel\n")) ;
    status = UMF_kernel (Ap, Ai, Ax,
#ifdef COMPLEX
	Az,
#endif
	Numeric, Work, Symbolic) ;

    Info [UMFPACK_STATUS] = status ;
    if (status < UMFPACK_OK)
    {
	/* out of memory, or pattern has changed */
	error (&Numeric, Work) ;
	return (status) ;
    }

    Info [UMFPACK_FORCED_UPDATES] = Work->nforced ;
    Info [UMFPACK_VARIABLE_INIT] = Numeric->init_usage ;
    if (Symbolic->prefer_diagonal)
    {
	Info [UMFPACK_NOFF_DIAG] = Work->noff_diagonal ;
    }

    DEBUG0 (("malloc: init_count "ID" UMF_malloc_count "ID"\n",
	init_count, UMF_malloc_count)) ;

    npiv = Numeric->npiv ;	/* = n_inner for nonsingular matrices */
    ulen = Numeric->ulen ;	/* = 0 for square nonsingular matrices */

    /* ---------------------------------------------------------------------- */
    /* free Work object */
    /* ---------------------------------------------------------------------- */

    /* (4) After numerical factorization all of the objects allocated in step
     * (1) are freed via UMF_free, except that one object of size n_col+1 is
     * kept if there are off-diagonal nonzeros in the last pivot row (can only
     * occur for singular or rectangular matrices).  This is Work->Upattern,
     * which is transfered to Numeric->Upattern if ulen > 0.
     */

    DEBUG0 (("malloc: init_count "ID" UMF_malloc_count "ID"\n",
	init_count, UMF_malloc_count)) ;

    free_work (Work) ;

    DEBUG0 (("malloc: init_count "ID" UMF_malloc_count "ID"\n",
	init_count, UMF_malloc_count)) ;
    DEBUG0 (("Numeric->ulen: "ID" scale: "ID"\n", ulen, scale)) ;
    ASSERT (UMF_malloc_count == init_count + (ulen > 0) +
	(11 + (scale != UMFPACK_SCALE_NONE))) ;

    /* ---------------------------------------------------------------------- */
    /* reduce Lpos, Lilen, Lip, Upos, Uilen and Uip to size npiv+1 */
    /* ---------------------------------------------------------------------- */

    /* (5) Six components of the Numeric object are reduced in size if the
     * matrix is singular or rectangular.   The original size is 3*(n_row+1) +
     * 3*(n_col+1) integers.  The new size is 6*(npiv+1) integers.  For
     * square non-singular matrices, these two sizes are the same.
     */

    if (npiv < n_row)
    {
	/* reduce Lpos, Uilen, and Uip from size n_row+1 to size npiv */
	inew = (Int *) UMF_realloc (Numeric->Lpos, npiv+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Lpos = inew ;
	}
	inew = (Int *) UMF_realloc (Numeric->Uilen, npiv+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Uilen = inew ;
	}
	inew = (Int *) UMF_realloc (Numeric->Uip, npiv+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Uip = inew ;
	}
    }

    if (npiv < n_col)
    {
	/* reduce Upos, Lilen, and Lip from size n_col+1 to size npiv */
	inew = (Int *) UMF_realloc (Numeric->Upos, npiv+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Upos = inew ;
	}
	inew = (Int *) UMF_realloc (Numeric->Lilen, npiv+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Lilen = inew ;
	}
	inew = (Int *) UMF_realloc (Numeric->Lip, npiv+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Lip = inew ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* reduce Numeric->Upattern from size n_col+1 to size ulen+1 */
    /* ---------------------------------------------------------------------- */

    /* (6) The size of Numeric->Upattern (formerly Work->Upattern) is reduced
     * from size n_col+1 to size ulen + 1.  If ulen is zero, the object does
     * not exist. */

    DEBUG4 (("ulen: "ID" Upattern "ID"\n", ulen, (Int) Numeric->Upattern)) ;
    ASSERT (IMPLIES (ulen == 0, Numeric->Upattern == (Int *) NULL)) ;
    if (ulen > 0 && ulen < n_col)
    {
	inew = (Int *) UMF_realloc (Numeric->Upattern, ulen+1, sizeof (Int)) ;
	if (inew)
	{
	    Numeric->Upattern = inew ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* reduce Numeric->Memory to hold just the LU factors at the head */
    /* ---------------------------------------------------------------------- */

    /* (7) The variable-sized block (Numeric->Memory) is reduced to hold just L
     * and U, via a call to UMF_realloc, since the frontal matrices are no
     * longer needed.
     */

    newsize = Numeric->ihead ;
    if (newsize < Numeric->size)
    {
	mnew = (Unit *) UMF_realloc (Numeric->Memory, newsize, sizeof (Unit)) ;
	if (mnew)
	{
	    /* realloc succeeded (how can it fail since the size is reduced?) */
	    Numeric->Memory = mnew ;
	    Numeric->size = newsize ;
	}
    }
    Numeric->ihead = Numeric->size ;
    Numeric->itail = Numeric->ihead ;
    Numeric->tail_usage = 0 ;
    Numeric->ibig = EMPTY ;
    /* UMF_mem_alloc_tail_block can no longer be called (no tail marker) */

    /* ---------------------------------------------------------------------- */
    /* report the results and return the Numeric object */
    /* ---------------------------------------------------------------------- */

    UMF_set_stats (
	Info,
	Symbolic,
	(double) Numeric->max_usage,	/* actual peak Numeric->Memory */
	(double) Numeric->size,		/* actual final Numeric->Memory */
	Numeric->flops,			/* actual "true flops" */
	(double) Numeric->lnz + n_inner,		/* actual nz in L */
	(double) Numeric->unz + Numeric->nnzpiv,	/* actual nz in U */
	(double) Numeric->maxfrsize,	/* actual largest front size */
	(double) ulen,			/* actual Numeric->Upattern size */
	(double) npiv,			/* actual # pivots found */
	(double) Numeric->maxnrows,	/* actual largest #rows in front */
	(double) Numeric->maxncols,	/* actual largest #cols in front */
	scale != UMFPACK_SCALE_NONE,
	Symbolic->prefer_diagonal,
	ACTUAL) ;

    Info [UMFPACK_ALLOC_INIT_USED] = Numeric->alloc_init ;
    Info [UMFPACK_NUMERIC_DEFRAG] = Numeric->ngarbage ;
    Info [UMFPACK_NUMERIC_REALLOC] = Numeric->nrealloc ;
    Info [UMFPACK_NUMERIC_COSTLY_REALLOC] = Numeric->ncostly ;
    Info [UMFPACK_COMPRESSED_PATTERN] = Numeric->isize ;
    Info [UMFPACK_LU_ENTRIES] = Numeric->nLentries + Numeric->nUentries +
	    Numeric->npiv ;
    Info [UMFPACK_UDIAG_NZ] = Numeric->nnzpiv ;
    Info [UMFPACK_RSMIN] = Numeric->rsmin ;
    Info [UMFPACK_RSMAX] = Numeric->rsmax ;
    Info [UMFPACK_WAS_SCALED] = Numeric->scale ;

    /* estimate of the reciprocal of the condition number. */
    if (SCALAR_IS_ZERO (Numeric->min_udiag)
     || SCALAR_IS_ZERO (Numeric->max_udiag)
     ||	SCALAR_IS_NAN (Numeric->min_udiag)
     ||	SCALAR_IS_NAN (Numeric->max_udiag))
    {
	/* rcond is zero if there is any zero or NaN on the diagonal */
	Numeric->rcond = 0.0 ;
    }
    else
    {
	/* estimate of the recipricol of the condition number. */
	/* This is NaN if diagonal is zero-free, but has one or more NaN's. */
	Numeric->rcond = Numeric->min_udiag / Numeric->max_udiag ;
    }
    Info [UMFPACK_UMIN]  = Numeric->min_udiag ;
    Info [UMFPACK_UMAX]  = Numeric->max_udiag ;
    Info [UMFPACK_RCOND] = Numeric->rcond ;

    if (Numeric->nnzpiv < n_inner
    || SCALAR_IS_ZERO (Numeric->rcond) || SCALAR_IS_NAN (Numeric->rcond))
    {
	/* there are zeros and/or NaN's on the diagonal of U */
	DEBUG0 (("Warning, matrix is singular in umfpack_numeric\n")) ;
	DEBUG0 (("nnzpiv "ID" n_inner "ID" rcond %g\n", Numeric->nnzpiv,
	    n_inner, Numeric->rcond)) ;
	status = UMFPACK_WARNING_singular_matrix ;
	Info [UMFPACK_STATUS] = status ;
    }

    Numeric->valid = NUMERIC_VALID ;
    *NumericHandle = (void *) Numeric ;

    /* Numeric has 11 to 13 objects */
    ASSERT (UMF_malloc_count == init_count + 11 +
	+ (ulen > 0)			    /* Numeric->Upattern */
	+ (scale != UMFPACK_SCALE_NONE)) ;  /* Numeric->Rs */

    /* ---------------------------------------------------------------------- */
    /* get the time used by UMFPACK_numeric */
    /* ---------------------------------------------------------------------- */

    umfpack_toc (stats) ;
    Info [UMFPACK_NUMERIC_WALLTIME] = stats [0] ;
    Info [UMFPACK_NUMERIC_TIME] = stats [1] ;

    /* return UMFPACK_OK or UMFPACK_WARNING_singular_matrix */
    return (status) ;

}


/* ========================================================================== */
/* === numeric_alloc ======================================================== */
/* ========================================================================== */

/* Allocate the Numeric object */

PRIVATE Int numeric_alloc
(
    NumericType **NumericHandle,
    SymbolicType *Symbolic,
    double alloc_init,
    Int scale
)
{
    Int n_row, n_col, n_inner, min_usage, trying ;
    NumericType *Numeric ;
    double nsize, bsize ;

    DEBUG0 (("numeric alloc:\n")) ;

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;
    n_inner = MIN (n_row, n_col) ;
    *NumericHandle = (NumericType *) NULL ;

    /* 1 allocation:  accounted for in UMF_set_stats (num_On_size1),
     * free'd in umfpack_free_numeric */
    Numeric = (NumericType *) UMF_malloc (1, sizeof (NumericType)) ;

    if (!Numeric)
    {
	return (FALSE) ;	/* out of memory */
    }
    Numeric->valid = 0 ;
    *NumericHandle = Numeric ;

    /* 9 allocations:  accounted for in UMF_set_stats (num_On_size1),
     * free'd in umfpack_free_numeric */
    Numeric->D = (Entry *) UMF_malloc (n_inner+1, sizeof (Entry)) ;
    Numeric->Rperm = (Int *) UMF_malloc (n_row+1, sizeof (Int)) ;
    Numeric->Cperm = (Int *) UMF_malloc (n_col+1, sizeof (Int)) ;
    Numeric->Lpos = (Int *) UMF_malloc (n_row+1, sizeof (Int)) ;
    Numeric->Lilen = (Int *) UMF_malloc (n_col+1, sizeof (Int)) ;
    Numeric->Lip = (Int *) UMF_malloc (n_col+1, sizeof (Int)) ;
    Numeric->Upos = (Int *) UMF_malloc (n_col+1, sizeof (Int)) ;
    Numeric->Uilen = (Int *) UMF_malloc (n_row+1, sizeof (Int)) ;
    Numeric->Uip = (Int *) UMF_malloc (n_row+1, sizeof (Int)) ;

    /* 1 allocation if scaling:  in UMF_set_stats (num_On_size1),
     * free'd in umfpack_free_numeric */
    if (scale != UMFPACK_SCALE_NONE)
    {
	DEBUG0 (("Allocating scale factors\n")) ;
	Numeric->Rs = (double *) UMF_malloc (n_row, sizeof (double)) ;
    }
    else
    {
	DEBUG0 (("No scale factors allocated (R = I)\n")) ;
	Numeric->Rs = (double *) NULL ;
    }

    Numeric->Memory = (Unit *) NULL ;

    /* Upattern has already been allocated as part of the Work object.  If
     * the matrix is singular or rectangular, and there are off-diagonal
     * nonzeros in the last pivot row, then Work->Upattern is not free'd.
     * Instead it is transfered to Numeric->Upattern.  If it exists,
     * Numeric->Upattern is free'd in umfpack_free_numeric. */
    Numeric->Upattern = (Int *) NULL ;	/* used for singular matrices only */

    if (!Numeric->D || !Numeric->Rperm || !Numeric->Cperm || !Numeric->Upos ||
	!Numeric->Lpos || !Numeric->Lilen || !Numeric->Uilen || !Numeric->Lip ||
	!Numeric->Uip || (scale != UMFPACK_SCALE_NONE && !Numeric->Rs))
    {
	return (FALSE) ;	/* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* allocate initial Numeric->Memory for LU factors and elements */
    /* ---------------------------------------------------------------------- */

    if (alloc_init < 0)
    {
	/* -alloc_init is the exact size to initially allocate */
	nsize = -alloc_init ;
    }
    else
    {
	/* alloc_init is a ratio of the upper bound memory usage */
	nsize = (alloc_init * Symbolic->num_mem_usage_est) + 1 ;
    }
    min_usage = Symbolic->num_mem_init_usage ;

    /* Numeric->Memory must be large enough for UMF_kernel_init */
    nsize = MAX (min_usage, nsize) ;

    /* Numeric->Memory cannot be larger in size than Int_MAX / sizeof(Unit) */
    /* For ILP32 mode:  2GB (nsize cannot be bigger than 256 Mwords) */
    bsize = ((double) Int_MAX) / sizeof (Unit) - 1 ;
    DEBUG0 (("bsize %g\n", bsize)) ;
    nsize = MIN (nsize, bsize) ;

    Numeric->size = (Int) nsize ;

    DEBUG0 (("Num init %g usage_est %g numsize "ID" minusage "ID"\n",
	alloc_init, Symbolic->num_mem_usage_est, Numeric->size, min_usage)) ;

    /* allocates 1 object: */
    /* keep trying until successful, or memory request is too small */
    trying = TRUE ;
    while (trying)
    {
	Numeric->Memory = (Unit *) UMF_malloc (Numeric->size, sizeof (Unit)) ;
	if (Numeric->Memory)
	{
	    DEBUG0 (("Successful Numeric->size: "ID"\n", Numeric->size)) ;
	    return (TRUE) ;
	}
	/* too much, reduce the request (but not below the minimum) */
	/* and try again */
	trying = Numeric->size > min_usage ;
	Numeric->size = (Int)
	    (UMF_REALLOC_REDUCTION * ((double) Numeric->size)) ;
	Numeric->size = MAX (min_usage, Numeric->size) ;
    }

    return (FALSE) ;	/* we failed to allocate Numeric->Memory */
}


/* ========================================================================== */
/* === work_alloc =========================================================== */
/* ========================================================================== */

/* Allocate the Work object.  Return TRUE if successful. */

PRIVATE Int work_alloc
(
    WorkType *Work,
    SymbolicType *Symbolic
)
{
    Int n_row, n_col, nn, maxnrows, maxncols, nfr, ok, maxnrc, n1 ;

    n_row = Work->n_row ;
    n_col = Work->n_col ;
    nn = MAX (n_row, n_col) ;
    nfr = Work->nfr ;
    n1 = Symbolic->n1 ;
    ASSERT (n1 <= n_row && n1 <= n_col) ;

    maxnrows = Symbolic->maxnrows + Symbolic->nb ;
    maxnrows = MIN (n_row, maxnrows) ;
    maxncols = Symbolic->maxncols + Symbolic->nb ;
    maxncols = MIN (n_col, maxncols) ;
    maxnrc = MAX (maxnrows, maxncols) ;

    DEBUG0 (("work alloc:  maxnrows+nb "ID" maxncols+nb "ID"\n",
	maxnrows, maxncols)) ;

    /* 15 allocations, freed in free_work: */
    /* accounted for in UMF_set_stats (work_usage) */
    Work->Wx = (Entry *) UMF_malloc (maxnrows + 1, sizeof (Entry)) ;
    Work->Wy = (Entry *) UMF_malloc (maxnrows + 1, sizeof (Entry)) ;
    Work->Frpos    = (Int *) UMF_malloc (n_row + 1, sizeof (Int)) ;
    Work->Lpattern = (Int *) UMF_malloc (n_row + 1, sizeof (Int)) ;
    Work->Fcpos = (Int *) UMF_malloc (n_col + 1, sizeof (Int)) ;
    Work->Wp = (Int *) UMF_malloc (nn + 1, sizeof (Int)) ;
    Work->Wrp = (Int *) UMF_malloc (MAX (n_col,maxnrows) + 1, sizeof (Int)) ;
    Work->Frows = (Int *) UMF_malloc (maxnrows + 1, sizeof (Int)) ;
    Work->Wm    = (Int *) UMF_malloc (maxnrows + 1, sizeof (Int)) ;
    Work->Fcols = (Int *) UMF_malloc (maxncols + 1, sizeof (Int)) ;
    Work->Wio   = (Int *) UMF_malloc (maxncols + 1, sizeof (Int)) ;
    Work->Woi   = (Int *) UMF_malloc (maxncols + 1, sizeof (Int)) ;
    Work->Woo = (Int *) UMF_malloc (maxnrc + 1, sizeof (Int));
    Work->elen = (n_col - n1) + (n_row - n1) + MIN (n_col-n1, n_row-n1) + 1 ;
    Work->E = (Int *) UMF_malloc (Work->elen, sizeof (Int)) ;
    Work->Front_new1strow = (Int *) UMF_malloc (nfr + 1, sizeof (Int)) ;

    ok = (Work->Frpos && Work->Fcpos && Work->Lpattern
	&& Work->Wp && Work->Wrp && Work->Frows && Work->Fcols
	&& Work->Wio && Work->Woi && Work->Woo && Work->Wm
	&& Work->E && Work->Front_new1strow && Work->Wx && Work->Wy) ;

    /* 2 allocations: accounted for in UMF_set_stats (work_usage) */
    if (Symbolic->prefer_diagonal)
    {
	Work->Diagonal_map  = (Int *) UMF_malloc (nn, sizeof (Int)) ;
	Work->Diagonal_imap = (Int *) UMF_malloc (nn, sizeof (Int)) ;
	ok = ok && Work->Diagonal_map && Work->Diagonal_imap ;
    }
    else
    {
	/* no diagonal map needed for rectangular matrices */
	Work->Diagonal_map = (Int *) NULL ;
	Work->Diagonal_imap = (Int *) NULL ;
    }

    /* 1 allocation, may become part of Numeric (if singular or rectangular): */
    Work->Upattern = (Int *) UMF_malloc (n_col + 1, sizeof (Int)) ;
    ok = ok && Work->Upattern ;

    /* current frontal matrix does not yet exist */
    Work->Flublock = (Entry *) NULL ;
    Work->Flblock  = (Entry *) NULL ;
    Work->Fublock  = (Entry *) NULL ;
    Work->Fcblock  = (Entry *) NULL ;

    DEBUG0 (("work alloc done.\n")) ;
    return (ok) ;
}


/* ========================================================================== */
/* === free_work ============================================================ */
/* ========================================================================== */

PRIVATE void free_work
(
    WorkType *Work
)
{
    DEBUG0 (("work free:\n")) ;
    if (Work)
    {
	/* these 16 objects do exist */
	Work->Wx = (Entry *) UMF_free ((void *) Work->Wx) ;
	Work->Wy = (Entry *) UMF_free ((void *) Work->Wy) ;
	Work->Frpos = (Int *) UMF_free ((void *) Work->Frpos) ;
	Work->Fcpos = (Int *) UMF_free ((void *) Work->Fcpos) ;
	Work->Lpattern = (Int *) UMF_free ((void *) Work->Lpattern) ;
	Work->Upattern = (Int *) UMF_free ((void *) Work->Upattern) ;
	Work->Wp = (Int *) UMF_free ((void *) Work->Wp) ;
	Work->Wrp = (Int *) UMF_free ((void *) Work->Wrp) ;
	Work->Frows = (Int *) UMF_free ((void *) Work->Frows) ;
	Work->Fcols = (Int *) UMF_free ((void *) Work->Fcols) ;
	Work->Wio = (Int *) UMF_free ((void *) Work->Wio) ;
	Work->Woi = (Int *) UMF_free ((void *) Work->Woi) ;
	Work->Woo = (Int *) UMF_free ((void *) Work->Woo) ;
	Work->Wm = (Int *) UMF_free ((void *) Work->Wm) ;
	Work->E = (Int *) UMF_free ((void *) Work->E) ;
	Work->Front_new1strow =
	    (Int *) UMF_free ((void *) Work->Front_new1strow) ;

	/* these objects might not exist */
	Work->Diagonal_map = (Int *) UMF_free ((void *) Work->Diagonal_map) ;
	Work->Diagonal_imap = (Int *) UMF_free ((void *) Work->Diagonal_imap) ;
    }
    DEBUG0 (("work free done.\n")) ;
}


/* ========================================================================== */
/* === error ================================================================ */
/* ========================================================================== */

/* Error return from UMFPACK_numeric.  Free all allocated memory. */

PRIVATE void error
(
    NumericType **Numeric,
    WorkType *Work
)
{
    free_work (Work) ;
    UMFPACK_free_numeric ((void **) Numeric) ;
    ASSERT (UMF_malloc_count == init_count) ;
}
