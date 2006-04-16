/* ========================================================================== */
/* === UMF_kernel_init ====================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Initialize the kernel: scale the matrix, load the initial elements, and
    build the tuple lists.

    Returns TRUE if successful, FALSE if out of memory or if the pattern has
    changed since UMFPACK_*symbolic.  UMFPACK_numeric allocates at least enough
    space for UMF_kernel_init to succeed; otherwise it does not call
    UMF_kernel_init.  So an out-of-memory condition means that the pattern must
    have gotten larger.
*/

#include "umf_internal.h"
#include "umf_tuple_lengths.h"
#include "umf_build_tuples.h"
#include "umf_mem_init_memoryspace.h"
#include "umf_mem_alloc_element.h"
#include "umf_mem_alloc_head_block.h"
#include "umf_mem_alloc_tail_block.h"
#include "umf_mem_free_tail_block.h"
#include "umf_free.h"
#include "umf_scale.h"

/* ========================================================================== */
/* === packsp =============================================================== */
/* ========================================================================== */

/* remove zero entries from a singleton column of L or a row of U */

PRIVATE Int packsp	/* returns new value of pnew */
(
    Int pnew,		/* index into Memory of next free space */
    Int *p_p,		/* ptr to index of old pattern in Memory on input,
			   new pattern on output */
    Int *p_len,		/* ptr to length of old pattern on input,
			   new pattern on output */
    Unit *Memory	/* contains the sparse vector on input and output */
)
{
    Entry x, *Bx, *Bx2 ;
    Int p, i, len, len_new, *Bi, *Bi2 ;
    double s ;

    /* get the pointers to the sparse vector, and its length */
    p = *p_p ;
    len = *p_len ;
    Bi = (Int   *) (Memory + p) ; p += UNITS (Int,   len) ;
    Bx = (Entry *) (Memory + p) ; p += UNITS (Entry, len) ;
    DEBUGm4 (("  p "ID" len "ID" pnew "ID"\n", p, len, pnew)) ;

    /* the vector resides in Bi [0..len-1] and Bx [0..len-1] */

    /* first, compact the vector in place */
    len_new = 0 ;
    for (p = 0 ; p < len ; p++)
    {
	i = Bi [p] ;
	x = Bx [p] ;
	DEBUGm4 (("    old vector: i "ID" value: ", i)) ;
	EDEBUGk (-4, x) ;
	DEBUGm4 (("\n")) ;
	ASSERT (i >= 0) ;
	/* skip if zero */
	if (IS_ZERO (x)) continue ;
	/* store the value back into the vector */
	if (len_new != p)
	{
	    Bi [len_new] = i ;
	    Bx [len_new] = x ;
	}
	len_new++ ;
    }
    ASSERT (len_new <= len) ;

    /* the vector is now in Bi [0..len_new-1] and Bx [0..len_new-1] */

#ifndef NDEBUG
    for (p = 0 ; p < len_new ; p++)
    {
	DEBUGm4 (("    new vector: i "ID" value: ", Bi [p])) ;
	EDEBUGk (-4, Bx [p]) ;
	DEBUGm4 (("\n")) ;
	ASSERT (Bi [p] >= 0) ;
    }
#endif

    /* allocate new space for the compacted vector */
    *p_p = pnew ;
    *p_len = len_new ;
    Bi2 = (Int   *) (Memory + pnew) ; pnew += UNITS (Int,   len_new) ;
    Bx2 = (Entry *) (Memory + pnew) ; pnew += UNITS (Entry, len_new) ;
    DEBUGm4 (("  pnew "ID" len_new "ID"\n", pnew, len_new)) ;

    /* shift the vector upwards, into its new space */
    for (p = 0 ; p < len_new ; p++)
    {
	Bi2 [p] = Bi [p] ;
    }
    for (p = 0 ; p < len_new ; p++)
    {
	Bx2 [p] = Bx [p] ;
    }

#ifndef NDEBUG
    for (p = 0 ; p < len_new ; p++)
    {
	DEBUGm4 (("    packed vec: i "ID" value: ", Bi2 [p])) ;
	EDEBUGk (-4, Bx2 [p]) ;
	DEBUGm4 (("\n")) ;
	ASSERT (Bi2 [p] >= 0) ;
    }
#endif

    /* return the pointer to the space just after the new vector */
    return (pnew) ;
}


/* ========================================================================== */
/* === UMF_kernel_init ====================================================== */
/* ========================================================================== */

GLOBAL Int UMF_kernel_init
(
    const Int Ap [ ],		/* user's input matrix (not modified) */
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int row, k, oldcol, size, e, p1, p2, p, nz, *Rows, *Cols, *E, i, *Upos,
	*Lpos, n_row, n_col, *Wp, *Cperm_init, *Frpos, *Fcpos, *Row_degree, nn,
	*Row_tlen, *Col_degree, *Col_tlen, oldrow, newrow, ilast, *Wrp,
	*Rperm_init, col, n_inner, prefer_diagonal, *Diagonal_map, nempty,
	*Diagonal_imap, fixQ, rdeg, cdeg, nempty_col, *Esize, esize, pnew,
	*Lip, *Uip, *Lilen, *Uilen, llen, pa, *Cdeg, *Rdeg, n1, clen, do_scale,
	lnz, unz, lip, uip, k1, *Rperm, *Cperm, pivcol, *Li, lilen,
	**Rpi, nempty_row, dense_row_threshold, empty_elements, rpi, rpx ;
    double unused = 0, *Rs, rsmin, rsmax, rs ;
    Entry *D, *C, x, *Lval, pivot_value, **Rpx ;
    Element *ep ;
    Unit *Memory ;

#ifndef NRECIPROCAL
    Int do_recip = FALSE ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    DEBUG0 (("KERNEL INIT\n")) ;

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;
    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;
    nempty_col = Symbolic->nempty_col ;
    nempty_row = Symbolic->nempty_row ;
    nempty = MIN (nempty_row, nempty_col) ;
    ASSERT (n_row > 0 && n_col > 0) ;
    Cperm_init = Symbolic->Cperm_init ;
    Rperm_init = Symbolic->Rperm_init ;
    Cdeg = Symbolic->Cdeg ;
    Rdeg = Symbolic->Rdeg ;
    n1 = Symbolic->n1 ;
    dense_row_threshold = Symbolic->dense_row_threshold ;
    DEBUG0 (("Singletons: "ID"\n", n1)) ;
    Work->nforced = 0 ;
    Work->ndiscard = 0 ;
    Work->noff_diagonal = 0 ;

    nz = Ap [n_col] ;
    if (nz < 0 || Ap [0] != 0 || nz != Symbolic->nz)
    {
	DEBUGm4 (("nz or Ap [0] bad\n")) ;
	return (FALSE) ;	/* pattern changed */
    }

    prefer_diagonal = Symbolic->prefer_diagonal ;
    Diagonal_map = Work->Diagonal_map ;
    Diagonal_imap = Work->Diagonal_imap ;

    /* ---------------------------------------------------------------------- */
    /* initialize the Numeric->Memory space for LU, elements, and tuples */
    /* ---------------------------------------------------------------------- */

    UMF_mem_init_memoryspace (Numeric) ;
    DEBUG1 (("Kernel init head usage, before allocs: "ID"\n", Numeric->ihead)) ;

    /* ---------------------------------------------------------------------- */
    /* initialize the Work and Numeric objects */
    /* ---------------------------------------------------------------------- */

    /* current front is empty */
    Work->fnpiv = 0 ;
    Work->fncols = 0 ;
    Work->fnrows = 0 ;
    Work->fncols_max = 0 ;
    Work->fnrows_max = 0 ;
    Work->fnzeros = 0 ;
    Work->fcurr_size = 0 ;
    Work->fnr_curr = 0 ;
    Work->fnc_curr = 0 ;

    Work->nz = nz ;
    Work->prior_element = EMPTY ;
    Work->ulen = 0 ;
    Work->llen = 0 ;
    Work->npiv = n1 ;
    Work->frontid = 0 ;
    Work->nextcand = n1 ;

    Memory = Numeric->Memory ;
    Rperm = Numeric->Rperm ;
    Cperm = Numeric->Cperm ;
    Row_degree = Numeric->Rperm ;
    Col_degree = Numeric->Cperm ;
    /* Row_tuples = Numeric->Uip ; not needed */
    Row_tlen   = Numeric->Uilen ;
    /* Col_tuples = Numeric->Lip ; not needed */
    Col_tlen   = Numeric->Lilen ;

    Lip = Numeric->Lip ;
    Uip = Numeric->Uip ;
    Lilen = Numeric->Lilen ;
    Uilen = Numeric->Uilen ;

    Frpos = Work->Frpos ;
    Fcpos = Work->Fcpos ;
    Wp = Work->Wp ;
    Wrp = Work->Wrp ;

    D = Numeric->D ;
    Upos = Numeric->Upos ;
    Lpos = Numeric->Lpos ;
    for (k = 0 ; k < n_inner ; k++)
    {
	CLEAR (D [k]) ;
    }

    Rs = Numeric->Rs ;

    for (row = 0 ; row <= n_row ; row++)
    {
	Lpos [row] = EMPTY ;
	/* Row_tuples [row] = 0 ; set in UMF_build_tuples */
	/* Row_degree [row] = 0 ; initialized below */
	Row_tlen [row] = 0 ;
	/* Frpos [row] = EMPTY ;  do this later */
    }

    for (col = 0 ; col <= n_col ; col++)
    {
	Upos [col] = EMPTY ;
	/* Col_tuples [col] = 0 ; set in UMF_build_tuples */
	/* Col_degree [col] = 0 ; initialized below */
	Col_tlen [col] = 0 ;
	Fcpos [col] = EMPTY ;
	Wrp [col] = 0 ;
    }
    Work->Wrpflag = 1 ;

    /* When cleared, Wp [0..nn] is < 0 */
    for (i = 0 ; i <= nn ; i++)
    {
	Wp [i] = EMPTY ;
    }
    /* In col search, Wp [row] is set to a position, which is >= 0. */

    /* When cleared, Wrp [0..n_col] is < Wrpflag */
    /* In row search, Wrp [col] is set to Wrpflag. */

    /* no need to initialize Wm, Wio, Woi, and Woo */

    /* clear the external degree counters */
    Work->cdeg0 = 1 ;
    Work->rdeg0 = 1 ;

    fixQ = Symbolic->fixQ ;

    E = Work->E ;

    Numeric->n_row = n_row ;
    Numeric->n_col = n_col ;
    Numeric->npiv = 0 ;
    Numeric->nnzpiv = 0 ;
    Numeric->min_udiag = 0.0 ;
    Numeric->max_udiag = 0.0 ;
    Numeric->rcond = 0.0 ;
    Numeric->isize = 0 ;
    Numeric->nLentries = 0 ;
    Numeric->nUentries = 0 ;
    Numeric->lnz = 0 ;
    Numeric->unz = 0 ;
    Numeric->maxfrsize = 0 ;
    Numeric->maxnrows = 0 ;
    Numeric->maxncols = 0 ;
    Numeric->flops = 0. ;
    Numeric->n1 = n1 ;

    /* ---------------------------------------------------------------------- */
    /* compute the scale factors, if requested, and check the input matrix */
    /* ---------------------------------------------------------------------- */

    /* UMFPACK_SCALE_SUM: Rs [i] = sum of the absolute values in row i.
     * UMFPACK_SCALE_MAX: Rs [i] = max of the absolute values in row i.
     *
     * If A is complex, an approximate abs is used (|xreal| + |ximag|).
     *
     * If min (Rs [0..n_row]) >= RECIPROCAL_TOLERANCE, then the scale
     * factors are inverted, and the rows of A are multiplied by the scale
     * factors.  Otherwise, the rows are divided by the scale factors.  If
     * NRECIPROCAL is defined, then the rows are always divided by the scale
     * factors.
     *
     * For MATLAB (either built-in routine or mexFunction), or for gcc,
     * the rows are always divided by the scale factors.
     */

    do_scale = (Numeric->scale != UMFPACK_SCALE_NONE) ;

    if (do_scale)
    {
	int do_max = Numeric->scale == UMFPACK_SCALE_MAX ;
	for (row = 0 ; row < n_row ; row++)
	{
	    Rs [row] = 0.0 ;
	}
	for (col = 0 ; col < n_col ; col++)
	{
	    ilast = EMPTY ;
	    p1 = Ap [col] ;
	    p2 = Ap [col+1] ;
	    if (p1 > p2)
	    {
		/* invalid matrix */
		DEBUGm4 (("invalid matrix (Ap)\n")) ;
		return (FALSE) ;
	    }
	    for (p = p1 ; p < p2 ; p++)
	    {
		double value ;
		Entry aij ;
		row = Ai [p] ;
		if (row <= ilast || row >= n_row)
		{
		    /* invalid matrix, columns must be sorted, no duplicates */
		    DEBUGm4 (("invalid matrix (Ai)\n")) ;
		    return (FALSE) ;
		}
		ASSIGN (aij, Ax [p], Az [p]) ;
		APPROX_ABS (value, aij) ;
		rs = Rs [row] ;
		if (!SCALAR_IS_NAN (rs))
		{
		    if (SCALAR_IS_NAN (value))
		    {
			/* if any entry in the row is NaN, then the scale factor
			 * is NaN too (for now) and then set to 1.0 below */
			Rs [row] = value ;
		    }
		    else if (do_max)
		    {
			Rs [row] = MAX (rs, value) ;
		    }
		    else
		    {
			Rs [row] += value ;
		    }
		}
		DEBUG4 (("i "ID" j "ID" value %g,  Rs[i]: %g\n",
		    row, col, value, Rs[row])) ;
		ilast = row ;
	    }
	}
	DEBUG2 (("Rs[0] = %30.20e\n", Rs [0])) ;
	for (row = 0 ; row < n_row ; row++)
	{
	    rs = Rs [row] ;
	    if (SCALAR_IS_ZERO (rs) || SCALAR_IS_NAN (rs))
	    {
		/* don't scale a completely zero row, or one with NaN's */
		Rs [row] = 1.0 ;
	    }
	}
	rsmin = Rs [0] ;
	rsmax = Rs [0] ;
	for (row = 0 ; row < n_row ; row++)
	{
	    DEBUG2 (("sum %30.20e ", Rs [row])) ;
	    rsmin = MIN (rsmin, Rs [row]) ;
	    rsmax = MAX (rsmax, Rs [row]) ;
	    DEBUG2 (("Rs["ID"] = %30.20e\n", row, Rs [row])) ;
	}
#ifndef NRECIPROCAL
	/* multiply by the reciprocal if Rs is not too small */
	do_recip = (rsmin >= RECIPROCAL_TOLERANCE) ;
	if (do_recip)
	{
	    /* invert the scale factors */
	    for (row = 0 ; row < n_row ; row++)
	    {
		Rs [row] = 1.0 / Rs [row] ;
	    }
	}
#endif
    }
    else
    {
	/* no scaling, rsmin and rsmax not computed */
	rsmin = -1 ;
	rsmax = -1 ;
#ifndef NRECIPROCAL
	do_recip = FALSE ;
#endif
	/* check the input matrix */
	if (!AMD_valid (n_row, n_col, Ap, Ai))
	{
	    /* matrix is invalid */
	    return (FALSE) ;
	}
    }

    Numeric->rsmin = rsmin ;
    Numeric->rsmax = rsmax ;
#ifndef NRECIPROCAL
    Numeric->do_recip = do_recip ;
#else
    Numeric->do_recip = FALSE ;
#endif

    /* ---------------------------------------------------------------------- */
    /* construct the inverse row Rperm_init permutation (use Frpos as temp) */
    /* ---------------------------------------------------------------------- */

    DEBUG3 (("\n\n===================LOAD_MATRIX:\n")) ;

    for (newrow = 0 ; newrow < n_row ; newrow++)
    {
	oldrow = Rperm_init [newrow] ;
	ASSERT (oldrow >= 0 && oldrow < n_row) ;
	Frpos [oldrow] = newrow ;
    }

    /* ---------------------------------------------------------------------- */
    /* construct the diagonal imap if doing symmetric pivoting */
    /* ---------------------------------------------------------------------- */

    if (prefer_diagonal)
    {
	ASSERT (n_row == n_col) ;
	ASSERT (nempty_col == Symbolic->nempty_row) ;
	ASSERT (nempty_col == nempty) ;
	for (i = 0 ; i < nn ; i++)
	{
	    Diagonal_map [i] = EMPTY ;
	    Diagonal_imap [i] = EMPTY ;
	}
	for (k = n1 ; k < nn - nempty ; k++)
	{
	    newrow = Symbolic->Diagonal_map [k] ;
	    Diagonal_map [k] = newrow ;
	    Diagonal_imap [newrow] = k ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* allocate O (n_row) workspace at the tail end of Memory */
    /* ---------------------------------------------------------------------- */

    rpi = UMF_mem_alloc_tail_block (Numeric, UNITS (Int *, n_row+1)) ;
    rpx = UMF_mem_alloc_tail_block (Numeric, UNITS (Entry *, n_row+1)) ;
    if (!rpi || !rpx)
    {
	/* :: pattern change (out of memory for Rpx, Rpx) :: */
	/* out of memory, which can only mean that the pattern has changed */
	return (FALSE) ;	/* pattern changed */
    }
    Rpi = (Int   **) (Memory + rpx) ;
    Rpx = (Entry **) (Memory + rpi) ;

    /* ---------------------------------------------------------------------- */
    /* allocate the LU factors for the columns of the singletons */
    /* ---------------------------------------------------------------------- */

    DEBUG1 (("Allocating singletons:\n")) ;
    for (k = 0 ; k < n1 ; k++)
    {
	lnz = Cdeg [k] - 1 ;
	unz = Rdeg [k] - 1 ;

	DEBUG1 (("Singleton k "ID" pivrow "ID" pivcol "ID" cdeg "ID" rdeg "
	    ID"\n", k, Rperm_init [k], Cperm_init [k], Cdeg [k], Rdeg [k])) ;
	ASSERT (unz >= 0 && lnz >= 0 && (lnz == 0 || unz == 0)) ;
	DEBUG1 (("   lnz "ID" unz "ID"\n", lnz, unz)) ;

	size = UNITS (Int, lnz) + UNITS (Entry, lnz)
	     + UNITS (Int, unz) + UNITS (Entry, unz) ;
	p = UMF_mem_alloc_head_block (Numeric, size) ;
	DEBUG1 (("Kernel init head usage: "ID"\n", Numeric->ihead)) ;
	if (!p)
	{
	    /* :: pattern change (out of memory for singletons) :: */
	    DEBUG0 (("Pattern has gotten larger - kernel init failed\n")) ;
	    return (FALSE) ;	/* pattern changed */
	}

	/* allocate the column of L */
	lip = p ;
	p += UNITS (Int, lnz) ;
	p += UNITS (Entry, lnz) ;

	/* allocate the row of U */
	uip = p ;
	Rpi [k] = (Int *) (Memory + p) ;
	p += UNITS (Int, unz) ;
	Rpx [k] = (Entry *) (Memory + p) ;
	/* p += UNITS (Entry, unz) ; (not needed) */

	/* a single column of L (no Lchains) */
	Lip [k] = lip ;
	Lilen [k] = lnz ;

	/* a single row of L (no Uchains) */
	Uip [k] = uip ;
	Uilen [k] = unz ;

	Wp [k] = unz ;

	/* save row and column inverse permutation */
	k1 = ONES_COMPLEMENT (k) ;
	Rperm [k] = k1 ;			/* aliased with Row_degree */
	Cperm [k] = k1 ;			/* aliased with Col_degree */
    }

    /* ---------------------------------------------------------------------- */
    /* current frontal matrix is empty */
    /* ---------------------------------------------------------------------- */

    e = 0 ;
    E [e] = 0 ;
    Work->Flublock = (Entry *) NULL ;
    Work->Flblock  = (Entry *) NULL ;
    Work->Fublock  = (Entry *) NULL ;
    Work->Fcblock  = (Entry *) NULL ;

    /* ---------------------------------------------------------------------- */
    /* allocate the column elements */
    /* ---------------------------------------------------------------------- */

    Esize = Symbolic->Esize ;
    empty_elements = FALSE  ;
    for (k = n1 ; k < n_col - nempty_col ; k++)
    {
	e = k - n1 + 1 ;
	ASSERT (e < Work->elen) ;
	esize = Esize ? Esize [k-n1] : Cdeg [k] ;
	if (esize > 0)
	{
	    /* allocate an element for this column */
	    E [e] = UMF_mem_alloc_element (Numeric, esize, 1, &Rows, &Cols, &C,
		&size, &ep) ;
	    if (E [e] <= 0)
	    {
		/* :: pattern change (out of memory for column elements) :: */
		return (FALSE) ;	/* pattern has changed */
	    }
	    Cols [0] = k ;
	    DEBUG0 (("Got column element e "ID" esize "ID"\n", e, esize)) ;
	}
	else
	{
	    /* all rows in this column are dense, or empty */
	    E [e] = 0 ;
	    empty_elements = TRUE  ;
	    DEBUG0 (("column element e is empty "ID"\n", e)) ;
	}
    }
    DEBUG0 (("e "ID" n_col "ID" nempty_col "ID" n1 "ID"\n", e, n_col,
		nempty_col, n1)) ;
    ASSERT (e == n_col - nempty_col - n1) ;

    /* ---------------------------------------------------------------------- */
    /* allocate the row elements for dense rows of A (if any) */
    /* ---------------------------------------------------------------------- */

    if (Esize)
    {
	for (k = n1 ; k < n_row - nempty_row ; k++)
	{
	    rdeg = Rdeg [k] ;
	    if (rdeg > dense_row_threshold)
	    {
		/* allocate an element for this dense row */
		e++ ;
		ASSERT (e < Work->elen) ;
		E [e] = UMF_mem_alloc_element (Numeric, 1, rdeg, &Rows, &Cols,
		    &C, &size, &ep) ;
		if (E [e] <= 0)
		{
		    /* :: pattern change (out of memory for row elements) :: */
		    return (FALSE) ;	/* pattern has changed */
		}
		Rows [0] = k ;
		Rpi [k] = Cols ;
		Rpx [k] = C ;
		Wp [k] = rdeg ;
		DEBUG0 (("Got row element e "ID" rdeg "ID"\n", e, rdeg)) ;
	    }
	}
    }

    /* elements are currently in the range 0 to e */
    Work->nel = e ;

    /* ---------------------------------------------------------------------- */
    /* create the first n1 columns of L and U */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < n1 ; k++)
    {
	pivcol = Cperm_init [k] ;
	p2 = Ap [pivcol+1] ;

	/* get the kth column of L */
	p = Lip [k] ;
	Li = (Int *) (Memory + p) ;
	lilen = Lilen [k] ;
	p += UNITS (Int, lilen) ;
	Lval = (Entry *) (Memory + p) ;

	llen = 0 ;

	for (pa = Ap [pivcol] ; pa < p2 ; pa++)
	{
	    oldrow = Ai [pa] ;
	    newrow = Frpos [oldrow] ;
	    ASSIGN (x, Ax [pa], Az [pa]) ;

	    /* scale the value using the scale factors, Rs */
	    if (do_scale)
	    {
#ifndef NRECIPROCAL
		if (do_recip)
		{
		    SCALE_RECIP (x, Rs [oldrow]) ;
		}
		else
#endif
		{
		    SCALE_DIV (x, Rs [oldrow]) ;
		}
	    }

	    if (newrow == k)
	    {
		/* this is the pivot entry itself */
		ASSERT (oldrow == Rperm_init [k]) ;
		D [k] = x ;
	    }
	    else if (newrow < k)
	    {
		/* this entry goes in a row of U */
		DEBUG1 (("Singleton row of U: k "ID" newrow "ID"\n",
		    k, newrow)) ;
		if (--(Wp [newrow]) < 0)
		{
		    /* :: pattern change (singleton row too long) :: */
		    DEBUGm4 (("bad U singleton row (too long)\n")) ;
		    return (FALSE) ;	/* pattern changed */
		}
		*(Rpi [newrow]++) = k ;
		*(Rpx [newrow]++) = x ;
	    }
	    else
	    {
		/* this entry goes in a column of L */
		DEBUG1 (("Singleton col of L: k "ID" newrow "ID"\n",
		    k, newrow)) ;
		if (llen >= lilen)
		{
		    DEBUGm4 (("bad L singleton col (too long)\n")) ;
		    return (FALSE) ;	/* pattern changed */
		}
		Li   [llen] = newrow ;
		Lval [llen] = x ;
		llen++ ;
	    }
	}

	if (llen != lilen)
	{
	    /* :: pattern change (singleton column too long) :: */
	    DEBUGm4 (("bad L singleton col (too short)\n")) ;
	    return (FALSE) ;	/* pattern changed */
	}

	/* scale the column of L */
	if (llen > 0)
	{
	    pivot_value = D [k] ;
	    UMF_scale (llen, pivot_value, Lval) ;
	}

    }

    /* ---------------------------------------------------------------------- */
    /* allocate the elements and copy the columns of A */
    /* ---------------------------------------------------------------------- */

    /* also apply the row and column pre-ordering.  */
    for (k = n1 ; k < n_col ; k++)
    {
	/* The newcol is k, which is what the name of the column is in the
	 * UMFPACK kernel.  The user's name for the column is oldcol. */
	oldcol = Cperm_init [k] ;

	ASSERT (oldcol >= 0 && oldcol < n_col) ;

	p2 = Ap [oldcol+1] ;

	cdeg = Cdeg [k] ;
	ASSERT (cdeg >= 0) ;
	ASSERT (IMPLIES (
	    (Symbolic->ordering != UMFPACK_ORDERING_GIVEN) && n1 > 0,
	    cdeg > 1 || cdeg == 0)) ;

	/* if fixQ: set Col_degree to 0 for the NON_PIVOTAL_COL macro */
	Col_degree [k] = fixQ ? 0 : cdeg ;

	/* get the element for this column (if any) */
	e = k - n1 + 1 ;
	if (k < n_col - nempty_col)
	{
	    esize = Esize ? Esize [k-n1] : cdeg ;
	    if (E [e])
	    {
		Int ncols, nrows ;
		Unit *pp ;
		pp = Memory + E [e] ;
		GET_ELEMENT (ep, pp, Cols, Rows, ncols, nrows, C) ;
		ASSERT (ncols == 1) ;
		ASSERT (nrows == esize) ;
		ASSERT (Cols [0] == k) ;
	    }
	}
	else
	{
	    ASSERT (cdeg == 0) ;
	    esize = 0 ;
	}

	clen = 0 ;

	for (pa = Ap [oldcol] ; pa < p2 ; pa++)
	{
	    oldrow = Ai [pa] ;
	    newrow = Frpos [oldrow] ;
	    ASSIGN (x, Ax [pa], Az [pa]) ;

	    /* scale the value using the scale factors, Rs */
	    if (do_scale)
	    {
#ifndef NRECIPROCAL
		if (do_recip)
		{
		    /* multiply by the reciprocal */
		    SCALE_RECIP (x, Rs [oldrow]) ;
		}
		else
#endif
		{
		    /* divide instead */
		    SCALE_DIV (x, Rs [oldrow]) ;
		}
	    }

	    rdeg = Rdeg [newrow] ;
	    if (newrow < n1 || rdeg > dense_row_threshold)
	    {
		/* this entry goes in a row of U or into a dense row */
		DEBUG1 (("Singleton/dense row of U: k "ID" newrow "ID"\n",
		    k, newrow)) ;
		if (--(Wp [newrow]) < 0)
		{
		    DEBUGm4 (("bad row of U or A (too long)\n")) ;
		    return (FALSE) ;	/* pattern changed */
		}
		*(Rpi [newrow]++) = k ;
		*(Rpx [newrow]++) = x ;
	    }
	    else
	    {
		/* this entry goes in an initial element */
		DEBUG1 (("In element k "ID" e "ID" newrow "ID"\n",
		    k, e, newrow)) ;
		if (clen >= esize)
		{
		    DEBUGm4 (("bad A column (too long)\n")) ;
		    return (FALSE) ;	/* pattern changed */
		}
		ASSERT (E [e]) ;
		ASSERT (k < n_col - nempty_col) ;
		Rows [clen] = newrow ;
		C    [clen] = x ;
		clen++ ;
#ifndef NDEBUG
		if (Diagonal_map && (newrow == Diagonal_map [k]))
		{
		    DEBUG0 (("Diagonal: old: row "ID" col "ID" : "
			"new: row "ID" col "ID" : ",
			oldrow, oldcol, newrow, k)) ;
		    EDEBUGk (0, x) ;
		}
#endif
	    }
	}

	if (clen != esize)
	{
	    /* :: pattern change (singleton column too short) :: */
	    DEBUGm4 (("bad A column (too short)\n")) ;
	    return (FALSE) ;	/* pattern changed */
	}
    }

    /* ---------------------------------------------------------------------- */
    /* free the Rpi and Rpx workspace at the tail end of memory */
    /* ---------------------------------------------------------------------- */

    UMF_mem_free_tail_block (Numeric, rpi) ;
    UMF_mem_free_tail_block (Numeric, rpx) ;

    /* ---------------------------------------------------------------------- */
    /* prune zeros from the singleton rows and columns */
    /* ---------------------------------------------------------------------- */

    if (n1 > 0)
    {
	pnew = Lip [0] ;
	ASSERT (pnew == 1) ;
	for (k = 0 ; k < n1 ; k++)
	{
	    DEBUGm4 (("\nPrune singleton L col "ID"\n", k)) ;
	    pnew = packsp (pnew, &Lip [k], &Lilen [k], Memory) ;
	    Numeric->lnz += Lilen [k] ;
	    DEBUGm4 (("\nPrune singleton U row "ID"\n", k)) ;
	    pnew = packsp (pnew, &Uip [k], &Uilen [k], Memory) ;
	    Numeric->unz += Uilen [k] ;
	}
	/* free the unused space at the head of memory */
	Numeric->ihead = pnew ;
    }

    /* ---------------------------------------------------------------------- */
    /* initialize row degrees */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < n1 ; k++)
    {
	if (Wp [k] != 0)
	{
	    /* :: pattern change (singleton row too short) :: */
	    DEBUGm4 (("bad U singleton row (too short)\n")) ;
	    return (FALSE) ;	/* pattern changed */
	}
    }

    for (k = n1 ; k < n_row ; k++)
    {
	DEBUG1 (("Initial row degree k "ID" oldrow "ID" Rdeg "ID"\n",
	    k, Rperm_init [k], Rdeg [k])) ;
	rdeg = Rdeg [k] ;
	Row_degree [k] = rdeg ;
	if (rdeg > dense_row_threshold && Wp [k] != 0)
	{
	    /* :: pattern change (dense row too short) :: */
	    DEBUGm4 (("bad dense row (too short)\n")) ;
	    return (FALSE) ;
	}
    }

#ifndef NDEBUG
    if (prefer_diagonal)
    {
	Int *InvCperm, newcol ;
	Entry aij ;
	UMF_dump_diagonal_map (Diagonal_map, Diagonal_imap, n1, nn, nempty) ;
	InvCperm = (Int *) malloc (n_col * sizeof (Int)) ;
	ASSERT (InvCperm != (Int *) NULL) ;
	for (newcol = 0 ; newcol < n_col ; newcol++)
	{
	    oldcol = Cperm_init [newcol] ;
	    InvCperm [oldcol] = newcol ;
	}
	DEBUGm3 (("Diagonal of P2*A:\n")) ;
	for (oldcol = 0 ; oldcol < n_col ; oldcol++)
	{
	    newcol = InvCperm [oldcol] ;
	    for (p = Ap [oldcol] ; p < Ap [oldcol+1] ; p++)
	    {
		oldrow = Ai [p] ;
		newrow = Frpos [oldrow] ;
		ASSIGN (aij, Ax [p], Az [p]) ;
		if (newrow == Diagonal_map [newcol])
		{
		    DEBUG0 (("old row "ID" col "ID" new row "ID" col "ID,
			oldrow, oldcol, newrow, newcol)) ;
		    EDEBUGk (0, aij) ;
		    DEBUG0 ((" scaled ")) ;
		    if (do_scale)
		    {
#ifndef NRECIPROCAL
			if (do_recip)
			{
			    SCALE_RECIP (aij, Rs [oldrow]) ;
			}
			else
#endif
			{
			    SCALE_DIV (aij, Rs [oldrow]) ;
			}
		    }
		    EDEBUGk (0, aij) ;
		    DEBUG0 (("\n")) ;
		}
	    }
	}
	free (InvCperm) ;
    }
#endif

    Col_degree [n_col] = 0 ;

    /* ---------------------------------------------------------------------- */
    /* pack the element name space */
    /* ---------------------------------------------------------------------- */

    if (empty_elements)
    {
	Int e2 = 0 ;
	DEBUG0 (("\n\n============= Packing element space\n")) ;
	for (e = 1 ; e <= Work->nel ; e++)
	{
	    if (E [e])
	    {
		e2++ ;
		E [e2] = E [e] ;
	    }
	}
	Work->nel = e2 ;
    }

#ifndef NDEBUG
    DEBUG0 (("Number of initial elements: "ID"\n", Work->nel)) ;
    for (e = 0 ; e <= Work->nel ; e++) UMF_dump_element (Numeric, Work,e,TRUE) ;
#endif

    for (e = Work->nel + 1 ; e < Work->elen ; e++)
    {
	E [e] = 0 ;
    }

    /* Frpos no longer needed */
    for (row = 0 ; row <= n_row ; row++)
    {
	Frpos [row] = EMPTY ;
    }

    /* clear Wp */
    for (i = 0 ; i <= nn ; i++)
    {
	Wp [i] = EMPTY ;
    }

    DEBUG1 (("Kernel init head usage: "ID"\n", Numeric->ihead)) ;

    /* ---------------------------------------------------------------------- */
    /* build the tuple lists */
    /* ---------------------------------------------------------------------- */

    /* if the memory usage changes, then the pattern has changed */

    (void) UMF_tuple_lengths (Numeric, Work, &unused) ;
    if (!UMF_build_tuples (Numeric, Work))
    {
	/* :: pattern change (out of memory in umf_build_tuples) :: */
	/* We ran out of memory, which can only mean that */
	/* the pattern (Ap and or Ai) has changed (gotten larger). */
	DEBUG0 (("Pattern has gotten larger - build tuples failed\n")) ;
	return (FALSE) ;	/* pattern changed */
    }

    Numeric->init_usage = Numeric->max_usage ;

    /* ---------------------------------------------------------------------- */
    /* construct the row merge sets */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i <= Symbolic->nfr ; i++)
    {
	Work->Front_new1strow [i] = Symbolic->Front_1strow [i] ;
    }

#ifndef NDEBUG
    UMF_dump_rowmerge (Numeric, Symbolic, Work) ;
    DEBUG6 (("Column form of original matrix:\n")) ;
    UMF_dump_col_matrix (Ax,
#ifdef COMPLEX
	Az,
#endif
	Ai, Ap, n_row, n_col, nz) ;
    UMF_dump_memory (Numeric) ;
    UMF_dump_matrix (Numeric, Work, FALSE) ;
    DEBUG0 (("kernel init done...\n")) ;
#endif

    return (TRUE) ;
}
