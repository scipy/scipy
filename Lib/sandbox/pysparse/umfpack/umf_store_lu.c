/* ========================================================================== */
/* === UMF_store_lu ========================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Store the LU factors.  Called by the kernel.
    Returns TRUE if successful, FALSE if out of memory.
*/

#include "umf_internal.h"
#include "umf_mem_alloc_head_block.h"
#include "umf_mem_free_tail_block.h"
#include "umf_get_memory.h"

/* ========================================================================== */

GLOBAL Int UMF_store_lu
(
    NumericType *Numeric,
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int i, k, fnr_curr, fnrows, fncols, row, col, pivrow, pivcol, *Frows,
	*Fcols, *Lpattern, *Upattern, *Lpos, *Upos, llen, ulen, fnc_curr, fnpiv,
	uilen, lnz, unz, nb, *Lilen,
	*Uilen, *Lip, *Uip, *Li, *Ui, pivcol_position, newLchain, newUchain,
	pivrow_position, p, size, lip, uip, lnzi, lnzx, unzx, lnz2i, lnz2x,
	unz2i, unz2x, zero_pivot, *Pivrow, *Pivcol, kk,
	Lnz [MAXNB] ;
    Entry *D, pivot_value, *Lval, *Uval, *Fl1, *Fl2, *Fu1, *Fu2,
	*Flublock, *Flblock, *Fublock ;

#ifndef NDEBUG
    Int *Col_degree, *Row_degree ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    fnrows = Work->fnrows ;
    fncols = Work->fncols ;
    fnpiv = Work->fnpiv ;

    Lpos = Numeric->Lpos ;
    Upos = Numeric->Upos ;
    Lilen = Numeric->Lilen ;
    Uilen = Numeric->Uilen ;

    Lip = Numeric->Lip ;
    Uip = Numeric->Uip ;
    D = Numeric->D ;

    Flublock = Work->Flublock ;
    Flblock  = Work->Flblock ;
    Fublock  = Work->Fublock ;

    fnr_curr = Work->fnr_curr ;
    fnc_curr = Work->fnc_curr ;
    Frows = Work->Frows ;
    Fcols = Work->Fcols ;

#ifndef NDEBUG
    Col_degree = Numeric->Cperm ;	/* for NON_PIVOTAL_COL macro */
    Row_degree = Numeric->Rperm ;	/* for NON_PIVOTAL_ROW macro */
#endif

    Lpattern = Work->Lpattern ;
    llen = Work->llen ;
    Upattern = Work->Upattern ;
    ulen = Work->ulen ;

    nb = Work->nb ;

#ifndef NDEBUG
    DEBUG1 (("\n##################################### STORE LU: fnrows "ID
	" fncols "ID"\n", fnrows, fncols)) ;

    DEBUG2 (("\nFrontal matrix, including all space:\n"
		"fnr_curr "ID" fnc_curr "ID" nb    "ID"\n"
		"fnrows   "ID" fncols   "ID" fnpiv "ID"\n",
		fnr_curr, fnc_curr, nb, fnrows, fncols, fnpiv)) ;

    DEBUG2 (("\nJust the active part:\n")) ;
    DEBUG7 (("C  block: ")) ;
    UMF_dump_dense (Work->Fcblock,  fnr_curr, fnrows, fncols) ;
    DEBUG7 (("L  block: ")) ;
    UMF_dump_dense (Work->Flblock,  fnr_curr, fnrows, fnpiv);
    DEBUG7 (("U' block: ")) ;
    UMF_dump_dense (Work->Fublock,  fnc_curr, fncols, fnpiv) ;
    DEBUG7 (("LU block: ")) ;
    UMF_dump_dense (Work->Flublock, nb, fnpiv, fnpiv) ;
    DEBUG7 (("Current frontal matrix: (prior to store LU)\n")) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif

    Pivrow = Work->Pivrow ;
    Pivcol = Work->Pivcol ;

    /* ---------------------------------------------------------------------- */
    /* store the columns of L */
    /* ---------------------------------------------------------------------- */

    for (kk = 0 ; kk < fnpiv ; kk++)
    {

	/* ------------------------------------------------------------------ */
	/* one more pivot row and column is being stored into L and U */
	/* ------------------------------------------------------------------ */

	k = Work->npiv + kk ;

	/* ------------------------------------------------------------------ */
	/* find the kth pivot row and pivot column */
	/* ------------------------------------------------------------------ */

	pivrow = Pivrow [kk] ;
	pivcol = Pivcol [kk] ;

#ifndef NDEBUG
	ASSERT (pivrow >= 0 && pivrow < Work->n_row) ;
	ASSERT (pivcol >= 0 && pivcol < Work->n_col) ;

	DEBUGm4 ((
	"\n -------------------------------------------------------------"
	"Store LU: step " ID"\n", k))  ;
	ASSERT (k < MIN (Work->n_row, Work->n_col)) ;
	DEBUG2 (("Store column of L, k = "ID", llen "ID"\n", k, llen)) ;
	for (i = 0 ; i < llen ; i++)
	{
	    row = Lpattern [i] ;
	    ASSERT (row >= 0 && row < Work->n_row) ;
	    DEBUG2 (("    Lpattern["ID"] "ID" Lpos "ID, i, row, Lpos [row])) ;
	    if (row == pivrow) DEBUG2 ((" <- pivot row")) ;
	    DEBUG2 (("\n")) ;
	    ASSERT (i == Lpos [row]) ;
	}
#endif

	/* ------------------------------------------------------------------ */
	/* remove pivot row from L */
	/* ------------------------------------------------------------------ */

	/* remove pivot row index from current column of L */
	/* if a new Lchain starts, then all entries are removed later */
	DEBUG2 (("Removing pivrow from Lpattern, k = "ID"\n", k)) ;
	ASSERT (!NON_PIVOTAL_ROW (pivrow)) ;
	pivrow_position = Lpos [pivrow] ;
	if (pivrow_position != EMPTY)
	{
	    /* place the last entry in the column in the */
	    /* position of the pivot row index */
	    ASSERT (pivrow == Lpattern [pivrow_position]) ;
	    row = Lpattern [--llen] ;
	    /* ASSERT (NON_PIVOTAL_ROW (row)) ; */
	    Lpattern [pivrow_position] = row ;
	    Lpos [row] = pivrow_position ;
	    Lpos [pivrow] = EMPTY ;
	}

	/* ------------------------------------------------------------------ */
	/* store the pivot value, for the diagonal matrix D */
	/* ------------------------------------------------------------------ */

	/* kk-th column of LU block */
	Fl1 = Flublock + kk * nb ;

	/* kk-th column of L in the L block */
	Fl2 = Flblock + kk * fnr_curr ;

	/* kk-th pivot in frontal matrix located in Flublock [kk, kk] */
	pivot_value = Fl1 [kk] ;

	D [k] = pivot_value ;
	zero_pivot = IS_ZERO (pivot_value) ;

	DEBUG4 (("Pivot D["ID"]=", k)) ;
	EDEBUG4 (pivot_value) ;
	DEBUG4 (("\n")) ;

	/* ------------------------------------------------------------------ */
	/* count nonzeros in kth column of L */
	/* ------------------------------------------------------------------ */

	lnz = 0 ;
	lnz2i = 0 ;
	lnz2x = llen ;

	for (i = kk + 1 ; i < fnpiv ; i++)
	{
	    if (IS_ZERO (Fl1 [i])) continue ;
	    lnz++ ;
	    if (Lpos [Pivrow [i]] == EMPTY) lnz2i++ ;
	}

	for (i = 0 ; i < fnrows ; i++)
	{
	    if (IS_ZERO (Fl2 [i])) continue ;
	    lnz++ ;
	    if (Lpos [Frows [i]] == EMPTY) lnz2i++ ;
	}

	lnz2x += lnz2i ;

	/* determine if we start a new Lchain or continue the old one */
	if (llen == 0 || zero_pivot)
	{
	    /* llen == 0 means there is no prior Lchain */
	    /* D [k] == 0 means the pivot column is empty */
	    newLchain = TRUE ;
	}
	else
	{
	    newLchain =
		    /* storage for starting a new Lchain */
		    UNITS (Entry, lnz) + UNITS (Int, lnz)
		<=
		    /* storage for continuing a prior Lchain */
		    UNITS (Entry, lnz2x) + UNITS (Int, lnz2i) ;
	}

	if (newLchain)
	{
	    /* start a new chain for column k of L */
	    DEBUG2 (("Start new Lchain, k = "ID"\n", k)) ;

	    pivrow_position = EMPTY ;

	    /* clear the prior Lpattern */
	    for (i = 0 ; i < llen ; i++)
	    {
		row = Lpattern [i] ;
		Lpos [row] = EMPTY ;
	    }
	    llen = 0 ;

	    lnzi = lnz ;
	    lnzx = lnz ;
	}
	else
	{
	    /* continue the prior Lchain */
	    DEBUG2 (("Continue  Lchain, k = "ID"\n", k)) ;
	    lnzi = lnz2i ;
	    lnzx = lnz2x ;
	}

	/* ------------------------------------------------------------------ */
	/* allocate space for the column of L */
	/* ------------------------------------------------------------------ */

	size = UNITS (Int, lnzi) + UNITS (Entry, lnzx) ;

#ifndef NDEBUG
	UMF_allocfail = FALSE ;
	if (UMF_gprob > 0)
	{
	    double rrr = ((double) (rand ( ))) / (((double) RAND_MAX) + 1) ;
	    DEBUG4 (("Check random %e %e\n", rrr, UMF_gprob)) ;
	    UMF_allocfail = rrr < UMF_gprob ;
	    if (UMF_allocfail) DEBUGm2 (("Random garbage coll. (store LU)\n"));
	}
#endif

	p = UMF_mem_alloc_head_block (Numeric, size) ;
	if (!p)
	{
	    Int r2, c2 ;
	    /* Do garbage collection, realloc, and try again. */
	    /* Note that there are pivot rows/columns in current front. */
	    if (Work->do_grow)
	    {
		/* full compaction of current frontal matrix, since
		 * UMF_grow_front will be called next anyway. */
		r2 = fnrows ;
		c2 = fncols ;
	    }
	    else
	    {
		/* partial compaction. */
		r2 = MAX (fnrows, Work->fnrows_new + 1) ;
		c2 = MAX (fncols, Work->fncols_new + 1) ;
	    }
	    DEBUGm3 (("get_memory from umf_store_lu:\n")) ;
	    if (!UMF_get_memory (Numeric, Work, size, r2, c2, TRUE))
	    {
		DEBUGm4 (("out of memory: store LU (1)\n")) ;
		return (FALSE) ;	/* out of memory */
	    }
	    p = UMF_mem_alloc_head_block (Numeric, size) ;
	    if (!p)
	    {
		DEBUGm4 (("out of memory: store LU (2)\n")) ;
		return (FALSE) ;	/* out of memory */
	    }
	    /* garbage collection may have moved the current front */
	    fnc_curr = Work->fnc_curr ;
	    fnr_curr = Work->fnr_curr ;
	    Flublock = Work->Flublock ;
	    Flblock  = Work->Flblock ;
	    Fublock  = Work->Fublock ;
	    Fl1 = Flublock + kk * nb ;
	    Fl2 = Flblock  + kk * fnr_curr ;
	}

	/* ------------------------------------------------------------------ */
	/* store the column of L */
	/* ------------------------------------------------------------------ */

	lip = p ;

	Li = (Int *) (Numeric->Memory + p) ;
	p += UNITS (Int, lnzi) ;
	Lval = (Entry *) (Numeric->Memory + p) ;
	p += UNITS (Entry, lnzx) ;

	for (i = 0 ; i < lnzx ; i++)
	{
	    CLEAR (Lval [i]) ;
	}

	/* store the numerical entries */

	if (newLchain)
	{
	    /* flag the first column in the Lchain by negating Lip [k] */
	    lip = -lip ;

	    ASSERT (llen == 0) ;

	    for (i = kk + 1 ; i < fnpiv ; i++)
	    {
		Int row2, pos ;
		Entry x = Fl1 [i] ;
		if (IS_ZERO (x)) continue ;
		row2 = Pivrow [i] ;
		pos = llen++ ;
		Lpattern [pos] = row2 ;
		Lpos [row2] = pos ;
		Li [pos] = row2 ;
		Lval [pos] = x ;
	    }

	    for (i = 0 ; i < fnrows ; i++)
	    {
		Int row2, pos ;
		Entry x = Fl2 [i] ;
		if (IS_ZERO (x)) continue ;
		row2 = Frows [i] ;
		pos = llen++ ;
		Lpattern [pos] = row2 ;
		Lpos [row2] = pos ;
		Li [pos] = row2 ;
		Lval [pos] = x ;
	    }

	}
	else
	{
	    ASSERT (llen > 0) ;

	    for (i = kk + 1 ; i < fnpiv ; i++)
	    {
		Int row2, pos ;
		Entry x = Fl1 [i] ;
		if (IS_ZERO (x)) continue ;
		row2 = Pivrow [i] ;
		pos = Lpos [row2] ;
		if (pos == EMPTY)
		{
		    pos = llen++ ;
		    Lpattern [pos] = row2 ;
		    Lpos [row2] = pos ;
		    *Li++ = row2 ;
		}
		Lval [pos] = x ;
	    }

	    for (i = 0 ; i < fnrows ; i++)
	    {
		Int row2, pos ;
		Entry x = Fl2 [i] ;
		if (IS_ZERO (x)) continue ;
		row2 = Frows [i] ;
		pos = Lpos [row2] ;
		if (pos == EMPTY)
		{
		    pos = llen++ ;
		    Lpattern [pos] = row2 ;
		    Lpos [row2] = pos ;
		    *Li++ = row2 ;
		}
		Lval [pos] = x ;
	    }

	}
	DEBUG4 (("llen "ID" lnzx "ID"\n", llen, lnzx)) ;
	ASSERT (llen == lnzx) ;
	ASSERT (lnz <= llen) ;
	DEBUG4 (("lnz "ID" \n", lnz)) ;

	Numeric->lnz += lnz ;
	Lnz [kk] = lnz ;

	Numeric->nLentries += lnzx ;
	Work->llen = llen ;
	Numeric->isize += lnzi ;

	/* ------------------------------------------------------------------ */
	/* the pivot column is fully assembled and scaled, and is now the */
	/* k-th column of L */
	/* ------------------------------------------------------------------ */

	Lpos [pivrow] = pivrow_position ;	/* not aliased */
	Lip [pivcol] = lip ;			/* aliased with Col_tuples */
	Lilen [pivcol] = lnzi ;			/* aliased with Col_tlen */

    }

    /* ---------------------------------------------------------------------- */
    /* store the rows of U */
    /* ---------------------------------------------------------------------- */

    for (kk = 0 ; kk < fnpiv ; kk++)
    {

	/* ------------------------------------------------------------------ */
	/* one more pivot row and column is being stored into L and U */
	/* ------------------------------------------------------------------ */

	k = Work->npiv + kk ;

	/* ------------------------------------------------------------------ */
	/* find the kth pivot row and pivot column */
	/* ------------------------------------------------------------------ */

	pivrow = Pivrow [kk] ;
	pivcol = Pivcol [kk] ;

#ifndef NDEBUG
	ASSERT (pivrow >= 0 && pivrow < Work->n_row) ;
	ASSERT (pivcol >= 0 && pivcol < Work->n_col) ;

	DEBUG2 (("Store row of U, k = "ID", ulen "ID"\n", k, ulen)) ;
	for (i = 0 ; i < ulen ; i++)
	{
	    col = Upattern [i] ;
	    DEBUG2 (("    Upattern["ID"] "ID, i, col)) ;
	    if (col == pivcol) DEBUG2 ((" <- pivot col")) ;
	    DEBUG2 (("\n")) ;
	    ASSERT (col >= 0 && col < Work->n_col) ;
	    ASSERT (i == Upos [col]) ;
	}
#endif

	/* ------------------------------------------------------------------ */
	/* get the pivot value, for the diagonal matrix D */
	/* ------------------------------------------------------------------ */

	zero_pivot = IS_ZERO (D [k]) ;

	/* ------------------------------------------------------------------ */
	/* count the nonzeros in the row of U */
	/* ------------------------------------------------------------------ */

	/* kk-th row of U in the LU block */
	Fu1 = Flublock + kk ;

	/* kk-th row of U in the U block */
	Fu2 = Fublock + kk * fnc_curr ;

	unz = 0 ;
	unz2i = 0 ;
	unz2x = ulen ;
	DEBUG2 (("unz2x is "ID", lnzx "ID"\n", unz2x, lnzx)) ;

	/* if row k does not end a Uchain, pivcol not included in ulen */
	ASSERT (!NON_PIVOTAL_COL (pivcol)) ;
	pivcol_position = Upos [pivcol] ;
	if (pivcol_position != EMPTY)
	{
	    unz2x-- ;
	    DEBUG2 (("(exclude pivcol) unz2x is now "ID"\n", unz2x)) ;
	}

	ASSERT (unz2x >= 0) ;

	for (i = kk + 1 ; i < fnpiv ; i++)
	{
	    if (IS_ZERO (Fu1 [i*nb])) continue ;
	    unz++ ;
	    if (Upos [Pivcol [i]] == EMPTY) unz2i++ ;
	}

	for (i = 0 ; i < fncols ; i++)
	{
	    if (IS_ZERO (Fu2 [i])) continue ;
	    unz++ ;
	    if (Upos [Fcols [i]] == EMPTY) unz2i++ ;
	}

	unz2x += unz2i ;

	ASSERT (IMPLIES (k == 0, ulen == 0)) ;

	/* determine if we start a new Uchain or continue the old one */
	if (ulen == 0 || zero_pivot)
	{
	    /* ulen == 0 means there is no prior Uchain */
	    /* D [k] == 0 means the matrix is singular (pivot row might */
	    /* not be empty, however, but start a new Uchain to prune zero */
	    /* entries for the deg > 0 test in UMF_u*solve) */
	    newUchain = TRUE ;
	}
	else
	{
	    newUchain =
		    /* approximate storage for starting a new Uchain */
		    UNITS (Entry, unz) + UNITS (Int, unz)
		<=
		    /* approximate storage for continuing a prior Uchain */
		    UNITS (Entry, unz2x) + UNITS (Int, unz2i) ;

	    /* this would be exact, except for the Int to Unit rounding, */
	    /* because the Upattern is stored only at the end of the Uchain */
	}

	/* ------------------------------------------------------------------ */
	/* allocate space for the row of U */
	/* ------------------------------------------------------------------ */

	size = 0 ;
	if (newUchain)
	{
	    /* store the pattern of the last row in the prior Uchain */
	    size += UNITS (Int, ulen) ;
	    unzx = unz ;
	}
	else
	{
	    unzx = unz2x ;
	}
	size += UNITS (Entry, unzx) ;

#ifndef NDEBUG
	UMF_allocfail = FALSE ;
	if (UMF_gprob > 0)
	{
	    double rrr = ((double) (rand ( ))) / (((double) RAND_MAX) + 1) ;
	    DEBUG4 (("Check random %e %e\n", rrr, UMF_gprob)) ;
	    UMF_allocfail = rrr < UMF_gprob ;
	    if (UMF_allocfail) DEBUGm2 (("Random garbage coll. (store LU)\n"));
	}
#endif

	p = UMF_mem_alloc_head_block (Numeric, size) ;
	if (!p)
	{
	    Int r2, c2 ;
	    /* Do garbage collection, realloc, and try again. */
	    /* Note that there are pivot rows/columns in current front. */
	    if (Work->do_grow)
	    {
		/* full compaction of current frontal matrix, since
		 * UMF_grow_front will be called next anyway. */
		r2 = fnrows ;
		c2 = fncols ;
	    }
	    else
	    {
		/* partial compaction. */
		r2 = MAX (fnrows, Work->fnrows_new + 1) ;
		c2 = MAX (fncols, Work->fncols_new + 1) ;
	    }
	    DEBUGm3 (("get_memory from umf_store_lu:\n")) ;
	    if (!UMF_get_memory (Numeric, Work, size, r2, c2, TRUE))
	    {
		/* :: get memory, column of L :: */
		DEBUGm4 (("out of memory: store LU (1)\n")) ;
		return (FALSE) ;	/* out of memory */
	    }
	    p = UMF_mem_alloc_head_block (Numeric, size) ;
	    if (!p)
	    {
		/* :: out of memory, column of U :: */
		DEBUGm4 (("out of memory: store LU (2)\n")) ;
		return (FALSE) ;	/* out of memory */
	    }
	    /* garbage collection may have moved the current front */
	    fnc_curr = Work->fnc_curr ;
	    fnr_curr = Work->fnr_curr ;
	    Flublock = Work->Flublock ;
	    Flblock  = Work->Flblock ;
	    Fublock  = Work->Fublock ;
	    Fu1 = Flublock + kk ;
	    Fu2 = Fublock  + kk * fnc_curr ;
	}

	/* ------------------------------------------------------------------ */
	/* store the row of U */
	/* ------------------------------------------------------------------ */

	uip = p ;

	if (newUchain)
	{
	    /* starting a new Uchain - flag this by negating Uip [k] */
	    uip = -uip ;
	    DEBUG2 (("Start new Uchain, k = "ID"\n", k)) ;

	    pivcol_position = EMPTY ;

	    /* end the prior Uchain */
	    /* save the current Upattern, and then */
	    /* clear it and start a new Upattern */
	    DEBUG2 (("Ending prior chain, k-1 = "ID"\n", k-1)) ;
	    uilen = ulen ;
	    Ui = (Int *) (Numeric->Memory + p) ;
	    Numeric->isize += ulen ;
	    p += UNITS (Int, ulen) ;
	    for (i = 0 ; i < ulen ; i++)
	    {
		col = Upattern [i] ;
		ASSERT (col >= 0 && col < Work->n_col) ;
		Upos [col] = EMPTY ;
		Ui [i] = col ;
	    }

	    ulen = 0 ;

	}
	else
	{
	    /* continue the prior Uchain */
	    DEBUG2 (("Continue  Uchain, k = "ID"\n", k)) ;
	    ASSERT (k > 0) ;

	    /* remove pivot col index from current row of U */
	    /* if a new Uchain starts, then all entries are removed later */
	    DEBUG2 (("Removing pivcol from Upattern, k = "ID"\n", k)) ;

	    if (pivcol_position != EMPTY)
	    {
		/* place the last entry in the row in the */
		/* position of the pivot col index */
		ASSERT (pivcol == Upattern [pivcol_position]) ;
		col = Upattern [--ulen] ;
		ASSERT (col >= 0 && col < Work->n_col) ;
		Upattern [pivcol_position] = col ;
		Upos [col] = pivcol_position ;
		Upos [pivcol] = EMPTY ;
	    }

	    /* this row continues the Uchain.  Keep track of how much */
	    /* to trim from the k-th length to get the length of the */
	    /* (k-1)st row of U */
	    uilen = unz2i ;

	}

	Uval = (Entry *) (Numeric->Memory + p) ;
	/* p += UNITS (Entry, unzx), no need to increment p */

	for (i = 0 ; i < unzx ; i++)
	{
	    CLEAR (Uval [i]) ;
	}

	if (newUchain)
	{
	    ASSERT (ulen == 0) ;

	    for (i = kk + 1 ; i < fnpiv ; i++)
	    {
		Int col2, pos ;
		Entry x = Fu1 [i*nb] ;
		if (IS_ZERO (x)) continue ;
		col2 = Pivcol [i] ;
		pos = ulen++ ;
		Upattern [pos] = col2 ;
		Upos [col2] = pos ;
		Uval [pos] = x ;
	    }

	    for (i = 0 ; i < fncols ; i++)
	    {
		Int col2, pos ;
		Entry x = Fu2 [i] ;
		if (IS_ZERO (x)) continue ;
		col2 = Fcols [i] ;
		pos = ulen++ ;
		Upattern [pos] = col2 ;
		Upos [col2] = pos ;
		Uval [pos] = x ;
	    }

	}
	else
	{

	    ASSERT (ulen > 0) ;

	    /* store the numerical entries and find new nonzeros */

	    for (i = kk + 1 ; i < fnpiv ; i++)
	    {
		Int col2, pos ;
		Entry x = Fu1 [i*nb] ;
		if (IS_ZERO (x)) continue ;
		col2 = Pivcol [i] ;
		pos = Upos [col2] ;
		if (pos == EMPTY)
		{
		    pos = ulen++ ;
		    Upattern [pos] = col2 ;
		    Upos [col2] = pos ;
		}
		Uval [pos] = x ;
	    }

	    for (i = 0 ; i < fncols ; i++)
	    {
		Int col2, pos ;
		Entry x = Fu2 [i] ;
		if (IS_ZERO (x)) continue ;
		col2 = Fcols [i] ;
		pos = Upos [col2] ;
		if (pos == EMPTY)
		{
		    pos = ulen++ ;
		    Upattern [pos] = col2 ;
		    Upos [col2] = pos ;
		}
		Uval [pos] = x ;
	    }

	}

	ASSERT (ulen == unzx) ;
	ASSERT (unz <= ulen) ;
	DEBUG4 (("unz "ID" \n", unz)) ;

	Numeric->unz += unz ;
	/* count the "true" flops, based on LU pattern only */
	Numeric->flops += DIV_FLOPS * Lnz [kk]	/* scale pivot column */
	    + MULTSUB_FLOPS * (Lnz [kk] * unz) ;    /* outer product */

	Numeric->nUentries += unzx ;
	Work->ulen = ulen ;
	DEBUG1 (("Work->ulen = "ID" at end of pivot step, k: "ID"\n", ulen, k));

	/* ------------------------------------------------------------------ */
	/* the pivot row is the k-th row of U */
	/* ------------------------------------------------------------------ */

	Upos [pivcol] = pivcol_position ;	/* not aliased */
	Uip [pivrow] = uip ;			/* aliased with Row_tuples */
	Uilen [pivrow] = uilen ;		/* aliased with Row_tlen */

    }

    /* ---------------------------------------------------------------------- */
    /* no more pivots in frontal working array */
    /* ---------------------------------------------------------------------- */

    Work->npiv += fnpiv ;
    Work->fnpiv = 0 ;
    Work->fnzeros = 0 ;
    return (TRUE) ;
}
