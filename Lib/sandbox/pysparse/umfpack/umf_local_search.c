/* ========================================================================== */
/* === UMF_local_search ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Perform pivot search to find pivot row and pivot column.
    The pivot column is selected from the candidate set.  The candidate set
    corresponds to a supercolumn from colamd or UMF_analyze.  The pivot column
    is then removed from that set.  Constructs the pivot column pattern and
    values.  Called by umf_kernel.  Returns UMFPACK_OK if successful, or
    UMFPACK_WARNING_singular_matrix or UMFPACK_ERROR_different_pattern if not.
*/

#include "umf_internal.h"
#include "umf_row_search.h"
#include "umf_mem_free_tail_block.h"

/* Version 4.1:  relaxed amalgamation control parameters are now fixed, and
 * cannot be changed via Control [..] settings, as they could in Version 4.0. */
#define RELAX1 0.25	    /* this was UMFPACK_DEFAULT_RELAXED_AMALGAMATION */
#define SYM_RELAX1 0.0	    /* this is new to Version 4.1 */
#define RELAX2 0.1	    /* this was UMFPACK_DEFAULT_RELAXED2_AMALGAMATION */
#define RELAX3 0.125	    /* this was UMFPACK_DEFAULT_RELAXED3_AMALGAMATION */

/* ========================================================================== */
/* === remove_candidate ===================================================== */
/* ========================================================================== */

/* Remove a column from the set of candidate pivot columns. */

PRIVATE void remove_candidate (Int jj, WorkType *Work, SymbolicType *Symbolic)
{

#ifndef NDEBUG
    Int j ;
    DEBUGm2 (("pivot column Candidates before remove: nCand "ID" ncand "ID
	" lo "ID" hi "ID" jj "ID"\n", Work->nCandidates, Work->ncand,
	Work->lo, Work->hi, jj)) ;
    for (j = 0 ; j < Work->nCandidates ; j++)
    {
	Int col = Work->Candidates [j] ;
	DEBUGm2 ((ID" ", col));
	ASSERT (col >= 0 && col < Work->n_col) ;
	/* ASSERT (NON_PIVOTAL_COL (col)) ; */
	ASSERT (col >= Work->lo && col <= Work->hi) ;
    }
    DEBUGm2 (("\n")) ;
#endif

    if (Symbolic->fixQ)
    {
	DEBUGm2 (("FixQ\n")) ;
	/* do not modify the column ordering */
	ASSERT (Work->nCandidates == 1) ;
	ASSERT (jj == 0) ;
	if (Work->ncand > 1)
	{
	    Work->Candidates [0] = Work->nextcand++ ;
	}
	else
	{
	    Work->nCandidates = 0 ;
	}
    }
    else
    {
	/* place the next candidate in the set */
	if (Work->ncand > MAX_CANDIDATES)
	{
	    Work->Candidates [jj] = Work->nextcand++ ;
	}
	else
	{
	    ASSERT (Work->nCandidates == Work->ncand) ;
	    Work->Candidates [jj] = Work->Candidates [Work->ncand - 1] ;
	    Work->Candidates [Work->ncand - 1] = EMPTY ;
	    Work->nCandidates-- ;
	}
    }
    Work->ncand-- ;

#ifndef NDEBUG
    DEBUGm2 (("pivot column Candidates after remove: nCand "ID" ncand "ID
	" lo "ID" hi "ID" jj "ID"\n", Work->nCandidates, Work->ncand, Work->lo,
	Work->hi, jj)) ;
    for (j = 0 ; j < Work->nCandidates ; j++)
    {
	Int col = Work->Candidates [j] ;
	DEBUGm2 ((ID" ", col));
	ASSERT (col >= 0 && col < Work->n_col) ;
	/* ASSERT (NON_PIVOTAL_COL (col)) ; */
	ASSERT (col >= Work->lo && col <= Work->hi) ;
    }
    DEBUGm2 (("\n")) ;
    ASSERT (Work->ncand >= 0) ;
    ASSERT (Work->nCandidates <= Work->ncand) ;
#endif
}

/* ========================================================================== */
/* === UMF_local_search ===================================================== */
/* ========================================================================== */

GLOBAL Int UMF_local_search
(
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Entry *Flblock, *Fublock, *Fs, *Fcblock, *C, *Wx, *Wy, *Fu, *Flublock,
	*Flu ;
    double relax1 ;
    Int pos, nrows, *Cols, *Rows, e, f, status, max_cdeg, fnzeros, nb, j, col,
	i, row, cdeg_in, rdeg [2][2], fnpiv, nothing [2], new_LUsize,
	pivrow [2][2], pivcol [2], *Wp, *Fcpos, *Frpos, new_fnzeros, cdeg_out,
	*Wm, *Wio, *Woi, *Woo, *Frows, *Fcols, fnrows, fncols, *E, deg, nr_in,
	nc, thiscost, bestcost, nr_out, do_update, extra_cols, extra_rows,
	extra_zeros, relaxed_front, do_extend, fnr_curr, fnc_curr, tpi,
	*Col_tuples, *Col_degree, *Col_tlen, jj, jcand [2], freebie [2],
	did_rowmerge, fnrows_new [2][2], fncols_new [2][2], search_pivcol_out,
	*Diagonal_map, *Diagonal_imap, row2, col2 ;
    Unit *Memory, *p ;
    Tuple *tp, *tpend, *tp1, *tp2 ;
    Element *ep ;

#ifndef NDEBUG
    Int debug_ok, n_row, n_col, *Row_degree ;
    Row_degree = Numeric->Rperm ;	/* for NON_PIVOTAL_ROW macro only */
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    Memory = Numeric->Memory ;
    E = Work->E ;
    Col_degree = Numeric->Cperm ;

    Col_tuples = Numeric->Lip ;
    Col_tlen   = Numeric->Lilen ;

    Wx = Work->Wx ;
    Wy = Work->Wy ;
    Wp = Work->Wp ;
    Wm = Work->Wm ;
    Woi = Work->Woi ;
    Wio = Work->Wio ;
    Woo = Work->Woo ;
    Fcpos = Work->Fcpos ;
    Frpos = Work->Frpos ;
    Frows = Work->Frows ;
    Fcols = Work->Fcols ;
    fnrows = Work->fnrows ;
    fncols = Work->fncols ;
    nb = Work->nb ;
    fnr_curr = Work->fnr_curr ;
    fnc_curr = Work->fnc_curr ;
    fnpiv = Work->fnpiv ;
    nothing [0] = EMPTY ;
    nothing [1] = EMPTY ;
    relax1 = (Symbolic->prefer_diagonal) ? SYM_RELAX1 : RELAX1 ;
    fnzeros = Work->fnzeros ;
    new_fnzeros = fnzeros ;
    jj = EMPTY ;

    Fcblock = Work->Fcblock ;	    /* current contribution block */
    Flblock = Work->Flblock ;	    /* current L block */
    Fublock = Work->Fublock ;	    /* current U block */
    Flublock = Work->Flublock ;	    /* current LU block */

    /* The pivot column degree cannot exceed max_cdeg */
    max_cdeg = Work->fnrows_max ;
    ASSERT (Work->fnrows_max <= Symbolic->maxnrows) ;
    ASSERT (Work->fncols_max <= Symbolic->maxncols) ;

    if (fnrows == 0 && fncols == 0)
    {
	/* frontal matrix is empty */
	Work->firstsuper = Work->ksuper ;
    }

#ifndef NDEBUG
    n_row = Work->n_row ;
    n_col = Work->n_col ;
    DEBUG2 (("\n========LOCAL SEARCH:  current frontal matrix: ========= \n")) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
    if (UMF_debug > 0 || MAX (n_row, n_col) < 1000)
    {
	for (i = 0 ; i < MAX (n_row, n_col) ; i++)
	{
	    ASSERT (Wp [i] < 0) ;
	}
    }

    DEBUGm2 ((ID" pivot column Candidates: lo "ID" hi "ID"\n",
	Work->nCandidates, Work->lo, Work->hi)) ;
    for (j = 0 ; j < Work->nCandidates ; j++)
    {
	col = Work->Candidates [j] ;
	DEBUGm2 ((ID" ", col));
	ASSERT (col >= 0 && col < n_col) ;
	ASSERT (NON_PIVOTAL_COL (col)) ;
	ASSERT (col >= Work->lo && col <= Work->hi) ;
    }

    DEBUGm2 (("\n")) ;
    /* there are no 0-by-c or r-by-0 fronts, where c and r are > 0 */
    /* a front is either 0-by-0, or r-by-c */
    DEBUG2 (("\n\n::: "ID" : Npiv: "ID" + fnpiv "ID" = "ID". "
	"size "ID"-by-"ID"\n", Work->frontid,
	Work->npiv, Work->fnpiv, Work->npiv + Work->fnpiv, fnrows, fncols)) ;
    ASSERT ((fnrows == 0 && fncols == 0) ||(fnrows != 0 && fncols != 0)) ;
#endif

    /* ====================================================================== */
    /* === PIVOT SEARCH ===================================================== */
    /* ====================================================================== */

    /* initialize */

    pivcol [IN] = EMPTY ;
    pivcol [OUT] = EMPTY ;

    cdeg_in = Int_MAX ;
    cdeg_out = Int_MAX ;

    pivrow [IN][IN] = EMPTY ;
    pivrow [IN][OUT] = EMPTY ;
    pivrow [OUT][IN] = EMPTY ;
    pivrow [OUT][OUT] = EMPTY ;

    rdeg [IN][IN] = Int_MAX ;
    rdeg [IN][OUT] = Int_MAX ;
    rdeg [OUT][IN] = Int_MAX ;
    rdeg [OUT][OUT] = Int_MAX ;

    freebie [IN] = FALSE ;
    freebie [OUT] = FALSE ;

    Work->pivot_case = EMPTY ;
    bestcost = EMPTY ;

    nr_out = EMPTY ;
    nr_in = EMPTY ;

    jcand [IN] = EMPTY ;
    jcand [OUT] = EMPTY ;

    fnrows_new [IN][IN] = EMPTY ;
    fnrows_new [IN][OUT] = EMPTY ;
    fnrows_new [OUT][IN] = EMPTY ;
    fnrows_new [OUT][OUT] = EMPTY ;

    fncols_new [IN][IN] = EMPTY ;
    fncols_new [IN][OUT] = EMPTY ;
    fncols_new [OUT][IN] = EMPTY ;
    fncols_new [OUT][OUT] = EMPTY ;

#ifndef NDEBUG
	/* check Frpos */
	DEBUG4 (("Check Frpos : fnrows "ID" col "ID" maxcdeg "ID"\n",
		fnrows, pivcol [IN], max_cdeg)) ;
	for (i = 0 ; i < fnrows ; i++)
	{
	    row = Frows [i] ;
	    DEBUG4 (("  row: "ID"\n", row)) ;
	    ASSERT (row >= 0 && row < n_row) ;
	    ASSERT (Frpos [row] == i) ;
	}
	DEBUG4 (("All:\n")) ;
	if (UMF_debug > 0 || n_row < 1000)
	{
	    Int cnt = fnrows ;
	    for (row = 0 ; row < n_row ; row++)
	    {
		if (Frpos [row] == EMPTY)
		{
		    cnt++ ;
		}
		else
		{
		    DEBUG4 (("  row: "ID" pos "ID"\n", row, Frpos [row])) ;
		}
	    }
	    ASSERT (cnt == n_row) ;
	}
#endif

    /* ---------------------------------------------------------------------- */
    /* find shortest column in the front, and shortest column not in the */
    /* front, from the candidate pivot column set */
    /* ---------------------------------------------------------------------- */

    /* If there are too many candidates, then only look at the first */
    /* MAX_CANDIDATES of them.   Otherwise, if there are O(n) candidates, */
    /* this code could take O(n^2) time. */

    /* ------------------------------------------------------------------ */
    /* look in the candidate set for the best column */
    /* ------------------------------------------------------------------ */

    DEBUG2 (("Max candidates %d, Work->ncand "ID" jmax "ID"\n",
	MAX_CANDIDATES, Work->ncand, Work->nCandidates)) ;
    col = Work->Candidates [0] ;
    ASSERT (Work->nCandidates > 0) ;
    DEBUG3 (("Pivot column candidate: "ID" j = "ID"\n", col, j)) ;
    ASSERT (col >= 0 && col < n_col) ;

    /* there is no Col_degree if fixQ is true */
    deg = Symbolic->fixQ ? EMPTY : Col_degree [col] ;

#ifndef NDEBUG
    DEBUG3 (("Pivot column candidate: "ID" cost: "ID"  Fcpos[col] "ID"\n",
	col, deg, Fcpos [col])) ;
    UMF_dump_rowcol (1, Numeric, Work, col, !Symbolic->fixQ) ;
    if (Symbolic->fixQ)
    {
	DEBUG1 (("FIXQ: Candidates "ID" pivcol "ID" npiv "ID" fnpiv "ID
	    " ndiscard "ID "\n", Work->nCandidates, col, Work->npiv,
	    Work->fnpiv, Work->ndiscard)) ;
	ASSERT (Work->nCandidates == 1) ;
	ASSERT (col == Work->npiv + Work->fnpiv + Work->ndiscard) ;
    }
#endif

    if (Fcpos [col] >= 0)
    {
	/* best column in front, so far */
	pivcol [IN] = col ;
	cdeg_in = deg ;		/* ignored, if fixQ is true */
	jcand [IN] = 0 ;
    }
    else
    {
	/* best column not in front, so far */
	pivcol [OUT] = col ;
	cdeg_out = deg ;	/* ignored, if fixQ is true */
	jcand [OUT] = 0 ;
    }

    /* look at the rest of the candidates */
    for (j = 1 ; j < Work->nCandidates ; j++)
    {
	col = Work->Candidates [j] ;

	DEBUG3 (("Pivot col candidate: "ID" j = "ID"\n", col, j)) ;
	ASSERT (col >= 0 && col < n_col) ;
	ASSERT (!Symbolic->fixQ) ;
	deg = Col_degree [col] ;
#ifndef NDEBUG
	DEBUG3 (("Pivot col candidate: "ID" cost: "ID" Fcpos[col] "ID"\n",
		col, deg, Fcpos [col])) ;
	UMF_dump_rowcol (1, Numeric, Work, col, !Symbolic->fixQ) ;
#endif
	if (Fcpos [col] >= 0)
	{
#ifndef NDEBUG
	    Int fs ;
	    fs = Fcpos [col] / fnr_curr ;
	    ASSERT (fs >= 0 && fs < fncols) ;
#endif
	    if (deg < cdeg_in || (deg == cdeg_in && col < pivcol [IN]))
	    {
		/* best column in front, so far */
		pivcol [IN] = col ;
		cdeg_in = deg ;
		jcand [IN] = j ;
	    }
	}
	else
	{
	    if (deg < cdeg_out || (deg == cdeg_out && col < pivcol [OUT]))
	    {
		/* best column not in front, so far */
		pivcol [OUT] = col ;
		cdeg_out = deg ;
		jcand [OUT] = j ;
	    }
	}
    }

    DEBUG2 (("Pivcol in "ID" out "ID"\n", pivcol [IN], pivcol [OUT])) ;
    ASSERT ((pivcol [IN] >= 0 && pivcol [IN] < n_col)
	|| (pivcol [OUT] >= 0 && pivcol [OUT] < n_col)) ;

    cdeg_in = EMPTY ;
    cdeg_out = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* construct candidate column in front, and search for pivot rows */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    /* check Frpos */
    DEBUG4 (("Prior to col update: fnrows "ID" col "ID" maxcdeg "ID"\n",
	    fnrows, pivcol [IN], max_cdeg)) ;
    for (i = 0 ; i < fnrows ; i++)
    {
	row = Frows [i] ;
	DEBUG4 (("  row: "ID"\n", row)) ;
	ASSERT (row >= 0 && row < n_row) ;
	ASSERT (Frpos [row] == i) ;
    }
    DEBUG4 (("All:\n")) ;
    if (UMF_debug > 0 || n_row < 1000)
    {
	Int cnt = fnrows ;
	for (row = 0 ; row < n_row ; row++)
	{
	    if (Frpos [row] == EMPTY)
	    {
		cnt++ ;
	    }
	    else
	    {
		DEBUG4 (("  row: "ID" pos "ID"\n", row, Frpos [row])) ;
	    }
	}
	ASSERT (cnt == n_row) ;
    }
#endif

    if (pivcol [IN] != EMPTY)
    {

#ifndef NDEBUG
	DEBUG2 (("col[IN] column "ID" in front at position = "ID"\n",
		pivcol [IN], Fcpos [pivcol [IN]])) ;
	UMF_dump_rowcol (1, Numeric, Work, pivcol [IN], !Symbolic->fixQ) ;
#endif

	/* the only way we can have a pivcol[IN] is if the front is not empty */
	ASSERT (fnrows > 0 && fncols > 0) ;

	DEBUG4 (("Update pivot column:\n")) ;
	Fs  = Fcblock  +  Fcpos [pivcol [IN]] ;
	Fu  = Fublock  + (Fcpos [pivcol [IN]] / fnr_curr) ;
	Flu = Flublock + fnpiv * nb ;

	/* ------------------------------------------------------------------ */
	/* copy the pivot column from the U block into the LU block */
	/* ------------------------------------------------------------------ */

	/* This copy is permanent if the pivcol [IN] is chosen. */
	for (i = 0 ; i < fnpiv ; i++)
	{
	    Flu [i] = Fu [i*fnc_curr] ;
	}

	/* ------------------------------------------------------------------ */
	/* update the pivot column in the LU block using a triangular solve */
	/* ------------------------------------------------------------------ */

	/* This work will be discarded if the pivcol [OUT] is chosen instead.
	 * It is permanent if the pivcol [IN] is chosen. */

	if (fnpiv > 1)
	{
	    /* solve Lx=b, where b = U (:,k), stored in the LU block */

#ifdef USE_NO_BLAS

	    /* no BLAS available - use plain C code instead */
	    Entry *Flub = Flublock ;
	    for (j = 0 ; j < fnpiv ; j++)
	    {
		Entry Fuj = Flu [j] ;
#pragma ivdep
		for (i = j+1 ; i < fnpiv ; i++)
		{
		    /* Flu [i] -= Flublock [i + j*nb] * Flu [j] ; */
		    MULT_SUB (Flu [i], Flub [i], Fuj) ;
		}
		Flub += nb ;
	    }

#else
	    BLAS_TRSV (fnpiv, Flublock, Flu, nb) ;
#endif

	}

	/* ------------------------------------------------------------------ */
	/* copy the pivot column from the C block into Wy */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < fnrows ; i++)
	{
	    Wy [i] = Fs [i] ;
	}

	/* ------------------------------------------------------------------ */
	/* update the pivot column of L using a matrix-vector multiply */
	/* ------------------------------------------------------------------ */

	/* this work will be discarded if the pivcol [OUT] is chosen instead */

#ifdef USE_NO_BLAS
	/* no BLAS available - use plain C code instead */
	for (j = 0 ; j < fnpiv ; j++)
	{
	    Entry Fuj, *Flub = Flblock + j * fnr_curr ;
	    Fuj = Flu [j] ;
	    if (IS_NONZERO (Fuj))
	    {
#pragma ivdep
		for (i = 0 ; i < fnrows ; i++)
		{
		    /* Wy [i] -= Flblock [i+j*fnr_curr] * Fuj ; */
		    MULT_SUB (Wy [i], Flub [i], Fuj) ;
		}
	    }
	    /* Flblock += fnr_curr ; */
	}
#else
	/* Using 1-based notation:
	 * Wy (1:fnrows) -= Flblock (1:fnrows,1:fnpiv) * Flu (1:fnpiv) */
	BLAS_GEMV (fnrows, fnpiv, Flblock, Flu, Wy, fnr_curr) ;
#endif

	/* ------------------------------------------------------------------ */

#ifndef NDEBUG
	DEBUG2 (("Wy after update: fnrows="ID"\n", fnrows)) ;
	DEBUG4 ((" fnpiv="ID" \n", fnpiv)) ;
	for (i = 0 ; i < fnrows ; i++)
	{
	    DEBUG4 ((ID" "ID" "ID, i, Frows [i], Frpos [Frows [i]])) ;
	    EDEBUG4 (Wy [i]) ;
	    DEBUG4 (("\n")) ;
	}
#endif

	/* ------------------------------------------------------------------ */
	/* construct the candidate column */
	/* ------------------------------------------------------------------ */

	cdeg_in = fnrows ;

#ifndef NDEBUG
	/* check Frpos */
	DEBUG4 (("After col update: fnrows "ID" col "ID" maxcdeg "ID"\n",
		fnrows, pivcol [IN], max_cdeg)) ;
	for (i = 0 ; i < fnrows ; i++)
	{
	    row = Frows [i] ;
	    DEBUG4 (("  row: "ID"\n", row)) ;
	    ASSERT (row >= 0 && row < n_row) ;
	    ASSERT (Frpos [row] == i) ;
	}
	DEBUG4 (("All:\n")) ;
	if (UMF_debug > 0 || n_row < 1000)
	{
	    Int cnt = fnrows ;
	    for (row = 0 ; row < n_row ; row++)
	    {
		if (Frpos [row] == EMPTY)
		{
		    cnt++ ;
		}
		else
		{
		    DEBUG4 (("  row: "ID" pos "ID"\n", row, Frpos [row])) ;
		}
	    }
	    ASSERT (cnt == n_row) ;
	}
#endif

#ifndef NDEBUG
	/* check Frpos */
	DEBUG4 (("COL ASSEMBLE: cdeg "ID"\nREDUCE COL in "ID" max_cdeg "ID"\n",
		cdeg_in, pivcol [IN], max_cdeg)) ;
	for (i = 0 ; i < cdeg_in ; i++)
	{
	    row = Frows [i] ;
	    ASSERT (row >= 0 && row < n_row) ;
	    ASSERT (Frpos [row] == i) ;
	}
	if (UMF_debug > 0 || n_row < 1000)
	{
	    Int cnt = cdeg_in ;
	    for (row = 0 ; row < n_row ; row++)
	    {
		if (Frpos [row] == EMPTY) cnt++ ;
	    }
	    ASSERT (cnt == n_row) ;
	}
#endif

	/* assemble column into Wy */

	ASSERT (pivcol [IN] >= 0 && pivcol [IN] < n_col) ;
	ASSERT (NON_PIVOTAL_COL (pivcol [IN])) ;

	tpi = Col_tuples [pivcol [IN]] ;
	if (tpi)
	{
	    tp = (Tuple *) (Memory + tpi) ;
	    tp1 = tp ;
	    tp2 = tp ;
	    tpend = tp + Col_tlen [pivcol [IN]] ;
	    for ( ; tp < tpend ; tp++)
	    {
		e = tp->e ;
		ASSERT (e > 0 && e <= Work->nel) ;
		if (!E [e]) continue ;	/* element already deallocated */
		f = tp->f ;
		p = Memory + E [e] ;
		ep = (Element *) p ;
		p += UNITS (Element, 1) ;
		Cols = (Int *) p ;
		if (Cols [f] == EMPTY) continue ; /* column already assembled */
		ASSERT (pivcol [IN] == Cols [f]) ;

		Rows = Cols + ep->ncols ;
		nrows = ep->nrows ;
		p += UNITS (Int, ep->ncols + nrows) ;
		C = ((Entry *) p) + f * nrows ;

		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    if (row >= 0) /* skip this if already gone from element */
		    {
			ASSERT (row < n_row) ;
			pos = Frpos [row] ;
			if (pos < 0)
			{
			    /* new entry in the pattern - save Frpos */
			    ASSERT (cdeg_in < n_row) ;
			    if (cdeg_in >= max_cdeg)
			    {
				/* :: pattern change (cdeg in failure) :: */
				DEBUGm4 (("cdeg_in failure\n")) ;
				return (UMFPACK_ERROR_different_pattern) ;
			    }
			    Frpos [row] = cdeg_in ;
			    Frows [cdeg_in] = row ;
			    Wy [cdeg_in++] = C [i] ;
			}
			else
			{
			    /* entry already in pattern - sum values in Wy */
			    /* Wy [pos] += C [i] ; */
			    ASSERT (pos < max_cdeg) ;
			    ASSEMBLE (Wy [pos], C [i]) ;
			}
		    }
		}
		*tp2++ = *tp ;	/* leave the tuple in the list */
	    }
	    Col_tlen [pivcol [IN]] = tp2 - tp1 ;
	}

	/* ------------------------------------------------------------------ */

#ifndef NDEBUG
	/* check Frpos again */
	DEBUG4 (("COL DONE: cdeg "ID"\nREDUCE COL in "ID" max_cdeg "ID"\n",
		cdeg_in, pivcol [IN], max_cdeg)) ;
	for (i = 0 ; i < cdeg_in ; i++)
	{
	    row = Frows [i] ;
	    ASSERT (row >= 0 && row < n_row) ;
	    ASSERT (Frpos [row] == i) ;
	}
	if (UMF_debug > 0 || n_row < 1000)
	{
	    Int cnt = cdeg_in ;
	    for (row = 0 ; row < n_row ; row++)
	    {
		if (Frpos [row] == EMPTY) cnt++ ;
	    }
	    ASSERT (cnt == n_row) ;
	}
#endif

#ifndef NDEBUG
	DEBUG4 (("Reduced column: cdeg in  "ID" fnrows_max "ID"\n",
	    cdeg_in, Work->fnrows_max)) ;
	for (i = 0 ; i < cdeg_in ; i++)
	{
	    DEBUG4 ((" "ID" "ID" "ID, i, Frows [i], Frpos [Frows [i]])) ;
	    EDEBUG4 (Wy [i]) ;
	    DEBUG4 (("\n")) ;
	    ASSERT (i == Frpos [Frows [i]]) ;
	}
	ASSERT (cdeg_in <= Work->fnrows_max) ;
#endif

	/* ------------------------------------------------------------------ */
	/* cdeg_in is now the exact degree of this column */
	/* ------------------------------------------------------------------ */

	nr_in = cdeg_in - fnrows ;

	/* since there are no 0-by-x fronts, if there is a pivcol [IN] the */
	/* front must have at least one row. */
	ASSERT (cdeg_in > 0) ;

	/* new degree of pivcol [IN], excluding current front is nr_in */
	/* column expands by nr_in rows */

	/* ------------------------------------------------------------------ */
	/* search for two candidate pivot rows */
	/* ------------------------------------------------------------------ */

	/* for the IN_IN pivot row (if any), */
	/* extend the pattern in place, using Fcols */
	status = UMF_row_search (Numeric, Work, Symbolic,
	    fnrows, cdeg_in, Frows, Frpos,   /* pattern of column to search */
	    pivrow [IN], rdeg [IN], Fcols, Wio, nothing, Wy,
	    pivcol [IN], freebie) ;
	ASSERT (!freebie [IN] && !freebie [OUT]) ;

	/* ------------------------------------------------------------------ */
	/* fatal error if matrix pattern has changed since symbolic analysis */
	/* ------------------------------------------------------------------ */

	if (status == UMFPACK_ERROR_different_pattern)
	{
	    /* :: pattern change (row search IN failure) :: */
	    DEBUGm4 (("row search IN failure\n")) ;
	    return (UMFPACK_ERROR_different_pattern) ;
	}

	/* ------------------------------------------------------------------ */
	/* we now must have a structural pivot */
	/* ------------------------------------------------------------------ */

	/* Since the pivcol[IN] exists, there must be at least one row in the */
	/* current frontal matrix, and so we must have found a structural */
	/* pivot.  The numerical value might be zero, of course. */

	ASSERT (status != UMFPACK_WARNING_singular_matrix) ;

	/* ------------------------------------------------------------------ */
	/* evaluate IN_IN option */
	/* ------------------------------------------------------------------ */

	if (pivrow [IN][IN] != EMPTY)
	{
	    /* The current front would become an (implicit) LUson.
	     * Both candidate pivot row and column are in the current front.
	     * Cost is how much the current front would expand */

	    /* pivrow[IN][IN] candidates are not found via row merge search */

	    ASSERT (fnrows >= 0 && fncols >= 0) ;

	    ASSERT (cdeg_in > 0) ;
	    nc = rdeg [IN][IN] - fncols ;

	    thiscost =
		/* each column in front (except pivot column) grows by nr_in: */
		(nr_in * (fncols - 1)) +
		/* new columns not in old front: */
		(nc * (cdeg_in - 1)) ;

	    /* no extra cost to relaxed amalgamation */

	    ASSERT (fnrows + nr_in == cdeg_in) ;
	    ASSERT (fncols + nc == rdeg [IN][IN]) ;

	    /* size of relaxed front (after pivot row column removed): */
	    fnrows_new [IN][IN] = (fnrows-1) + nr_in ;
	    fncols_new [IN][IN] = (fncols-1) + nc ;
	    /* relaxed_front = fnrows_new [IN][IN] * fncols_new [IN][IN] ; */

	    do_extend = TRUE ;

	    DEBUG2 (("Evaluating option IN-IN:\n")) ;
	    DEBUG2 (("Work->fnzeros "ID" fnpiv "ID" nr_in "ID" nc "ID"\n",
		Work->fnzeros, fnpiv, nr_in, nc)) ;
	    DEBUG2 (("fncols "ID" fnrows "ID"\n", fncols, fnrows)) ;

	    /* determine if BLAS-3 update should be applied before extending. */
	    /* update if too many zero entries accumulate in the LU block */
	    fnzeros = Work->fnzeros + fnpiv * (nr_in + nc) ;

	    DEBUG2 (("fnzeros "ID"\n", fnzeros)) ;

	    new_LUsize = (fnpiv+1) * (fnrows + nr_in + fncols + nc) ;

	    DEBUG2 (("new_LUsize "ID"\n", new_LUsize)) ;

	    /* There are fnpiv pivots currently in the front.  This one
	     * will be the (fnpiv+1)st pivot, if it is extended. */

	    /* RELAX2 parameter uses a double relop, but ignore NaN case: */
	    do_update = fnpiv > 0 &&
		(((double) fnzeros) / ((double) new_LUsize)) > RELAX2 ;

	    DEBUG2 (("do_update "ID"\n", do_update))

	    DEBUG2 (("option IN  IN : nr "ID" nc "ID" cost "ID"(0) relax "ID
		"\n", nr_in, nc, thiscost, do_extend)) ;

	    /* this is the best option seen so far */
	    Work->pivot_case = IN_IN ;
	    bestcost = thiscost ;

	    /* do the amalgamation and extend the front */
	    Work->do_extend = do_extend ;
	    Work->do_update = do_update ;
	    new_fnzeros = fnzeros ;

	}

	/* ------------------------------------------------------------------ */
	/* evaluate IN_OUT option */
	/* ------------------------------------------------------------------ */

	if (pivrow [IN][OUT] != EMPTY)
	{
	    /* The current front would become a Uson of the new front.
	     * Candidate pivot column is in the current front, but the
	     * candidate pivot row is not. */

	    ASSERT (fnrows >= 0 && fncols > 0) ;
	    ASSERT (cdeg_in > 0) ;

	    /* must be at least one row outside the front */
	    /* (the pivrow [IN][OUT] itself) */
	    ASSERT (nr_in >= 1) ;

	    /* count columns not in current front */
	    nc = 0 ;
#ifndef NDEBUG
	    debug_ok = FALSE ;
#endif
	    for (i = 0 ; i < rdeg [IN][OUT] ; i++)
	    {
		col = Wio [i] ;
		DEBUG4 (("counting col "ID" Fcpos[] = "ID"\n", col,
		    Fcpos [col])) ;
		ASSERT (col >= 0 && col < n_col && NON_PIVOTAL_COL (col)) ;
		if (Fcpos [col] < 0) nc++ ;
#ifndef NDEBUG
		/* we must see the pivot column somewhere */
		if (col == pivcol [IN])
		{
		    ASSERT (Fcpos [col] >= 0) ;
		    debug_ok = TRUE ;
		}
#endif
	    }
	    ASSERT (debug_ok) ;

	    thiscost =
		/* each row in front grows by nc: */
		(nc * fnrows) +
		/* new rows not affected by front: */
		((nr_in-1) * (rdeg [IN][OUT]-1)) ;

	    /* check the cost of relaxed IN_OUT amalgamation */

	    extra_cols = ((fncols-1) + nc ) - (rdeg [IN][OUT] - 1) ;
	    ASSERT (extra_cols >= 0) ;
	    ASSERT (fncols + nc == extra_cols + rdeg [IN][OUT]) ;
	    extra_zeros = (nr_in-1) * extra_cols ;	/* symbolic fill-in */

	    ASSERT (fnrows + nr_in == cdeg_in) ;
	    ASSERT (fncols + nc == rdeg [IN][OUT] + extra_cols) ;

	    /* size of relaxed front (after pivot column removed): */
	    fnrows_new [IN][OUT] = fnrows + (nr_in-1) ;
	    fncols_new [IN][OUT] = (fncols-1) + nc ;
	    relaxed_front = fnrows_new [IN][OUT] * fncols_new [IN][OUT] ;

	    /* do relaxed amalgamation if the extra zeros are no more */
	    /* than a fraction (default 0.25) of the relaxed front */
	    /* if relax = 0: no extra zeros allowed */
	    /* if relax = +inf: always amalgamate */

	    /* relax parameter uses a double relop, but ignore NaN case: */
	    if (extra_zeros == 0)
	    {
		do_extend = TRUE ;
	    }
	    else
	    {
		do_extend = ((double) extra_zeros) <
		   (relax1 * (double) relaxed_front) ;
	    }

	    if (do_extend)
	    {
		/* count the cost of relaxed amalgamation */
		thiscost += extra_zeros ;

		DEBUG2 (("Evaluating option IN-OUT:\n")) ;
		DEBUG2 (("Work->fnzeros "ID" fnpiv "ID" nr_in "ID" nc "ID"\n",
		    Work->fnzeros, fnpiv, nr_in, nc)) ;
		DEBUG2 (("fncols "ID" fnrows "ID"\n", fncols, fnrows)) ;

		/* determine if BLAS-3 update to be applied before extending. */
		/* update if too many zero entries accumulate in the LU block */
		fnzeros = Work->fnzeros + fnpiv * (nr_in + nc) ;

		DEBUG2 (("fnzeros "ID"\n", fnzeros)) ;

		new_LUsize = (fnpiv+1) * (fnrows + nr_in + fncols + nc) ;

		DEBUG2 (("new_LUsize "ID"\n", new_LUsize)) ;

		/* RELAX3 parameter uses a double relop, ignore NaN case: */
		do_update = fnpiv > 0 &&
		    (((double) fnzeros) / ((double) new_LUsize)) > RELAX3 ;
		DEBUG2 (("do_update "ID"\n", do_update))

	    }
	    else
	    {
		/* the current front would not be extended */
		do_update = fnpiv > 0 ;
		fnzeros = 0 ;
		DEBUG2 (("IN-OUT do_update forced true: "ID"\n", do_update)) ;

		/* The new front would be just big enough to hold the new
		 * pivot row and column. */
		fnrows_new [IN][OUT] = cdeg_in - 1 ;
		fncols_new [IN][OUT] = rdeg [IN][OUT] - 1 ;

	    }

	    DEBUG2 (("option IN  OUT: nr "ID" nc "ID" cost "ID"("ID") relax "ID
		"\n", nr_in, nc, thiscost, extra_zeros, do_extend)) ;

	    if (bestcost == EMPTY || thiscost < bestcost)
	    {
		/* this is the best option seen so far */
		Work->pivot_case = IN_OUT ;
		bestcost = thiscost ;
		Work->do_extend = do_extend ;
		Work->do_update = do_update ;
		new_fnzeros = fnzeros ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* construct candidate column not in front, and search for pivot rows */
    /* ---------------------------------------------------------------------- */

    search_pivcol_out = (bestcost != 0 && pivcol [OUT] != EMPTY) ;
    if (Symbolic->prefer_diagonal)
    {
	search_pivcol_out = search_pivcol_out && (pivrow [IN][IN] == EMPTY) ;
    }

    if (search_pivcol_out)
    {

#ifndef NDEBUG
	DEBUG2 (("out_col column "ID" NOT in front at position = "ID"\n",
		pivcol [OUT], Fcpos [pivcol [OUT]])) ;
	UMF_dump_rowcol (1, Numeric, Work, pivcol [OUT], !Symbolic->fixQ) ;
	DEBUG2 (("fncols "ID" fncols_max "ID"\n", fncols, Work->fncols_max)) ;
	ASSERT (fncols < Work->fncols_max) ;
#endif

	/* Use Wx as temporary workspace to construct the pivcol [OUT] */


	/* ------------------------------------------------------------------ */
	/* construct the candidate column (currently not in the front) */
	/* ------------------------------------------------------------------ */

	/* Construct the column in Wx, Wm, using Wp for the positions: */
	/* Wm [0..cdeg_out-1]	list of row indices in the column */
	/* Wx [0..cdeg_out-1]	list of corresponding numerical values */
	/* Wp [0..n-1] starts as all negative, and ends that way too. */

	cdeg_out = 0 ;

#ifndef NDEBUG
	/* check Wp */
	DEBUG4 (("COL ASSEMBLE: cdeg 0\nREDUCE COL out "ID"\n", pivcol [OUT])) ;
	if (UMF_debug > 0 || MAX (n_row, n_col) < 1000)
	{
	    for (i = 0 ; i < MAX (n_row, n_col) ; i++)
	    {
		ASSERT (Wp [i] < 0) ;
	    }
	}
	DEBUG4 (("max_cdeg: "ID"\n", max_cdeg)) ;
#endif

	ASSERT (pivcol [OUT] >= 0 && pivcol [OUT] < n_col) ;
	ASSERT (NON_PIVOTAL_COL (pivcol [OUT])) ;

	tpi = Col_tuples [pivcol [OUT]] ;
	if (tpi)
	{
	    tp = (Tuple *) (Memory + tpi) ;
	    tp1 = tp ;
	    tp2 = tp ;
	    tpend = tp + Col_tlen [pivcol [OUT]] ;
	    for ( ; tp < tpend ; tp++)
	    {
		e = tp->e ;
		ASSERT (e > 0 && e <= Work->nel) ;
		if (!E [e]) continue ;	/* element already deallocated */
		f = tp->f ;
		p = Memory + E [e] ;
		ep = (Element *) p ;
		p += UNITS (Element, 1) ;
		Cols = (Int *) p ;
		if (Cols [f] == EMPTY) continue ; /* column already assembled */
		ASSERT (pivcol [OUT] == Cols [f]) ;

		Rows = Cols + ep->ncols ;
		nrows = ep->nrows ;
		p += UNITS (Int, ep->ncols + nrows) ;
		C = ((Entry *) p) + f * nrows ;

		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    if (row >= 0) /* skip this if already gone from element */
		    {
			ASSERT (row < n_row) ;
			pos = Wp [row] ;
			if (pos < 0)
			{
			    /* new entry in the pattern - save Wp */
			    ASSERT (cdeg_out < n_row) ;
			    if (cdeg_out >= max_cdeg)
			    {
				/* :: pattern change (cdeg out failure) :: */
				DEBUGm4 (("cdeg out failure\n")) ;
				return (UMFPACK_ERROR_different_pattern) ;
			    }
			    Wp [row] = cdeg_out ;
			    Wm [cdeg_out] = row ;
			    Wx [cdeg_out++] = C [i] ;
			}
			else
			{
			    /* entry already in pattern - sum the values */
			    /* Wx [pos] += C [i] ; */
			    ASSEMBLE (Wx [pos], C [i]) ;
			}
		    }
		}
		*tp2++ = *tp ;	/* leave the tuple in the list */
	    }
	    Col_tlen [pivcol [OUT]] = tp2 - tp1 ;
	}

	/* ------------------------------------------------------------------ */

#ifndef NDEBUG
	DEBUG4 (("Reduced column: cdeg out "ID"\n", cdeg_out)) ;
	for (i = 0 ; i < cdeg_out ; i++)
	{
	    DEBUG4 ((" "ID" "ID" "ID, i, Wm [i], Wp [Wm [i]])) ;
	    EDEBUG4 (Wx [i]) ;
	    DEBUG4 (("\n")) ;
	    ASSERT (i == Wp [Wm [i]]) ;
	}
#endif

	/* ------------------------------------------------------------------ */
	/* new degree of pivcol [OUT] is cdeg_out */
	/* ------------------------------------------------------------------ */

	/* search for two candidate pivot rows */
	status = UMF_row_search (Numeric, Work, Symbolic,
	    0, cdeg_out, Wm, Wp, /* pattern of column to search */
	    pivrow [OUT], rdeg [OUT], Woi, Woo, pivrow [IN], Wx,
	    pivcol [OUT], freebie) ;

	/* ------------------------------------------------------------------ */
	/* fatal error if matrix pattern has changed since symbolic analysis */
	/* ------------------------------------------------------------------ */

	if (status == UMFPACK_ERROR_different_pattern)
	{
	    /* :: pattern change detected in umf_local_search :: */
	    return (UMFPACK_ERROR_different_pattern) ;
	}

	/* ------------------------------------------------------------------ */
	/* Clear Wp */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < cdeg_out ; i++)
	{
	    Wp [Wm [i]] = EMPTY ;	/* clear Wp */
	}

	/* ------------------------------------------------------------------ */
	/* check for rectangular, singular matrix */
	/* ------------------------------------------------------------------ */

	if (status == UMFPACK_WARNING_singular_matrix)
	{
	    /* Pivot column is empty, and row-merge set is empty too.  The
	     * matrix is structurally singular.  The current frontal matrix must
	     * be empty, too.  It it weren't, and pivcol [OUT] exists, then
	     * there would be at least one row that could be selected.  Since
	     * the current front is empty, pivcol [IN] must also be EMPTY.
	     */

	    DEBUGm4 (("Note: pivcol [OUT]: "ID" discard\n", pivcol [OUT])) ;
	    ASSERT ((Work->fnrows == 0 && Work->fncols == 0)) ;
	    ASSERT (pivcol [IN] == EMPTY) ;

	    /* remove the failed pivcol [OUT] from candidate set */
	    ASSERT (pivcol [OUT] == Work->Candidates [jcand [OUT]]) ;
	    remove_candidate (jcand [OUT], Work, Symbolic) ;
	    Work->ndiscard++ ;

	    /* delete all of the tuples, and all contributions to this column */
	    DEBUG1 (("Prune tuples of dead outcol: "ID"\n", pivcol [OUT])) ;
	    Col_tlen [pivcol [OUT]] = 0 ;
	    UMF_mem_free_tail_block (Numeric, Col_tuples [pivcol [OUT]]) ;
	    Col_tuples [pivcol [OUT]] = 0 ;

	    /* no pivot found at all */
	    return (UMFPACK_WARNING_singular_matrix) ;
	}

	/* ------------------------------------------------------------------ */

	if (freebie [IN])
	{
	    /* the "in" row is the same as the "in" row for the "in" column */
	    Woi = Fcols ;
	    rdeg [OUT][IN] = rdeg [IN][IN] ;
	    DEBUG4 (("Freebie in, row "ID"\n", pivrow [IN][IN])) ;
	}

	if (freebie [OUT])
	{
	    /* the "out" row is the same as the "out" row for the "in" column */
	    Woo = Wio ;
	    rdeg [OUT][OUT] = rdeg [IN][OUT] ;
	    DEBUG4 (("Freebie out, row "ID"\n", pivrow [IN][OUT])) ;
	}

	/* ------------------------------------------------------------------ */
	/* evaluate OUT_IN option */
	/* ------------------------------------------------------------------ */

	if (pivrow [OUT][IN] != EMPTY)
	{
	    /* The current front would become an Lson of the new front.
	     * The candidate pivot row is in the current front, but the
	     * candidate pivot column is not. */

	    ASSERT (fnrows > 0 && fncols >= 0) ;

	    did_rowmerge = (cdeg_out == 0) ;
	    if (did_rowmerge)
	    {
		/* pivrow [OUT][IN] was found via row merge search */
		/* it is not (yet) in the pivot column pattern (add it now) */
		DEBUGm4 (("did row merge OUT col, IN row\n")) ;
		Wm [0] = pivrow [OUT][IN] ;
		CLEAR (Wx [0]) ;
		cdeg_out = 1 ;
		ASSERT (nr_out == EMPTY) ;
	    }

	    nc = rdeg [OUT][IN] - fncols ;
	    ASSERT (nc >= 1) ;

	    /* count rows not in current front */
	    nr_out = 0 ;
#ifndef NDEBUG
	    debug_ok = FALSE ;
#endif
	    for (i = 0 ; i < cdeg_out ; i++)
	    {
		row = Wm [i] ;
		ASSERT (row >= 0 && row < n_row && NON_PIVOTAL_ROW (row)) ;
		if (Frpos [row] < 0 || Frpos [row] >= fnrows) nr_out++ ;
#ifndef NDEBUG
		/* we must see the pivot row somewhere */
		if (row == pivrow [OUT][IN])
		{
		    ASSERT (Frpos [row] >= 0) ;
		    debug_ok = TRUE ;
		}
#endif
	    }
	    ASSERT (debug_ok) ;

	    thiscost =
		/* each column in front grows by nr_out: */
		(nr_out * fncols) +
		/* new cols not affected by front: */
		((nc-1) * (cdeg_out-1)) ;

	    /* check the cost of relaxed OUT_IN amalgamation */

	    extra_rows = ((fnrows-1) + nr_out) - (cdeg_out - 1) ;
	    ASSERT (extra_rows >= 0) ;
	    ASSERT (fnrows + nr_out == extra_rows + cdeg_out) ;
	    extra_zeros = (nc-1) * extra_rows ;	/* symbolic fill-in */

	    ASSERT (fnrows + nr_out == cdeg_out + extra_rows) ;
	    ASSERT (fncols + nc == rdeg [OUT][IN]) ;

	    /* size of relaxed front (after pivot row removed): */
	    fnrows_new [OUT][IN] = (fnrows-1) + nr_out ;
	    fncols_new [OUT][IN] = fncols + (nc-1) ;
	    relaxed_front = fnrows_new [OUT][IN] * fncols_new [OUT][IN] ;

	    /* do relaxed amalgamation if the extra zeros are no more */
	    /* than a fraction (default 0.25) of the relaxed front */
	    /* if relax = 0: no extra zeros allowed */
	    /* if relax = +inf: always amalgamate */
	    if (did_rowmerge)
	    {
		do_extend = FALSE ;
	    }
	    else
	    {
		/* relax parameter uses a double relop, but ignore NaN case: */
		if (extra_zeros == 0)
		{
		    do_extend = TRUE ;
		}
		else
		{
		    do_extend = ((double) extra_zeros) <
		       (relax1 * (double) relaxed_front) ;
		}
	    }

	    if (do_extend)
	    {
		/* count the cost of relaxed amalgamation */
		thiscost += extra_zeros ;

		DEBUG2 (("Evaluating option OUT-IN:\n")) ;
		DEBUG2 ((" Work->fnzeros "ID" fnpiv "ID" nr_out "ID" nc "ID"\n",
		Work->fnzeros, fnpiv, nr_out, nc)) ;
		DEBUG2 (("fncols "ID" fnrows "ID"\n", fncols, fnrows)) ;

		/* determine if BLAS-3 update to be applied before extending. */
		/* update if too many zero entries accumulate in the LU block */
		fnzeros = Work->fnzeros + fnpiv * (nr_out + nc) ;

		DEBUG2 (("fnzeros "ID"\n", fnzeros)) ;

		new_LUsize = (fnpiv+1) * (fnrows + nr_out + fncols + nc) ;

		DEBUG2 (("new_LUsize "ID"\n", new_LUsize)) ;

		/* RELAX3 parameter uses a double relop, ignore NaN case: */
		do_update = fnpiv > 0 &&
		    (((double) fnzeros) / ((double) new_LUsize)) > RELAX3 ;
		DEBUG2 (("do_update "ID"\n", do_update))
	    }
	    else
	    {
		/* the current front would not be extended */
		do_update = fnpiv > 0 ;
		fnzeros = 0 ;
		DEBUG2 (("OUT-IN do_update forced true: "ID"\n", do_update)) ;

		/* The new front would be just big enough to hold the new
		 * pivot row and column. */
		fnrows_new [OUT][IN] = cdeg_out - 1 ;
		fncols_new [OUT][IN] = rdeg [OUT][IN] - 1 ;
	    }

	    DEBUG2 (("option OUT IN : nr "ID" nc "ID" cost "ID"("ID") relax "ID
		"\n", nr_out, nc, thiscost, extra_zeros, do_extend)) ;

	    if (bestcost == EMPTY || thiscost < bestcost)
	    {
		/* this is the best option seen so far */
		Work->pivot_case = OUT_IN ;
		bestcost = thiscost ;
		Work->do_extend = do_extend ;
		Work->do_update = do_update ;
		new_fnzeros = fnzeros ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* evaluate OUT_OUT option */
	/* ------------------------------------------------------------------ */

	if (pivrow [OUT][OUT] != EMPTY)
	{
	    /* Neither the candidate pivot row nor the candidate pivot column
	     * are in the current front. */

	    ASSERT (fnrows >= 0 && fncols >= 0) ;

	    did_rowmerge = (cdeg_out == 0) ;
	    if (did_rowmerge)
	    {
		/* pivrow [OUT][OUT] was found via row merge search */
		/* it is not (yet) in the pivot column pattern (add it now) */
		DEBUGm4 (("did row merge OUT col, OUT row\n")) ;
		Wm [0] = pivrow [OUT][OUT] ;
		CLEAR (Wx [0]) ;
		cdeg_out = 1 ;
		ASSERT (nr_out == EMPTY) ;
		nr_out = 1 ;
	    }

	    if (fnrows == 0 && fncols == 0)
	    {
		/* the current front is completely empty */
		ASSERT (fnpiv == 0) ;
		nc = rdeg [OUT][OUT] ;
		extra_cols = 0 ;
		nr_out = cdeg_out ;
		extra_rows = 0 ;
		extra_zeros = 0 ;

		thiscost = (nc-1) * (cdeg_out-1) ; /* new columns only */

		/* size of new front: */
		fnrows_new [OUT][OUT] = nr_out-1 ;
		fncols_new [OUT][OUT] = nc-1 ;
		relaxed_front = fnrows_new [OUT][OUT] * fncols_new [OUT][OUT] ;
	    }
	    else
	    {

		/* count rows not in current front */
		if (nr_out == EMPTY)
		{
		    nr_out = 0 ;
#ifndef NDEBUG
		    debug_ok = FALSE ;
#endif
		    for (i = 0 ; i < cdeg_out ; i++)
		    {
			row = Wm [i] ;
			ASSERT (row >= 0 && row < n_row) ;
			ASSERT (NON_PIVOTAL_ROW (row)) ;
			if (Frpos [row] < 0 || Frpos [row] >= fnrows) nr_out++ ;
#ifndef NDEBUG
			/* we must see the pivot row somewhere */
			if (row == pivrow [OUT][OUT])
			{
			    ASSERT (Frpos [row] < 0 || Frpos [row] >= fnrows) ;
			    debug_ok = TRUE ;
			}
#endif
		    }
		    ASSERT (debug_ok) ;
		}

		/* count columns not in current front */
		nc = 0 ;
#ifndef NDEBUG
		debug_ok = FALSE ;
#endif
		for (i = 0 ; i < rdeg [OUT][OUT] ; i++)
		{
		    col = Woo [i] ;
		    ASSERT (col >= 0 && col < n_col && NON_PIVOTAL_COL (col)) ;
		    if (Fcpos [col] < 0) nc++ ;
#ifndef NDEBUG
		    /* we must see the pivot column somewhere */
		    if (col == pivcol [OUT])
		    {
			ASSERT (Fcpos [col] < 0) ;
			debug_ok = TRUE ;
		    }
#endif
		}
		ASSERT (debug_ok) ;

		extra_cols = (fncols + (nc-1)) - (rdeg [OUT][OUT] - 1) ;
		extra_rows = (fnrows + (nr_out-1)) - (cdeg_out - 1) ;
		ASSERT (extra_rows >= 0) ;
		ASSERT (extra_cols >= 0) ;
		extra_zeros = ((nc-1) * extra_rows) + ((nr_out-1) * extra_cols);

		ASSERT (fnrows + nr_out == cdeg_out + extra_rows) ;
		ASSERT (fncols + nc == rdeg [OUT][OUT] + extra_cols) ;

		thiscost =
		    /* new columns: */
		    ((nc-1) * (cdeg_out-1)) +
		    /* old columns in front grow by nr_out-1: */
		    ((nr_out-1) * (fncols - extra_cols)) ;

		/* size of relaxed front: */
		fnrows_new [OUT][OUT] = fnrows + (nr_out-1) ;
		fncols_new [OUT][OUT] = fncols + (nc-1) ;
		relaxed_front = fnrows_new [OUT][OUT] * fncols_new [OUT][OUT] ;

	    }

	    /* do relaxed amalgamation if the extra zeros are no more */
	    /* than a fraction (default 0.25) of the relaxed front */
	    /* if relax = 0: no extra zeros allowed */
	    /* if relax = +inf: always amalgamate */
	    if (did_rowmerge)
	    {
		do_extend = FALSE ;
	    }
	    else
	    {
		/* relax parameter uses a double relop, but ignore NaN case: */
		if (extra_zeros == 0)
		{
		    do_extend = TRUE ;
		}
		else
		{
		    do_extend = ((double) extra_zeros) <
		       (relax1 * (double) relaxed_front) ;
		}
	    }

	    if (do_extend)
	    {
		/* count the cost of relaxed amalgamation */
		thiscost += extra_zeros ;

		DEBUG2 (("Evaluating option OUT-OUT:\n")) ;
		DEBUG2 (("Work->fnzeros "ID" fnpiv "ID" nr_out "ID" nc "ID"\n",
		    Work->fnzeros, fnpiv, nr_out, nc)) ;
		DEBUG2 (("fncols "ID" fnrows "ID"\n", fncols, fnrows)) ;

		/* determine if BLAS-3 update to be applied before extending. */
		/* update if too many zero entries accumulate in the LU block */
		fnzeros = Work->fnzeros + fnpiv * (nr_out + nc) ;

		DEBUG2 (("fnzeros "ID"\n", fnzeros)) ;

		new_LUsize = (fnpiv+1) * (fnrows + nr_out + fncols + nc) ;

		DEBUG2 (("new_LUsize "ID"\n", new_LUsize)) ;

		/* RELAX3 parameter uses a double relop, ignore NaN case: */
		do_update = fnpiv > 0 &&
		    (((double) fnzeros) / ((double) new_LUsize)) > RELAX3 ;
		DEBUG2 (("do_update "ID"\n", do_update))
	    }
	    else
	    {
		/* the current front would not be extended */
		do_update = fnpiv > 0 ;
		fnzeros = 0 ;
		DEBUG2 (("OUT-OUT do_update forced true: "ID"\n", do_update)) ;

		/* The new front would be just big enough to hold the new
		 * pivot row and column. */
		fnrows_new [OUT][OUT] = cdeg_out - 1 ;
		fncols_new [OUT][OUT] = rdeg [OUT][OUT] - 1 ;
	    }

	    DEBUG2 (("option OUT OUT: nr "ID" nc "ID" cost "ID"\n",
		rdeg [OUT][OUT], cdeg_out, thiscost)) ;

	    if (bestcost == EMPTY || thiscost < bestcost)
	    {
		/* this is the best option seen so far */
		Work->pivot_case = OUT_OUT ;
		bestcost = thiscost ;
		Work->do_extend = do_extend ;
		Work->do_update = do_update ;
		new_fnzeros = fnzeros ;
	    }
	}
    }

    /* At this point, a structural pivot has been found. */
    /* It may be numerically zero, however. */
    ASSERT (Work->pivot_case != EMPTY) ;
    DEBUG2 (("local search, best option "ID", best cost "ID"\n",
	Work->pivot_case, bestcost)) ;

    /* ====================================================================== */
    /* Pivot row and column, and extension, now determined */
    /* ====================================================================== */

    Work->fnzeros = new_fnzeros ;

    /* ---------------------------------------------------------------------- */
    /* finalize the pivot row and column */
    /* ---------------------------------------------------------------------- */

    switch (Work->pivot_case)
    {
	case IN_IN:
	    DEBUG2 (("IN-IN option selected\n")) ;
	    ASSERT (fnrows > 0 && fncols > 0) ;
	    Work->pivcol_in_front = TRUE ;
	    Work->pivrow_in_front = TRUE ;
	    Work->pivcol = pivcol [IN] ;
	    Work->pivrow = pivrow [IN][IN] ;
	    Work->ccdeg = nr_in ;
	    Work->Wrow = Fcols ;
	    Work->rrdeg = rdeg [IN][IN] ;
	    jj = jcand [IN] ;
	    Work->fnrows_new = fnrows_new [IN][IN] ;
	    Work->fncols_new = fncols_new [IN][IN] ;
	    break ;

	case IN_OUT:
	    DEBUG2 (("IN-OUT option selected\n")) ;
	    ASSERT (fnrows >= 0 && fncols > 0) ;
	    Work->pivcol_in_front = TRUE ;
	    Work->pivrow_in_front = FALSE ;
	    Work->pivcol = pivcol [IN] ;
	    Work->pivrow = pivrow [IN][OUT] ;
	    Work->ccdeg = nr_in ;
	    Work->Wrow = Wio ;
	    Work->rrdeg = rdeg [IN][OUT] ;
	    jj = jcand [IN] ;
	    Work->fnrows_new = fnrows_new [IN][OUT] ;
	    Work->fncols_new = fncols_new [IN][OUT] ;
	    break ;

	case OUT_IN:
	    DEBUG2 (("OUT-IN option selected\n")) ;
	    ASSERT (fnrows > 0 && fncols >= 0) ;
	    Work->pivcol_in_front = FALSE ;
	    Work->pivrow_in_front = TRUE ;
	    Work->pivcol = pivcol [OUT] ;
	    Work->pivrow = pivrow [OUT][IN] ;
	    Work->ccdeg = cdeg_out ;
	    /* Wrow might be equivalenced to Fcols (Freebie in): */
	    Work->Wrow = Woi ;
	    Work->rrdeg = rdeg [OUT][IN] ;
	    /* Work->Wrow[0..fncols-1] is not there.  See Fcols instead */
	    jj = jcand [OUT] ;
	    Work->fnrows_new = fnrows_new [OUT][IN] ;
	    Work->fncols_new = fncols_new [OUT][IN] ;
	    break ;

	case OUT_OUT:
	    DEBUG2 (("OUT-OUT option selected\n")) ;
	    ASSERT (fnrows >= 0 && fncols >= 0) ;
	    Work->pivcol_in_front = FALSE ;
	    Work->pivrow_in_front = FALSE ;
	    Work->pivcol = pivcol [OUT] ;
	    Work->pivrow = pivrow [OUT][OUT] ;
	    Work->ccdeg = cdeg_out ;
	    /* Wrow might be equivalenced to Wio (Freebie out): */
	    Work->Wrow = Woo ;
	    Work->rrdeg = rdeg [OUT][OUT] ;
	    jj = jcand [OUT] ;
	    Work->fnrows_new = fnrows_new [OUT][OUT] ;
	    Work->fncols_new = fncols_new [OUT][OUT] ;
	    break ;

    }

    ASSERT (IMPLIES (fnrows == 0 && fncols == 0, Work->pivot_case == OUT_OUT)) ;

    if (!Work->pivcol_in_front && pivcol [IN] != EMPTY)
    {
	/* clear Frpos if pivcol [IN] was searched, but not selected */
	for (i = fnrows ; i < cdeg_in ; i++)
	{
	    Frpos [Frows [i]] = EMPTY;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* Pivot row and column have been found */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* remove pivot column from candidate pivot column set */
    /* ---------------------------------------------------------------------- */

    ASSERT (jj >= 0 && jj < Work->nCandidates) ;
    ASSERT (Work->pivcol == Work->Candidates [jj]) ;
    remove_candidate (jj, Work, Symbolic) ;

    /* ---------------------------------------------------------------------- */
    /* check for frontal matrix growth */
    /* ---------------------------------------------------------------------- */

    DEBUG1 (("Check frontal growth:\n")) ;
    DEBUG1 (("fnrows_new "ID" + 1 = "ID",  fnr_curr "ID"\n",
	    Work->fnrows_new, Work->fnrows_new + 1, fnr_curr)) ;
    DEBUG1 (("fncols_new "ID" + 1 = "ID",  fnc_curr "ID"\n",
	    Work->fncols_new, Work->fncols_new + 1, fnc_curr)) ;

    Work->do_grow = (Work->fnrows_new + 1 > fnr_curr
		  || Work->fncols_new + 1 > fnc_curr) ;
    if (Work->do_grow)
    {
	DEBUG0 (("\nNeed to grow frontal matrix, force do_update true\n")) ;
	/* If the front must grow, then apply the pending updates and remove
	 * the current pivot rows/columns from the front prior to growing the
	 * front.  This frees up as much space as possible for the new front. */
	if (!Work->do_update && fnpiv > 0)
	{
	    /* This update would not have to be done if the current front
	     * was big enough. */
	    Work->nforced++ ;
	    Work->do_update = TRUE ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* current pivot column */
    /* ---------------------------------------------------------------------- */

    /*
	c1) If pivot column index is in the current front:

	    The pivot column pattern is in Frows [0 .. fnrows-1] and
	    the extension is in Frows [fnrows ... fnrows+ccdeg-1].

	    Frpos [Frows [0 .. fnrows+ccdeg-1]] is
	    equal to 0 .. fnrows+ccdeg-1.  Wm is not needed.

	    The values are in Wy [0 .. fnrows+ccdeg-1].

	c2) Otherwise, if the pivot column index is not in the current front:

	    c2a) If the front is being extended, old row indices in the the
		pivot column pattern are in Frows [0 .. fnrows-1].

		All entries are in Wm [0 ... ccdeg-1], with values in
		Wx [0 .. ccdeg-1].  These may include entries already in
		Frows [0 .. fnrows-1].

		Frpos [Frows [0 .. fnrows-1]] is equal to 0 .. fnrows-1.
		Frpos [Wm [0 .. ccdeg-1]] for new entries is < 0.

	    c2b) If the front is not being extended, then the entire pivot
		column pattern is in Wm [0 .. ccdeg-1].  It includes
		the pivot row index.  It is does not contain the pattern
		Frows [0..fnrows-1].  The intersection of these two
		sets may or may not be empty.  The values are in Wx [0..ccdeg-1]

	In both cases c1 and c2, Frpos [Frows [0 .. fnrows-1]] is equal
	to 0 .. fnrows-1, which is the pattern of the current front.
	Any entry of Frpos that is not specified above is < 0.
    */


#ifndef NDEBUG
    DEBUG2 (("\n\nSEARCH DONE: Pivot col "ID" in: ("ID") pivot row "ID" in: ("ID
	") extend: "ID"\n\n", Work->pivcol, Work->pivcol_in_front,
	Work->pivrow, Work->pivrow_in_front, Work->do_extend)) ;
    UMF_dump_rowcol (1, Numeric, Work, Work->pivcol, !Symbolic->fixQ) ;
    DEBUG2 (("Pivot col "ID": fnrows "ID" ccdeg "ID"\n", Work->pivcol, fnrows,
	Work->ccdeg)) ;
    if (Work->pivcol_in_front)	/* case c1 */
    {
	Int found = FALSE ;
	DEBUG3 (("Pivcol in front\n")) ;
	for (i = 0 ; i < fnrows ; i++)
	{
	    row = Frows [i] ;
	    DEBUG3 ((ID":   row:: "ID" in front ", i, row)) ;
	    ASSERT (row >= 0 && row < n_row && NON_PIVOTAL_ROW (row)) ;
	    ASSERT (Frpos [row] == i) ;
	    EDEBUG3 (Wy [i]) ;
	    if (row == Work->pivrow)
	    {
		DEBUG3 ((" <- pivrow")) ;
		found = TRUE ;
	    }
	    DEBUG3 (("\n")) ;
	}
	ASSERT (found == Work->pivrow_in_front) ;
	found = FALSE ;
	for (i = fnrows ; i < fnrows + Work->ccdeg ; i++)
	{
	    row = Frows [i] ;
	    DEBUG3 ((ID":   row:: "ID" (new)", i, row)) ;
	    ASSERT (row >= 0 && row < n_row && NON_PIVOTAL_ROW (row)) ;
	    ASSERT (Frpos [row] == i) ;
	    EDEBUG3 (Wy [i]) ;
	    if (row == Work->pivrow)
	    {
		DEBUG3 ((" <- pivrow")) ;
		found = TRUE ;
	    }
	    DEBUG3 (("\n")) ;
	}
	ASSERT (found == !Work->pivrow_in_front) ;
    }
    else
    {
	if (Work->do_extend)
	{
	    Int found = FALSE ;
	    DEBUG3 (("Pivcol not in front (extend)\n")) ;
	    for (i = 0 ; i < fnrows ; i++)
	    {
		row = Frows [i] ;
		DEBUG3 ((ID":   row:: "ID" in front ", i, row)) ;
		ASSERT (row >= 0 && row < n_row && NON_PIVOTAL_ROW (row)) ;
		ASSERT (Frpos [row] == i) ;
		if (row == Work->pivrow)
		{
		    DEBUG3 ((" <- pivrow")) ;
		    found = TRUE ;
		}
		DEBUG3 (("\n")) ;
	    }
	    ASSERT (found == Work->pivrow_in_front) ;
	    found = FALSE ;
	    DEBUG3 (("----\n")) ;
	    for (i = 0 ; i < Work->ccdeg ; i++)
	    {
		row = Wm [i] ;
		ASSERT (row >= 0 && row < n_row && NON_PIVOTAL_ROW (row)) ;
		DEBUG3 ((ID":   row:: "ID" ", i, row)) ;
		EDEBUG3 (Wx [i]) ;
		if (Frpos [row] < 0)
		{
		    DEBUG3 ((" (new) ")) ;
		}
		if (row == Work->pivrow)
		{
		    DEBUG3 ((" <- pivrow")) ;
		    found = TRUE ;
		    /* ... */
		    if (Work->pivrow_in_front) ASSERT (Frpos [row] >= 0) ;
		    else ASSERT (Frpos [row] < 0) ;
		}
		DEBUG3 (("\n")) ;
	    }
	    ASSERT (found) ;
	}
	else
	{
	    Int found = FALSE ;
	    DEBUG3 (("Pivcol not in front (no extend)\n")) ;
	    for (i = 0 ; i < Work->ccdeg ; i++)
	    {
		row = Wm [i] ;
		ASSERT (row >= 0 && row < n_row && NON_PIVOTAL_ROW (row)) ;
		DEBUG3 ((ID":   row:: "ID" ", i, row)) ;
		EDEBUG3 (Wx [i]) ;
		DEBUG3 ((" (new) ")) ;
		if (row == Work->pivrow)
		{
		    DEBUG3 ((" <- pivrow")) ;
		    found = TRUE ;
		}
		DEBUG3 (("\n")) ;
	    }
	    ASSERT (found) ;
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* current pivot row */
    /* ---------------------------------------------------------------------- */

    /*
	r1) If the pivot row index is in the current front:

	    The pivot row pattern is in Fcols [0..fncols-1] and the extenson is
	    in Wrow [fncols .. rrdeg-1].  If the pivot column is in the current
	    front, then Fcols and Wrow are equivalenced.

	r2)  If the pivot row index is not in the current front:

	    r2a) If the front is being extended, the pivot row pattern is in
		Fcols [0 .. fncols-1].  New entries are in Wrow [0 .. rrdeg-1],
		but these may include entries already in Fcols [0 .. fncols-1].

	    r2b) Otherwise, the pivot row pattern is Wrow [0 .. rrdeg-1].

	Fcpos [Fcols [0..fncols-1]] is (0..fncols-1) * fnr_curr.
	All other entries in Fcpos are < 0.

	These conditions are asserted below.

	------------------------------------------------------------------------
	Other items in Work structure that are relevant:

	pivcol: the pivot column index
	pivrow: the pivot column index

	rrdeg:
	ccdeg:

	fnrows: the number of rows in the currnt contribution block
	fncols: the number of columns in the current contribution block

	fnrows_new: the number of rows in the new contribution block
	fncols_new: the number of rows in the new contribution block

	------------------------------------------------------------------------
    */


#ifndef NDEBUG
    UMF_dump_rowcol (0, Numeric, Work, Work->pivrow, TRUE) ;
    DEBUG2 (("Pivot row "ID":\n", Work->pivrow)) ;
    if (Work->pivrow_in_front)
    {
	Int found = FALSE ;
	for (i = 0 ; i < fncols ; i++)
	{
	    col = Fcols [i] ;
	    DEBUG3 (("   col:: "ID" in front\n", col)) ;
	    ASSERT (col >= 0 && col < n_col && NON_PIVOTAL_COL (col)) ;
	    ASSERT (Fcpos [col] == i * fnr_curr) ;
	    if (col == Work->pivcol) found = TRUE ;
	}
	ASSERT (found == Work->pivcol_in_front) ;
	found = FALSE ;
	ASSERT (IMPLIES (Work->pivcol_in_front, Fcols == Work->Wrow)) ;
	for (i = fncols ; i < Work->rrdeg ; i++)
	{
	    col = Work->Wrow [i] ;
	    ASSERT (col >= 0 && col < n_col && NON_PIVOTAL_COL (col)) ;
	    ASSERT (Fcpos [col] < 0) ;
	    if (col == Work->pivcol) found = TRUE ;
	    else DEBUG3 (("   col:: "ID" (new)\n", col)) ;
	}
	ASSERT (found == !Work->pivcol_in_front) ;
    }
    else
    {
	if (Work->do_extend)
	{
	    Int found = FALSE ;
	    for (i = 0 ; i < fncols ; i++)
	    {
		col = Fcols [i] ;
		DEBUG3 (("   col:: "ID" in front\n", col)) ;
		ASSERT (col >= 0 && col < n_col && NON_PIVOTAL_COL (col)) ;
		ASSERT (Fcpos [col] == i * fnr_curr) ;
		if (col == Work->pivcol) found = TRUE ;
	    }
	    ASSERT (found == Work->pivcol_in_front) ;
	    found = FALSE ;
	    for (i = 0 ; i < Work->rrdeg ; i++)
	    {
		col = Work->Wrow [i] ;
		ASSERT (col >= 0 && col < n_col && NON_PIVOTAL_COL (col)) ;
		if (Fcpos [col] >= 0) continue ;
		if (col == Work->pivcol) found = TRUE ;
		else DEBUG3 (("   col:: "ID" (new, extend)\n", col)) ;
	    }
	    ASSERT (found == !Work->pivcol_in_front) ;
	}
	else
	{
	    Int found = FALSE ;
	    for (i = 0 ; i < Work->rrdeg ; i++)
	    {
		col = Work->Wrow [i] ;
		ASSERT (col >= 0 && col < n_col && NON_PIVOTAL_COL (col)) ;
		if (col == Work->pivcol) found = TRUE ;
		else DEBUG3 (("   col:: "ID" (all new)\n", col)) ;
	    }
	    ASSERT (found) ;
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* determine whether to do scan2-row and scan2-col */
    /* ---------------------------------------------------------------------- */

    if (Work->do_extend)
    {
	Work->do_scan2row = (fncols > 0) ;
	Work->do_scan2col = (fnrows > 0) ;
    }
    else
    {
	Work->do_scan2row = (fncols > 0) && Work->pivrow_in_front ;
	Work->do_scan2col = (fnrows > 0) && Work->pivcol_in_front ;
    }

    /* ---------------------------------------------------------------------- */

    DEBUG2 (("LOCAL SEARCH DONE: pivot column "ID" pivot row: "ID"\n",
	Work->pivcol, Work->pivrow)) ;
    DEBUG2 (("do_extend: "ID"\n", Work->do_extend)) ;
    DEBUG2 (("do_update: "ID"\n", Work->do_update)) ;
    DEBUG2 (("do_grow:   "ID"\n", Work->do_grow)) ;

    /* ---------------------------------------------------------------------- */
    /* keep track of the diagonal */
    /* ---------------------------------------------------------------------- */

    if (Symbolic->prefer_diagonal
	&& Work->pivcol < Work->n_col - Symbolic->nempty_col)
    {
	Diagonal_map = Work->Diagonal_map ;
	Diagonal_imap = Work->Diagonal_imap ;
	ASSERT (Diagonal_map != (Int *) NULL) ;
	ASSERT (Diagonal_imap != (Int *) NULL) ;

	row2 = Diagonal_map  [Work->pivcol] ;
	col2 = Diagonal_imap [Work->pivrow] ;

	if (row2 < 0)
	{
	    /* this was an off-diagonal pivot row */
	    Work->noff_diagonal++ ;
	    row2 = UNFLIP (row2) ;
	}

	ASSERT (Diagonal_imap [row2] == Work->pivcol) ;
	ASSERT (UNFLIP (Diagonal_map [col2]) == Work->pivrow) ;

	if (row2 != Work->pivrow)
	{
	    /* swap the diagonal map to attempt to maintain symmetry later on.
	     * Also mark the map for col2 (via FLIP) to denote that the entry
	     * now on the diagonal is not the original entry on the diagonal. */

	    DEBUG0 (("Unsymmetric pivot\n")) ;
	    Diagonal_map  [Work->pivcol] = FLIP (Work->pivrow) ;
	    Diagonal_imap [Work->pivrow] = Work->pivcol ;

	    Diagonal_map  [col2] = FLIP (row2) ;
	    Diagonal_imap [row2] = col2 ;

	}
	ASSERT (n_row == n_col) ;
#ifndef NDEBUG
	UMF_dump_diagonal_map (Diagonal_map, Diagonal_imap, Symbolic->n1,
	    Symbolic->n_col, Symbolic->nempty_col) ;
#endif
    }

    return (UMFPACK_OK) ;
}
