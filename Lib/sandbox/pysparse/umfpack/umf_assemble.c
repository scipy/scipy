/* ========================================================================== */
/* === UMF_assemble ========================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*  Degree update and numerical assembly.  This is compiled twice (with and
 *  without FIXQ) for each real/complex int/long version, for a total of 8
 *  versions.*/

#include "umf_internal.h"
#include "umf_mem_free_tail_block.h"

/* ========================================================================== */
/* === row_assemble ========================================================= */
/* ========================================================================== */

PRIVATE void row_assemble
(
    Int row,
    NumericType *Numeric,
    WorkType *Work
)
{

    Int tpi, e, *E, *Fcpos, *Frpos, *Row_degree, *Row_tuples, *Row_tlen, rdeg0,
	f, nrows, ncols, *Rows, *Cols, col, ncolsleft, j ;
    Tuple *tp, *tp1, *tp2, *tpend ;
    Unit *Memory, *p ;
    Element *ep ;
    Entry *S, *Fcblock, *Frow ;

#ifndef FIXQ
    Int *Col_degree ;
    Col_degree = Numeric->Cperm ;
#endif

    Row_tuples = Numeric->Uip ;
    tpi = Row_tuples [row] ;
    if (!tpi) return ;

    Memory = Numeric->Memory ;
    E = Work->E ;
    Fcpos = Work->Fcpos ;
    Frpos = Work->Frpos ;
    Row_degree = Numeric->Rperm ;
    Row_tlen   = Numeric->Uilen ;
    E = Work->E ;
    Memory = Numeric->Memory ;
    rdeg0 = Work->rdeg0 ;
    Fcblock = Work->Fcblock ;

#ifndef NDEBUG
    DEBUG6 (("SCAN2-row: "ID"\n", row)) ;
    UMF_dump_rowcol (0, Numeric, Work, row, FALSE) ;
#endif

    ASSERT (NON_PIVOTAL_ROW (row)) ;

    tp = (Tuple *) (Memory + tpi) ;
    tp1 = tp ;
    tp2 = tp ;
    tpend = tp + Row_tlen [row] ;
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
	Rows = Cols + ep->ncols ;
	if (Rows [f] == EMPTY) continue ;   /* row already assembled */
	ASSERT (row == Rows [f] && row >= 0 && row < Work->n_row) ;

	if (ep->rdeg == rdeg0)
	{
	    /* ------------------------------------------------------ */
	    /* this is an old Lson - assemble just one row */
	    /* ------------------------------------------------------ */

	    /* flag the row as assembled from the Lson */
	    Rows [f] = EMPTY ;

	    nrows = ep->nrows ;
	    ncols = ep->ncols ;

	    p += UNITS (Int, ncols + nrows) ;
	    S = ((Entry *) p) + f ;

	    DEBUG6 (("Old LSON: "ID"\n", e)) ;
#ifndef NDEBUG
	    UMF_dump_element (Numeric, Work, e, FALSE) ;
#endif

	    ncolsleft = ep->ncolsleft ;

	    Frow = Fcblock + Frpos [row] ;
	    DEBUG6 (("LSON found (in scan2-row): "ID"\n", e)) ;

	    Row_degree [row] -= ncolsleft ;

	    if (ncols == ncolsleft)
	    {
		/* -------------------------------------------------- */
		/* no columns assembled out this Lson yet */
		/* -------------------------------------------------- */

#pragma ivdep
		for (j = 0 ; j < ncols ; j++)
		{
		    col = Cols [j] ;
		    ASSERT (col >= 0 && col < Work->n_col) ;
#ifndef FIXQ
		    Col_degree [col] -- ;
#endif
		    /* Frow [Fcpos [col]] += *S ; */
		    ASSEMBLE (Frow [Fcpos [col]], *S) ;
		    S += nrows ;
		}

	    }
	    else
	    {
		/* -------------------------------------------------- */
		/* some columns have been assembled out of this Lson */
		/* -------------------------------------------------- */

#pragma ivdep
		for (j = 0 ; j < ncols ; j++)
		{
		    col = Cols [j] ;
		    if (col >= 0)
		    {
			ASSERT (col < Work->n_col) ;
#ifndef FIXQ
			Col_degree [col] -- ;
#endif
			/* Frow [Fcpos [col]] += *S ; */
			ASSEMBLE (Frow [Fcpos [col]], *S) ;
		    }
		    S += nrows ;
		}

	    }
	    ep->nrowsleft-- ;
	    ASSERT (ep->nrowsleft > 0) ;
	}
	else
	{
	    *tp2++ = *tp ;	/* leave the tuple in the list */
	}
    }
    Row_tlen [row] = tp2 - tp1 ;

#ifndef NDEBUG
    DEBUG7 (("row assembled in scan2-row: "ID"\n", row)) ;
    UMF_dump_rowcol (0, Numeric, Work, row, FALSE) ;
    DEBUG7 (("Current frontal matrix: (scan 1b)\n")) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif
}

/* ========================================================================== */
/* === col_assemble ========================================================= */
/* ========================================================================== */

PRIVATE void col_assemble
(
    Int col,
    NumericType *Numeric,
    WorkType *Work
)
{

    Int tpi, e, *E, *Fcpos, *Frpos, *Row_degree, *Col_tuples, *Col_tlen, cdeg0,
	f, nrows, ncols, *Rows, *Cols, row, nrowsleft, i ;
    Tuple *tp, *tp1, *tp2, *tpend ;
    Unit *Memory, *p ;
    Element *ep ;
    Entry *S, *Fcblock, *Fcol ;

#if !defined (FIXQ) || !defined (NDEBUG)
    Int *Col_degree ;
    Col_degree = Numeric->Cperm ;
#endif

    Col_tuples = Numeric->Lip ;
    tpi = Col_tuples [col] ;
    if (!tpi) return ;

    Memory = Numeric->Memory ;
    E = Work->E ;
    Fcpos = Work->Fcpos ;
    Frpos = Work->Frpos ;
    Row_degree = Numeric->Rperm ;
    Col_tlen   = Numeric->Lilen ;
    E = Work->E ;
    Memory = Numeric->Memory ;
    cdeg0 = Work->cdeg0 ;
    Fcblock = Work->Fcblock ;

    DEBUG6 (("SCAN2-col: "ID"\n", col)) ;
#ifndef NDEBUG
    UMF_dump_rowcol (1, Numeric, Work, col, FALSE) ;
#endif

    ASSERT (NON_PIVOTAL_COL (col)) ;
    tp = (Tuple *) (Memory + tpi) ;
    tp1 = tp ;
    tp2 = tp ;
    tpend = tp + Col_tlen [col] ;
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

	if (Cols [f] == EMPTY) continue ;   /* col already assembled */
	ASSERT (col == Cols [f] && col >= 0 && col < Work->n_col) ;

	if (ep->cdeg == cdeg0)
	{
	    /* ------------------------------------------------------ */
	    /* this is an old Uson - assemble just one col */
	    /* ------------------------------------------------------ */

	    /* flag the col as assembled from the Uson */
	    Cols [f] = EMPTY ;

	    nrows = ep->nrows ;
	    ncols = ep->ncols ;
	    Rows = Cols + ncols ;
	    p += UNITS (Int, ncols + nrows) ;
	    S = ((Entry *) p) + f * nrows ;

	    DEBUG6 (("Old USON: "ID"\n", e)) ;
#ifndef NDEBUG
	    UMF_dump_element (Numeric, Work, e, FALSE) ;
#endif

	    nrowsleft = ep->nrowsleft ;

	    Fcol = Fcblock + Fcpos [col] ;
	    DEBUG6 (("USON found (in scan2-col): "ID"\n", e)) ;
#ifndef FIXQ
	    Col_degree [col] -= nrowsleft ;
#endif
	    if (nrows == nrowsleft)
	    {
		/* -------------------------------------------------- */
		/* no rows assembled out of this Uson yet */
		/* -------------------------------------------------- */

#pragma ivdep
		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    ASSERT (row >= 0 && row < Work->n_row) ;
		    Row_degree [row]-- ;
		    /* Fcol [Frpos [row]] += S [i] ; */
		    ASSEMBLE (Fcol [Frpos [row]], S [i]) ;
		}
	    }
	    else
	    {
		/* -------------------------------------------------- */
		/* some rows have been assembled out of this Uson */
		/* -------------------------------------------------- */

#pragma ivdep
		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    if (row >= 0)
		    {
			ASSERT (row < Work->n_row) ;
			Row_degree [row]-- ;
			/* Fcol [Frpos [row]] += S [i] ; */
			ASSEMBLE (Fcol [Frpos [row]], S [i]) ;
		    }
		}
	    }
	    ep->ncolsleft-- ;
	    ASSERT (ep->ncolsleft > 0) ;
	}
	else
	{
	    *tp2++ = *tp ;	/* leave the tuple in the list */
	}
    }
    Col_tlen [col] = tp2 - tp1 ;

#ifndef NDEBUG
    DEBUG7 (("Column assembled in scan2-col: "ID"\n", col)) ;
    UMF_dump_rowcol (1, Numeric, Work, col, FALSE) ;
    DEBUG7 (("Current frontal matrix: after scan2-col\n")) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif

}


/* ========================================================================== */
/* === UMF_assemble / UMF_assemble_fixq ===================================== */
/* ========================================================================== */

#ifndef FIXQ
GLOBAL void UMF_assemble
#else
GLOBAL void UMF_assemble_fixq
#endif
(
    NumericType *Numeric,
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int e, i, row, col, i2, nrows, ncols, f, tpi, extcdeg, extrdeg, rdeg0,
	cdeg0, son_list, next, nrows_to_assemble,
	ncols_to_assemble, ngetrows, j, j2,
	nrowsleft,	/* number of rows remaining in S */
	ncolsleft,	/* number of columns remaining in S */
	prior_Lson, prior_Uson, *E, *Cols, *Rows, *Wm, *Woo,
	*Row_tuples, *Row_degree, *Row_tlen,
	*Col_tuples, *Col_tlen ;
    Unit *Memory, *p ;
    Element *ep ;
    Tuple *tp, *tp1, *tp2, *tpend ;
    Entry
	*S,		/* a pointer into the contribution block of a son */
	*Fcblock,	/* current contribution block */
	*Fcol ;		/* a column of FC */
    Int *Frpos,
	*Fcpos,
	fnrows,		/* number of rows in contribution block in F */
	fncols ;	/* number of columns in contribution block in F */

#if !defined (FIXQ) || !defined (NDEBUG)
    Int *Col_degree ;
#endif

#ifndef NDEBUG
    Int n_row, n_col ;
    n_row = Work->n_row ;
    n_col = Work->n_col ;
    DEBUG3 (("::Assemble SCANS 1-4\n")) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif

#if !defined (FIXQ) || !defined (NDEBUG)
    Col_degree = Numeric->Cperm ;   /* not updated if FIXQ is true */
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    fncols = Work->fncols ;
    fnrows = Work->fnrows ;
    Fcpos = Work->Fcpos ;
    Frpos = Work->Frpos ;
    Row_degree = Numeric->Rperm ;
    Row_tuples = Numeric->Uip ;
    Row_tlen   = Numeric->Uilen ;
    Col_tuples = Numeric->Lip ;
    Col_tlen   = Numeric->Lilen ;
    E = Work->E ;
    Memory = Numeric->Memory ;
    Wm = Work->Wm ;
    Woo = Work->Woo ;
    rdeg0 = Work->rdeg0 ;
    cdeg0 = Work->cdeg0 ;

#ifndef NDEBUG
    DEBUG6 (("============================================\n")) ;
    DEBUG6 (("Degree update, assembly.\n")) ;
    DEBUG6 (("pivot row pattern:  fncols="ID"\n", fncols)) ;
    for (j = 0 ; j < fncols ; j++)
    {
	col = Work->Fcols [j] ;
	DEBUG6 ((ID" ", col)) ;
	ASSERT (Fcpos [col] == j * Work->fnr_curr) ;
	ASSERT (NON_PIVOTAL_COL (col)) ;
    }
    ASSERT (Fcpos [Work->pivcol] >= 0) ;
    DEBUG6 (("pivcol: "ID" pos "ID" fnr_curr "ID" fncols "ID"\n",
	Work->pivcol, Fcpos [Work->pivcol], Work->fnr_curr, fncols)) ;
    ASSERT (Fcpos [Work->pivcol] <  fncols * Work->fnr_curr) ;
    DEBUG6 (("\npivot col pattern:  fnrows="ID"\n", fnrows)) ;
    for (i = 0 ; i < fnrows ; i++)
    {
	row = Work->Frows [i] ;
	DEBUG6 ((ID" ", row)) ;
	ASSERT (Frpos [row] == i) ;
	ASSERT (NON_PIVOTAL_ROW (row)) ;
    }
    DEBUG6 (("\n")) ;
    ASSERT (Frpos [Work->pivrow] >= 0) ;
    ASSERT (Frpos [Work->pivrow] < fnrows) ;
    ASSERT (Work->Flublock == (Entry *) (Numeric->Memory + E [0])) ;
    ASSERT (Work->Fcblock == Work->Flublock + Work->nb *
	(Work->nb + Work->fnr_curr + Work->fnc_curr)) ;
#endif

    Fcblock = Work->Fcblock ;

    /* ---------------------------------------------------------------------- */
    /* determine the largest actual frontal matrix size (for Info only) */
    /* ---------------------------------------------------------------------- */

    ASSERT (fnrows == Work->fnrows_new + 1) ;
    ASSERT (fncols == Work->fncols_new + 1) ;

    Numeric->maxnrows = MAX (Numeric->maxnrows, fnrows) ;
    Numeric->maxncols = MAX (Numeric->maxncols, fncols) ;

    /* this is safe from integer overflow, since the current frontal matrix
     * is already allocated. */
    Numeric->maxfrsize = MAX (Numeric->maxfrsize, fnrows * fncols) ;

    /* ---------------------------------------------------------------------- */
    /* assemble from prior elements into the current frontal matrix */
    /* ---------------------------------------------------------------------- */

    DEBUG2 (("New assemble start [prior_element:"ID"\n", Work->prior_element)) ;

    /* Currently no rows or columns are marked.  No elements are scanned, */
    /* that is, (ep->next == EMPTY) is true for all elements */

    son_list = 0 ;	/* start creating son_list [ */

    /* ---------------------------------------------------------------------- */
    /* determine if most recent element is Lson or Uson of current front */
    /* ---------------------------------------------------------------------- */

    if (!Work->do_extend)
    {
	prior_Uson = ( Work->pivcol_in_front && !Work->pivrow_in_front) ;
	prior_Lson = (!Work->pivcol_in_front &&  Work->pivrow_in_front) ;
	if (prior_Uson || prior_Lson)
	{
	    e = Work->prior_element ;
	    if (e != EMPTY)
	    {
		ASSERT (E [e]) ;
		p = Memory + E [e] ;
		ep = (Element *) p ;
		ep->next = son_list ;
		son_list = e ;
#ifndef NDEBUG
		DEBUG2 (("e "ID" is Prior son "ID" "ID"\n",
		    e, prior_Uson, prior_Lson)) ;
		UMF_dump_element (Numeric, Work, e, FALSE) ;
#endif
		ASSERT (E [e]) ;
	    }
	}
    }
    Work->prior_element = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* SCAN1-row:  scan the element lists of each new row in the pivot col */
    /* and compute the external column degree for each frontal */
    /* ---------------------------------------------------------------------- */

    for (i2 = Work->fscan_row ; i2 < fnrows ; i2++)
    {
	/* Get a row */
	row = Work->NewRows [i2] ;
	if (row < 0) row = FLIP (row) ;
	ASSERT (row >= 0 && row < n_row) ;

	DEBUG6 (("SCAN1-row: "ID"\n", row)) ;
#ifndef NDEBUG
	UMF_dump_rowcol (0, Numeric, Work, row, FALSE) ;
#endif

	ASSERT (NON_PIVOTAL_ROW (row)) ;
	tpi = Row_tuples [row] ;
	if (!tpi) continue ;
	tp = (Tuple *) (Memory + tpi) ;
	tp1 = tp ;
	tp2 = tp ;
	tpend = tp + Row_tlen [row] ;
	for ( ; tp < tpend ; tp++)
	{
	    e = tp->e ;
	    ASSERT (e > 0 && e <= Work->nel) ;
	    if (!E [e]) continue ;	/* element already deallocated */
	    f = tp->f ;
	    p = Memory + E [e] ;
	    ep = (Element *) p ;
	    p += UNITS (Element, 1) ;
	    Rows = ((Int *) p) + ep->ncols ;
	    if (Rows [f] == EMPTY) continue ;	/* row already assembled */
	    ASSERT (row == Rows [f]) ;

	    if (ep->cdeg < cdeg0)
	    {
		/* first time seen in scan1-row */
		ep->cdeg = ep->nrowsleft + cdeg0 ;
		DEBUG6 (("e "ID" First seen: cdeg: "ID" ", e, ep->cdeg-cdeg0)) ;
		ASSERT (ep->ncolsleft > 0 && ep->nrowsleft > 0) ;
	    }

	    ep->cdeg-- ;	/* decrement external column degree */
	    DEBUG6 (("e "ID" New ext col deg: "ID"\n", e, ep->cdeg - cdeg0)) ;

	    /* this element is not yet in the new son list */
	    if (ep->cdeg == cdeg0 && ep->next == EMPTY)
	    {
		/* A new LUson or Uson has been found */
		ep->next = son_list ;
		son_list = e ;
	    }

	    ASSERT (ep->cdeg >= cdeg0) ;
	    *tp2++ = *tp ;	/* leave the tuple in the list */
	}
	Row_tlen [row] = tp2 - tp1 ;
    }

    /* ---------------------------------------------------------------------- */
    /* SCAN1-col:  scan the element lists of each new col in the pivot row */
    /*	 and compute the external row degree for each frontal */
    /* ---------------------------------------------------------------------- */

    for (j2 = Work->fscan_col ; j2 < fncols ; j2++)
    {
	/* Get a column */
	col = Work->NewCols [j2] ;
	if (col < 0) col = FLIP (col) ;
	ASSERT (col >= 0 && col < n_col) ;

	DEBUG6 (("SCAN 1-col: "ID"\n", col)) ;
#ifndef NDEBUG
	UMF_dump_rowcol (1, Numeric, Work, col, FALSE) ;
#endif

	ASSERT (NON_PIVOTAL_COL (col)) ;
	tpi = Col_tuples [col] ;
	if (!tpi) continue ;
	tp = (Tuple *) (Memory + tpi) ;
	tp1 = tp ;
	tp2 = tp ;
	tpend = tp + Col_tlen [col] ;
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
	    if (Cols [f] == EMPTY) continue ;	/* column already assembled */
	    ASSERT (col == Cols [f]) ;

	    if (ep->rdeg < rdeg0)
	    {
		/* first time seen in scan1-col */
		ep->rdeg = ep->ncolsleft + rdeg0 ;
		DEBUG6 (("e "ID" First seen: rdeg: "ID" ", e, ep->rdeg-rdeg0)) ;
		ASSERT (ep->ncolsleft > 0 && ep->nrowsleft > 0) ;
	    }

	    ep->rdeg-- ;	/* decrement external row degree */
	    DEBUG6 (("e "ID" New ext row degree: "ID"\n", e, ep->rdeg-rdeg0)) ;

	    if (ep->rdeg == rdeg0 && ep->next == EMPTY)
	    {
		/* A new LUson or Lson has been found */
		ep->next = son_list ;
		son_list = e ;
	    }

	    ASSERT (ep->rdeg >= rdeg0) ;
	    *tp2++ = *tp ;	/* leave the tuple in the list */
	}
	Col_tlen [col] = tp2 - tp1 ;
    }

    /* ---------------------------------------------------------------------- */
    /* assemble new sons via full scans */
    /* ---------------------------------------------------------------------- */

    next = EMPTY ;

    for (e = son_list ; e > 0 ; e = next)
    {
	ASSERT (e > 0 && e <= Work->nel && E [e]) ;
	p = Memory + E [e] ;
	DEBUG2 (("New son: "ID"\n", e)) ;
#ifndef NDEBUG
	UMF_dump_element (Numeric, Work, e, FALSE) ;
#endif
	GET_ELEMENT (ep, p, Cols, Rows, ncols, nrows, S) ;
	nrowsleft = ep->nrowsleft ;
	ncolsleft = ep->ncolsleft ;
	next = ep->next ;
	ep->next = EMPTY ;

	extrdeg = (ep->rdeg < rdeg0) ? ncolsleft : (ep->rdeg - rdeg0) ;
	extcdeg = (ep->cdeg < cdeg0) ? nrowsleft : (ep->cdeg - cdeg0) ;
	ncols_to_assemble = ncolsleft - extrdeg ;
	nrows_to_assemble = nrowsleft - extcdeg ;
	DEBUG2 (("extrdeg "ID" extcdeg "ID"\n", extrdeg, extcdeg)) ;

	if (extrdeg == 0 && extcdeg == 0)
	{

	    /* -------------------------------------------------------------- */
	    /* this is an LUson - assemble an entire contribution block */
	    /* -------------------------------------------------------------- */

	    DEBUG6 (("LUson found: "ID"\n", e)) ;

	    if (nrows == nrowsleft)
	    {
		/* ---------------------------------------------------------- */
		/* no rows assembled out of this LUson yet */
		/* ---------------------------------------------------------- */

		/* compute the compressed column offset vector*/
		/* [ use Wm [0..nrows-1] for offsets */
#pragma ivdep
		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    Row_degree [row] -= ncolsleft ;
		    Wm [i] = Frpos [row] ;
		}

		if (ncols == ncolsleft)
		{
		    /* ------------------------------------------------------ */
		    /* no rows or cols assembled out of LUson yet */
		    /* ------------------------------------------------------ */

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
#ifndef FIXQ
			Col_degree [col] -= nrowsleft ;
#endif
			Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
			for (i = 0 ; i < nrows ; i++)
			{
			    /* Fcol [Wm [i]] += S [i] ; */
			    ASSEMBLE (Fcol [Wm [i]], S [i]) ;
			}
			S += nrows ;
		    }


		}
		else
		{
		    /* ------------------------------------------------------ */
		    /* only cols have been assembled out of LUson */
		    /* ------------------------------------------------------ */

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			if (col >= 0)
			{
#ifndef FIXQ
			    Col_degree [col] -= nrowsleft ;
#endif
			    Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
			    for (i = 0 ; i < nrows ; i++)
			    {
				/* Fcol [Wm [i]] += S [i] ; */
				ASSEMBLE (Fcol [Wm [i]], S [i]) ;
			    }
			}
			S += nrows ;
		    }

		}
		/* ] done using Wm [0..nrows-1] for offsets */
	    }
	    else
	    {
		/* ---------------------------------------------------------- */
		/* some rows have been assembled out of this LUson */
		/* ---------------------------------------------------------- */

		/* compute the compressed column offset vector*/
		/* [ use Woo,Wm [0..nrowsleft-1] for offsets */
		ngetrows = 0 ;
		for (i = 0 ; i < nrows ; i++)
		{
		    row = Rows [i] ;
		    if (row >= 0)
		    {
			Row_degree [row] -= ncolsleft ;
			Woo [ngetrows] = i ;
			Wm [ngetrows++] = Frpos [row] ;
		    }
		}
		ASSERT (ngetrows == nrowsleft) ;

		if (ncols == ncolsleft)
		{
		    /* ------------------------------------------------------ */
		    /* only rows have been assembled out of this LUson */
		    /* ------------------------------------------------------ */

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
#ifndef FIXQ
			Col_degree [col] -= nrowsleft ;
#endif
			Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
			for (i = 0 ; i < nrowsleft ; i++)
			{
			    /* Fcol [Wm [i]] += S [Woo [i]] ; */
			    ASSEMBLE (Fcol [Wm [i]], S [Woo [i]]) ;
			}
			S += nrows ;
		    }

		}
		else
		{
		    /* ------------------------------------------------------ */
		    /* both rows and columns have been assembled out of LUson */
		    /* ------------------------------------------------------ */

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			if (col >= 0)
			{
#ifndef FIXQ
			    Col_degree [col] -= nrowsleft ;
#endif
			    Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
			    for (i = 0 ; i < nrowsleft ; i++)
			    {
				/* Fcol [Wm [i]] += S [Woo [i]] ; */
				ASSEMBLE (Fcol [Wm [i]], S [Woo [i]]) ;
			    }
			}
			S += nrows ;
		    }

		}
		/* ] done using Woo,Wm [0..nrowsleft-1] */
	    }

	    /* deallocate the element: remove from ordered list */
	    UMF_mem_free_tail_block (Numeric, E [e]) ;
	    E [e] = 0 ;

	}
	else if (extcdeg == 0)
	{

	    /* -------------------------------------------------------------- */
	    /* this is a Uson - assemble all possible columns */
	    /* -------------------------------------------------------------- */

	    DEBUG6 (("New USON: "ID"\n", e)) ;
	    ASSERT (extrdeg > 0) ;

	    DEBUG6 (("New uson "ID" cols to do "ID"\n", e, ncols_to_assemble)) ;

	    if (ncols_to_assemble > 0)
	    {

		Int skip = FALSE ;
		if (ncols_to_assemble * 16 < ncols && nrows == 1)
		{
		    /* this is a tall and thin frontal matrix consisting of
		     * only one column (most likely an original column). Do
		     * not assemble it.   It cannot be the pivot column, since
		     * the pivot column element would be an LU son, not an Lson,
		     * of the current frontal matrix. */
		    ASSERT (nrowsleft == 1) ;
		    ASSERT (Rows [0] >= 0 && Rows [0] < Work->n_row) ;
		    skip = TRUE ;
		    Work->any_skip = TRUE ;
		}

		if (!skip)
		{

		    if (nrows == nrowsleft)
		    {
			/* -------------------------------------------------- */
			/* no rows have been assembled out of this Uson yet */
			/* -------------------------------------------------- */

			/* compute the compressed column offset vector */
			/* [ use Wm [0..nrows-1] for offsets */
#pragma ivdep
			for (i = 0 ; i < nrows ; i++)
			{
			    row = Rows [i] ;
			    ASSERT (row >= 0 && row < n_row) ;
			    Row_degree [row] -= ncols_to_assemble ;
			    Wm [i] = Frpos [row] ;
			}

			for (j = 0 ; j < ncols ; j++)
			{
			    col = Cols [j] ;
			    if ((col >= 0) && (Fcpos [col] >= 0))
			    {
#ifndef FIXQ
				Col_degree [col] -= nrowsleft ;
#endif
				Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
				for (i = 0 ; i < nrows ; i++)
				{
				    /* Fcol [Wm [i]] += S [i] ; */
				    ASSEMBLE (Fcol [Wm [i]], S [i]) ;
				}
				/* flag the column as assembled from Uson */
				Cols [j] = EMPTY ;
			    }
			    S += nrows ;
			}


			/* ] done using Wm [0..nrows-1] for offsets */
		    }
		    else
		    {
			/* -------------------------------------------------- */
			/* some rows have been assembled out of this Uson */
			/* -------------------------------------------------- */

			/* compute the compressed column offset vector*/
			/* [ use Woo,Wm [0..nrows-1] for offsets */
			ngetrows = 0 ;
			for (i = 0 ; i < nrows ; i++)
			{
			    row = Rows [i] ;
			    if (row >= 0)
			    {
				Row_degree [row] -= ncols_to_assemble ;
				ASSERT (row < n_row && Frpos [row] >= 0) ;
				Woo [ngetrows] = i ;
				Wm [ngetrows++] = Frpos [row] ;
			    }
			}
			ASSERT (ngetrows == nrowsleft) ;

			for (j = 0 ; j < ncols ; j++)
			{
			    col = Cols [j] ;
			    if ((col >= 0) && (Fcpos [col] >= 0))
			    {
#ifndef FIXQ
				Col_degree [col] -= nrowsleft ;
#endif
				Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
				for (i = 0 ; i < nrowsleft ; i++)
				{
				    /* Fcol [Wm [i]] += S [Woo [i]] ; */
				    ASSEMBLE (Fcol [Wm [i]], S [Woo [i]]) ;
				}
				/* flag the column as assembled from Uson */
				Cols [j] = EMPTY ;
			    }
			    S += nrows ;
			}

			/* ] done using Woo,Wm */
		    }
		    ep->ncolsleft = extrdeg ;
		}
	    }

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* this is an Lson - assemble all possible rows */
	    /* -------------------------------------------------------------- */

	    DEBUG6 (("New LSON: "ID"\n", e)) ;
	    ASSERT (extrdeg == 0 && extcdeg > 0) ;

	    DEBUG6 (("New lson "ID" rows to do "ID"\n", e, nrows_to_assemble)) ;

	    if (nrows_to_assemble > 0)
	    {

		Int skip = FALSE ;
		if (nrows_to_assemble * 16 < nrows && ncols == 1)
		{
		    /* this is a tall and thin frontal matrix consisting of
		     * only one column (most likely an original column). Do
		     * not assemble it.   It cannot be the pivot column, since
		     * the pivot column element would be an LU son, not an Lson,
		     * of the current frontal matrix. */
		    ASSERT (ncolsleft == 1) ;
		    ASSERT (Cols [0] >= 0 && Cols [0] < Work->n_col) ;
		    Work->any_skip = TRUE ;
		    skip = TRUE ;
		}

		if (!skip)
		{

		    /* compute the compressed column offset vector */
		    /* [ use Woo,Wm [0..nrows-1] for offsets */
		    ngetrows = 0 ;
		    for (i = 0 ; i < nrows ; i++)
		    {
			row = Rows [i] ;
			if ((row >= 0) && (Frpos [row] >= 0))
			{
			    ASSERT (row < n_row) ;
			    Row_degree [row] -= ncolsleft ;
			    Woo [ngetrows] = i ;
			    Wm [ngetrows++] = Frpos [row] ;
			    /* flag the row as assembled from the Lson */
			    Rows [i] = EMPTY ;
			}
		    }
		    ASSERT (nrowsleft - ngetrows == extcdeg) ;
		    ASSERT (ngetrows == nrows_to_assemble) ;

		    if (ncols == ncolsleft)
		    {
			/* -------------------------------------------------- */
			/* no columns assembled out this Lson yet */
			/* -------------------------------------------------- */

			for (j = 0 ; j < ncols ; j++)
			{
			    col = Cols [j] ;
			    ASSERT (col >= 0 && col < n_col) ;
#ifndef FIXQ
			    Col_degree [col] -= nrows_to_assemble ;
#endif
			    Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
			    for (i = 0 ; i < nrows_to_assemble ; i++)
			    {
				/* Fcol [Wm [i]] += S [Woo [i]] ; */
				ASSEMBLE (Fcol [Wm [i]], S [Woo [i]]) ;
			    }
			    S += nrows ;
			}


		    }
		    else
		    {
			/* -------------------------------------------------- */
			/* some columns have been assembled out of this Lson */
			/* -------------------------------------------------- */

			for (j = 0 ; j < ncols ; j++)
			{
			    col = Cols [j] ;
			    ASSERT (col < n_col) ;
			    if (col >= 0)
			    {
#ifndef FIXQ
				Col_degree [col] -= nrows_to_assemble ;
#endif
				Fcol = Fcblock + Fcpos [col] ;
#pragma ivdep
				for (i = 0 ; i < nrows_to_assemble ; i++)
				{
				    /* Fcol [Wm [i]] += S [Woo [i]] ; */
				    ASSEMBLE (Fcol [Wm [i]], S [Woo [i]]) ;
				}
			    }
			    S += nrows ;
			}

		    }

		    /* ] done using Woo,Wm */

		    ep->nrowsleft = extcdeg ;
		}
	    }
	}
    }

    /* Note that garbage collection, and build tuples */
    /* both destroy the son list. */

    /* ] son_list now empty */

    /* ---------------------------------------------------------------------- */
    /* If frontal matrix extended, assemble old L/Usons from new rows/cols */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* SCAN2-row:  assemble rows of old Lsons from the new rows */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    DEBUG7 (("Current frontal matrix: (prior to scan2-row)\n")) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif

    /* rescan the pivot row */
    if (Work->any_skip)
    {
	row_assemble (Work->pivrow, Numeric, Work) ;
    }

    if (Work->do_scan2row)
    {
	for (i2 = Work->fscan_row ; i2 < fnrows ; i2++)
	{
	    /* Get a row */
	    row = Work->NewRows [i2] ;
	    if (row < 0) row = FLIP (row) ;
	    ASSERT (row >= 0 && row < n_row) ;
	    if (!(row == Work->pivrow && Work->any_skip))
	    {
		/* assemble it */
		row_assemble (row, Numeric, Work) ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* SCAN2-col:  assemble columns of old Usons from the new columns */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    DEBUG7 (("Current frontal matrix: (prior to scan2-col)\n")) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif

    /* rescan the pivot col */
    if (Work->any_skip)
    {
	col_assemble (Work->pivcol, Numeric, Work) ;
    }

    if (Work->do_scan2col)
    {

	for (j2 = Work->fscan_col ; j2 < fncols ; j2++)
	{
	    /* Get a column */
	    col = Work->NewCols [j2] ;
	    if (col < 0) col = FLIP (col) ;
	    ASSERT (col >= 0 && col < n_col) ;
	    if (!(col == Work->pivcol && Work->any_skip))
	    {
		/* assemble it */
		col_assemble (col, Numeric, Work) ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* done.  the remainder of this routine is used only when in debug mode */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG

    /* ---------------------------------------------------------------------- */
    /* when debugging: make sure the assembly did everything that it could */
    /* ---------------------------------------------------------------------- */

    DEBUG3 (("::Assemble done\n")) ;

    for (i2 = 0 ; i2 < fnrows ; i2++)
    {
	/* Get a row */
	row = Work->Frows [i2] ;
	ASSERT (row >= 0 && row < n_row) ;

	DEBUG6 (("DEBUG SCAN 1: "ID"\n", row)) ;
	UMF_dump_rowcol (0, Numeric, Work, row, TRUE) ;

	ASSERT (NON_PIVOTAL_ROW (row)) ;
	tpi = Row_tuples [row] ;
	if (!tpi) continue ;
	tp = (Tuple *) (Memory + tpi) ;
	tpend = tp + Row_tlen [row] ;
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
	    Rows = ((Int *) p) + ep->ncols ;
	    if (Rows [f] == EMPTY) continue ;	/* row already assembled */
	    ASSERT (row == Rows [f]) ;
	    extrdeg = (ep->rdeg < rdeg0) ? ep->ncolsleft : (ep->rdeg - rdeg0) ;
	    extcdeg = (ep->cdeg < cdeg0) ? ep->nrowsleft : (ep->cdeg - cdeg0) ;
	    DEBUG6 ((
		"e "ID" After assembly ext row deg: "ID" ext col degree "ID"\n",
		e, extrdeg, extcdeg)) ;

	    if (Work->any_skip)
	    {
		/* no Lsons in any row, except for very tall and thin ones */
		ASSERT (extrdeg >= 0) ;
		if (extrdeg == 0)
		{
		    /* this is an unassemble Lson */
		    ASSERT (ep->ncols == 1) ;
		    ASSERT (ep->ncolsleft == 1) ;
		    col = Cols [0] ;
		    ASSERT (col != Work->pivcol) ;
		}
	    }
	    else
	    {
		/* no Lsons in any row */
		ASSERT (extrdeg > 0) ;
		/* Uson external row degree is = number of cols left */
		ASSERT (IMPLIES (extcdeg == 0, extrdeg == ep->ncolsleft)) ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */

    for (j2 = 0 ; j2 < fncols ; j2++)
    {
	/* Get a column */
	col = Work->Fcols [j2] ;
	ASSERT (col >= 0 && col < n_col) ;

	DEBUG6 (("DEBUG SCAN 2: "ID"\n", col)) ;
#ifndef FIXQ
	UMF_dump_rowcol (1, Numeric, Work, col, TRUE) ;
#else
	UMF_dump_rowcol (1, Numeric, Work, col, FALSE) ;
#endif

	ASSERT (NON_PIVOTAL_COL (col)) ;
	tpi = Col_tuples [col] ;
	if (!tpi) continue ;
	tp = (Tuple *) (Memory + tpi) ;
	tpend = tp + Col_tlen [col] ;
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
	    Rows = ((Int *) p) + ep->ncols ;
	    if (Cols [f] == EMPTY) continue ;	/* column already assembled */
	    ASSERT (col == Cols [f]) ;
	    extrdeg = (ep->rdeg < rdeg0) ? ep->ncolsleft : (ep->rdeg - rdeg0) ;
	    extcdeg = (ep->cdeg < cdeg0) ? ep->nrowsleft : (ep->cdeg - cdeg0) ;
	    DEBUG6 (("e "ID" After assembly ext col deg: "ID"\n", e, extcdeg)) ;

	    if (Work->any_skip)
	    {
		/* no Usons in any column, except for very tall and thin ones */
		ASSERT (extcdeg >= 0) ;
		if (extcdeg == 0)
		{
		    /* this is an unassemble Uson */
		    ASSERT (ep->nrows == 1) ;
		    ASSERT (ep->nrowsleft == 1) ;
		    row = Rows [0] ;
		    ASSERT (row != Work->pivrow) ;
		}
	    }
	    else
	    {
		/* no Usons in any column */
		ASSERT (extcdeg > 0) ;
		/* Lson external column degree is = number of rows left */
		ASSERT (IMPLIES (extrdeg == 0, extcdeg == ep->nrowsleft)) ;
	    }
	}
    }
#endif /* NDEBUG */
}
