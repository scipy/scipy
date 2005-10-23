/* ========================================================================== */
/* === UMF_dump ============================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* These routines, and external variables, are used only when debugging. */
/* If debugging is disabled (for normal operation) then this entire file */
/* becomes empty */

#include "umf_internal.h"

#ifndef NDEBUG

/* These global debugging variables and arrays do not exist if debugging */
/* is disabled at compile time (which is the default). */
GLOBAL Int UMF_debug = -999 ;
GLOBAL Int UMF_allocfail = FALSE ;
GLOBAL double UMF_gprob = -1.0 ;

/* static debugging arrays used only in UMF_dump_rowcol */
PRIVATE Int UMF_DBflag = 0 ;
PRIVATE Int UMF_DBpacked [UMF_DBMAX+1] ;
PRIVATE Int UMF_DBscatter [UMF_DBMAX+1] ;

/* ========================================================================== */
/* === UMF_DBinit =========================================================== */
/* ========================================================================== */

/* clear the debugging arrays */

PRIVATE void UMF_DBinit
(
    void
)
{
    Int i ;

    /* Int_MAX is defined in umfpack.h */
    if (UMF_DBflag < 1 || UMF_DBflag == Int_MAX)
    {
	/* clear the debugging arrays */
	UMF_DBflag = 0 ;
	for (i = 0 ; i <= UMF_DBMAX ; i++)
	{
	    UMF_DBscatter [i] = 0 ;
	    UMF_DBpacked  [i] = 0 ;
	}
    }

    UMF_DBflag++ ;

    /* UMF_DBflag > UMF_DBscatter [0...UMF_DBmax] is now true */
}

/* ========================================================================== */
/* === UMF_dump_dense ======================================================= */
/* ========================================================================== */

GLOBAL void UMF_dump_dense
(
    Entry *C,
    Int dim,
    Int m,
    Int n
)
{

    /* dump C [1..m,1..n], with column dimenstion dim */
    Int i, j;

    if (UMF_debug < 7) return ;
    if (C == (Entry *) NULL)
    {
	DEBUG7 (("No dense matrix allocated\n")) ;
	return ;
    }
    DEBUG8 ((" dimension= "ID" rows= "ID" cols= "ID"\n", dim, m, n)) ;

    for (i = 0 ; i < m ; i++)
    {
	DEBUG9 ((ID": ", i)) ;
	for (j = 0 ; j < n ; j++)
	{
	    EDEBUG9 (C [i+j*dim]) ;
	    if (j % 6 == 5) DEBUG9 (("\n     ")) ;
	}
	DEBUG9 (("\n")) ;
    }

    for (i = 0 ; i < m ; i++)
    {
	for (j = 0 ; j < n ; j++)
	{
	    if (IS_ZERO (C [i+j*dim]))
	    {
		DEBUG8 ((".")) ;
	    }
	    else
	    {
		DEBUG8 (("X")) ;
	    }
	}
	DEBUG8 (("\n")) ;
    }
}

/* ========================================================================== */
/* === UMF_dump_element ===================================================== */
/* ========================================================================== */

GLOBAL void UMF_dump_element
(
    NumericType *Numeric,
    WorkType *Work,
    Int e,
    Int clean
)
{

    Int i, j, k, *Rows, *Cols, nrows, ncols, *E, row, col,
	*Row_degree, *Col_degree ;
    Entry *C ;
    Element *ep ;
    Unit *p ;

    if (UMF_debug < 7) return ;

    if (e == 0)
    {
	UMF_dump_current_front (Numeric, Work, FALSE) ;
	return ;
    }

    DEBUG7 (("\n====================ELEMENT: "ID" ", e)) ;
    if (!Numeric || !Work || !Numeric->Memory)
    {
	DEBUG7 ((" No Numeric, Work\n")) ;
	return ;
    }
    DEBUG7 ((" nel: "ID" of "ID, e, Work->nel)) ;
    E = Work->E ;
    if (!E)
    {
	DEBUG7 ((" No elements\n")) ;
	return ;
    }
    if (e < 0 || e > Work->nel)
    {
	DEBUG7 (("e out of range!\n")) ;
	return ;
    }
    if (!E [e])
    {
	DEBUG7 ((" deallocated\n")) ;
	return ;
    }
    DEBUG7 (("\n")) ;
    Col_degree = Numeric->Cperm ;
    Row_degree = Numeric->Rperm ;

    p = Numeric->Memory + E [e] ;
    DEBUG7 (("ep "ID"\n", (Int) (p-Numeric->Memory))) ;
    GET_ELEMENT (ep, p, Cols, Rows, ncols, nrows, C) ;
    DEBUG7 (("nrows "ID" nrowsleft "ID"\n", nrows, ep->nrowsleft)) ;
    DEBUG7 (("ncols "ID" ncolsleft "ID"\n", ncols, ep->ncolsleft)) ;
    DEBUG7 (("cdeg-cdeg0 "ID" rdeg-rdeg0 "ID" next "ID"\n",
    ep->cdeg - Work->cdeg0, ep->rdeg - Work->rdeg0, ep->next)) ;

    DEBUG8 (("rows: ")) ;
    k = 0 ;
    for (i = 0 ; i < ep->nrows ; i++)
    {
	row = Rows [i] ;
	if (row >= 0)
	{
	    DEBUG8 ((" "ID, row)) ;
	    ASSERT (row < Work->n_row) ;
	    if ((k++ % 10) == 9) DEBUG8 (("\n")) ;
	    ASSERT (IMPLIES (clean, NON_PIVOTAL_ROW (row))) ;
	}
    }

    DEBUG8 (("\ncols: ")) ;
    k = 0 ;
    for (j = 0 ; j < ep->ncols ; j++)
    {
	col = Cols [j] ;
	if (col >= 0)
	{
	    DEBUG8 ((" "ID, col)) ;
	    ASSERT (col < Work->n_col) ;
	    if ((k++ % 10) == 9) DEBUG8 (("\n")) ;
	    ASSERT (IMPLIES (clean, NON_PIVOTAL_COL (col))) ;
	}
    }

    DEBUG8 (("\nvalues:\n")) ;
    if (UMF_debug >= 9)
    {
	for (i = 0 ; i < ep->nrows ; i++)
	{
	    row = Rows [i] ;
	    if (row >= 0)
	    {
		DEBUG9 ((ID": ", row)) ;
		k = 0 ;
		for (j = 0 ; j < ep->ncols ; j++)
		{
		    col = Cols [j] ;
		    if (col >= 0)
		    {
			EDEBUG9 (C [i+j*ep->nrows]) ;
			if (k++ % 6 == 5) DEBUG9 (("\n     ")) ;
		    }
		}
		DEBUG9 (("\n")) ;
	    }
	}
    }

    DEBUG7 (("====================\n")) ;
}


/* ========================================================================== */
/* === UMF_dump_rowcol ====================================================== */
/* ========================================================================== */

/* dump a row or a column, from one or more memory spaces */
/* return exact degree */

GLOBAL void UMF_dump_rowcol
(
    Int dumpwhich,		/* 0 for row, 1 for column */
    NumericType *Numeric,
    WorkType *Work,
    Int dumpindex,		/* row or column index to dump */
    Int check_degree	/* true if degree is to be checked */
)
{
    Int f, nrows, j, jj, len, e, deg, index, n_row, n_col, *Cols, *Rows, nn,
	dumpdeg, ncols, preve, *E, tpi, *Pattern, approx_deg, not_in_use ;
    Tuple *tp, *tend ;
    Element *ep ;
    Int *Row_tuples, *Row_degree, *Row_tlen ;
    Int *Col_tuples, *Col_degree, *Col_tlen ;
    Entry value, *C ;
    Unit *p ;
    Int is_there ;

    /* clear the debugging arrays */
    UMF_DBinit () ;

    if (dumpwhich == 0)
    {
	DEBUG7 (("\n====================ROW: "ID, dumpindex)) ;
    }
    else
    {
	DEBUG7 (("\n====================COL: "ID, dumpindex)) ;
    }

    if (dumpindex == EMPTY)
    {
	DEBUG7 ((" (EMPTY)\n")) ;
	return ;
    }

    deg = 0 ;
    approx_deg = 0 ;

    if (!Numeric || !Work)
    {
	DEBUG7 ((" No Numeric, Work\n")) ;
	return ;
    }

    n_row = Work->n_row ;
    n_col = Work->n_col ;
    nn = MAX (n_row, n_col) ;
    E = Work->E ;

    Col_degree = Numeric->Cperm ;
    Row_degree = Numeric->Rperm ;

    Row_tuples = Numeric->Uip ;
    Row_tlen   = Numeric->Uilen ;
    Col_tuples = Numeric->Lip ;
    Col_tlen   = Numeric->Lilen ;

	if (!E
	|| !Row_tuples || !Row_degree || !Row_tlen
	|| !Col_tuples || !Col_degree || !Col_tlen)
	{
	    DEBUG7 ((" No E, Rows, Cols\n")) ;
	    return ;
	}

	if (dumpwhich == 0)
	{
	    /* dump a row */
	    ASSERT (dumpindex >= 0 && dumpindex < n_row) ;
	    if (!NON_PIVOTAL_ROW (dumpindex))
	    {
		DEBUG7 ((" Pivotal\n")) ;
		return ;
	    }
	    len = Row_tlen [dumpindex] ;
	    dumpdeg = Row_degree [dumpindex] ;
	    tpi = Row_tuples [dumpindex] ;
	}
	else
	{
	    /* dump a column */
	    ASSERT (dumpindex >= 0 && dumpindex < n_col) ;
	    if (!NON_PIVOTAL_COL (dumpindex))
	    {
		DEBUG7 ((" Pivotal\n")) ;
		return ;
	    }
	    len = Col_tlen [dumpindex] ;
	    dumpdeg = Col_degree [dumpindex] ;
	    tpi = Col_tuples [dumpindex] ;
	}

	p = Numeric->Memory + tpi ;
	tp = (Tuple *) p ;
	if (!tpi)
	{
	    DEBUG7 ((" Nonpivotal, No tuple list tuples "ID" tlen "ID"\n",
		tpi, len)) ;
	    return ;
	}
	ASSERT (p >= Numeric->Memory + Numeric->itail) ;
	ASSERT (p <  Numeric->Memory + Numeric->size) ;

	DEBUG7 ((" degree: "ID" len: "ID"\n", dumpdeg, len)) ;
	not_in_use = (p-1)->header.size - UNITS (Tuple, len) ;
	DEBUG7 ((" Tuple list: p+1: "ID" size: "ID" units, "ID" not in use\n",
		(Int) (p-Numeric->Memory), (p-1)->header.size, not_in_use)) ;
	ASSERT (not_in_use >= 0) ;
	tend = tp + len ;
	preve = 0 ;
	for ( ; tp < tend ; tp++)
	{
	    /* row/col of element e, offset is f: */
	    /* DEBUG8 (("    (tp="ID")\n", tp)) ; */
	    e = tp->e ;
	    f = tp->f ;
	    DEBUG8 (("    (e="ID", f="ID")\n", e, f)) ;
	    ASSERT (e > 0 && e <= Work->nel) ;
	    /* dump the pattern and values */
	    if (E [e])
	    {
		p = Numeric->Memory + E [e] ;
		GET_ELEMENT (ep, p, Cols, Rows, ncols, nrows, C) ;
		if (dumpwhich == 0)
		{
		    Pattern = Cols ;
		    jj = ep->ncols ;
		    is_there = Rows [f] >= 0 ;
		    if (is_there) approx_deg += ep->ncolsleft ;
		}
		else
		{
		    Pattern = Rows ;
		    jj = ep->nrows ;
		    is_there = Cols [f] >= 0 ;
		    if (is_there) approx_deg += ep->nrowsleft ;
		}
		if (!is_there)
		{
			DEBUG8 (("\t\tnot present\n")) ;
		}
		else
		{
		    for (j = 0 ; j < jj ; j++)
		    {
			index = Pattern [j] ;
			value =
			    C [ (dumpwhich == 0) ? (f+nrows*j) : (j+nrows*f) ] ;
			if (index >= 0)
			{
			    DEBUG8 (("\t\t"ID":", index)) ;
			    EDEBUG8 (value) ;
			    DEBUG8 (("\n")) ;
			    if (dumpwhich == 0)
			    {
				/* col must be in the range 0..n_col-1 */
				ASSERT (index < n_col) ;
			    }
			    else
			    {
				/* row must be in the range 0..n_row-1 */
				ASSERT (index < n_row) ;
			    }

			    if (nn <= UMF_DBMAX)
			    {
				if (UMF_DBscatter [index] != UMF_DBflag)
				{
				    UMF_DBpacked [deg++] = index ;
				    UMF_DBscatter [index] = UMF_DBflag ;
				}
			    }
			}
		    }
		}
		/* the (e,f) tuples should be in order of their creation */
		/* this means that garbage collection will not jumble them */
		ASSERT (preve < e) ;
		preve = e ;
	    }
	    else
	    {
		DEBUG8 (("\t\tdeallocated\n")) ;
	    }
	}

    if (nn <= UMF_DBMAX)
    {
	if (deg > 0)
	{
	    DEBUG7 ((" Assembled, actual deg: "ID" : ", deg)) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		index = UMF_DBpacked [j] ;
		DEBUG8 ((ID" ", index)) ;
		if (j % 20 == 19) DEBUG8 (("\n ")) ;
		ASSERT (UMF_DBscatter [index] == UMF_DBflag) ;
	    }
	    DEBUG7 (("\n")) ;
	}
    }

    /* Col_degree is not maintained when fixQ is true */
    if (check_degree)
    {
	DEBUG8 (("  approx_deg "ID"  dumpdeg "ID"\n", approx_deg, dumpdeg)) ;
	ASSERT (approx_deg == dumpdeg) ;
    }

    DEBUG7 (("====================\n")) ;

    /* deg is now the exact degree */
    /* if nn <= UMF_DBMAX, then UMF_DBscatter [i] == UMF_DBflag for every i */
    /* in the row or col, and != UMF_DBflag if not */

    return ;
}


/* ========================================================================== */
/* === UMF_dump_matrix ====================================================== */
/* ========================================================================== */

GLOBAL void UMF_dump_matrix
(
    NumericType *Numeric,
    WorkType *Work,
    Int check_degree
)
{

    Int e, row, col, intfrag, frag, n_row, n_col, *E, fullsize, actualsize ;
    Element *ep ;
    Unit *p ;

    DEBUG6 (("=================================================== MATRIX:\n")) ;
    if (!Numeric || !Work)
    {
	DEBUG6 (("No Numeric or Work allocated\n")) ;
	return ;
    }
    if (!Numeric->Memory)
    {
	DEBUG6 (("No Numeric->Memory\n")) ;
	return ;
    }

	n_row = Work->n_row ;
	n_col = Work->n_col ;
	DEBUG6 (("n_row "ID" n_col "ID" nz "ID"\n", n_row, n_col, Work->nz)) ;
	DEBUG6 (("============================ ELEMENTS: "ID" \n", Work->nel)) ;
	intfrag = 0 ;
	E = Work->E ;
	if (!E)
	{
	    DEBUG6 (("No elements allocated\n")) ;
	}
	else
	{
	    for (e = 0 ; e <= Work->nel ; e++)
	    {
		UMF_dump_element (Numeric, Work, e, FALSE) ;
		if (e > 0 && E [e])
		{
		    p = Numeric->Memory + E [e] ;
		    ep = (Element *) p ;
		    ASSERT (ep->nrowsleft > 0 || ep->ncolsleft > 0) ;
		    fullsize = GET_BLOCK_SIZE (p) ;
		    actualsize = GET_ELEMENT_SIZE (ep->nrowsleft,ep->ncolsleft);
		    frag =  fullsize - actualsize ;
		    intfrag += frag ;
		    DEBUG7 (("dump el: "ID", full "ID" actual "ID" frag: "ID
			" intfrag: "ID"\n", e, fullsize, actualsize, frag,
			intfrag)) ;
		}
	    }
	}

	DEBUG6 (("CURRENT INTERNAL FRAG in elements: "ID" \n", intfrag)) ;



    DEBUG6 (("======================================== ROWS: "ID"\n", n_row)) ;
    UMF_debug -= 2 ;
    for (row = 0 ; row < n_row ; row++)
    {
	UMF_dump_rowcol (0, Numeric, Work, row, check_degree) ;
    }
    UMF_debug += 2 ;
    DEBUG6 (("======================================== COLS: "ID"\n", n_col)) ;
    UMF_debug -= 2 ;
    for (col = 0 ; col < n_col ; col++)
    {
	UMF_dump_rowcol (1, Numeric, Work, col, FALSE) ;
    }
    UMF_debug += 2 ;
    DEBUG6 (("============================================= END OF MATRIX:\n"));
}


/* ========================================================================== */
/* === UMF_dump_current_front =============================================== */
/* ========================================================================== */

GLOBAL void UMF_dump_current_front
(
    NumericType *Numeric,
    WorkType *Work,
    Int check
)
{

    Entry *Flublock, *Flblock, *Fublock, *Fcblock ;
    Int fnrows_max, fncols_max, fnrows, fncols, fnpiv, *Frows, *Fcols,
	i, j, *Fcpos, *Frpos, fnr_curr, fnc_curr, *E ;
    if (!Work) return ;
    DEBUG7 (("\n\n========CURRENT FRONTAL MATRIX:\n")) ;

    Flublock = Work->Flublock ;
    Flblock = Work->Flblock ;
    Fublock = Work->Fublock ;
    Fcblock = Work->Fcblock ;

    Frows = Work->Frows ;
    Fcols = Work->Fcols ;
    Frpos = Work->Frpos ;
    Fcpos = Work->Fcpos ;
    fnrows_max = Work->fnrows_max ;
    fncols_max = Work->fncols_max ;
    fnr_curr = Work->fnr_curr ;
    fnc_curr = Work->fnc_curr ;
    fnrows = Work->fnrows ;
    fncols = Work->fncols ;
    fnpiv = Work->fnpiv ;
    E = Work->E ;

    DEBUG6 (("=== fnpiv= "ID"\n", fnpiv)) ;
    DEBUG6 (("fnrows_max      fncols_max "ID" "ID"\n",fnrows_max, fncols_max)) ;
    DEBUG6 (("fnr_curr        fnc_curr   "ID" "ID"\n",fnr_curr,   fnc_curr)) ;
    DEBUG6 (("fnrows          fncols     "ID" "ID"\n",fnrows,     fncols)) ;
    ASSERT ((fnr_curr % 2 == 1) || fnr_curr == 0) ;
    DEBUG6 (("Pivot row pattern:\n")) ;
    for (j = 0 ; j < fncols ; j++)
    {
	DEBUG7 ((ID" "ID" "ID" %d\n", j, Fcols [j], Fcpos [Fcols [j]],
	    j < fncols)) ;
	if (check)
	{
	    ASSERT (Fcols [j] >= 0 && Fcols [j] < Work->n_col) ;
	    ASSERT (Fcpos [Fcols [j]] == j * fnr_curr) ;
	}
    }
    DEBUG6 (("Pivot col pattern:\n")) ;
    for (i = 0 ; i < fnrows ; i++)
    {
	DEBUG7 ((ID" "ID" "ID" %d\n", i, Frows [i], Frpos [Frows [i]],
	    i < fnrows)) ;
	if (check)
	{
	    ASSERT (Frows [i] >= 0 && Frows [i] < Work->n_row) ;
	    ASSERT (Frpos [Frows [i]] == i) ;
	}
    }
    if (UMF_debug < 7) return ;

    if (!E [0])
    {
	DEBUG6 (("current front not allocated\n")) ;
	ASSERT (!Work->Flublock) ;
	return ;
    }

    ASSERT (Work->Flublock == (Entry *) (Numeric->Memory + E [0])) ;
    DEBUG7 (("C  block: ")) ;
    UMF_dump_dense (Fcblock,  fnr_curr, fnrows, fncols) ;
    DEBUG7 (("L  block: ")) ;
    UMF_dump_dense (Flblock,  fnr_curr, fnrows, fnpiv) ;
    DEBUG7 (("U' block: ")) ;
    UMF_dump_dense (Fublock,  fnc_curr, fncols, fnpiv) ;
    DEBUG7 (("LU block: ")) ;
    UMF_dump_dense (Flublock, Work->nb, fnpiv, fnpiv) ;
    if (fnpiv > 0)
    {
	DEBUG7 (("Pivot entry: ")) ;
	EDEBUG7 (Flublock [(fnpiv-1)+(fnpiv-1)*Work->nb]) ;
	DEBUG7 (("\n")) ;
    }
}

/* ========================================================================== */
/* === UMF_dump_lu ========================================================== */
/* ========================================================================== */

GLOBAL void UMF_dump_lu
(
    NumericType *Numeric
)
{
    Int i, n_row, n_col, *Cperm, *Rperm ;

    DEBUG6 (("=============================================== LU factors:\n")) ;
    if (!Numeric)
    {
	DEBUG6 (("No LU factors allocated\n")) ;
	return ;
    }
    n_row = Numeric->n_row ;
    n_col = Numeric->n_col ;
    DEBUG6 (("n_row: "ID" n_col: "ID"\n", n_row, n_col)) ;
    DEBUG6 (("nLentries: "ID" nUentries: "ID"\n",
	Numeric->nLentries, Numeric->nUentries)) ;

    if (Numeric->Cperm)
    {
	Cperm = Numeric->Cperm ;
	DEBUG7 (("Column permutations: (new: old)\n")) ;
	for (i = 0 ; i < n_col ; i++)
	{
	    if (Cperm [i] != EMPTY)
	    {
		DEBUG7 ((ID": "ID"\n", i, Cperm [i])) ;
	    }
	}
    }
    else
    {
	DEBUG7 (("No Numeric->Cperm allocatated\n")) ;
    }

    if (Numeric->Rperm)
    {
	Rperm = Numeric->Rperm ;
	DEBUG7 (("row permutations: (new: old)\n")) ;
	for (i = 0 ; i < n_row ; i++)
	{
	    if (Rperm [i] != EMPTY)
	    {
		DEBUG7 ((ID": "ID"\n", i, Rperm [i])) ;
	    }
	}
    }
    else
    {
	DEBUG7 (("No Numeric->Rperm allocatated\n")) ;
    }

    DEBUG6 (("========================================= END OF LU factors:\n"));
}


/* ========================================================================== */
/* === UMF_dump_memory ====================================================== */
/* ========================================================================== */

GLOBAL void UMF_dump_memory
(
    NumericType *Numeric
)
{

    Unit *p ;
    Int prevsize, s ;
    Int found ;

    if (!Numeric)
    {
	DEBUG6 (("No memory space S allocated\n")) ;
	return ;
    }

    DEBUG6 (("\n ============================================== MEMORY:\n")) ;
    if (!Numeric || !Numeric->Memory)
    {
	DEBUG6 (("No memory space Numeric allocated\n")) ;
	return ;
    }

    DEBUG6 (("S: "ID"\n", (Int) Numeric)) ;
    DEBUG6 (("S->ihead           : "ID"\n", Numeric->ihead)) ;
    DEBUG6 (("S->itail           : "ID"\n", Numeric->itail)) ;
    DEBUG6 (("S->size            : "ID"\n", Numeric->size)) ;
    DEBUG6 (("S->ngarbage        : "ID"\n", Numeric->ngarbage)) ;
    DEBUG6 (("S->nrealloc        : "ID"\n", Numeric->nrealloc)) ;
    DEBUG6 (("   in use at head           : "ID"\n", Numeric->ihead)) ;
    DEBUG6 (("   free space               : "ID"\n",
	Numeric->itail - Numeric->ihead)) ;
    DEBUG6 (("   blocks in use at tail    : "ID"\n",
	Numeric->size - Numeric->itail)) ;
    DEBUG6 (("   total in use             : "ID"\n",
	Numeric->size - (Numeric->itail - Numeric->ihead))) ;

    prevsize = 0 ;
    found = FALSE ;

    ASSERT (0 <= Numeric->ihead) ;
    ASSERT (Numeric->ihead <= Numeric->itail) ;
    ASSERT (Numeric->itail <= Numeric->size) ;

    p = Numeric->Memory + Numeric->itail ;

    while (p < Numeric->Memory + Numeric->size)
    {
	DEBUG8 (("p: "ID" p+1: "ID" prevsize: "ID" size: "ID,
	    (Int) (p-Numeric->Memory), (Int) (p+1-Numeric->Memory),
	    p->header.prevsize, p->header.size)) ;
	if (p->header.size < 0)
	{
	    DEBUG8 ((" free")) ;
	}

	if (p == Numeric->Memory + Numeric->itail)
	{
	    ASSERT (p->header.prevsize == 0) ;
	}
	else
	{
	    ASSERT (p->header.prevsize > 0) ;
	}

	ASSERT (p->header.size != 0) ;
	s = prevsize >= 0 ? prevsize : -prevsize ;
	ASSERT (p->header.prevsize == s) ;
	/* no adjacent free blocks */
	ASSERT (p->header.size > 0 || prevsize > 0) ;
	if (Numeric->ibig != EMPTY)
	{
	    if (p == Numeric->Memory + Numeric->ibig)
	    {
		ASSERT (p->header.size < 0) ;
		DEBUG8 ((" <===== Numeric->ibig")) ;
		found = TRUE ;
	    }
	}
	s = p->header.size ;
	prevsize = s ;
	s = s >= 0 ? s : -s ;
	p = p + 1 + s ;
	DEBUG8 (("\n")) ;
    }

    ASSERT (p == Numeric->Memory + Numeric->size) ;
    ASSERT (IMPLIES (Numeric->ibig != EMPTY, found)) ;
    DEBUG6 (("============================================= END OF MEMORY:\n"));

}


/* ========================================================================== */
/* === UMF_dump_packed_memory =============================================== */
/* ========================================================================== */

GLOBAL void UMF_dump_packed_memory
(
    NumericType *Numeric,
    WorkType *Work
)
{
    Unit *p, *p3 ;
    Int prevsize, col, row, *Rows, *Cols, ncols, nrows, k, esize,
	*Row_tuples, *Row_degree, *Col_tuples, *Col_degree ;
    Entry *C ;
    Element *ep ;

    Col_degree = Numeric->Cperm ;	/* for NON_PIVOTAL_COL macro */
    Row_degree = Numeric->Rperm ;	/* for NON_PIVOTAL_ROW macro */
    Row_tuples = Numeric->Uip ;
    Col_tuples = Numeric->Lip ;

    DEBUG6 (("============================================ PACKED MEMORY:\n")) ;
    if (!Numeric || !Numeric->Memory)
    {
	DEBUG6 (("No memory space S allocated\n")) ;
	return ;
    }
    DEBUG6 (("S: "ID"\n", (Int) Numeric)) ;
    DEBUG6 (("S->ihead           : "ID"\n", Numeric->ihead)) ;
    DEBUG6 (("S->itail           : "ID"\n", Numeric->itail)) ;
    DEBUG6 (("S->size            : "ID"\n", Numeric->size)) ;
    DEBUG6 (("S->ngarbage        : "ID"\n", Numeric->ngarbage)) ;
    DEBUG6 (("S->nrealloc        : "ID"\n", Numeric->nrealloc)) ;
    DEBUG6 (("   in use at head           : "ID"\n", Numeric->ihead)) ;
    DEBUG6 (("   free space               : "ID"\n",
	Numeric->itail - Numeric->ihead)) ;
    DEBUG6 (("   blocks in use at tail    : "ID"\n",
	Numeric->size - Numeric->itail)) ;
    DEBUG6 (("   total in use             : "ID"\n",
	Numeric->size - (Numeric->itail - Numeric->ihead))) ;

    ASSERT (0 <= Numeric->ihead) ;
    ASSERT (Numeric->ihead <= Numeric->itail) ;
    ASSERT (Numeric->itail <= Numeric->size) ;

    for (row = 0 ; row < Work->n_row ; row++)
    {
	ASSERT (IMPLIES (NON_PIVOTAL_ROW (row), !Row_tuples [row])) ;
    }
    for (col = 0 ; col < Work->n_col ; col++)
    {
	ASSERT (IMPLIES (NON_PIVOTAL_COL (col), !Col_tuples [col])) ;
    }

    prevsize = 0 ;
    p = Numeric->Memory + Numeric->itail ;
    while (p < Numeric->Memory + Numeric->size)
    {
	DEBUG9 (("====================\n")) ;
	DEBUG7 (("p: "ID" p+1: "ID" prevsize: "ID" size: "ID"\n",
	    (Int) (p-Numeric->Memory), (Int) (p+1-Numeric->Memory),
	    p->header.prevsize, p->header.size)) ;
	ASSERT (p->header.size > 0) ;

	if (p == Numeric->Memory + Numeric->itail)
	{
	    ASSERT (p->header.prevsize == 0) ;
	}
	else
	{
	    ASSERT (p->header.prevsize > 0) ;
	}

	ASSERT (p->header.prevsize == prevsize) ;
	prevsize = p->header.size ;

	if (p != Numeric->Memory + Numeric->size - 2)
	{

	    p3 = p + 1 ;
	    if (p3 == Numeric->Memory + Work->E [0])
	    {
		/* this is the current frontal matrix */
		UMF_dump_current_front (Numeric, Work, FALSE) ;
	    }
	    else
	    {

		/* this is a packed element */
		GET_ELEMENT (ep, p3, Cols, Rows, ncols, nrows, C) ;
		DEBUG9 (("ep "ID"\n nrows "ID" ncols "ID"\n",
		    (Int) ((p+1)-Numeric->Memory), ep->nrows, ep->ncols)) ;
		DEBUG9 (("rows:")) ;
		for (k = 0 ; k < ep->nrows; k++)
		{
		    row = Rows [k] ;
		    DEBUG9 ((" "ID, row)) ;
		    ASSERT (row >= 0 && row <= Work->n_row) ;
		    if ((k % 10) == 9) DEBUG9 (("\n")) ;
		}
		DEBUG9 (("\ncols:")) ;
		for (k = 0 ; k < ep->ncols; k++)
		{
		    col = Cols [k] ;
		    DEBUG9 ((" "ID, col)) ;
		    ASSERT (col >= 0 && col <= Work->n_col) ;
		    if ((k % 10) == 9) DEBUG9 (("\n")) ;
		}
		DEBUG9 (("\nvalues: ")) ;
		if (UMF_debug >= 9)
		{
		    UMF_dump_dense (C, ep->nrows, ep->nrows, ep->ncols) ;
		}
		esize = GET_ELEMENT_SIZE (ep->nrows, ep->ncols) ;
		DEBUG9 (("esize: "ID"\n", esize)) ;
		ASSERT (esize <= p->header.size) ;
	    }

	}
	else
	{
	    /* this is the final marker block */
	    ASSERT (p->header.size == 1) ;
	}
	p = p + 1 + p->header.size ;
    }

    ASSERT (Numeric->ibig == EMPTY) ;
    ASSERT (p == Numeric->Memory + Numeric->size) ;
    DEBUG6 (("======================================END OF PACKED MEMORY:\n")) ;

}

/* ========================================================================== */
/* === UMF_dump_col_matrix ================================================== */
/* ========================================================================== */

/* This code is the same for real or complex matrices. */

GLOBAL void UMF_dump_col_matrix
(
    const double Ax [ ],	/* Ax [0..nz-1]: real values, in column order */
#ifdef COMPLEX
    const double Az [ ],	/* Az [0..nz-1]: imag values, in column order */
#endif
    const Int Ai [ ],		/* Ai [0..nz-1]: row indices, in column order */
    const Int Ap [ ],		/* Ap [0..n_col]: column pointers */
    Int n_row,			/* number of rows of A */
    Int n_col,			/* number of columns of A */
    Int nz			/* number of entries */
)
{
    Int col, p, p1, p2, row ;
    if (!Ai || !Ap) return ;
    DEBUG6 (("============================================ COLUMN FORM:\n")) ;


    ASSERT (n_col >= 0) ;
    nz = Ap [n_col] ;
    DEBUG2 (("UMF_dump_col:  nz "ID"\n", nz)) ;
    DEBUG2 (("n_row "ID"  \n", n_row)) ;
    DEBUG2 (("n_col "ID"  \n", n_col)) ;

    DEBUG6 ((" n_row = "ID", n_col ="ID" nz = "ID" Ap [0] "ID", Ap [n] "ID"\n",
	n_row, n_col, nz, Ap [0], Ap [n_col])) ;
    ASSERT (Ap [0] == 0) ;
    ASSERT (Ap [n_col] == nz) ;
    for (col = 0 ; col < n_col ; col++)
    {
	p1 = Ap [col] ;
	p2 = Ap [col+1] ;
	DEBUG6 (("col: "ID", length "ID"\n", col, p2 - p1)) ;
	ASSERT (p2 >= p1) ;
	for (p = p1 ; p < p2 ; p++)
	{
	    row = Ai [p] ;
	    ASSERT (row >= 0 && row < n_row) ;
	    DEBUG6 (("\t"ID" ", row)) ;
	    if (Ax != (double *) NULL)
	    {
#ifdef COMPLEX
		if (Az != (double *) NULL)
		{
		    DEBUG6 ((" (%e+%ei) ", Ax [p], Az [p])) ;
		}
		else
		{
		    DEBUG6 ((" %e", Ax [p])) ;
		}
#else
		DEBUG6 ((" %e", Ax [p])) ;
#endif
	    }
	    DEBUG6 (("\n")) ;
	}
    }
    DEBUG6 (("========================================== COLUMN FORM done\n")) ;
}


/* ========================================================================== */
/* === UMF_dump_chain ======================================================= */
/* ========================================================================== */

GLOBAL void UMF_dump_chain
(
    Int frontid,
    Int Front_parent [ ],
    Int Front_npivcol [ ],
    Int Front_nrows [ ],
    Int Front_ncols [ ],
    Int nfr
)
{
    Int i, len = 0 ;

    /* print a list of contiguous parents */
    i = frontid ;
    ASSERT (Front_parent [i] == EMPTY ||
	(Front_parent [i] > i && Front_parent [i] < nfr)) ;

    len++ ;
    DEBUG3 (("Chain:\n	"ID" ["ID","ID"]("ID"-by-"ID")\n", i,
		Front_npivcol [i],
		MIN (Front_npivcol [i], Front_nrows [i]),
		Front_nrows [i],
		Front_ncols [i])) ;

    for (i = frontid ; i < nfr ; i++)
    {
	ASSERT (Front_parent [i] == EMPTY ||
	(Front_parent [i] > i && Front_parent [i] < nfr)) ;
	if (Front_parent [i] == i+1)
	{
	    len++ ;
	    DEBUG3 (("\t"ID" ["ID","ID"]("ID"-by-"ID")\n", i+1,
		Front_npivcol [i+1],
		MIN (Front_npivcol [i+1], Front_nrows [i+1]),
		Front_nrows [i+1],
		Front_ncols [i+1])) ;
	}
	else
	{
	    DEBUG2 (("Length of chain: "ID"\n", len)) ;
	    return ;
	}
    }
}


/* ========================================================================== */
/* === UMF_dump_start ======================================================= */
/* ========================================================================== */

GLOBAL void UMF_dump_start
(
    void
)
{
    FILE *ff ;

    /* get the debug print level from the "debug.umf" file, if it exists */
    UMF_debug = -999 ;
    ff = fopen ("debug.umf", "r") ;
    if (ff)
    {
	(void) fscanf (ff, ID, &UMF_debug) ;
	(void) fclose (ff) ;
    }

    DEBUG0 (("umfpack: debug version (SLOW!) ")) ;

    DEBUG0 ((" BLAS: ")) ;

#if defined (USE_NO_BLAS)
    DEBUG0 (("none.")) ;
#elif defined (USE_C_BLAS)
    DEBUG0 (("C-BLAS.")) ;
#elif defined (USE_MATLAB_BLAS)
    DEBUG0 (("built-in MATLAB BLAS.")) ;
#elif defined (USE_SUNPERF_BLAS)
    DEBUG0 (("Sun Performance Library BLAS.")) ;
#elif defined (USE_SCSL_BLAS)
    DEBUG0 (("SGI SCSL BLAS.")) ;
#elif defined (USE_FORTRAN_BLAS)
    DEBUG0 (("Fortran BLAS.")) ;
#endif

    DEBUG0 ((" MATLAB: ")) ;
#ifdef MATLAB_MEX_FILE
    DEBUG0 (("mexFunction.\n")) ;
#else
#ifdef MATHWORKS
    DEBUG0 (("yes (uses MathWorks internal ut* routines).\n")) ;
#else
    DEBUG0 (("no.\n")) ;
#endif
#endif

    UMF_gprob = -1.0 ;
    ff = fopen ("gprob.umf", "r") ;
    if (ff)
    {
	(void) fscanf (ff, "%lg", &UMF_gprob) ;
	(void) fclose (ff) ;
	srand (1) ;	/* restart the random number generator */
    }

    if (UMF_gprob > 1.0) UMF_gprob = 1.0 ;
    DEBUG1 (("factor: UMF_gprob: %e UMF_debug "ID"\n", UMF_gprob, UMF_debug)) ;

    DEBUG2 (("sizeof: (bytes / int / Units) \n")) ;
    DEBUG2 (("sizeof (Int)           %u %u %u\n",
    sizeof (Int), sizeof (Int) / sizeof (int), UNITS (Int, 1) )) ;
    DEBUG2 (("sizeof (int)           %u %u %u\n",
    sizeof (int), sizeof (int) / sizeof (int), UNITS (int, 1) )) ;
    DEBUG2 (("sizeof (size_t)        %u %u %u\n",
    sizeof (size_t), sizeof (size_t) / sizeof (size_t), UNITS (size_t, 1) )) ;
    DEBUG2 (("sizeof (long)          %u %u %u\n",
    sizeof (long), sizeof (long) / sizeof (long), UNITS (long, 1) )) ;
    DEBUG2 (("sizeof (double)        %u %u %u\n",
    sizeof (double), sizeof (double) / sizeof (int), UNITS (double, 1) )) ;
    DEBUG2 (("sizeof (Unit)          %u %u %u\n",
    sizeof (Unit), sizeof (Unit) / sizeof (int), UNITS (Unit, 1) )) ;
    DEBUG2 (("sizeof (Entry)         %u %u %u\n",
    sizeof (Entry), sizeof (Entry) / sizeof (int), UNITS (Entry, 1) )) ;
    DEBUG2 (("sizeof (Tuple)         %u %u %u\n",
    sizeof (Tuple), sizeof (Tuple) / sizeof (int), UNITS (Tuple, 1) )) ;
    DEBUG2 (("sizeof (Tuple *)       %u %u %u\n",
    sizeof (Tuple *), sizeof (Tuple *) / sizeof (int), UNITS (Tuple *, 1) )) ;
    DEBUG2 (("sizeof (Element)       %u %u %u\n",
    sizeof (Element), sizeof (Element) / sizeof (int), UNITS (Element, 1) )) ;
    DEBUG2 (("sizeof (Element *)     %u %u %u\n",
    sizeof (Element *), sizeof (Element *) / sizeof (int),
    UNITS (Element *, 1) )) ;
    DEBUG2 (("sizeof (WorkType)      %u %u %u\n",
    sizeof (WorkType), sizeof (WorkType) / sizeof (int),
    UNITS (WorkType, 1) )) ;
    DEBUG2 (("sizeof (NumericType)   %u %u %u\n",
    sizeof (NumericType), sizeof (NumericType) / sizeof (int),
    UNITS (NumericType, 1) )) ;
    DEBUG2 (("sizeof (SymbolicType)  %u %u %u\n",
    sizeof (SymbolicType), sizeof (SymbolicType) / sizeof (int),
    UNITS (SymbolicType, 1) )) ;

}


/* ========================================================================== */
/* === UMF_dump_rowmerge ==================================================== */
/* ========================================================================== */

GLOBAL void UMF_dump_rowmerge
(
    NumericType *Numeric,
    SymbolicType *Symbolic,
    WorkType *Work
)
{
    Int *Front_leftmostdesc, *Front_1strow, *Front_new1strow, row1, row2,
	fleftmost, nfr, n_row, *Row_degree, i, frontid, row ;

    nfr = Symbolic->nfr ;
    DEBUG3 (("\n================== Row merge sets: nfr "ID"\n", nfr)) ;
    Front_leftmostdesc = Symbolic->Front_leftmostdesc ;
    Front_1strow = Symbolic->Front_1strow ;
    Front_new1strow = Work->Front_new1strow ;
    n_row = Symbolic->n_row ;
    Row_degree = Numeric->Rperm ;
    frontid = Work->frontid ;

    for (i = frontid ; i <= nfr ; i++)
    {
	DEBUG3 (("----------------------\n")) ;
	if (i == nfr) DEBUG3 (("Dummy: ")) ;
	DEBUG3 (("Front "ID" 1strow "ID" new1strow "ID" leftmostdesc "ID,
	    i, Front_1strow [i], Front_new1strow [i], Front_leftmostdesc [i])) ;
	DEBUG3 ((" parent "ID" pivcol "ID"\n", Symbolic->Front_parent [i],
	    Symbolic->Front_npivcol [i])) ;

	if (i == nfr)
	{
	    fleftmost = -1 ;
	    row1 = Front_new1strow [i] ;
	    row2 = n_row-1 ;
	}
	else
	{
	    fleftmost = Front_leftmostdesc [i] ;
	    row1 = Front_new1strow [fleftmost] ;
	    row2 = Front_1strow [i+1] - 1 ;
	}
	DEBUG3 (("Leftmost: "ID"  Rows ["ID" to "ID"], search ["ID" to "ID"]\n",
	    fleftmost, Front_1strow [i], row2, row1, row2)) ;

	for (row = row1 ; row <= row2 ; row++)
	{
	    ASSERT (row >= 0 && row < n_row) ;
	    DEBUG3 (("   Row "ID" live: %d\n", row, NON_PIVOTAL_ROW (row))) ;
	}
    }
}

/* ========================================================================== */
/* === UMF_dump_diagonal_map ================================================ */
/* ========================================================================== */

GLOBAL void UMF_dump_diagonal_map
(
    Int Diagonal_map [ ],
    Int Diagonal_imap [ ],
    Int n1,
    Int nn,
    Int nempty
)
{
    Int row, col ;
    if (Diagonal_map != (Int *) NULL)
    {
	DEBUG2 (("\nDump the Diagonal_map: n1 "ID" nn "ID" nempty "ID"\n",
	    n1, nn, nempty)) ;
	for (col = n1 ; col < nn - nempty ; col++)
	{
	    row = Diagonal_map [col] ;
	    DEBUG2 (("     Diagonal_map [col = "ID"] gives "ID": ",
		col, row)) ;
	    row = UNFLIP (row) ;
	    DEBUG2 ((" row "ID"\n", row)) ;
	    ASSERT (Diagonal_imap [row] == col) ;
	}
    }
}

#endif /* NDEBUG */
