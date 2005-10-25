/* ========================================================================== */
/* === UMF_extend_front ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Called by kernel. */

#include "umf_internal.h"
#include "umf_grow_front.h"

/* ========================================================================== */
/* === zero_front =========================================================== */
/* ========================================================================== */

PRIVATE void zero_front (
    Entry *Flblock, Entry *Fublock, Entry *Fcblock,
    Int fnrows, Int fncols, Int fnr_curr, Int fnc_curr,
    Int fnpiv, Int fnrows_extended, Int fncols_extended)
{
    Int j, i ;
    Entry *F, *Fj, *Fi ;

    Fj = Fcblock + fnrows ;
    for (j = 0 ; j < fncols ; j++)
    {
	/* zero the new rows in the contribution block: */
	F = Fj ;
	Fj += fnr_curr ;
#pragma ivdep
	for (i = fnrows ; i < fnrows_extended ; i++)
	{
	    /* CLEAR (Fcblock [i + j*fnr_curr]) ; */
	    CLEAR_AND_INCREMENT (F) ;
	}
    }

    Fj -= fnrows ;
    for (j = fncols ; j < fncols_extended ; j++)
    {
	/* zero the new columns in the contribution block: */
	F = Fj ;
	Fj += fnr_curr ;
#pragma ivdep
	for (i = 0 ; i < fnrows_extended ; i++)
	{
	    /* CLEAR (Fcblock [i + j*fnr_curr]) ; */
	    CLEAR_AND_INCREMENT (F) ;
	}
    }

    Fj = Flblock + fnrows ;
    for (j = 0 ; j < fnpiv ; j++)
    {
	/* zero the new rows in L block: */
	F = Fj ;
	Fj += fnr_curr ;
#pragma ivdep
	for (i = fnrows ; i < fnrows_extended ; i++)
	{
	    /* CLEAR (Flblock [i + j*fnr_curr]) ; */
	    CLEAR_AND_INCREMENT (F) ;
	}
    }

    Fi = Fublock + fncols ;
    for (i = 0 ; i < fnpiv ; i++)
    {
	/* zero the new columns in U block: */
	F = Fi ;
	Fi += fnc_curr ;
#pragma ivdep
	for (j = fncols ; j < fncols_extended ; j++)
	{
	    /* CLEAR (Fublock [i*fnc_curr + j]) ; */
	    CLEAR_AND_INCREMENT (F) ;
	}
    }

}

/* ========================================================================== */
/* === UMF_extend_front ===================================================== */
/* ========================================================================== */

GLOBAL Int UMF_extend_front
(
    NumericType *Numeric,
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int j, i, *Frows, row, col, *Wrow, fnr2, fnc2, *Frpos, *Fcpos, *Fcols,
	fnrows_extended, rrdeg, ccdeg, fncols_extended, fnr_curr, fnc_curr,
	fnrows, fncols, pos, fnpiv, *Wm ;
    Entry *Wx, *Wy, *Fu, *Fl ;

    /* ---------------------------------------------------------------------- */
    /* get current frontal matrix and check for frontal growth */
    /* ---------------------------------------------------------------------- */

    fnpiv = Work->fnpiv ;

#ifndef NDEBUG
    DEBUG2 (("EXTEND FRONT\n")) ;
    DEBUG2 (("Work->fnpiv "ID"\n", fnpiv)) ;
    ASSERT (Work->Flblock  == Work->Flublock + Work->nb*Work->nb) ;
    ASSERT (Work->Fublock  == Work->Flblock  + Work->fnr_curr*Work->nb) ;
    ASSERT (Work->Fcblock  == Work->Fublock  + Work->nb*Work->fnc_curr) ;
    DEBUG7 (("C  block: ")) ;
    UMF_dump_dense (Work->Fcblock,  Work->fnr_curr, Work->fnrows, Work->fncols) ;
    DEBUG7 (("L  block: ")) ;
    UMF_dump_dense (Work->Flblock,  Work->fnr_curr, Work->fnrows, fnpiv);
    DEBUG7 (("U' block: ")) ;
    UMF_dump_dense (Work->Fublock,  Work->fnc_curr, Work->fncols, fnpiv) ;
    DEBUG7 (("LU block: ")) ;
    UMF_dump_dense (Work->Flublock, Work->nb, fnpiv, fnpiv) ;
#endif

    if (Work->do_grow)
    {
	fnr2 = UMF_FRONTAL_GROWTH * Work->fnrows_new + 2 ;
	fnc2 = UMF_FRONTAL_GROWTH * Work->fncols_new + 2 ;
	if (!UMF_grow_front (Numeric, fnr2, fnc2, Work, 1))
	{
	    DEBUGm4 (("out of memory: extend front\n")) ;
	    return (FALSE) ;
	}
    }

    fnr_curr = Work->fnr_curr ;
    fnc_curr = Work->fnc_curr ;
    ASSERT (Work->fnrows_new + 1 <= fnr_curr) ;
    ASSERT (Work->fncols_new + 1 <= fnc_curr) ;
    ASSERT (fnr_curr >= 0 && fnr_curr % 2 == 1) ;

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    Frows = Work->Frows ;
    Frpos = Work->Frpos ;
    Fcols = Work->Fcols ;
    Fcpos = Work->Fcpos ;
    fnrows = Work->fnrows ;
    fncols = Work->fncols ;
    rrdeg = Work->rrdeg ;
    ccdeg = Work->ccdeg ;

    /* scan starts at the first new column in Fcols */
    /* also scan the pivot column if it was not in the front */
    Work->fscan_col = fncols ;
    Work->NewCols = Fcols ;

    /* scan1 starts at the first new row in Frows */
    /* also scan the pivot row if it was not in the front */
    Work->fscan_row = fnrows ;
    Work->NewRows = Frows ;

    /* ---------------------------------------------------------------------- */
    /* extend row pattern of the front with the new pivot column */
    /* ---------------------------------------------------------------------- */

    fnrows_extended = fnrows ;
    fncols_extended = fncols ;

#ifndef NDEBUG
    DEBUG2 (("Pivot col, before extension: "ID"\n", fnrows)) ;
    for (i = 0 ; i < fnrows ; i++)
    {
	DEBUG2 ((" "ID": row "ID"\n", i, Frows [i])) ;
	ASSERT (Frpos [Frows [i]] == i) ;
    }
    DEBUG2 (("Extending pivot column: pivcol_in_front: "ID"\n",
	Work->pivcol_in_front)) ;
#endif

    Fl = Work->Flblock + fnpiv * fnr_curr ;

    if (Work->pivcol_in_front)
    {
	/* extended pattern and position already in Frows, Frpos.  Values above
	 * the diagonal are already in LU block.  Values on and below the
	 * diagonal are in Wy [0 .. fnrows_extended-1].  Copy into the L
	 * block. */
	fnrows_extended += ccdeg ;
	Wy = Work->Wy ;

	for (i = 0 ; i < fnrows_extended ; i++)
	{
	    Fl [i] = Wy [i] ;
#ifndef NDEBUG
	    row = Frows [i] ;
	    DEBUG2 ((" "ID": row "ID" ", i, row)) ;
	    EDEBUG2 (Fl [i]) ;
	    if (row == Work->pivrow) DEBUG2 ((" <- pivrow")) ;
	    DEBUG2 (("\n")) ;
	    if (i == fnrows - 1) DEBUG2 ((" :::::::\n")) ;
	    ASSERT (row >= 0 && row < Work->n_row) ;
	    ASSERT (Frpos [row] == i) ;
#endif
	}

    }
    else
    {
	/* extended pattern,values is in (Wm,Wx), not yet in the front */
	Entry *F ;
	Fu = Work->Flublock + fnpiv * Work->nb ;
	Wm = Work->Wm ;
	Wx = Work->Wx ;
	F = Fu ;
	for (i = 0 ; i < fnpiv ; i++)
	{
	    CLEAR_AND_INCREMENT (F) ;
	}
	F = Fl ;
	for (i = 0 ; i < fnrows ; i++)
	{
	    CLEAR_AND_INCREMENT (F) ;
	}
	for (i = 0 ; i < ccdeg ; i++)
	{
	    row = Wm [i] ;
#ifndef NDEBUG
	    DEBUG2 ((" "ID": row "ID" (ext) ", fnrows_extended, row)) ;
	    EDEBUG2 (Wx [i]) ;
	    if (row == Work->pivrow) DEBUG2 ((" <- pivrow")) ;
	    DEBUG2 (("\n")) ;
	    ASSERT (row >= 0 && row < Work->n_row) ;
#endif
	    pos = Frpos [row] ;
	    if (pos < 0)
	    {
		pos = fnrows_extended++ ;
		Frows [pos] = row ;
		Frpos [row] = pos ;
	    }
	    Fl [pos] = Wx [i] ;
	}
    }

    ASSERT (fnrows_extended <= fnr_curr) ;

    /* ---------------------------------------------------------------------- */
    /* extend the column pattern of the front with the new pivot row */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    DEBUG6 (("Pivot row, before extension: "ID"\n", fncols)) ;
    for (j = 0 ; j < fncols ; j++)
    {
	DEBUG7 ((" "ID": col "ID"\n", j, Fcols [j])) ;
	ASSERT (Fcpos [Fcols [j]] == j * fnr_curr) ;
    }
    DEBUG6 (("Extending pivot row:\n")) ;
#endif

    if (Work->pivrow_in_front)
    {
	if (Work->pivcol_in_front)
	{
	    ASSERT (Fcols == Work->Wrow) ;
	    for (j = fncols ; j < rrdeg ; j++)
	    {
#ifndef NDEBUG
		col = Fcols [j] ;
		DEBUG2 ((" "ID": col "ID" (ext)\n", j, col)) ;
		ASSERT (col != Work->pivcol) ;
		ASSERT (col >= 0 && col < Work->n_col) ;
		ASSERT (Fcpos [col] < 0) ;
#endif
		Fcpos [Fcols [j]] = j * fnr_curr ;
	    }
	}
	else
	{
	    /* OUT-IN option: pivcol not in front, but pivrow is in front */
	    Wrow = Work->Wrow ;
	    ASSERT (IMPLIES (Work->pivcol_in_front, Wrow == Fcols)) ;
	    if (Wrow == Fcols)
	    {
		/* Wrow and Fcols are equivalenced */
		for (j = fncols ; j < rrdeg ; j++)
		{
		    col = Wrow [j] ;
		    DEBUG2 ((" "ID": col "ID" (ext)\n", j, col)) ;
		    ASSERT (Fcpos [col] < 0) ;
		    /* Fcols [j] = col ;  not needed */
		    Fcpos [col] = j * fnr_curr ;
		}
	    }
	    else
	    {
		for (j = fncols ; j < rrdeg ; j++)
		{
		    col = Wrow [j] ;
		    DEBUG2 ((" "ID": col "ID" (ext)\n", j, col)) ;
		    ASSERT (Fcpos [col] < 0) ;
		    Fcols [j] = col ;
		    Fcpos [col] = j * fnr_curr ;
		}
	    }
	}
	fncols_extended = rrdeg ;
    }
    else
    {
	ASSERT (Fcols != Work->Wrow) ;
	Wrow = Work->Wrow ;
	for (j = 0 ; j < rrdeg ; j++)
	{
	    col = Wrow [j] ;
	    ASSERT (col >= 0 && col < Work->n_col) ;
	    if (Fcpos [col] < 0)
	    {
		DEBUG2 ((" col:: "ID" (ext)\n", col)) ;
		Fcols [fncols_extended] = col ;
		Fcpos [col] = fncols_extended * fnr_curr ;
		fncols_extended++ ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* pivot row and column have been extended */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    ASSERT (fncols_extended <= fnc_curr) ;
    ASSERT (fnrows_extended <= fnr_curr) ;

    DEBUG6 (("Pivot col, after ext: "ID" "ID"\n", fnrows,fnrows_extended)) ;
    for (i = 0 ; i < fnrows_extended ; i++)
    {
	row = Frows [i] ;
	DEBUG7 ((" "ID": row "ID" pos "ID" old: %d", i, row, Frpos [row],
	    i < fnrows)) ;
	if (row == Work->pivrow ) DEBUG7 (("  <-- pivrow")) ;
	DEBUG7 (("\n")) ;
	ASSERT (Frpos [Frows [i]] == i) ;
    }

    DEBUG6 (("Pivot row position: "ID"\n", Frpos [Work->pivrow])) ;
    ASSERT (Frpos [Work->pivrow] >= 0) ;
    ASSERT (Frpos [Work->pivrow] < fnrows_extended) ;

    DEBUG6 (("Pivot row, after ext: "ID" "ID"\n", fncols,fncols_extended)) ;
    for (j = 0 ; j < fncols_extended ; j++)
    {
	col = Fcols [j] ;
	DEBUG7 ((" "ID": col "ID" pos "ID" old: %d", j, col, Fcpos [col],
	    j < fncols)) ;
	if (col == Work->pivcol ) DEBUG7 (("  <-- pivcol")) ;
	DEBUG7 (("\n")) ;
	ASSERT (Fcpos [Fcols [j]] == j * fnr_curr) ;
    }

    DEBUG6 (("Pivot col position: "ID"\n", Fcpos [Work->pivcol])) ;
    ASSERT (Fcpos [Work->pivcol] >= 0) ;
    ASSERT (Fcpos [Work->pivcol] < fncols_extended * fnr_curr) ;

#endif

    /* ---------------------------------------------------------------------- */
    /* Zero the newly extended frontal matrix */
    /* ---------------------------------------------------------------------- */

    zero_front (Work->Flblock, Work->Fublock, Work->Fcblock,
	fnrows, fncols, fnr_curr, fnc_curr,
	fnpiv, fnrows_extended, fncols_extended) ;

    /* ---------------------------------------------------------------------- */
    /* finalize extended row and column pattern of the frontal matrix */
    /* ---------------------------------------------------------------------- */

    Work->fnrows = fnrows_extended ;
    Work->fncols = fncols_extended ;

    ASSERT (fnrows_extended == Work->fnrows_new + 1) ;
    ASSERT (fncols_extended == Work->fncols_new + 1) ;

    return (TRUE) ;

}
