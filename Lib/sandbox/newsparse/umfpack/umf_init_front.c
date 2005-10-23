/* ========================================================================== */
/* === UMF_init_front ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

#include "umf_internal.h"
#include "umf_grow_front.h"

/* ========================================================================== */
/* === zero_init_front ====================================================== */
/* ========================================================================== */

/* Set the initial frontal matrix to zero. */

PRIVATE void zero_init_front (Int m, Int n, Entry *Fcblock, Int d)
{
    Int i, j ;
    Entry *F, *Fj = Fcblock ;
    for (j = 0 ; j < m ; j++)
    {
	F = Fj ;
	Fj += d ;
	for (i = 0 ; i < n ; i++)
	{
	    /* CLEAR (Fcblock [i + j*d]) ; */
	    CLEAR (*F) ;
	    F++ ;
	}
    }
}

/* ========================================================================== */
/* === UMF_init_front ======================================================= */
/* ========================================================================== */

GLOBAL Int UMF_init_front
(
    NumericType *Numeric,
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int i, j, fnr_curr, row, col, *Frows, *Fcols,
	*Fcpos, *Frpos, fncols, fnrows, *Wrow, fnr2, fnc2, rrdeg, ccdeg, *Wm,
	fnrows_extended ;
    Entry *Fcblock, *Fl, *Wy, *Wx ;

    /* ---------------------------------------------------------------------- */
    /* get current frontal matrix and check for frontal growth */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    DEBUG0 (("INIT FRONT\n")) ;
    DEBUG1 (("CURR before init:\n")) ;
    UMF_dump_current_front (Numeric, Work, FALSE) ;
#endif
    if (Work->do_grow)
    {
	fnr2 = UMF_FRONTAL_GROWTH * Work->fnrows_new + 2 ;
	fnc2 = UMF_FRONTAL_GROWTH * Work->fncols_new + 2 ;
	if (!UMF_grow_front (Numeric, fnr2, fnc2, Work,
	    Work->pivrow_in_front ? 2 : 0))
	{
	    /* :: out of memory in umf_init_front :: */
	    DEBUGm4 (("out of memory: init front\n")) ;
	    return (FALSE) ;
	}
    }
#ifndef NDEBUG
    DEBUG1 (("CURR after grow:\n")) ;
    UMF_dump_current_front (Numeric, Work, FALSE) ;
    DEBUG1 (("fnrows new "ID" fncols new "ID"\n",
	Work->fnrows_new, Work->fncols_new)) ;
#endif
    ASSERT (Work->fnpiv == 0) ;
    fnr_curr = Work->fnr_curr ;
    ASSERT (Work->fnrows_new + 1 <= fnr_curr) ;
    ASSERT (Work->fncols_new + 1 <= Work->fnc_curr) ;
    ASSERT (fnr_curr >= 0) ;
    ASSERT (fnr_curr % 2 == 1) ;

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    /* current front is defined by pivot row and column */

    Frows = Work->Frows ;
    Fcols = Work->Fcols ;
    Frpos = Work->Frpos ;
    Fcpos = Work->Fcpos ;

    Work->fnzeros = 0 ;

    ccdeg = Work->ccdeg ;
    rrdeg = Work->rrdeg ;

    fnrows = Work->fnrows ;
    fncols = Work->fncols ;

    /* if both pivrow and pivcol are in front, then we extend the old one */
    /* in UMF_extend_front, rather than starting a new one here. */
    ASSERT (! (Work->pivrow_in_front && Work->pivcol_in_front)) ;

    /* ---------------------------------------------------------------------- */
    /* place pivot column pattern in frontal matrix */
    /* ---------------------------------------------------------------------- */

    Fl = Work->Flblock ;

    if (Work->pivcol_in_front)
    {
	/* Append the pivot column extension.
	 * Note that all we need to do is increment the size, since the
	 * candidate pivot column pattern is already in place in
	 * Frows [0 ... fnrows-1] (the old pattern), and
	 * Frows [fnrows ... fnrows + Work->ccdeg - 1] (the new
	 * pattern).  Frpos is also properly defined. */
	/* make a list of the new rows to scan */
	Work->fscan_row = fnrows ;	/* only scan the new rows */
	Work->NewRows = Work->Wrp ;
	Wy = Work->Wy ;
	for (i = 0 ; i < fnrows ; i++)
	{
	    Fl [i] = Wy [i] ;
	}
	fnrows_extended = fnrows + ccdeg ;
	for (i = fnrows ; i < fnrows_extended ; i++)
	{
	    Fl [i] = Wy [i] ;
	    /* flip the row index, since Wrp must be < 0 */
	    row = Frows [i] ;
	    Work->NewRows [i] = FLIP (row) ;
	}
	fnrows = fnrows_extended ;
    }
    else
    {
	/* this is a completely new column */
	Work->fscan_row = 0 ;			/* scan all the rows */
	Work->NewRows = Frows ;
	Wm = Work->Wm ;
	Wx = Work->Wx ;
	for (i = 0 ; i < ccdeg ; i++)
	{
	    Fl [i] = Wx [i] ;
	    row = Wm [i] ;
	    Frows [i] = row ;
	    Frpos [row] = i ;
	}
	fnrows = ccdeg ;
    }

    Work->fnrows = fnrows ;

#ifndef NDEBUG
    DEBUG3 (("New Pivot col "ID" now in front, length "ID"\n",
	Work->pivcol, fnrows)) ;
    for (i = 0 ; i < fnrows ; i++)
    {
	DEBUG4 ((" "ID": row "ID"\n", i, Frows [i])) ;
	ASSERT (Frpos [Frows [i]] == i) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* place pivot row pattern in frontal matrix */
    /* ---------------------------------------------------------------------- */

    Wrow = Work->Wrow ;
    if (Work->pivrow_in_front)
    {
	/* append the pivot row extension */
	Work->fscan_col = fncols ;	/* only scan the new columns */
	Work->NewCols = Work->Wp ;
#ifndef NDEBUG
	for (j = 0 ; j < fncols ; j++)
	{
	    col = Fcols [j] ;
	    ASSERT (col >= 0 && col < Work->n_col) ;
	    ASSERT (Fcpos [col] == j * fnr_curr) ;
	}
#endif
	/* Wrow == Fcol for the IN_IN case, and for the OUT_IN case when
	 * the pivrow [IN][IN] happens to be the same as pivrow [OUT][IN].
	 * See UMF_local_search for more details. */
	ASSERT (IMPLIES (Work->pivcol_in_front, Wrow == Fcols)) ;
	if (Wrow == Fcols)
	{
	    for (j = fncols ; j < rrdeg ; j++)
	    {
		col = Wrow [j] ;
		/* Fcols [j] = col ; not needed */
		/* flip the col index, since Wp must be < 0 */
		Work->NewCols [j] = FLIP (col) ;
		Fcpos [col] = j * fnr_curr ;
	    }
	}
	else
	{
	    for (j = fncols ; j < rrdeg ; j++)
	    {
		col = Wrow [j] ;
		Fcols [j] = col ;
		/* flip the col index, since Wp must be < 0 */
		Work->NewCols [j] = FLIP (col) ;
		Fcpos [col] = j * fnr_curr ;
	    }
	}
    }
    else
    {
	/* this is a completely new row */
	Work->fscan_col = 0 ;			/* scan all the columns */
	Work->NewCols = Fcols ;
	for (j = 0 ; j < rrdeg ; j++)
	{
	    col = Wrow [j] ;
	    Fcols [j] = col ;
	    Fcpos [col] = j * fnr_curr ;
	}
    }

    DEBUGm1 (("rrdeg "ID" fncols "ID"\n", rrdeg, fncols)) ;
    fncols = rrdeg ;
    Work->fncols = fncols ;

    /* ---------------------------------------------------------------------- */
    /* clear the frontal matrix */
    /* ---------------------------------------------------------------------- */

    ASSERT (fnrows == Work->fnrows_new + 1) ;
    ASSERT (fncols == Work->fncols_new + 1) ;

    Fcblock = Work->Fcblock ;
    ASSERT (Fcblock != (Entry *) NULL) ;

    zero_init_front (fncols, fnrows, Fcblock, fnr_curr) ;

#ifndef NDEBUG
    DEBUG3 (("New Pivot row "ID" now in front, length "ID" fnr_curr "ID"\n",
		Work->pivrow, fncols, fnr_curr)) ;
    for (j = 0 ; j < fncols ; j++)
    {
	DEBUG4 (("col "ID" position "ID"\n", j, Fcols [j])) ;
	ASSERT (Fcpos [Fcols [j]] == j * fnr_curr) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* current workspace usage: */
    /* ---------------------------------------------------------------------- */

    /* Fcblock [0..fnr_curr-1, 0..fnc_curr-1]: space for the new frontal
     * matrix.  Fcblock (i,j) is located at Fcblock [i+j*fnr_curr] */

    return (TRUE) ;

}
