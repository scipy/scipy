/* ========================================================================== */
/* === UMF_row_search ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Find two candidate pivot rows in a column: the best one in the front,
    and the best one not in the front.  Return the two pivot row patterns and
    their exact degrees.  Called by UMF_local_search.

    Returns UMFPACK_OK if successful, or UMFPACK_WARNING_singular_matrix or
    UMFPACK_ERROR_different_pattern if not.

*/

#include "umf_internal.h"
#include "umf_row_search.h"

GLOBAL Int UMF_row_search
(
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic,
    Int cdeg0,			/* length of column in Front */
    Int cdeg1,			/* length of column outside Front */
    const Int Pattern [ ],	/* pattern of column, Pattern [0..cdeg1 -1] */
    const Int Pos [ ],		/* Pos [Pattern [0..cdeg1 -1]] = 0..cdeg1 -1 */
    Int pivrow [2],		/* pivrow [IN] and pivrow [OUT] */
    Int rdeg [2],		/* rdeg [IN] and rdeg [OUT] */
    Int W_i [ ],		/* pattern of pivrow [IN], */
				/* either Fcols or Woi */
    Int W_o [ ],		/* pattern of pivrow [OUT], */
				/* either Wio or Woo */
    Int prior_pivrow [2],	/* the two other rows just scanned, if any */
    const Entry Wxy [ ],	/* numerical values Wxy [0..cdeg1-1],
				   either Wx or Wy */

    Int pivcol,			/* the candidate column being searched */
    Int freebie [ ]
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    double maxval, toler, toler2, value, pivot [2] ;
    Int i, row, deg, col, *Frpos, fnrows, *E, j, ncols, *Cols, *Rows,
	e, f, Wrpflag, *Fcpos, fncols, tpi, max_rdeg, nans_in_col, was_offdiag,
	diag_row, prefer_diagonal, *Wrp, found, *Diagonal_map ;
    Tuple *tp, *tpend, *tp1, *tp2 ;
    Unit *Memory, *p ;
    Element *ep ;
    Int *Row_tuples, *Row_degree, *Row_tlen ;

#ifndef NDEBUG
    Int *Col_degree ;
    DEBUG2 (("Row_search:\n")) ;
    for (i = 0 ; i < cdeg1 ; i++)
    {
	row = Pattern [i] ;
	DEBUG4 (("   row: "ID"\n", row)) ;
	ASSERT (row >= 0 && row < Numeric->n_row) ;
	ASSERT (i == Pos [row]) ;
    }
    /* If row is not in Pattern [0..cdeg1-1], then Pos [row] == EMPTY */
    if (UMF_debug > 0 || Numeric->n_row < 1000)
    {
	Int cnt = cdeg1 ;
	DEBUG4 (("Scan all rows:\n")) ;
	for (row = 0 ; row < Numeric->n_row ; row++)
	{
	    if (Pos [row] < 0)
	    {
		cnt++ ;
	    }
	    else
	    {
		DEBUG4 (("   row: "ID" pos "ID"\n", row, Pos [row])) ;
	    }
	}
	ASSERT (cnt == Numeric->n_row) ;
    }
    Col_degree = Numeric->Cperm ;   /* for NON_PIVOTAL_COL macro only */
    ASSERT (pivcol >= 0 && pivcol < Work->n_col) ;
    ASSERT (NON_PIVOTAL_COL (pivcol)) ;
#endif

    pivot [IN] = 0. ;
    pivot [OUT] = 0. ;

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    Row_degree = Numeric->Rperm ;
    Row_tuples = Numeric->Uip ;
    Row_tlen   = Numeric->Uilen ;
    Wrp = Work->Wrp ;
    Frpos = Work->Frpos ;
    E = Work->E ;
    Memory = Numeric->Memory ;
    fnrows = Work->fnrows ;

    prefer_diagonal = Symbolic->prefer_diagonal ;
    Diagonal_map = Work->Diagonal_map ;

    if (Diagonal_map)
    {
	diag_row = Diagonal_map [pivcol] ;
	was_offdiag = diag_row < 0 ;
	if (was_offdiag)
	{
	    /* the "diagonal" entry in this column was permuted here by an
	     * earlier pivot choice.  The tighter off-diagonal tolerance will
	     * be used instead of the symmetric tolerance. */
	    diag_row = FLIP (diag_row) ;
	}
	ASSERT (diag_row >= 0 && diag_row < Numeric->n_row) ;
    }
    else
    {
	diag_row = EMPTY ;	/* unused */
	was_offdiag = EMPTY ;	/* unused */
    }

    /* pivot row degree cannot exceed max_rdeg */
    max_rdeg = Work->fncols_max ;

    /* ---------------------------------------------------------------------- */
    /* scan pivot column for candidate rows */
    /* ---------------------------------------------------------------------- */

    maxval = 0.0 ;
    nans_in_col = FALSE ;

    for (i = 0 ; i < cdeg1 ; i++)
    {
	APPROX_ABS (value, Wxy [i]) ;
	if (SCALAR_IS_NAN (value))
	{
	    nans_in_col = TRUE ;
	    maxval = value ;
	    break ;
	}
	/* This test can now ignore the NaN case: */
	maxval = MAX (maxval, value) ;
    }

    /* if maxval is zero, the matrix is numerically singular */

    toler = Numeric->relpt * maxval ;
    toler2 = Numeric->relpt2 * maxval ;
    toler2 = was_offdiag ? toler : toler2 ;

    DEBUG5 (("Row_search begins [ maxval %g toler %g %g\n",
	maxval, toler, toler2)) ;
    if (SCALAR_IS_NAN (toler) || SCALAR_IS_NAN (toler2))
    {
	nans_in_col = TRUE ;
    }

    if (!nans_in_col)
    {

	/* look for the diagonal entry, if it exists */
	found = FALSE ;
	ASSERT (!SCALAR_IS_NAN (toler)) ;

	if (prefer_diagonal)
	{
	    ASSERT (diag_row != EMPTY) ;
	    i = Pos [diag_row] ;
	    if (i >= 0)
	    {
		double a ;
		ASSERT (i < cdeg1) ;
		ASSERT (diag_row == Pattern [i]) ;

		APPROX_ABS (a, Wxy [i]) ;

		ASSERT (!SCALAR_IS_NAN (a)) ;
		ASSERT (!SCALAR_IS_NAN (toler2)) ;

		if (SCALAR_IS_NONZERO (a) && a >= toler2)
		{
		    /* found it! */
		    DEBUG3 (("Symmetric pivot: "ID" "ID"\n", pivcol, diag_row));
		    found = TRUE ;
		    if (Frpos [diag_row] >= 0 && Frpos [diag_row] < fnrows)
		    {
			pivrow [IN] = diag_row ;
			pivrow [OUT] = EMPTY ;
		    }
		    else
		    {
			pivrow [IN] = EMPTY ;
			pivrow [OUT] = diag_row ;
		    }
		}
	    }
	}

	/* either no diagonal found, or we didn't look for it */
	if (!found)
	{
	    if (cdeg0 > 0)
	    {

		/* this is a column in the front */
		for (i = 0 ; i < cdeg0 ; i++)
		{
		    double a ;
		    APPROX_ABS (a, Wxy [i]) ;
		    ASSERT (!SCALAR_IS_NAN (a)) ;
		    ASSERT (!SCALAR_IS_NAN (toler)) ;
		    if (SCALAR_IS_NONZERO (a) && a >= toler)
		    {
			row = Pattern [i] ;
			deg = Row_degree [row] ;
#ifndef NDEBUG
			DEBUG6 ((ID" candidate row "ID" deg "ID" absval %g\n",
			    i, row, deg, a)) ;
			UMF_dump_rowcol (0, Numeric, Work, row, TRUE) ;
#endif
			ASSERT (Frpos [row] >= 0 && Frpos [row] < fnrows) ;
			ASSERT (Frpos [row] == i) ;
			/* row is in the current front */
			DEBUG4 ((" in front\n")) ;
			if (deg < rdeg [IN]
			    /* break ties by picking the largest entry: */
			       || (deg == rdeg [IN] && a > pivot [IN])
			    /* break ties by picking the diagonal entry: */
			    /* || (deg == rdeg [IN] && row == diag_row) */
			   )
			{
			    /* best row in front, so far */
			    pivrow [IN] = row ;
			    rdeg [IN] = deg ;
			    pivot [IN] = a ;
			}
		    }
		}
		for ( ; i < cdeg1 ; i++)
		{
		    double a ;
		    APPROX_ABS (a, Wxy [i]) ;
		    ASSERT (!SCALAR_IS_NAN (a)) ;
		    ASSERT (!SCALAR_IS_NAN (toler)) ;
		    if (SCALAR_IS_NONZERO (a) && a >= toler)
		    {
			row = Pattern [i] ;
			deg = Row_degree [row] ;
#ifndef NDEBUG
			DEBUG6 ((ID" candidate row "ID" deg "ID" absval %g\n",
			    i, row, deg, a)) ;
			UMF_dump_rowcol (0, Numeric, Work, row, TRUE) ;
#endif
			ASSERT (Frpos [row] == i) ;
			/* row is not in the current front */
			DEBUG4 ((" NOT in front\n")) ;
			if (deg < rdeg [OUT]
			    /* break ties by picking the largest entry: */
			       || (deg == rdeg [OUT] && a > pivot [OUT])
			    /* break ties by picking the diagonal entry: */
			    /* || (deg == rdeg [OUT] && row == diag_row) */
			   )
			{
			    /* best row not in front, so far */
			    pivrow [OUT] = row ;
			    rdeg [OUT] = deg ;
			    pivot [OUT] = a ;
			}
		    }
		}

	    }
	    else
	    {

		/* this column is not in the front */
		for (i = 0 ; i < cdeg1 ; i++)
		{
		    double a ;
		    APPROX_ABS (a, Wxy [i]) ;
		    ASSERT (!SCALAR_IS_NAN (a)) ;
		    ASSERT (!SCALAR_IS_NAN (toler)) ;
		    if (SCALAR_IS_NONZERO (a) && a >= toler)
		    {
			row = Pattern [i] ;
			deg = Row_degree [row] ;
#ifndef NDEBUG
			DEBUG6 ((ID" candidate row "ID" deg "ID" absval %g\n",
			    i, row, deg, a)) ;
			UMF_dump_rowcol (0, Numeric, Work, row, TRUE) ;
#endif
			if (Frpos [row] >= 0 && Frpos [row] < fnrows)
			{
			    /* row is in the current front */
			    DEBUG4 ((" in front\n")) ;
			    if (deg < rdeg [IN]
			    /* break ties by picking the largest entry: */
			       || (deg == rdeg [IN] && a > pivot [IN])
			    /* break ties by picking the diagonal entry: */
			    /* || (deg == rdeg [IN] && row == diag_row) */
			       )
			    {
				/* best row in front, so far */
				pivrow [IN] = row ;
				rdeg [IN] = deg ;
				pivot [IN] = a ;
			    }
			}
			else
			{
			    /* row is not in the current front */
			    DEBUG4 ((" NOT in front\n")) ;
			    if (deg < rdeg [OUT]
			    /* break ties by picking the largest entry: */
			       || (deg == rdeg[OUT] && a > pivot [OUT])
			    /* break ties by picking the diagonal entry: */
			    /* || (deg == rdeg[OUT] && row == diag_row) */
			       )
			    {
				/* best row not in front, so far */
				pivrow [OUT] = row ;
				rdeg [OUT] = deg ;
				pivot [OUT] = a ;
			    }
			}
		    }
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* NaN handling */
    /* ---------------------------------------------------------------------- */

    /* if cdeg1 > 0 then we must have found a pivot row ... unless NaN's */
    /* exist.  Try with no numerical tests if no pivot found. */

    if (cdeg1 > 0 && pivrow [IN] == EMPTY && pivrow [OUT] == EMPTY)
    {
	/* cleanup for the NaN case */
	DEBUG0 (("Found a NaN in pivot column!\n")) ;

	/* grab the first entry in the pivot column, ignoring degree, */
	/* numerical stability, and symmetric preference */
	row = Pattern [0] ;
	deg = Row_degree [row] ;
	if (Frpos [row] >= 0 && Frpos [row] < fnrows)
	{
	    /* row is in the current front */
	    DEBUG4 ((" in front\n")) ;
	    pivrow [IN] = row ;
	    rdeg [IN] = deg ;
	}
	else
	{
	    /* row is not in the current front */
	    DEBUG4 ((" NOT in front\n")) ;
	    pivrow [OUT] = row ;
	    rdeg [OUT] = deg ;
	}

	/* We are now guaranteed to have a pivot, no matter how broken */
	/* (non-IEEE compliant) the underlying numerical operators are. */
	/* This is particularly a problem for Microsoft compilers (they do */
	/* not handle NaN's properly). Now try to find a sparser pivot, if */
	/* possible. */

	for (i = 1 ; i < cdeg1 ; i++)
	{
	    row = Pattern [i] ;
	    deg = Row_degree [row] ;

	    if (Frpos [row] >= 0 && Frpos [row] < fnrows)
	    {
		/* row is in the current front */
		DEBUG4 ((" in front\n")) ;
		if (deg < rdeg [IN] || (deg == rdeg [IN] && row == diag_row))
		{
		    /* best row in front, so far */
		    pivrow [IN] = row ;
		    rdeg [IN] = deg ;
		}
	    }
	    else
	    {
		/* row is not in the current front */
		DEBUG4 ((" NOT in front\n")) ;
		if (deg < rdeg [OUT] || (deg == rdeg [OUT] && row == diag_row))
		{
		    /* best row not in front, so far */
		    pivrow [OUT] = row ;
		    rdeg [OUT] = deg ;
		}
	    }
	}
    }

    /* We found a pivot if there are entries (even zero ones) in pivot col */
    ASSERT (IMPLIES (cdeg1 > 0, pivrow[IN] != EMPTY || pivrow[OUT] != EMPTY)) ;

    /* If there are no entries in the pivot column, then no pivot is found */
    ASSERT (IMPLIES (cdeg1 == 0, pivrow[IN] == EMPTY && pivrow[OUT] == EMPTY)) ;

    /* ---------------------------------------------------------------------- */
    /* check for singular matrix */
    /* ---------------------------------------------------------------------- */

    if (cdeg1  == 0)
    {
	if (fnrows > 0)
	{
	    /*
		Get the pivrow [OUT][IN] from the current front.
		The frontal matrix looks like this:

			pivcol[OUT]
			|
			v
		x x x x 0   <- so grab this row as the pivrow [OUT][IN].
		x x x x 0
		x x x x 0
		0 0 0 0 0

		The current frontal matrix has some rows in it.  The degree
		of the pivcol[OUT] is zero.  The column is empty, and the
		current front does not contribute to it.

	    */
	    pivrow [IN] = Work->Frows [0] ;
	    DEBUGm4 (("Got zero pivrow[OUT][IN] "ID" from current front\n",
		pivrow [IN])) ;
	}
	else
	{

	    /*
		Get a pivot row from the row-merge tree, use as
		pivrow [OUT][OUT].   pivrow [IN] remains EMPTY.
		This can only happen if the current front is 0-by-0.
	    */

	    Int *Front_leftmostdesc, *Front_1strow, *Front_new1strow, row1,
		row2, fleftmost, nfr, n_row, frontid ;

	    ASSERT (Work->fncols == 0) ;

	    Front_leftmostdesc = Symbolic->Front_leftmostdesc ;
	    Front_1strow = Symbolic->Front_1strow ;
	    Front_new1strow = Work->Front_new1strow ;
	    nfr = Symbolic->nfr ;
	    n_row = Numeric->n_row ;
	    frontid = Work->frontid ;

	    DEBUGm4 (("Note: pivcol: "ID" is empty front "ID"\n",
		pivcol, frontid)) ;
#ifndef NDEBUG
	    DEBUG1 (("Calling dump rowmerge\n")) ;
	    UMF_dump_rowmerge (Numeric, Symbolic, Work) ;
#endif

	    /* Row-merge set is the non-pivotal rows in the range */
	    /* Front_new1strow [Front_leftmostdesc [frontid]] to */
	    /* Front_1strow [frontid+1] - 1. */
	    /* If this is empty, then use the empty rows, in the range */
	    /* Front_new1strow [nfr] to n_row-1. */
	    /* If this too is empty, then pivrow [OUT] will be empty. */
	    /* In both cases, update Front_new1strow [...]. */

	    fleftmost = Front_leftmostdesc [frontid] ;
	    row1 = Front_new1strow [fleftmost] ;
	    row2 = Front_1strow [frontid+1] - 1 ;
	    DEBUG1 (("Leftmost: "ID" Rows ["ID" to "ID"] srch ["ID" to "ID"]\n",
		fleftmost, Front_1strow [frontid], row2, row1, row2)) ;

	    /* look in the range row1 ... row2 */
	    for (row = row1 ; row <= row2 ; row++)
	    {
		DEBUG3 (("   Row: "ID"\n", row)) ;
		if (NON_PIVOTAL_ROW (row))
		{
		    /* found it */
		    DEBUG3 (("   Row: "ID" found\n", row)) ;
		    ASSERT (Frpos [row] == EMPTY) ;
		    pivrow [OUT] = row ;
		    DEBUGm4 (("got row merge pivrow %d\n", pivrow [OUT])) ;
		    break ;
		}
	    }
	    Front_new1strow [fleftmost] = row ;

	    if (pivrow [OUT] == EMPTY)
	    {
		/* not found, look in empty row set in "dummy" front */
		row1 = Front_new1strow [nfr] ;
		row2 = n_row-1 ;
		DEBUG3 (("Empty: "ID" Rows ["ID" to "ID"] srch["ID" to "ID"]\n",
		    nfr, Front_1strow [nfr], row2, row1, row2)) ;

		/* look in the range row1 ... row2 */
		for (row = row1 ; row <= row2 ; row++)
		{
		    DEBUG3 (("   Empty Row: "ID"\n", row)) ;
		    if (NON_PIVOTAL_ROW (row))
		    {
			/* found it */
			DEBUG3 (("   Empty Row: "ID" found\n", row)) ;
			ASSERT (Frpos [row] == EMPTY) ;
			pivrow [OUT] = row ;
			DEBUGm4 (("got dummy row pivrow %d\n", pivrow [OUT])) ;
			break ;
		    }
		}
		Front_new1strow [nfr] = row ;
	    }

	    if (pivrow [OUT] == EMPTY)
	    {
		/* Row-merge set is empty.  We can just discard */
		/* the candidate pivot column. */
		DEBUG0 (("Note: row-merge set empty\n")) ;
		DEBUGm4 (("got no pivrow \n")) ;
		return (UMFPACK_WARNING_singular_matrix) ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* construct the candidate row in the front, if any */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    /* check Wrp */
    ASSERT (Work->Wrpflag > 0) ;
    if (UMF_debug > 0 || Work->n_col < 1000)
    {
	for (i = 0 ; i < Work->n_col ; i++)
	{
	    ASSERT (Wrp [i] < Work->Wrpflag) ;
	}
    }
#endif

#ifndef NDEBUG
    DEBUG4 (("pivrow [IN]: "ID"\n", pivrow [IN])) ;
    UMF_dump_rowcol (0, Numeric, Work, pivrow [IN], TRUE) ;
#endif

    if (pivrow [IN] != EMPTY)
    {

	/* the row merge candidate row is not pivrow [IN] */
	freebie [IN] = (pivrow [IN] == prior_pivrow [IN]) && (cdeg1  > 0) ;
	ASSERT (cdeg1  >= 0) ;

	if (!freebie [IN])
	{
	    /* include current front in the degree of this row */

	    Fcpos = Work->Fcpos ;
	    fncols = Work->fncols ;

	    Wrpflag = Work->Wrpflag ;

	    /* -------------------------------------------------------------- */
	    /* construct the pattern of the IN row */
	    /* -------------------------------------------------------------- */

#ifndef NDEBUG
	    /* check Fcols */
	    DEBUG5 (("ROW ASSEMBLE: rdeg "ID"\nREDUCE ROW "ID"\n",
		fncols, pivrow [IN])) ;
	    for (j = 0 ; j < fncols ; j++)
	    {
		col = Work->Fcols [j] ;
		ASSERT (col >= 0 && col < Work->n_col) ;
		ASSERT (Fcpos [col] >= 0) ;
	    }
	    if (UMF_debug > 0 || Work->n_col < 1000)
	    {
		Int cnt = fncols ;
		for (col = 0 ; col < Work->n_col ; col++)
		{
		    if (Fcpos [col] < 0) cnt++ ;
		}
		ASSERT (cnt == Work->n_col) ;
	    }
#endif

	    rdeg [IN] = fncols ;

	    ASSERT (pivrow [IN] >= 0 && pivrow [IN] < Work->n_row) ;
	    ASSERT (NON_PIVOTAL_ROW (pivrow [IN])) ;

	    /* add the pivot column itself */
	    ASSERT (Wrp [pivcol] != Wrpflag) ;
	    if (Fcpos [pivcol] < 0)
	    {
		DEBUG3 (("Adding pivot col to pivrow [IN] pattern\n")) ;
		if (rdeg [IN] >= max_rdeg)
		{
		    /* :: pattern change (in) :: */
		    return (UMFPACK_ERROR_different_pattern) ;
		}
		Wrp [pivcol] = Wrpflag ;
		W_i [rdeg [IN]++] = pivcol ;
	    }

	    tpi = Row_tuples [pivrow [IN]] ;
	    if (tpi)
	    {
		tp = (Tuple *) (Memory + tpi) ;
		tp1 = tp ;
		tp2 = tp ;
		tpend = tp + Row_tlen [pivrow [IN]] ;
		for ( ; tp < tpend ; tp++)
		{
		    e = tp->e ;
		    ASSERT (e > 0 && e <= Work->nel) ;
		    if (!E [e])
		    {
			continue ;	/* element already deallocated */
		    }
		    f = tp->f ;
		    p = Memory + E [e] ;
		    ep = (Element *) p ;
		    p += UNITS (Element, 1) ;
		    Cols = (Int *) p ;
		    ncols = ep->ncols ;
		    Rows = Cols + ncols ;
		    if (Rows [f] == EMPTY)
		    {
			continue ;	/* row already assembled */
		    }
		    ASSERT (pivrow [IN] == Rows [f]) ;

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			ASSERT (col >= EMPTY && col < Work->n_col) ;
			if ((col >= 0) && (Wrp [col] != Wrpflag)
			    && Fcpos [col] <0)
			{
			    ASSERT (NON_PIVOTAL_COL (col)) ;
			    if (rdeg [IN] >= max_rdeg)
			    {
				/* :: pattern change (rdeg in failure) :: */
				DEBUGm4 (("rdeg [IN] >= max_rdeg failure\n")) ;
				return (UMFPACK_ERROR_different_pattern) ;
			    }
			    Wrp [col] = Wrpflag ;
			    W_i [rdeg [IN]++] = col ;
			}
		    }

		    *tp2++ = *tp ;	/* leave the tuple in the list */
		}
		Row_tlen [pivrow [IN]] = tp2 - tp1 ;
	    }

#ifndef NDEBUG
	    DEBUG4 (("Reduced IN row:\n")) ;
	    for (j = 0 ; j < fncols ; j++)
	    {
		DEBUG6 ((" "ID" "ID" "ID"\n",
		    j, Work->Fcols [j], Fcpos [Work->Fcols [j]])) ;
		ASSERT (Fcpos [Work->Fcols [j]] >= 0) ;
	    }
	    for (j = fncols ; j < rdeg [IN] ; j++)
	    {
		DEBUG6 ((" "ID" "ID" "ID"\n", j, W_i [j], Wrp [W_i [j]]));
		ASSERT (W_i [j] >= 0 && W_i [j] < Work->n_col) ;
		ASSERT (Wrp [W_i [j]] == Wrpflag) ;
	    }
	    /* mark the end of the pattern in case we scan it by mistake */
	    /* Note that this means W_i must be of size >= fncols_max + 1 */
	    W_i [rdeg [IN]] = EMPTY ;
#endif

	    /* rdeg [IN] is now the exact degree of the IN row */

	    /* clear Work->Wrp. */
	    Work->Wrpflag++ ;
	    /* All Wrp [0..n_col] is now < Wrpflag */
	}
    }

#ifndef NDEBUG
    /* check Wrp */
    if (UMF_debug > 0 || Work->n_col < 1000)
    {
	for (i = 0 ; i < Work->n_col ; i++)
	{
	    ASSERT (Wrp [i] < Work->Wrpflag) ;
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* construct the candidate row not in the front, if any */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    DEBUG4 (("pivrow [OUT]: "ID"\n", pivrow [OUT])) ;
    UMF_dump_rowcol (0, Numeric, Work, pivrow [OUT], TRUE) ;
#endif

    /* If this is a candidate row from the row merge set, force it to be */
    /* scanned (ignore prior_pivrow [OUT]). */

    if (pivrow [OUT] != EMPTY)
    {
	freebie [OUT] = (pivrow [OUT] == prior_pivrow [OUT]) && cdeg1  > 0 ;
	ASSERT (cdeg1  >= 0) ;

	if (!freebie [OUT])
	{

	    Wrpflag = Work->Wrpflag ;

	    /* -------------------------------------------------------------- */
	    /* construct the pattern of the row */
	    /* -------------------------------------------------------------- */

	    rdeg [OUT] = 0 ;

	    ASSERT (pivrow [OUT] >= 0 && pivrow [OUT] < Work->n_row) ;
	    ASSERT (NON_PIVOTAL_ROW (pivrow [OUT])) ;

	    /* add the pivot column itself */
	    ASSERT (Wrp [pivcol] != Wrpflag) ;
	    DEBUG3 (("Adding pivot col to pivrow [OUT] pattern\n")) ;
	    if (rdeg [OUT] >= max_rdeg)
	    {
		/* :: pattern change (out) :: */
		return (UMFPACK_ERROR_different_pattern) ;
	    }
	    Wrp [pivcol] = Wrpflag ;
	    W_o [rdeg [OUT]++] = pivcol ;

	    tpi = Row_tuples [pivrow [OUT]] ;
	    if (tpi)
	    {
		tp = (Tuple *) (Memory + tpi) ;
		tp1 = tp ;
		tp2 = tp ;
		tpend = tp + Row_tlen [pivrow [OUT]] ;
		for ( ; tp < tpend ; tp++)
		{
		    e = tp->e ;
		    ASSERT (e > 0 && e <= Work->nel) ;
		    if (!E [e])
		    {
			continue ;	/* element already deallocated */
		    }
		    f = tp->f ;
		    p = Memory + E [e] ;
		    ep = (Element *) p ;
		    p += UNITS (Element, 1) ;
		    Cols = (Int *) p ;
		    ncols = ep->ncols ;
		    Rows = Cols + ncols ;
		    if (Rows [f] == EMPTY)
		    {
			continue ;	/* row already assembled */
		    }
		    ASSERT (pivrow [OUT] == Rows [f]) ;

		    for (j = 0 ; j < ncols ; j++)
		    {
			col = Cols [j] ;
			ASSERT (col >= EMPTY && col < Work->n_col) ;
			if ((col >= 0) && (Wrp [col] != Wrpflag))
			{
			    ASSERT (NON_PIVOTAL_COL (col)) ;
			    if (rdeg [OUT] >= max_rdeg)
			    {
				/* :: pattern change (rdeg out failure) :: */
				DEBUGm4 (("rdeg [OUT] failure\n")) ;
				return (UMFPACK_ERROR_different_pattern) ;
			    }
			    Wrp [col] = Wrpflag ;
			    W_o [rdeg [OUT]++] = col ;
			}
		    }
		    *tp2++ = *tp ;	/* leave the tuple in the list */
		}
		Row_tlen [pivrow [OUT]] = tp2 - tp1 ;
	    }

#ifndef NDEBUG
	    DEBUG4 (("Reduced row OUT:\n")) ;
	    for (j = 0 ; j < rdeg [OUT] ; j++)
	    {
		DEBUG6 ((" "ID" "ID" "ID"\n", j, W_o [j], Wrp [W_o [j]])) ;
		ASSERT (W_o [j] >= 0 && W_o [j] < Work->n_col) ;
		ASSERT (Wrp [W_o [j]] == Wrpflag) ;
	    }
	    /* mark the end of the pattern in case we scan it by mistake */
	    /* Note that this means W_o must be of size >= fncols_max + 1 */
	    W_o [rdeg [OUT]] = EMPTY ;
#endif

	    /* rdeg [OUT] is now the exact degree of the row */

	    /* clear Work->Wrp. */
	    Work->Wrpflag++ ;
	    /* All Wrp [0..n] is now < Wrpflag */

	}

    }
    DEBUG5 (("Row_search end ] \n")) ;

#ifndef NDEBUG
    /* check Wrp */
    if (UMF_debug > 0 || Work->n_col < 1000)
    {
	for (i = 0 ; i < Work->n_col ; i++)
	{
	    ASSERT (Wrp [i] < Work->Wrpflag) ;
	}
    }
#endif

    return (UMFPACK_OK) ;
}
