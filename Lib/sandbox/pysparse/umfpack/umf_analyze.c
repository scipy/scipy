/* ========================================================================== */
/* === UMF_analyze ========================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Symbolic LL' factorization of A'*A, to get upper bounds on the size of
    L and U for LU = PAQ, and to determine the frontal matrices and
    (supernodal) column elimination tree.  No fill-reducing column pre-ordering
    is used.

    Returns TRUE if successful, FALSE if out of memory.  UMF_analyze can only
    run out of memory if anzmax (which is Ap [n_row]) is too small.

    Uses workspace of size O(nonzeros in A).  On input, the matrix A is
    stored in row-form at the tail end of Ai.  It is destroyed on output.
    The rows of A must be sorted by increasing first column index.
    The matrix is assumed to be valid.

    Empty rows and columns have already been removed.

*/

#include "umf_internal.h"
#include "umf_apply_order.h"
#include "umf_malloc.h"
#include "umf_free.h"
#include "umf_fsize.h"

/* ========================================================================== */

GLOBAL Int UMF_analyze
(
    Int n_row,		/* A is n_row-by-n_col */
    Int n_col,
    Int Ai [ ],		/* Ai [Ap [0]..Ap[n_row]-1]: column indices */
			/* destroyed on output.  Note that this is NOT the */
			/* user's Ai that was passed to UMFPACK_*symbolic */
			/* size of Ai, Ap [n_row] = anzmax >= anz + n_col */
			/* Ap [0] must be => n_col.  The space to the */
			/* front of Ai is used as workspace. */

    Int Ap [ ],		/* of size MAX (n_row, n_col) + 1 */
			/* Ap [0..n_row]: row pointers */
			/* Row i is in Ai [Ap [i] ... Ap [i+1]-1] */

			/* rows must have smallest col index first, or be */
			/* in sorted form.  Used as workspace of size n_col */
			/* and destroyed. */

			/* Note that this is NOT the */
			/* user's Ap that was passed to UMFPACK_*symbolic */

    Int Up [ ],		/* workspace of size n_col, and output column perm.
			 * for column etree postorder. */

    Int fixQ,

    /* temporary workspaces: */
    Int W [ ],		/* W [0..n_col-1] */
    Int Link [ ],	/* Link [0..n_col-1] */

    /* output: information about each frontal matrix: */
    Int Front_ncols [ ],	/* size n_col */
    Int Front_nrows [ ],	/* of size n_col */
    Int Front_npivcol [ ],	/* of size n_col */
    Int Front_parent [ ],	/* of size n_col */
    Int *nfr_out,

    Int *p_ncompactions		/* number of compactions in UMF_analyze */
)
{
    /* ====================================================================== */
    /* ==== local variables ================================================= */
    /* ====================================================================== */

    Int j, j3, col, k, row, parent, j2, pdest, p, p2, thickness, npivots, nfr,
	i, *Winv, kk, npiv, jnext, krow, knext, pfirst, jlast, ncompactions,
	*Front_stack, *Front_order, *Front_child, *Front_sibling,
	Wflag, npivcol, fallrows, fallcols, fpiv, frows, fcols, *Front_size ;

    nfr = 0 ;
    DEBUG0 (("UMF_analyze: anzmax "ID" anrow "ID" ancol "ID"\n",
	Ap [n_row], n_row, n_col)) ;

    /* ====================================================================== */
    /* ==== initializations ================================================= */
    /* ====================================================================== */

#pragma ivdep
    for (j = 0 ; j < n_col ; j++)
    {
	Link [j] = EMPTY ;
	W [j] = EMPTY ;
	Up [j] = EMPTY ;

	/* Frontal matrix data structure: */
	Front_npivcol [j] = 0 ;		/* number of pivot columns */
	Front_nrows [j] = 0 ;		/* number of rows, incl. pivot rows */
	Front_ncols [j] = 0 ;		/* number of cols, incl. pivot cols */
	Front_parent [j] = EMPTY ;	/* parent front */
	/* Note that only non-pivotal columns are stored in a front (a "row" */
	/* of U) during elimination. */
    }

    /* the rows must be sorted by increasing min col */
    krow = 0 ;
    pfirst = Ap [0] ;
    jlast = EMPTY ;
    jnext = EMPTY ;
    Wflag = 0 ;

    /* this test requires the size of Ai to be >= n_col + nz */
    ASSERT (pfirst >= n_col) ;	/* Ai must be large enough */

    /* pdest points to the first free space in Ai */
    pdest = 0 ;
    ncompactions = 0 ;

    /* ====================================================================== */
    /* === compute symbolic LL' factorization (unsorted) ==================== */
    /* ====================================================================== */

    for (j = 0 ; j < n_col ; j = jnext)
    {
	DEBUG1 (("\n\n============Front "ID" starting. nfr = "ID"\n", j, nfr)) ;

	/* ================================================================== */
	/* === garbage collection =========================================== */
	/* ================================================================== */

	if (pdest + (n_col-j) > pfirst)
	{
	    /* we might run out ... compact the rows of U */

#ifndef NDEBUG
	    DEBUG0 (("UMF_analyze COMPACTION, j="ID" pfirst="ID"\n",
		j, pfirst)) ;
	    for (row = 0 ; row < j ; row++)
	    {
		if (Up [row] != EMPTY)
		{
		    /* this is a live row of U */
		    DEBUG1 (("Live row: "ID" cols: ", row)) ;
		    p = Up [row] ;
		    ASSERT (Front_ncols [row] > Front_npivcol [row]) ;
		    p2 = p + (Front_ncols [row] - Front_npivcol [row]) ;
		    for ( ; p < p2 ; p++)
		    {
			DEBUG1 ((ID, Ai [p])) ;
			ASSERT (p < pfirst) ;
			ASSERT (Ai [p] > row && Ai [p] < n_col) ;
		    }
		    DEBUG1 (("\n")) ;
		}
	    }
	    DEBUG1 (("\nStarting to compact:\n")) ;
#endif

	    pdest = 0 ;
	    ncompactions++ ;
	    for (row = 0 ; row < j ; row++)
	    {
		if (Up [row] != EMPTY)
		{
		    /* this is a live row of U */
		    DEBUG1 (("Live row: "ID" cols: ", row)) ;
		    ASSERT (row < n_col) ;
		    p = Up [row] ;
		    ASSERT (Front_ncols [row] > Front_npivcol [row]) ;
		    p2 = p + (Front_ncols [row] - Front_npivcol [row]) ;
		    Up [row] = pdest ;
		    for ( ; p < p2 ; p++)
		    {
			DEBUG1 ((ID, Ai [p])) ;
			ASSERT (p < pfirst) ;
			ASSERT (Ai [p] > row && Ai [p] < n_col) ;
			Ai [pdest++] = Ai [p] ;
			ASSERT (pdest <= pfirst) ;
		    }
		    DEBUG1 (("\n")) ;
		}
	    }

#ifndef NDEBUG
	    DEBUG1 (("\nAFTER COMPACTION, j="ID" pfirst="ID"\n", j, pfirst)) ;
	    for (row = 0 ; row < j ; row++)
	    {
		if (Up [row] != EMPTY)
		{
		    /* this is a live row of U */
		    DEBUG1 (("Live row: "ID" cols: ", row)) ;
		    p = Up [row] ;
		    ASSERT (Front_ncols [row] > Front_npivcol [row]) ;
		    p2 = p + (Front_ncols [row] - Front_npivcol [row]) ;
		    for ( ; p < p2 ; p++)
		    {
			DEBUG1 ((ID, Ai [p])) ;
			ASSERT (p < pfirst) ;
			ASSERT (Ai [p] > row && Ai [p] < n_col) ;
		    }
		    DEBUG1 (("\n")) ;
		}
	    }
#endif

	}

	if (pdest + (n_col-j) > pfirst)
	{
	    /* :: out of memory in umf_analyze :: */
	    /* it can't happen, if pfirst >= n_col */
	    return (FALSE) ;	/* internal error! */
	}

	/* ------------------------------------------------------------------ */
	/* is the last front a child of this one? */
	/* ------------------------------------------------------------------ */

	if (jlast != EMPTY && Link [j] == jlast)
	{
	    /* yes - create row j by appending to jlast */
	    DEBUG1 (("GOT:last front is child of this one: j "ID" jlast "ID"\n",
		j, jlast)) ;
	    ASSERT (jlast >= 0 && jlast < j) ;

	    Up [j] = Up [jlast] ;
	    Up [jlast] = EMPTY ;

	    /* find the parent, delete column j, and update W */
	    parent = n_col ;
	    for (p = Up [j] ; p < pdest ; )
	    {
		j3 = Ai [p] ;
		DEBUG1 (("Initial row of U: col "ID" ", j3)) ;
		ASSERT (j3 >= 0 && j3 < n_col) ;
		DEBUG1 (("W: "ID" \n", W [j3])) ;
		ASSERT (W [j3] == Wflag) ;
		if (j == j3)
		{
		    DEBUG1 (("Found column j at p = "ID"\n", p)) ;
		    Ai [p] = Ai [--pdest] ;
		}
		else
		{
		    if (j3 < parent)
		    {
			parent = j3 ;
		    }
		    p++ ;
		}
	    }

	    /* delete jlast from the link list of j */
	    Link [j] = Link [jlast] ;

	    ASSERT (Front_nrows [jlast] > Front_npivcol [jlast]) ;
	    thickness = (Front_nrows [jlast] - Front_npivcol [jlast]) ;
	    DEBUG1 (("initial thickness: "ID"\n", thickness)) ;

	}
	else
	{
	    Up [j] = pdest ;
	    parent = n_col ;
	    /* thickness: number of (nonpivotal) rows in frontal matrix j */
	    thickness = 0 ;
	    Wflag = j ;
	}

	/* ================================================================== */
	/* === compute row j of A*A' ======================================== */
	/* ================================================================== */

	/* ------------------------------------------------------------------ */
	/* flag the diagonal entry in row U, but do not add to pattern */
	/* ------------------------------------------------------------------ */

	ASSERT (pdest <= pfirst) ;
	W [j] = Wflag ;

	DEBUG1 (("\nComputing row "ID" of A'*A\n", j)) ;
	DEBUG2 (("	col: "ID" (diagonal)\n", j)) ;

	/* ------------------------------------------------------------------ */
	/* find the rows the contribute to this column j */
	/* ------------------------------------------------------------------ */

	jnext = n_col ;
	for (knext = krow ; knext < n_row ; knext++)
	{
	    ASSERT (Ap [knext] < Ap [knext+1]) ;
	    ASSERT (Ap [knext] >= pfirst && Ap [knext] <= Ap [n_row]) ;
	    jnext = Ai [Ap [knext]] ;
	    ASSERT (jnext >= j) ;
	    if (jnext != j)
	    {
		break ;
	    }
	}

	/* rows krow ... knext-1 all have first column index of j */
	/* (or are empty) */

	/* row knext has first column index of jnext */
	/* if knext = n_row, then jnext is n_col */
	if (knext == n_row)
	{
	    jnext = n_col ;
	}

	ASSERT (jnext > j) ;
	ASSERT (jnext <= n_col) ;

	/* ------------------------------------------------------------------ */
	/* for each nonzero A (k,j) in column j of A do: */
	/* ------------------------------------------------------------------ */

	for (k = krow ; k < knext ; k++)
	{
	    p = Ap [k] ;
	    p2 = Ap [k+1] ;
	    ASSERT (p < p2) ;

	    /* merge row k of A into W */
	    DEBUG2 (("	---- A row "ID" ", k)) ;
	    ASSERT (k >= 0 && k < n_row) ;
	    ASSERT (Ai [p] == j) ;
	    DEBUG2 (("  p "ID" p2 "ID"\n        cols:", p, p2)) ;
	    ASSERT (p  >= pfirst && p  < Ap [n_row]) ;
	    ASSERT (p2 >  pfirst && p2 <= Ap [n_row]) ;
	    for ( ; p < p2 ; p++)
	    {
		/* add to pattern if seen for the first time */
		col = Ai [p] ;
		ASSERT (col >= j && col < n_col) ;
		DEBUG3 ((" "ID, col)) ;
		if (W [col] != Wflag)
		{
		    Ai [pdest++] = col ;
		    ASSERT (pdest <= pfirst) ;
		    /* flag this column has having been seen for row j */
		    W [col] = Wflag ;
		    if (col < parent)
		    {
			parent = col ;
		    }
		}
	    }
	    DEBUG2 (("\n")) ;
	    thickness++ ;
	}

#ifndef NDEBUG
	DEBUG3 (("\nRow "ID" of A'A:\n", j)) ;
	for (p = Up [j] ; p < pdest ; p++)
	{
	    DEBUG3 ((" "ID, Ai [p])) ;
	}
	DEBUG3 (("\n")) ;
#endif

	/* ------------------------------------------------------------------ */
	/* delete rows up to but not including knext */
	/* ------------------------------------------------------------------ */

	krow = knext ;
	pfirst = Ap [knext] ;

	/* we can now use Ai [0..pfirst-1] as workspace for rows of U */

	/* ================================================================== */
	/* === compute jth row of U ========================================= */
	/* ================================================================== */

	/* for each nonzero U (k,j) in column j of U (1:j-1,:) do */
	for (k = Link [j] ; k != EMPTY ; k = Link [k])
	{
	    /* merge row k of U into W */
	    DEBUG2 (("	---- U row "ID, k)) ;
	    ASSERT (k >= 0 && k < n_col) ;
	    ASSERT (Up [k] != EMPTY) ;
	    p = Up [k] ;
	    ASSERT (Front_ncols [k] > Front_npivcol [k]) ;
	    p2 = p + (Front_ncols [k] - Front_npivcol [k]) ;
	    DEBUG2 (("  p "ID" p2 "ID"\n        cols:", p, p2)) ;
	    ASSERT (p <= pfirst) ;
	    ASSERT (p2 <= pfirst) ;
	    for ( ; p < p2 ; p++)
	    {
		/* add to pattern if seen for the first time */
		col = Ai [p] ;
		ASSERT (col >= j && col < n_col) ;
		DEBUG3 ((" "ID, col)) ;
		if (W [col] != Wflag)
		{
		    Ai [pdest++] = col ;
		    ASSERT (pdest <= pfirst) ;
		    /* flag this col has having been seen for row j */
		    W [col] = Wflag ;
		    if (col < parent)
		    {
			parent = col ;
		    }
		}
	    }
	    DEBUG2 (("\n")) ;

	    /* mark the row k as deleted */
	    Up [k] = EMPTY ;

	    ASSERT (Front_nrows [k] > Front_npivcol [k]) ;
	    thickness += (Front_nrows [k] - Front_npivcol [k]) ;
	    ASSERT (Front_parent [k] == j) ;
	}

#ifndef NDEBUG
	DEBUG3 (("\nRow "ID" of U prior to supercolumn detection:\n", j));
	for (p = Up [j] ; p < pdest ; p++)
	{
	    DEBUG3 ((" "ID, Ai [p])) ;
	}
	DEBUG3 (("\n")) ;
	DEBUG1 (("thickness, prior to supercol detect: "ID"\n", thickness)) ;
#endif

	/* ================================================================== */
	/* === quicky mass elimination ====================================== */
	/* ================================================================== */

	/* this code detects some supernodes, but it might miss */
	/* some because the elimination tree (created on the fly) */
	/* is not yet post-ordered, and because the pattern of A'*A */
	/* is also computed on the fly. */

	/* j2 is incremented because the pivot columns are not stored */

	for (j2 = j+1 ; j2 < jnext ; j2++)
	{
	    ASSERT (j2 >= 0 && j2 < n_col) ;
	    if (W [j2] != Wflag || Link [j2] != EMPTY)
	    {
		break ;
	    }
	}

	/* the loop above terminated with j2 at the first non-supernode */
	DEBUG1 (("jnext = "ID"\n", jnext)) ;
	ASSERT (j2 <= jnext) ;
	jnext = j2 ;
	j2-- ;
	DEBUG1 (("j2 = "ID"\n", j2)) ;
	ASSERT (j2 < n_col) ;

	npivots = j2-j+1 ;
	DEBUG1 (("Number of pivot columns: "ID"\n", npivots)) ;

	/* rows j:j2 have the same nonzero pattern, except for columns j:j2-1 */

	if (j2 > j)
	{
	    /* supernode detected, prune the pattern of new row j */
	    ASSERT (parent == j+1) ;
	    ASSERT (j2 < n_col) ;
	    DEBUG1 (("Supernode detected, j "ID" to j2 "ID"\n", j, j2)) ;

	    parent = n_col ;
	    p2 = pdest ;
	    pdest = Up [j] ;
	    for (p = Up [j] ; p < p2 ; p++)
	    {
		col = Ai [p] ;
		ASSERT (col >= 0 && col < n_col) ;
		ASSERT (W [col] == Wflag) ;
		if (col > j2)
		{
		    /* keep this col in the pattern of the new row j */
		    Ai [pdest++] = col ;
		    if (col < parent)
		    {
			parent = col ;
		    }
		}
	    }
	}

	DEBUG1 (("Parent ["ID"] = "ID"\n", j, parent)) ;
	ASSERT (parent > j2) ;

	if (parent == n_col)
	{
	    /* this front has no parent - it is the root of a subtree */
	    parent = EMPTY ;
	}

#ifndef NDEBUG
	DEBUG3 (("\nFinal row "ID" of U after supercolumn detection:\n", j)) ;
	for (p = Up [j] ; p < pdest ; p++)
	{
	    ASSERT (Ai [p] >= 0 && Ai [p] < n_col) ;
	    DEBUG3 ((" "ID" ("ID")", Ai [p], W [Ai [p]])) ;
	    ASSERT (W [Ai [p]] == Wflag) ;
	}
	DEBUG3 (("\n")) ;
#endif

	/* ================================================================== */
	/* === frontal matrix =============================================== */
	/* ================================================================== */

	/* front has Front_npivcol [j] pivot columns */
	/* entire front is Front_nrows [j] -by- Front_ncols [j] */
	/* j is first column in the front */

	npivcol = npivots ;
	fallrows = thickness ;
	fallcols = npivots + pdest - Up [j] ;

	/* number of pivots in the front (rows and columns) */
	fpiv = MIN (npivcol, fallrows) ;

	/* size of contribution block */
	frows = fallrows - fpiv ;
	fcols = fallcols - fpiv ;

	if (frows == 0 || fcols == 0)
	{
	    /* front has no contribution block and thus needs no parent */
	    DEBUG1 (("Frontal matrix evaporation\n")) ;
	    Up [j] = EMPTY ;
	    parent = EMPTY ;
	}

	Front_npivcol [j] = npivots ;
	Front_nrows [j] = fallrows ;
	Front_ncols [j] = fallcols ;
	Front_parent [j] = parent ;
	ASSERT (npivots > 0) ;

	/* Front_parent [j] is the first column of the parent frontal matrix */

	DEBUG1 (("\n\n==== Front "ID", nfr "ID" pivot columns "ID":"ID
	    " all front: "ID"-by-"ID" Parent: "ID"\n", j, nfr, j,j+npivots-1,
	    Front_nrows [j], Front_ncols [j], Front_parent [j])) ;
	nfr++ ;

	/* ================================================================== */
	/* === prepare this row for its parent ============================== */
	/* ================================================================== */

	if (parent != EMPTY)
	{
	    Link [j] = Link [parent] ;
	    Link [parent] = j ;
	}

	ASSERT (jnext > j) ;

	jlast = j ;
    }

    /* ====================================================================== */
    /* === postorder the fronts ============================================= */
    /* ====================================================================== */

    *nfr_out = nfr ;

    Front_order = W ;	/* use W for Front_order [ */

    if (fixQ)
    {
	/* do not postorder the fronts if Q is fixed */
	DEBUG1 (("\nNo postorder (Q is fixed)\n")) ;
	k = 0 ;
	for (j = 0 ; j < n_col ; j++)
	{
	    if (Front_npivcol [j] > 0)
	    {
		Front_order [j] = k++ ;
		DEBUG1 (("Front order of j: "ID" is:"ID"\n", j,
		    Front_order [j])) ;
	    }
	    else
	    {
		Front_order [j] = EMPTY ;
	    }
	}
    }
    else
    {

	/* use Ap for Front_child and use Link for Front_sibling [ */
	Front_child = Ap ;
	Front_sibling = Link ;

	/* use Ai for Front_stack, size of Ai is >= 2*n_col */
	Front_stack = Ai ;
	Front_size = Front_stack + n_col ;

	UMF_fsize (n_col, Front_size, Front_nrows, Front_ncols,
	    Front_parent, Front_npivcol) ;

	AMD_postorder (n_col, Front_parent, Front_npivcol, Front_size,
	    Front_order, Front_child, Front_sibling, Front_stack) ;

	/* done with Front_child, Front_sibling, Front_size, and Front_stack ]*/

	/* ------------------------------------------------------------------ */
	/* construct the column permutation (return in Up) */
	/* ------------------------------------------------------------------ */

	/* Front_order [i] = k means that front i is kth front in the new order.
	 * i is in the range 0 to n_col-1, and k is in the range 0 to nfr-1 */

	/* Use Ai as workspace for Winv [ */
	Winv = Ai ;
	for (k = 0 ; k < nfr ; k++)
	{
	    Winv [k] = EMPTY ;
	}

	/* compute the inverse of Front_order, so that Winv [k] = i */
	/* if Front_order [i] = k */

	DEBUG1 (("\n\nComputing output column permutation:\n")) ;
	for (i = 0 ; i < n_col ; i++)
	{
	    k = Front_order [i] ;
	    if (k != EMPTY)
	    {
		DEBUG1 (("Front "ID" new order: "ID"\n", i, k)) ;
		ASSERT (k >= 0 && k < nfr) ;
		ASSERT (Winv [k] == EMPTY) ;
		Winv [k] = i ;
	    }
	}

	/* Use Up as output permutation */
	kk = 0 ;
	for (k = 0 ; k < nfr ; k++)
	{
	    i = Winv [k] ;
	    DEBUG1 (("Old Front "ID" New Front "ID" npivots "ID" nrows "ID
		" ncols "ID"\n",
		i, k, Front_npivcol [i], Front_nrows [i], Front_ncols [i])) ;
	    ASSERT (i >= 0 && i < n_col) ;
	    ASSERT (Front_npivcol [i] > 0) ;
	    for (npiv = 0 ; npiv < Front_npivcol [i] ; npiv++)
	    {
		Up [kk] = i + npiv ;
		DEBUG1 (("    Cperm ["ID"] = "ID"\n", kk, Up [kk])) ;
		kk++ ;
	    }
	}
	ASSERT (kk == n_col) ;

	/* Winv no longer needed ] */
    }

    /* ---------------------------------------------------------------------- */
    /* apply the postorder traversal to renumber the frontal matrices */
    /* (or pack them in same order, if fixQ) */
    /* ---------------------------------------------------------------------- */

    /* use Ai as workspace */

    UMF_apply_order (Front_npivcol, Front_order, Ai, n_col, nfr) ;
    UMF_apply_order (Front_nrows,   Front_order, Ai, n_col, nfr) ;
    UMF_apply_order (Front_ncols,   Front_order, Ai, n_col, nfr) ;
    UMF_apply_order (Front_parent,  Front_order, Ai, n_col, nfr) ;

    /* fix the parent to refer to the new numbering */
    for (i = 0 ; i < nfr ; i++)
    {
	parent = Front_parent [i] ;
	if (parent != EMPTY)
	{
	    ASSERT (parent >= 0 && parent < n_col) ;
	    ASSERT (Front_order [parent] >= 0 && Front_order [parent] < nfr) ;
	    Front_parent [i] = Front_order [parent] ;
	}
    }

    /* Front_order longer needed ] */

#ifndef NDEBUG
    DEBUG1 (("\nFinal frontal matrices:\n")) ;
    for (i = 0 ; i < nfr ; i++)
    {
	DEBUG1 (("Final front "ID": npiv "ID" nrows "ID" ncols "ID" parent "
	    ID"\n", i, Front_npivcol [i], Front_nrows [i],
	    Front_ncols [i], Front_parent [i])) ;
    }
#endif

    *p_ncompactions = ncompactions ;
    return (TRUE) ;
}
