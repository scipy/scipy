/* ========================================================================== */
/* === UMF_singletons ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Find and order the row and column singletons of a matrix A.  If there are
 * row and column singletons, the output is a row and column permutation such
 * that the matrix is in the following form:
 *
 *	x x x x x x x x x
 *	0 x x x x x x x x
 *	0 0 x x x x x x x
 *	0 0 0 x 0 0 0 0 0
 *	0 0 0 x x 0 0 0 0
 *	0 0 0 x x s s s s
 *	0 0 0 x x s s s s
 *	0 0 0 x x s s s s
 *	0 0 0 x x s s s s
 *
 * The above example has 3 column singletons (the first three columns and
 * their corresponding pivot rows) and 2 row singletons.  The singletons are
 * ordered first, because they have zero Markowitz cost.  The LU factorization
 * for these first five rows and columns is free - there is no work to do
 * (except to scale the pivot columns for the 2 row singletons), and no
 * fill-in occurs.  * The remaining * submatrix (4-by-4 in the above example)
 * has no rows or columns with degree one.  It may have empty rows or columns.
 *
 * This algorithm does not perform a full permutation to block triangular
 * form.  If there are one or more singletons, then the matrix can be
 * permuted to block triangular form, but UMFPACK does not perform the full
 * BTF permutation (see also "dmperm" in MATLAB).
 */

#include "umf_internal.h"

#ifndef NDEBUG

/* ========================================================================== */
/* === debug routines ======================================================= */
/* ========================================================================== */

/* Dump the singleton queue */

PRIVATE void dump_singletons
(
    Int head,		/* head of the queue */
    Int tail,		/* tail of the queue */
    Int Next [ ],	/* Next [i] is the next object after i */
    char *name,		/* "row" or "col" */
    Int Deg [ ],	/* Deg [i] is the degree of object i */
    Int n		/* objects are in the range 0 to n-1 */
)
{
    Int i, next, cnt ;
    DEBUG6 (("%s Singleton list: head "ID" tail "ID"\n", name, head, tail)) ;
    i = head ;
    ASSERT (head >= EMPTY && head < n) ;
    ASSERT (tail >= EMPTY && tail < n) ;
    cnt = 0 ;
    while (i != EMPTY)
    {
	DEBUG7 ((" "ID": "ID" deg: "ID"\n", cnt, i, Deg [i])) ;
	ASSERT (i >= 0 && i < n) ;
	next = Next [i] ;
	if (i == tail) ASSERT (next == EMPTY) ;
	i = next ;
	cnt++ ;
	ASSERT (cnt <= n) ;
    }
}

PRIVATE void dump_mat
(
    char *xname,
    char *yname,
    Int nx,
    Int ny,
    const Int Xp [ ],
    const Int Xi [ ],
    Int Xdeg [ ],
    Int Ydeg [ ]
)
{
    Int x, y, p, p1, p2, xdeg, do_xdeg, ydeg ;
    DEBUG6 (("\n ==== Dump %s mat:\n", xname)) ;
    for (x = 0 ; x < nx ; x++)
    {
	p1 = Xp [x] ;
	p2 = Xp [x+1] ;
	xdeg = Xdeg [x] ;
	DEBUG6 (("Dump %s "ID" p1 "ID" p2 "ID" deg "ID"\n",
	    xname, x, p1, p2, xdeg)) ;
	do_xdeg = (xdeg >= 0) ;
	for (p = p1 ; p < p2 ; p++)
	{
	    y = Xi [p] ;
	    DEBUG7 (("    %s "ID" deg: ", yname, y)) ;
	    ASSERT (y >= 0 && y < ny) ;
	    ydeg = Ydeg [y] ;
	    DEBUG7 ((ID"\n", ydeg)) ;
	    if (do_xdeg && ydeg >= 0)
	    {
		xdeg-- ;
	    }
	}
	ASSERT (IMPLIES (do_xdeg, xdeg == 0)) ;
    }
}
#endif

/* ========================================================================== */
/* === create_row_form ====================================================== */
/* ========================================================================== */

/* Create the row-form R of the column-form input matrix A.  This could be done
 * by UMF_transpose, except that Rdeg has already been computed.
 */

PRIVATE void create_row_form
(
    /* input, not modified: */
    Int n_row,		    /* A is n_row-by-n_col, nz = Ap [n_col] */
    Int n_col,
    const Int Ap [ ],	    /* Ap [0..n_col]: column pointers for A */
    const Int Ai [ ],	    /* Ai [0..nz-1]:  row indices for A */
    Int Rdeg [ ],	    /* Rdeg [0..n_row-1]: row degrees */

    /* output, not defined on input: */
    Int Rp [ ],		    /* Rp [0..n_row]: row pointers for R */
    Int Ri [ ],		    /* Ri [0..nz-1]:  column indices for R */

    /* workspace, not defined on input or output */
    Int W [ ]		    /* size n_row */
)
{
    Int row, col, p, p2 ;

    /* create the row pointers */
    Rp [0] = 0 ;
    W [0] = 0 ;
    for (row = 0 ; row < n_row ; row++)
    {
	Rp [row+1] = Rp [row] + Rdeg [row] ;
	W [row] = Rp [row] ;
    }

    /* create the indices for the row-form */
    for (col = 0 ; col < n_col ; col++)
    {
	p2 = Ap [col+1] ;
	for (p = Ap [col] ; p < p2 ; p++)
	{
	    Ri [W [Ai [p]]++] = col ;
	}
    }
}

/* ========================================================================== */
/* === order_singletons ===================================================== */
/* ========================================================================== */

PRIVATE int order_singletons	/* return new number of singletons */
(
    Int k,	    /* the number of singletons so far */
    Int head,
    Int tail,
    Int Next [ ],
    Int Xdeg [ ], Int Xperm [ ], const Int Xp [ ], const Int Xi [ ],
    Int Ydeg [ ], Int Yperm [ ], const Int Yp [ ], const Int Yi [ ]
#ifndef NDEBUG
    , char *xname, char *yname, Int nx, Int ny
#endif
)
{
    Int xpivot, x, y, ypivot, p, p2, deg ;

#ifndef NDEBUG
    Int i, k1 = k ;
    dump_singletons (head, tail, Next, xname, Xdeg, nx) ;
    dump_mat (xname, yname, nx, ny, Xp, Xi, Xdeg, Ydeg) ;
    dump_mat (yname, xname, ny, nx, Yp, Yi, Ydeg, Xdeg) ;
#endif

    while (head != EMPTY)
    {
	/* remove the singleton at the head of the queue */
	xpivot = head ;
	DEBUG1 (("------ Order %s singleton: "ID"\n", xname, xpivot)) ;
	head = Next [xpivot] ;
	if (head == EMPTY) tail = EMPTY ;

#ifndef NDEBUG
	if (k % 100 == 0) dump_singletons (head, tail, Next, xname, Xdeg, nx) ;
#endif

	ASSERT (Xdeg [xpivot] >= 0) ;
	if (Xdeg [xpivot] != 1)
	{
	    /* This row/column x is empty.  The matrix is singular.
	     * x will be ordered last in Xperm. */
	    DEBUG1 (("empty %s, after singletons removed\n", xname)) ;
	    continue ;
	}

	/* find the ypivot to match with this xpivot */
#ifndef NDEBUG
	/* there can only be one ypivot, since the degree of x is 1 */
	deg = 0 ;
	p2 = Xp [xpivot+1] ;
	for (p = Xp [xpivot] ; p < p2 ; p++)
	{
	    y = Xi [p] ;
	    DEBUG1 (("%s: "ID"\n", yname, y)) ;
	    if (Ydeg [y] >= 0)
	    {
		/* this is a live index in this xpivot vector */
		deg++ ;
	    }
	}
	ASSERT (deg == 1) ;
#endif

	ypivot = EMPTY ;
	p2 = Xp [xpivot+1] ;
	for (p = Xp [xpivot] ; p < p2 ; p++)
	{
	    y = Xi [p] ;
	    DEBUG1 (("%s: "ID"\n", yname, y)) ;
	    if (Ydeg [y] >= 0)
	    {
		/* this is a live index in this xpivot vector */
		ypivot = y ;
		break ;
	    }
	}

	DEBUG1 (("Pivot %s: "ID"\n", yname, ypivot)) ;
	ASSERT (ypivot != EMPTY) ;
	DEBUG1 (("deg "ID"\n", Ydeg [ypivot])) ;
	ASSERT (Ydeg [ypivot] >= 0) ;

	/* decrement the degrees after removing this singleton */
	DEBUG1 (("p1 "ID"\n", Yp [ypivot])) ;
	DEBUG1 (("p2 "ID"\n", Yp [ypivot+1])) ;
	p2 = Yp [ypivot+1] ;
	for (p = Yp [ypivot] ; p < p2 ; p++)
	{
	    x = Yi [p] ;
	    DEBUG1 (("    %s: "ID" deg: "ID"\n", xname, x, Xdeg [x])) ;
	    if (Xdeg [x] < 0) continue ;
	    ASSERT (Xdeg [x] > 0) ;
	    if (x == xpivot) continue ;
	    deg = --(Xdeg [x]) ;
	    ASSERT (Xdeg [x] >= 0) ;
	    if (deg == 1)
	    {
		/* this is a new singleton, put at the end of the queue */
		Next [x] = EMPTY ;
		if (head == EMPTY)
		{
		    head = x ;
		}
		else
		{
		    ASSERT (tail != EMPTY) ;
		    Next [tail] = x ;
		}
		tail = x ;
		DEBUG1 ((" New %s singleton:  "ID"\n", xname, x)) ;
#ifndef NDEBUG
		if (k % 100 == 0)
		{
		    dump_singletons (head, tail, Next, xname, Xdeg, nx) ;
		}
#endif
	    }
	}

	/* flag the xpivot and ypivot by FLIP'ing the degrees */
	Xdeg [xpivot] = FLIP (1) ;
	Ydeg [ypivot] = FLIP (Ydeg [ypivot]) ;

	/* keep track of the pivot row and column */
	Xperm [k] = xpivot ;
	Yperm [k] = ypivot ;
	k++ ;

#ifndef NDEBUG
	if (k % 1000 == 0)
	{
	    dump_mat (xname, yname, nx, ny, Xp, Xi, Xdeg, Ydeg) ;
	    dump_mat (yname, xname, ny, nx, Yp, Yi, Ydeg, Xdeg) ;
	}
#endif
    }

#ifndef NDEBUG
    DEBUGm4 (("%s singletons: k = "ID"\n", xname, k)) ;
    for (i = k1 ; i < k ; i++)
    {
	DEBUG1 (("  %s: "ID" %s: "ID"\n", xname, Xperm [i], yname, Yperm [i])) ;
    }
    ASSERT (k > 0) ;
#endif

    return (k) ;
}

/* ========================================================================== */
/* === find_any_singletons ================================================== */
/* ========================================================================== */

PRIVATE Int find_any_singletons	    /* returns # of singletons found */
(
    /* input, not modified: */
    Int n_row,
    Int n_col,
    const Int Ap [ ],	    /* size n_col+1 */
    const Int Ai [ ],	    /* size nz = Ap [n_col] */

    /* input, modified on output: */
    Int Cdeg [ ],	    /* size n_col */
    Int Rdeg [ ],	    /* size n_row */

    /* output, not defined on input: */
    Int Cperm [ ],	    /* size n_col */
    Int Rperm [ ],	    /* size n_row */
    Int *p_n1r,		    /* # of row singletons */
    Int *p_n1c,		    /* # of col singletons */

    /* workspace, not defined on input or output */
    Int Rp [ ],		    /* size n_row+1 */
    Int Ri [ ],		    /* size nz */
    Int W [ ],		    /* size n_row */
    Int Next [ ]	    /* size MAX (n_row, n_col) */
)
{
    Int n1, col, row, row_form, head, tail, n1r, n1c ;

    /* ---------------------------------------------------------------------- */
    /* eliminate column singletons */
    /* ---------------------------------------------------------------------- */

    n1 = 0 ;
    n1r = 0 ;
    n1c = 0 ;
    row_form = FALSE ;

    head = EMPTY ;
    tail = EMPTY ;
    for (col = n_col-1 ; col >= 0 ; col--)
    {
	if (Cdeg [col] == 1)
	{
	    /* put the column singleton in the queue */
	    if (head == EMPTY) tail = col ;
	    Next [col] = head ;
	    head = col ;
	    DEBUG1 (("Column singleton: "ID"\n", col)) ;
	}
    }

    if (head != EMPTY)
    {

	/* ------------------------------------------------------------------ */
	/* create the row-form of A */
	/* ------------------------------------------------------------------ */

	create_row_form (n_row, n_col, Ap, Ai, Rdeg, Rp, Ri, W) ;
	row_form = TRUE ;

	/* ------------------------------------------------------------------ */
	/* find and order the column singletons */
	/* ------------------------------------------------------------------ */

	n1 = order_singletons (0, head, tail, Next,
		Cdeg, Cperm, Ap, Ai,
		Rdeg, Rperm, Rp, Ri
#ifndef NDEBUG
		, "col", "row", n_col, n_row
#endif
		) ;
	n1c = n1 ;
    }

    /* ---------------------------------------------------------------------- */
    /* eliminate row singletons */
    /* ---------------------------------------------------------------------- */

    head = EMPTY ;
    tail = EMPTY ;
    for (row = n_row-1 ; row >= 0 ; row--)
    {
	if (Rdeg [row] == 1)
	{
	    /* put the row singleton in the queue */
	    if (head == EMPTY) tail = row ;
	    Next [row] = head ;
	    head = row ;
	    DEBUG1 (("Row singleton: "ID"\n", row)) ;
	}
    }

    if (head != EMPTY)
    {

	/* ------------------------------------------------------------------ */
	/* create the row-form of A, if not already created */
	/* ------------------------------------------------------------------ */

	if (!row_form)
	{
	    create_row_form (n_row, n_col, Ap, Ai, Rdeg, Rp, Ri, W) ;
	}

	/* ------------------------------------------------------------------ */
	/* find and order the row singletons */
	/* ------------------------------------------------------------------ */

	n1 = order_singletons (n1, head, tail, Next,
		Rdeg, Rperm, Rp, Ri,
		Cdeg, Cperm, Ap, Ai
#ifndef NDEBUG
		, "row", "col", n_row, n_col
#endif
		) ;
	n1r = n1 - n1c ;
    }

    DEBUG0 (("n1 "ID"\n", n1)) ;
    *p_n1r = n1r ;
    *p_n1c = n1c ;
    return (n1) ;
}

/* ========================================================================== */
/* === find_user_singletons ================================================= */
/* ========================================================================== */

PRIVATE Int find_user_singletons	/* returns # singletons found */
(
    /* input, not modified: */
    Int n_row,
    Int n_col,
    const Int Ap [ ],	    /* size n_col+1 */
    const Int Ai [ ],	    /* size nz = Ap [n_col] */
    const Int Quser [ ],    /* size n_col if present */

    /* input, modified on output: */
    Int Cdeg [ ],	    /* size n_col */
    Int Rdeg [ ],	    /* size n_row */

    /* output, not defined on input */
    Int Cperm [ ],	    /* size n_col */
    Int Rperm [ ],	    /* size n_row */
    Int *p_n1r,		    /* # of row singletons */
    Int *p_n1c,		    /* # of col singletons */

    /* workspace, not defined on input or output */
    Int Rp [ ],		    /* size n_row+1 */
    Int Ri [ ],		    /* size nz */
    Int W [ ]		    /* size n_row */
)
{
    Int n1, col, row, p, p2, pivcol, pivrow, found, k, n1r, n1c ;

    n1 = 0 ;
    n1r = 0 ;
    n1c = 0 ;
    *p_n1r = 0 ;
    *p_n1c = 0 ;

    /* find singletons in the user column permutation, Quser */
    pivcol = Quser [0] ;
    found = (Cdeg [pivcol] == 1) ;
    DEBUG0 (("Is first col: "ID" a col singleton?: "ID"\n", pivcol, found)) ;
    if (!found)
    {
	/* the first column is not a column singleton, check for a row
	 * singleton in the first column. */
	for (p = Ap [pivcol] ; p < Ap [pivcol+1] ; p++)
	{
	    if (Rdeg [Ai [p]] == 1)
	    {
		DEBUG0 (("Row singleton in first col: "ID" row: "ID"\n",
		    pivcol, Ai [p])) ;
		found = TRUE ;
		break ;
	    }
	}
    }

    if (!found)
    {
	/* no singletons in the leading part of A (:,Quser) */
	return (0) ;
    }

    /* there is at least one row or column singleton.  Look for more. */
    create_row_form (n_row, n_col, Ap, Ai, Rdeg, Rp, Ri, W) ;

    n1 = 0 ;

    for (k = 0 ; k < n_col ; k++)
    {
	pivcol = Quser [k] ;
	pivrow = EMPTY ;

	/* ------------------------------------------------------------------ */
	/* check if col is a column singleton, or contains a row singleton */
	/* ------------------------------------------------------------------ */

	found = (Cdeg [pivcol] == 1) ;

	if (found)
	{

	    /* -------------------------------------------------------------- */
	    /* pivcol is a column singleton */
	    /* -------------------------------------------------------------- */

	    DEBUG0 (("Found a col singleton: k "ID" pivcol "ID"\n", k, pivcol));

	    /* find the pivrow to match with this pivcol */
#ifndef NDEBUG
	    /* there can only be one pivrow, since the degree of pivcol is 1 */
	    {
		Int deg = 0 ;
		p2 = Ap [pivcol+1] ;
		for (p = Ap [pivcol] ; p < p2 ; p++)
		{
		    row = Ai [p] ;
		    DEBUG1 (("row: "ID"\n", row)) ;
		    if (Rdeg [row] >= 0)
		    {
			/* this is a live index in this column vector */
			deg++ ;
		    }
		}
		ASSERT (deg == 1) ;
	    }
#endif

	    p2 = Ap [pivcol+1] ;
	    for (p = Ap [pivcol] ; p < p2 ; p++)
	    {
		row = Ai [p] ;
		DEBUG1 (("row: "ID"\n", row)) ;
		if (Rdeg [row] >= 0)
		{
		    /* this is a live index in this pivcol vector */
		    pivrow = row ;
		    break ;
		}
	    }

	    DEBUG1 (("Pivot row: "ID"\n", pivrow)) ;
	    ASSERT (pivrow != EMPTY) ;
	    DEBUG1 (("deg "ID"\n", Rdeg [pivrow])) ;
	    ASSERT (Rdeg [pivrow] >= 0) ;

	    /* decrement the degrees after removing this col singleton */
	    DEBUG1 (("p1 "ID"\n", Rp [pivrow])) ;
	    DEBUG1 (("p2 "ID"\n", Rp [pivrow+1])) ;
	    p2 = Rp [pivrow+1] ;
	    for (p = Rp [pivrow] ; p < p2 ; p++)
	    {
		col = Ri [p] ;
		DEBUG1 (("    col: "ID" deg: "ID"\n", col, Cdeg [col])) ;
		if (Cdeg [col] < 0) continue ;
		ASSERT (Cdeg [col] > 0) ;
		Cdeg [col]-- ;
		ASSERT (Cdeg [col] >= 0) ;
	    }

	    /* flag the pivcol and pivrow by FLIP'ing the degrees */
	    Cdeg [pivcol] = FLIP (1) ;
	    Rdeg [pivrow] = FLIP (Rdeg [pivrow]) ;
	    n1c++ ;

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* pivcol may contain a row singleton */
	    /* -------------------------------------------------------------- */

	    p2 = Ap [pivcol+1] ;
	    for (p = Ap [pivcol] ; p < p2 ; p++)
	    {
		pivrow = Ai [p] ;
		if (Rdeg [pivrow] == 1)
		{
		    DEBUG0 (("Row singleton in pivcol: "ID" row: "ID"\n",
			pivcol, pivrow)) ;
		    found = TRUE ;
		    break ;
		}
	    }

	    if (!found)
	    {
		DEBUG0 (("End of user singletons\n")) ;
		break ;
	    }

#ifndef NDEBUG
	    /* there can only be one pivrow, since the degree of pivcol is 1 */
	    {
		Int deg = 0 ;
		p2 = Rp [pivrow+1] ;
		for (p = Rp [pivrow] ; p < p2 ; p++)
		{
		    col = Ri [p] ;
		    DEBUG1 (("col: "ID" cdeg::: "ID"\n", col, Cdeg [col])) ;
		    if (Cdeg [col] >= 0)
		    {
			/* this is a live index in this column vector */
			ASSERT (col == pivcol) ;
			deg++ ;
		    }
		}
		ASSERT (deg == 1) ;
	    }
#endif

	    DEBUG1 (("Pivot row: "ID"\n", pivrow)) ;
	    DEBUG1 (("pivcol deg "ID"\n", Cdeg [pivcol])) ;
	    ASSERT (Cdeg [pivcol] > 1) ;

	    /* decrement the degrees after removing this row singleton */
	    DEBUG1 (("p1 "ID"\n", Ap [pivcol])) ;
	    DEBUG1 (("p2 "ID"\n", Ap [pivcol+1])) ;
	    p2 = Ap [pivcol+1] ;
	    for (p = Ap [pivcol] ; p < p2 ; p++)
	    {
		row = Ai [p] ;
		DEBUG1 (("    row: "ID" deg: "ID"\n", row, Rdeg [row])) ;
		if (Rdeg [row] < 0) continue ;
		ASSERT (Rdeg [row] > 0) ;
		Rdeg [row]-- ;
		ASSERT (Rdeg [row] >= 0) ;
	    }

	    /* flag the pivcol and pivrow by FLIP'ing the degrees */
	    Cdeg [pivcol] = FLIP (Cdeg [pivcol]) ;
	    Rdeg [pivrow] = FLIP (1) ;
	    n1r++ ;
	}

	/* keep track of the pivot row and column */
	Cperm [k] = pivcol ;
	Rperm [k] = pivrow ;
	n1++ ;

#ifndef NDEBUG
	dump_mat ("col", "row", n_col, n_row, Ap, Ai, Cdeg, Rdeg) ;
	dump_mat ("row", "col", n_row, n_col, Rp, Ri, Rdeg, Cdeg) ;
#endif

    }

    DEBUGm4 (("User singletons found: "ID"\n", n1)) ;
    ASSERT (n1 > 0) ;

    *p_n1r = n1r ;
    *p_n1c = n1c ;
    return (n1) ;
}

/* ========================================================================== */
/* === finish_permutation =================================================== */
/* ========================================================================== */

/* Complete the permutation for the pruned submatrix.  The singletons are
 * already ordered, but remove their flags.  Place rows/columns that are empty
 * in the pruned submatrix at the end of the output permutation.  This can only
 * occur if the matrix is singular.
 */

PRIVATE Int finish_permutation
(
    Int n1,
    Int nx,
    Int Xdeg [ ],
    const Int Xuser [ ],
    Int Xperm [ ],
    Int *p_max_deg
)
{
    Int nempty, x, deg, s, max_deg, k ;
    nempty = 0 ;
    s = n1 ;
    max_deg = 0 ;
    DEBUG0 (("n1 "ID" nempty "ID"\n", n1, nempty)) ;
    for (k = 0 ; k < nx ; k++)
    {
	x = (Xuser != (Int *) NULL) ? Xuser [k] : k ;
	DEBUG0 (("finish perm k "ID" x "ID" nx "ID"\n", k, x, nx)) ;
	deg = Xdeg [x] ;
	if (deg == 0)
	{
	    /* this row/col is empty in the pruned submatrix */
	    ASSERT (s < nx - nempty) ;
	    nempty++ ;
	    Xperm [nx - nempty] = x ;
	}
	else if (deg > 0)
	{
	    /* this row/col is nonempty in the pruned submatrix */
	    ASSERT (s < nx - nempty) ;
	    Xperm [s++] = x ;
	    max_deg = MAX (max_deg, deg) ;
	}
	else
	{
	    /* This is a singleton row/column - it is already ordered.
	     * Just clear the flag. */
	    Xdeg [x] = FLIP (deg) ;
	}
    }
    ASSERT (s == nx - nempty) ;
    *p_max_deg = max_deg ;
    return (nempty) ;
}

/* ========================================================================== */
/* === UMF_singletons ======================================================= */
/* ========================================================================== */

GLOBAL Int UMF_singletons
(

    /* input, not modified: */
    Int n_row,
    Int n_col,
    const Int Ap [ ],	    /* size n_col+1 */
    const Int Ai [ ],	    /* size nz = Ap [n_col] */
    const Int Quser [ ],    /* size n_col if present */

    /* output, not defined on input: */
    Int Cdeg [ ],	/* size n_col */
    Int Cperm [ ],	/* size n_col */
    Int Rdeg [ ],	/* size n_row */
    Int Rperm [ ],	/* size n_row */
    Int InvRperm [ ],	/* size n_row, the inverse of Rperm */
    Int *p_n1,		/* # of col and row singletons */
    Int *p_n1c,		/* # of col singletons */
    Int *p_n1r,		/* # of row singletons */
    Int *p_nempty_col,	/* # of empty columns in pruned submatrix */
    Int *p_nempty_row,	/* # of empty columns in pruned submatrix */
    Int *p_is_sym,	/* TRUE if pruned submatrix is square and has been
			 * symmetrically permuted by Cperm and Rperm */
    Int *p_max_rdeg,	/* maximum Rdeg in pruned submatrix */

    /* workspace, not defined on input or output */
    Int Rp [ ],		/* size n_row+1 */
    Int Ri [ ],		/* size nz */
    Int W [ ],		/* size n_row */
    Int Next [ ]	/* size MAX (n_row, n_col) */
)
{
    Int n1, s, col, row, p, p1, p2, cdeg, last_row, is_sym, k,
	nempty_row, nempty_col, max_cdeg, max_rdeg, n1c, n1r ;

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    UMF_dump_start ( ) ;
    DEBUGm4 (("Starting umf_singletons\n")) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* scan the columns, check for errors and count row degrees */
    /* ---------------------------------------------------------------------- */

    if (Ap [0] != 0 || Ap [n_col] < 0)
    {
	return (UMFPACK_ERROR_invalid_matrix) ;
    }
    for (row = 0 ; row < n_row ; row++)
    {
	Rdeg [row] = 0 ;
    }
    for (col = 0 ; col < n_col ; col++)
    {
	p1 = Ap [col] ;
	p2 = Ap [col+1] ;
	cdeg = p2 - p1 ;
	if (cdeg < 0)
	{
	    return (UMFPACK_ERROR_invalid_matrix) ;
	}
	last_row = EMPTY ;
	for (p = p1 ; p < p2 ; p++)
	{
	    row = Ai [p] ;
	    if (row <= last_row || row >= n_row)
	    {
		return (UMFPACK_ERROR_invalid_matrix) ;
	    }
	    Rdeg [row]++ ;
	    last_row = row ;
	}
	Cdeg [col] = cdeg ;
    }

    /* ---------------------------------------------------------------------- */
    /* find singletons */
    /* ---------------------------------------------------------------------- */

    if (Quser != (Int *) NULL)
    {
	/* look for singletons, but respect the user's input permutation */
	n1 = find_user_singletons (n_row, n_col, Ap, Ai, Quser,
		Cdeg, Rdeg, Cperm, Rperm, &n1r, &n1c, Rp, Ri, W) ;
    }
    else
    {
	/* look for singletons anywhere */
	n1 = find_any_singletons (n_row, n_col, Ap, Ai,
		Cdeg, Rdeg, Cperm, Rperm, &n1r, &n1c, Rp, Ri, W, Next) ;
    }

    /* ---------------------------------------------------------------------- */
    /* eliminate empty columns and complete the column permutation */
    /* ---------------------------------------------------------------------- */

    nempty_col = finish_permutation (n1, n_col, Cdeg, Quser, Cperm, &max_cdeg) ;

    /* ---------------------------------------------------------------------- */
    /* eliminate empty rows and complete the row permutation */
    /* ---------------------------------------------------------------------- */

    nempty_row = finish_permutation (n1, n_row, Rdeg, (Int *) NULL, Rperm,
	&max_rdeg) ;

    /* ---------------------------------------------------------------------- */
    /* compute the inverse of Rperm */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < n_row ; k++)
    {
	ASSERT (Rperm [k] >= 0 && Rperm [k] < n_row) ;
	InvRperm [Rperm [k]] = k ;
    }

    /* ---------------------------------------------------------------------- */
    /* see if pruned submatrix is square and has been symmetrically permuted */
    /* ---------------------------------------------------------------------- */

    if (n_row == n_col && nempty_row == nempty_col)
    {
	/* is_sym is true if the submatrix is square, and
	 * Rperm [n1..n_row-nempty_row-1] = Cperm [n1..n_col-nempty_col-1] */
	is_sym = TRUE ;
	for (s = n1 ; s < n_col - nempty_col ; s++)
	{
	    if (Cperm [s] != Rperm [s])
	    {
		is_sym = FALSE ;
		break ;
	    }
	}
    }
    else
    {
	is_sym = FALSE ;
    }
    DEBUGm4 (("Submatrix square and symmetrically permuted? "ID"\n", is_sym)) ;
    DEBUGm4 (("singletons "ID" row "ID" col "ID"\n", n1, n1r, n1c)) ;
    DEBUGm4 (("Empty cols "ID" rows "ID"\n", nempty_col, nempty_row)) ;
    *p_n1 = n1 ;
    *p_n1r = n1r ;
    *p_n1c = n1c ;
    *p_is_sym = is_sym ;
    *p_nempty_col = nempty_col ;
    *p_nempty_row = nempty_row ;
    *p_max_rdeg = max_rdeg ;
    return (UMFPACK_OK) ;
}
