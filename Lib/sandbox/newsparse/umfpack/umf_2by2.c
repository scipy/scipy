/* ========================================================================== */
/* === UMF_2by2 ============================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*  Not user-callable.  Computes a row permutation P so that A (P,:) has a
 *  mostly zero-free diagonal, with large entries on the diagonal.  It does this
 *  by swapping pairs of rows.  Once a row is swapped it is not swapped again.
 *  This is a "cheap" assignment, not a complete max. transversal or
 *  bi-partite matching.  It is only a partial matching.  For most matrices
 *  for which this algorithm is used, however, the matching is complete (in
 *  UMFPACK this algorithm is used for matrices with roughly symmetric pattern,
 *  and these matrices typically have a mostly-zero-free diagonal to begin with.
 *  This algorithm is not meant to be used on arbitrary unsymmetric matrices
 *  (for those matrices, UMFPACK uses its unsymmetric strategy and does not
 *  use this algorithm).
 *
 *  Even if incomplete, the matching is usually good enough for UMFPACK's
 *  symmetric strategy, which can easily pivot off the diagonal during numerical
 *  factorization if it finds a weak diagonal entry.
 *
 *  The algorithms works as follows.  First, row scaling factors are computed,
 *  and weak diagonal entries are found.  A weak entry is a value A(k,k) whose
 *  absolute value is < tol * max (abs (A (:,k))).  For each weak diagonal k in
 *  increasing order of degree in A+A', the algorithm finds an index j such
 *  that A (k,j) and A (j,k) are "large" (greater than or equal to tol times
 *  the largest magnitude in their columns).  Row j must also not have already
 *  been swapped.  Rows j and k are then swapped.  If we come to a diagonal k
 *  that has already been swapped, then it is not modified.  This case occurs
 *  for "oxo" pivots:
 *
 *    k j
 *  k o x
 *  j x o
 *
 *  which are swapped once to obtain
 *
 *    k j
 *  j x o
 *  k o x
 *
 *  These two rows are then not modified any further (A (j,j) was weak, but
 *  after one swap the permuted the jth diagonal entry is strong.
 *
 *  This algorithm only works on square matrices (real, complex, or pattern-
 *  only).  The numerical values are optional.  If not present, each entry is
 *  treated as numerically acceptable (tol is ignored), and the algorithm
 *  operates by just using the pattern, not the values.  Each column of the
 *  input matrix A must be sorted, with no duplicate entries.  The matrix A
 *  can be optionally scaled prior to the numerical test.  The matrix A (:,P)
 *  has the same diagonal entries as A (:,P), except in different order.  So
 *  the output permutation P can also be used to swap the columns of A.
 */

#include "umf_internal.h"

#ifndef NDEBUG
#include "umf_is_permutation.h"
#endif

/* x is "weak" if it is less than ctol.  If x or ctol are NaN, then define
 * x as not "weak".  This is a rather arbitrary choice, made to simplify the
 * computation.  On all but a PC with Microsoft C/C++, this test becomes
 * ((x) - ctol < 0). */
#define WEAK(x,ctol) (SCALAR_IS_LTZERO ((x)-(ctol)))

/* For flag value in Next [col] */
#define IS_WEAK -2

/* ========================================================================== */
/* === two_by_two =========================================================== */
/* ========================================================================== */

PRIVATE Int two_by_two	    /* returns # unmatched weak diagonals */
(
    /* input, not modified */
    Int n2,		/* C is n2-by-n2 */
    Int Cp [ ],		/* size n2+1, column pointers for C */
    Int Ci [ ],		/* size snz = Cp [n2], row indices for C */
    Int Degree [ ],	/* Degree [i] = degree of row i of C+C' */

    /* input, not defined on output */
    Int Next [ ],	/* Next [k] == IS_WEAK if k is a weak diagonal */
    Int Ri [ ],		/* Ri [i] is the length of row i in C */

    /* output, not defined on input */
    Int P [ ],

    /* workspace, not defined on input or output */
    Int Rp [ ],
    Int Head [ ]
)
{
    Int deg, newcol, row, col, p, p2, unmatched, k, j, j2, j_best, best, jdiff,
	jdiff_best, jdeg, jdeg_best, cp, cp1, cp2, rp, rp1, rp2, maxdeg,
	mindeg ;

    /* ---------------------------------------------------------------------- */
    /* place weak diagonals in the degree lists */
    /* ---------------------------------------------------------------------- */

    for (deg = 0 ; deg < n2 ; deg++)
    {
	Head [deg] = EMPTY ;
    }

    maxdeg = 0 ;
    mindeg = Int_MAX ;
    for (newcol = n2-1 ; newcol >= 0 ; newcol--)
    {
	if (Next [newcol] == IS_WEAK)
	{
	    /* add this column to the list of weak nodes */
	    DEBUGm1 (("    newcol "ID" has a weak diagonal deg "ID"\n",
		newcol, deg)) ;
	    deg = Degree [newcol] ;
	    ASSERT (deg >= 0 && deg < n2) ;
	    Next [newcol] = Head [deg] ;
	    Head [deg] = newcol ;
	    maxdeg = MAX (maxdeg, deg) ;
	    mindeg = MIN (mindeg, deg) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* construct R = C' (C = strong entries in pruned submatrix) */
    /* ---------------------------------------------------------------------- */

    /* Ri [0..n2-1] is the length of each row of R */
    /* use P as temporary pointer into the row form of R [ */
    Rp [0] = 0 ;
    for (row = 0 ; row < n2 ; row++)
    {
	Rp [row+1] = Rp [row] + Ri [row] ;
	P [row] = Rp [row] ;
    }
    /* Ri no longer needed for row counts */

    /* all entries in C are strong */
    for (col = 0 ; col < n2 ; col++)
    {
	p2 = Cp [col+1] ;
	for (p = Cp [col] ; p < p2 ; p++)
	{
	    /* place the column index in row = Ci [p] */
	    Ri [P [Ci [p]]++] = col ;
	}
    }

    /* contents of P no longer needed ] */

#ifndef NDEBUG
    DEBUG0 (("==================R: row form of strong entries in A:\n")) ;
    UMF_dump_col_matrix ((double *) NULL,
#ifdef COMPLEX
	    (double *) NULL,
#endif
	    Ri, Rp, n2, n2, Rp [n2]) ;
#endif
    ASSERT (AMD_valid (n2, n2, Rp, Ri)) ;

    /* ---------------------------------------------------------------------- */
    /* for each weak diagonal, find a pair of strong off-diagonal entries */
    /* ---------------------------------------------------------------------- */

    for (row = 0 ; row < n2 ; row++)
    {
	P [row] = EMPTY ;
    }

    unmatched = 0 ;
    best = EMPTY ;
    jdiff = EMPTY ;
    jdeg = EMPTY ;

    for (deg = mindeg ; deg <= maxdeg ; deg++)
    {
	/* find the next weak diagonal of lowest degree */
	DEBUGm2 (("---------------------------------- Deg: "ID"\n", deg)) ;
	for (k = Head [deg] ; k != EMPTY ; k = Next [k])
	{
	    DEBUGm2 (("k: "ID"\n", k)) ;
	    if (P [k] == EMPTY)
	    {
		/* C (k,k) is a weak diagonal entry.  Find an index j != k such
		 * that C (j,k) and C (k,j) are both strong, and also such
		 * that Degree [j] is minimized.  In case of a tie, pick
		 * the smallest index j.  C and R contain the pattern of
		 * strong entries only.
		 *
		 * Note that row k of R and column k of C are both sorted. */

		DEBUGm4 (("===== Weak diagonal k = "ID"\n", k)) ;
		DEBUG1 (("Column k of C:\n")) ;
		for (p = Cp [k] ; p < Cp [k+1] ; p++)
		{
		    DEBUG1 (("    "ID": deg "ID"\n", Ci [p], Degree [Ci [p]]));
		}
		DEBUG1 (("Row k of R (strong entries only):\n")) ;
		for (p = Rp [k] ; p < Rp [k+1] ; p++)
		{
		    DEBUG1 (("    "ID": deg "ID"\n", Ri [p], Degree [Ri [p]]));
		}

		/* no (C (k,j), C (j,k)) pair exists yet */
		j_best = EMPTY ;
		jdiff_best = Int_MAX ;
		jdeg_best = Int_MAX ;

		/* pointers into column k (including values) */
		cp1 = Cp [k] ;
		cp2 = Cp [k+1] ;
		cp = cp1 ;

		/* pointers into row k (strong entries only, no values) */
		rp1 = Rp [k] ;
		rp2 = Rp [k+1] ;
		rp = rp1 ;

		/* while entries searched in column k and row k */
		while (TRUE)
		{

		    if (cp >= cp2)
		    {
			/* no more entries in this column */
			break ;
		    }

		    /* get C (j,k), which is strong */
		    j = Ci [cp] ;

		    if (rp >= rp2)
		    {
			/* no more entries in this column */
			break ;
		    }

		    /* get R (k,j2), which is strong */
		    j2 = Ri [rp] ;

		    if (j < j2)
		    {
			/* C (j,k) is strong, but R (k,j) is not strong */
			cp++ ;
			continue ;
		    }

		    if (j2 < j)
		    {
			/* C (k,j2) is strong, but R (j2,k) is not strong */
			rp++ ;
			continue ;
		    }

		    /* j == j2: C (j,k) is strong and R (k,j) is strong */

		    best = FALSE ;

		    if (P [j] == EMPTY)
		    {
			/* j has not yet been matched */
			jdeg = Degree [j] ;
			jdiff = SCALAR_ABS (k-j) ;

			DEBUG1 (("Try candidate j "ID" deg "ID" diff "ID
				    "\n", j, jdeg, jdiff)) ;

			if (j_best == EMPTY)
			{
			    /* this is the first candidate seen */
			    DEBUG1 (("   first\n")) ;
			    best = TRUE ;
			}
			else
			{
			    if (jdeg < jdeg_best)
			    {
				/* the degree of j is best seen so far. */
				DEBUG1 (("   least degree\n")) ;
				best = TRUE ;
			    }
			    else if (jdeg == jdeg_best)
			    {
				/* degree of j and j_best are the same */
				/* tie break by nearest node number */
				if (jdiff < jdiff_best)
				{
				    DEBUG1 (("   tie degree, closer\n")) ;
				    best = TRUE ;
				}
				else if (jdiff == jdiff_best)
				{
				    /* |j-k| = |j_best-k|.  For any given k
				     * and j_best there is only one other j
				     * than can be just as close as j_best.
				     * Tie break by picking the smaller of
				     * j and j_best */
				    DEBUG1 (("   tie degree, as close\n"));
				    best = j < j_best ;
				}
			    }
			    else
			    {
				/* j has higher degree than best so far */
				best = FALSE ;
			    }
			}
		    }

		    if (best)
		    {
			/* j is best match for k */
			/* found a strong pair, A (j,k) and A (k,j) */
			DEBUG1 ((" --- Found pair k: "ID" j: " ID
			    " jdeg: "ID" jdiff: "ID"\n",
			    k, j, jdeg, jdiff)) ;
			ASSERT (jdiff != EMPTY) ;
			ASSERT (jdeg != EMPTY) ;
			j_best = j ;
			jdeg_best = jdeg ;
			jdiff_best = jdiff ;
		    }

		    /* get the next entries in column k and row k */
		    cp++ ;
		    rp++ ;
		}

		/* save the pair (j,k), if we found one */
		if (j_best != EMPTY)
		{
		    j = j_best ;
		    DEBUGm4 ((" --- best pair j: "ID" for k: "ID"\n", j, k)) ;
		    P [k] = j ;
		    P [j] = k ;
		}
		else
		{
		    /* no match was found for k */
		    unmatched++ ;
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* finalize the row permutation, P */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < n2 ; k++)
    {
	if (P [k] == EMPTY)
	{
	    P [k] = k ;
	}
    }
    ASSERT (UMF_is_permutation (P, Rp, n2, n2)) ;

    return (unmatched) ;
}


/* ========================================================================== */
/* === UMF_2by2 ============================================================= */
/* ========================================================================== */

GLOBAL void UMF_2by2
(
    /* input, not modified: */
    Int n,		    /* A is n-by-n */
    const Int Ap [ ],	    /* size n+1 */
    const Int Ai [ ],	    /* size nz = Ap [n] */
    const double Ax [ ],    /* size nz if present */
#ifdef COMPLEX
    const double Az [ ],    /* size nz if present */
#endif
    double tol,		/* tolerance for determining whether or not an
			 * entry is numerically acceptable.  If tol <= 0
			 * then all numerical values ignored. */
    Int scale,		/* scaling to perform (none, sum, or max) */
    Int Cperm1 [ ],	/* singleton permutations */
#ifndef NDEBUG
    Int Rperm1 [ ],	/* not needed, since Rperm1 = Cperm1 for submatrix S */
#endif
    Int InvRperm1 [ ],	/* inverse of Rperm1 */
    Int n1,		/* number of singletons */
    Int nempty,		/* number of empty rows/cols */

    /* input, contents undefined on output: */
    Int Degree [ ],	/* Degree [j] is the number of off-diagonal
			 * entries in row/column j of S+S', where
			 * where S = A (Cperm1 [n1..], Rperm1 [n1..]).
			 * Note that S is not used, nor formed. */

    /* output: */
    Int P [ ],		/* P [k] = i means original row i is kth row in S(P,:)
			 * where S = A (Cperm1 [n1..], Rperm1 [n1..]) */
    Int *p_nweak,
    Int *p_unmatched,

    /* workspace (not defined on input or output): */
    Int Ri [ ],		/* of size >= max (nz, n) */
    Int Rp [ ],		/* of size n+1 */
    double Rs [ ],	/* of size n if present.  Rs = sum (abs (A),2) or
			 * max (abs (A),2), the sum or max of each row.  Unused
			 * if scale is equal to UMFPACK_SCALE_NONE. */
    Int Head [ ],	/* of size n.  Head pointers for bucket sort */
    Int Next [ ],	/* of size n.  Next pointers for bucket sort */
    Int Ci [ ],		/* size nz */
    Int Cp [ ]		/* size n+1 */
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int k, p, row, col, do_values, do_sum, do_max, do_scale, nweak, weak,
	p1, p2, dfound, unmatched, n2, oldrow, newrow, oldcol, newcol, pp ;

    double cmax, value, rs, ctol, dvalue ;
    Entry aij ;

#ifndef NRECIPROCAL
    Int do_recip = FALSE ;
#endif

#ifndef NDEBUG
    /* UMF_debug += 99 ; */
    DEBUGm3 (("\n ==================================UMF_2by2: tol %g\n", tol)) ;
    ASSERT (AMD_valid (n, n, Ap, Ai)) ;
    for (k = n1 ; k < n - nempty ; k++)
    {
	ASSERT (Cperm1 [k] == Rperm1 [k]) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* determine scaling options */
    /* ---------------------------------------------------------------------- */

    /* use the values, but only if they are present */
    /* ignore the values if tol <= 0 */
    do_values = (tol > 0)
	&& (Ax != (double *) NULL)
#ifdef COMPLEX
	&& (Az != (double *) NULL)
#endif
    ;
    if (do_values && (Rs != (double *) NULL))
    {
	do_sum = (scale == UMFPACK_SCALE_SUM) ;
	do_max = (scale == UMFPACK_SCALE_MAX) ;
    }
    else
    {
	/* no scaling */
	do_sum = FALSE ;
	do_max = FALSE ;
    }
    do_scale = do_max || do_sum ;
    DEBUGm3 (("do_values "ID" do_sum "ID" do_max "ID" do_scale "ID"\n",
	do_values, do_sum, do_max, do_scale)) ;

    /* ---------------------------------------------------------------------- */
    /* compute the row scaling, if requested */
    /* ---------------------------------------------------------------------- */

    /* see also umf_kernel_init */

    if (do_scale)
    {
#ifndef NRECIPROCAL
	double rsmin ;
#endif
	for (row = 0 ; row < n ; row++)
	{
	    Rs [row] = 0.0 ;
	}
	for (col = 0 ; col < n ; col++)
	{
	    p2 = Ap [col+1] ;
	    for (p = Ap [col] ; p < p2 ; p++)
	    {
		row = Ai [p] ;
		ASSIGN (aij, Ax [p], Az [p]) ;
		APPROX_ABS (value, aij) ;
		rs = Rs [row] ;
		if (!SCALAR_IS_NAN (rs))
		{
		    if (SCALAR_IS_NAN (value))
		    {
			/* if any entry in a row is NaN, then the scale factor
			 * for the row is NaN.  It will be set to 1 later. */
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
	    }
	}
#ifndef NRECIPROCAL
	rsmin = Rs [0] ;
	if (SCALAR_IS_ZERO (rsmin) || SCALAR_IS_NAN (rsmin))
	{
	    rsmin = 1.0 ;
	}
#endif
	for (row = 0 ; row < n ; row++)
	{
	    /* do not scale an empty row, or a row with a NaN */
	    rs = Rs [row] ;
	    if (SCALAR_IS_ZERO (rs) || SCALAR_IS_NAN (rs))
	    {
		Rs [row] = 1.0 ;
	    }
#ifndef NRECIPROCAL
	    rsmin = MIN (rsmin, Rs [row]) ;
#endif
	}

#ifndef NRECIPROCAL
	/* multiply by the reciprocal if Rs is not too small */
	do_recip = (rsmin >= RECIPROCAL_TOLERANCE) ;
	if (do_recip)
	{
	    /* invert the scale factors */
	    for (row = 0 ; row < n ; row++)
	    {
		Rs [row] = 1.0 / Rs [row] ;
	    }
	}
#endif
    }

    /* ---------------------------------------------------------------------- */
    /* compute the max in each column and find diagonal */
    /* ---------------------------------------------------------------------- */

    nweak = 0 ;

#ifndef NDEBUG
    for (k = 0 ; k < n ; k++)
    {
	ASSERT (Rperm1 [k] >= 0 && Rperm1 [k] < n) ;
	ASSERT (InvRperm1 [Rperm1 [k]] == k) ;
    }
#endif

    n2 = n - n1 - nempty ;

    /* use Ri to count the number of strong entries in each row */
    for (row = 0 ; row < n2 ; row++)
    {
	Ri [row] = 0 ;
    }

    pp = 0 ;
    ctol = 0 ;
    dvalue = 1 ;

    /* construct C = pruned submatrix, strong values only, column form */

    for (k = n1 ; k < n - nempty ; k++)
    {
	oldcol = Cperm1 [k] ;
	newcol = k - n1 ;
	Next [newcol] = EMPTY ;
	DEBUGm1 (("Column "ID" newcol "ID" oldcol "ID"\n", k, newcol, oldcol)) ;

	Cp [newcol] = pp ;

	dfound = FALSE ;
	p1 = Ap [oldcol] ;
	p2 = Ap [oldcol+1] ;
	if (do_values)
	{
	    cmax = 0 ;
	    dvalue = 0 ;

	    if (!do_scale)
	    {
		/* no scaling */
		for (p = p1 ; p < p2 ; p++)
		{
		    oldrow = Ai [p] ;
		    ASSERT (oldrow >= 0 && oldrow < n) ;
		    newrow = InvRperm1 [oldrow] - n1 ;
		    ASSERT (newrow >= -n1 && newrow < n2) ;
		    if (newrow < 0) continue ;
		    ASSIGN (aij, Ax [p], Az [p]) ;
		    APPROX_ABS (value, aij) ;
		    /* if either cmax or value is NaN, define cmax as NaN */
		    if (!SCALAR_IS_NAN (cmax))
		    {
			if (SCALAR_IS_NAN (value))
			{
			    cmax = value ;
			}
			else
			{
			    cmax = MAX (cmax, value) ;
			}
		    }
		    if (oldrow == oldcol)
		    {
			/* we found the diagonal entry in this column */
			dvalue = value ;
			dfound = TRUE ;
			ASSERT (newrow == newcol) ;
		    }
		}
	    }
#ifndef NRECIPROCAL
	    else if (do_recip)
	    {
		/* multiply by the reciprocal */
		for (p = p1 ; p < p2 ; p++)
		{
		    oldrow = Ai [p] ;
		    ASSERT (oldrow >= 0 && oldrow < n) ;
		    newrow = InvRperm1 [oldrow] - n1 ;
		    ASSERT (newrow >= -n1 && newrow < n2) ;
		    if (newrow < 0) continue ;
		    ASSIGN (aij, Ax [p], Az [p]) ;
		    APPROX_ABS (value, aij) ;
		    value *= Rs [oldrow] ;
		    /* if either cmax or value is NaN, define cmax as NaN */
		    if (!SCALAR_IS_NAN (cmax))
		    {
			if (SCALAR_IS_NAN (value))
			{
			    cmax = value ;
			}
			else
			{
			    cmax = MAX (cmax, value) ;
			}
		    }
		    if (oldrow == oldcol)
		    {
			/* we found the diagonal entry in this column */
			dvalue = value ;
			dfound = TRUE ;
			ASSERT (newrow == newcol) ;
		    }
		}
	    }
#endif
	    else
	    {
		/* divide instead */
		for (p = p1 ; p < p2 ; p++)
		{
		    oldrow = Ai [p] ;
		    ASSERT (oldrow >= 0 && oldrow < n) ;
		    newrow = InvRperm1 [oldrow] - n1 ;
		    ASSERT (newrow >= -n1 && newrow < n2) ;
		    if (newrow < 0) continue ;
		    ASSIGN (aij, Ax [p], Az [p]) ;
		    APPROX_ABS (value, aij) ;
		    value /= Rs [oldrow] ;
		    /* if either cmax or value is NaN, define cmax as NaN */
		    if (!SCALAR_IS_NAN (cmax))
		    {
			if (SCALAR_IS_NAN (value))
			{
			    cmax = value ;
			}
			else
			{
			    cmax = MAX (cmax, value) ;
			}
		    }
		    if (oldrow == oldcol)
		    {
			/* we found the diagonal entry in this column */
			dvalue = value ;
			dfound = TRUE ;
			ASSERT (newrow == newcol) ;
		    }
		}
	    }

	    ctol = tol * cmax ;
	    DEBUGm1 (("    cmax col "ID" %g  ctol %g\n", oldcol, cmax, ctol)) ;
	}
	else
	{
	    for (p = p1 ; p < p2 ; p++)
	    {
		oldrow = Ai [p] ;
		ASSERT (oldrow >= 0 && oldrow < n) ;
		newrow = InvRperm1 [oldrow] - n1 ;
		ASSERT (newrow >= -n1 && newrow < n2) ;
		if (newrow < 0) continue ;
		Ci [pp++] = newrow ;
		if (oldrow == oldcol)
		{
		    /* we found the diagonal entry in this column */
		    ASSERT (newrow == newcol) ;
		    dfound = TRUE ;
		}
		/* count the entries in each column */
		Ri [newrow]++ ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* flag the weak diagonals */
	/* ------------------------------------------------------------------ */

	if (!dfound)
	{
	    /* no diagonal entry present */
	    weak = TRUE ;
	}
	else
	{
	    /* diagonal entry is present, check its value */
	    weak = (do_values) ?  WEAK (dvalue, ctol) : FALSE ;
	}
	if (weak)
	{
	    /* flag this column as weak */
	    DEBUG0 (("Weak!\n")) ;
	    Next [newcol] = IS_WEAK ;
	    nweak++ ;
	}

	/* ------------------------------------------------------------------ */
	/* count entries in each row that are not numerically weak */
	/* ------------------------------------------------------------------ */

	if (do_values)
	{
	    if (!do_scale)
	    {
		/* no scaling */
		for (p = p1 ; p < p2 ; p++)
		{
		    oldrow = Ai [p] ;
		    newrow = InvRperm1 [oldrow] - n1 ;
		    if (newrow < 0) continue ;
		    ASSIGN (aij, Ax [p], Az [p]) ;
		    APPROX_ABS (value, aij) ;
		    weak = WEAK (value, ctol) ;
		    if (!weak)
		    {
			DEBUG0 (("    strong: row "ID": %g\n", oldrow, value)) ;
			Ci [pp++] = newrow ;
			Ri [newrow]++ ;
		    }
		}
	    }
#ifndef NRECIPROCAL
	    else if (do_recip)
	    {
		/* multiply by the reciprocal */
		for (p = p1 ; p < p2 ; p++)
		{
		    oldrow = Ai [p] ;
		    newrow = InvRperm1 [oldrow] - n1 ;
		    if (newrow < 0) continue ;
		    ASSIGN (aij, Ax [p], Az [p]) ;
		    APPROX_ABS (value, aij) ;
		    value *= Rs [oldrow] ;
		    weak = WEAK (value, ctol) ;
		    if (!weak)
		    {
			DEBUG0 (("    strong: row "ID": %g\n", oldrow, value)) ;
			Ci [pp++] = newrow ;
			Ri [newrow]++ ;
		    }
		}
	    }
#endif
	    else
	    {
		/* divide instead */
		for (p = p1 ; p < p2 ; p++)
		{
		    oldrow = Ai [p] ;
		    newrow = InvRperm1 [oldrow] - n1 ;
		    if (newrow < 0) continue ;
		    ASSIGN (aij, Ax [p], Az [p]) ;
		    APPROX_ABS (value, aij) ;
		    value /= Rs [oldrow] ;
		    weak = WEAK (value, ctol) ;
		    if (!weak)
		    {
			DEBUG0 (("    strong: row "ID": %g\n", oldrow, value)) ;
			Ci [pp++] = newrow ;
			Ri [newrow]++ ;
		    }
		}
	    }
	}
    }
    Cp [n2] = pp ;
    ASSERT (AMD_valid (n2, n2, Cp, Ci)) ;

    if (nweak == 0)
    {
	/* nothing to do, quick return */
	DEBUGm2 (("\n =============================UMF_2by2: quick return\n")) ;
	for (k = 0 ; k < n ; k++)
	{
	    P [k] = k ;
	}
	*p_nweak = 0 ;
	*p_unmatched = 0 ;
	return ;
    }

#ifndef NDEBUG
    for (k = 0 ; k < n2 ; k++)
    {
	P [k] = EMPTY ;
    }
    for (k = 0 ; k < n2 ; k++)
    {
	ASSERT (Degree [k] >= 0 && Degree [k] < n2) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* find the 2-by-2 permutation */
    /* ---------------------------------------------------------------------- */

    /* The matrix S is now mapped to the index range 0 to n2-1.  We have
     * S = A (Rperm [n1 .. n-nempty-1], Cperm [n1 .. n-nempty-1]), and then
     * C = pattern of strong entries in S.  A weak diagonal k in S is marked
     * with Next [k] = IS_WEAK. */

    unmatched = two_by_two (n2, Cp, Ci, Degree, Next, Ri, P, Rp, Head) ;

    /* ---------------------------------------------------------------------- */

    *p_nweak = nweak ;
    *p_unmatched = unmatched ;

#ifndef NDEBUG
    DEBUGm4 (("UMF_2by2: weak "ID"  unmatched "ID"\n", nweak, unmatched)) ;
    for (row = 0 ; row < n ; row++)
    {
	DEBUGm2 (("P ["ID"] = "ID"\n", row, P [row])) ;
    }
    DEBUGm2 (("\n =============================UMF_2by2: done\n\n")) ;
#endif
}
