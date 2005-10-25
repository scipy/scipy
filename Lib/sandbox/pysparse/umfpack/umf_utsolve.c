/* ========================================================================== */
/* === UMF_utsolve ========================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*  solves U'x = b or U.'x=b, where U is the upper triangular factor of a */
/*  matrix.  B is overwritten with the solution X. */
/*  Returns the floating point operation count */

#include "umf_internal.h"

GLOBAL double
#ifdef CONJUGATE_SOLVE
UMF_uhsolve			/* solve U'x=b  (complex conjugate transpose) */
#else
UMF_utsolve			/* solve U.'x=b (array transpose) */
#endif
(
    NumericType *Numeric,
    Entry X [ ],		/* b on input, solution x on output */
    Int Pattern [ ]		/* a work array of size n */
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int k, deg, j, *ip, col, *Upos, *Uilen, kstart, kend, up,
	*Uip, n, uhead, ulen, pos, npiv, n1, *Ui ;
    Entry *xp, xk, *D, *Uval ;

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    if (Numeric->n_row != Numeric->n_col) return (0.) ;
    n = Numeric->n_row ;
    npiv = Numeric->npiv ;
    Upos = Numeric->Upos ;
    Uilen = Numeric->Uilen ;
    Uip = Numeric->Uip ;
    D = Numeric->D ;
    kend = 0 ;
    n1 = Numeric->n1 ;

#ifndef NDEBUG
    DEBUG4 (("Utsolve start: npiv "ID" n "ID"\n", npiv, n)) ;
    for (j = 0 ; j < n ; j++)
    {
	DEBUG4 (("Utsolve start "ID": ", j)) ;
	EDEBUG4 (X [j]) ;
	DEBUG4 (("\n")) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* singletons */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < n1 ; k++)
    {
	DEBUG4 (("Singleton k "ID"\n", k)) ;
	/* Go ahead and divide by zero if D [k] is zero. */
#ifdef CONJUGATE_SOLVE
	/* xk = X [k] / conjugate (D [k]) ; */
	DIV_CONJ (xk, X [k], D [k]) ;
#else
	/* xk = X [k] / D [k] ; */
	DIV (xk, X [k], D [k]) ;
#endif
	X [k] = xk ;
	deg = Uilen [k] ;
	if (deg > 0 && IS_NONZERO (xk))
	{
	    up = Uip [k] ;
	    Ui = (Int *) (Numeric->Memory + up) ;
	    up += UNITS (Int, deg) ;
	    Uval = (Entry *) (Numeric->Memory + up) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		DEBUG4 (("  k "ID" col "ID" value", k, Ui [j])) ;
		EDEBUG4 (Uval [j]) ;
		DEBUG4 (("\n")) ;
#ifdef CONJUGATE_SOLVE
		/* X [Ui [j]] -= xk * conjugate (Uval [j]) ; */
		MULT_SUB_CONJ (X [Ui [j]], xk, Uval [j]) ;
#else
		/* X [Ui [j]] -= xk * Uval [j] ; */
		MULT_SUB (X [Ui [j]], xk, Uval [j]) ;
#endif
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* nonsingletons */
    /* ---------------------------------------------------------------------- */

    for (kstart = n1 ; kstart < npiv ; kstart = kend + 1)
    {

	/* ------------------------------------------------------------------ */
	/* find the end of this Uchain */
	/* ------------------------------------------------------------------ */

	DEBUG4 (("kstart "ID" kend "ID"\n", kstart, kend)) ;
	/* for (kend = kstart ; kend < npiv && Uip [kend+1] > 0 ; kend++) ; */
	kend = kstart ;
	while (kend < npiv && Uip [kend+1] > 0)
	{
	    kend++ ;
	}

	/* ------------------------------------------------------------------ */
	/* scan the whole Uchain to find the pattern of the first row of U */
	/* ------------------------------------------------------------------ */

	k = kend+1 ;
	DEBUG4 (("\nKend "ID" K "ID"\n", kend, k)) ;

	/* ------------------------------------------------------------------ */
	/* start with last row in Uchain of U in Pattern [0..deg-1] */
	/* ------------------------------------------------------------------ */

	if (k == npiv)
	{
	    deg = Numeric->ulen ;
	    if (deg > 0)
	    {
		/* :: make last pivot row of U (singular matrices only) :: */
		for (j = 0 ; j < deg ; j++)
		{
		    Pattern [j] = Numeric->Upattern [j] ;
		}
	    }
	}
	else
	{
	    ASSERT (k >= 0 && k < npiv) ;
	    up = -Uip [k] ;
	    ASSERT (up > 0) ;
	    deg = Uilen [k] ;
	    DEBUG4 (("end of chain for row of U "ID" deg "ID"\n", k-1, deg)) ;
	    ip = (Int *) (Numeric->Memory + up) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		col = *ip++ ;
		DEBUG4 (("  k "ID" col "ID"\n", k-1, col)) ;
		ASSERT (k <= col) ;
		Pattern [j] = col ;
	    }
	}

	/* empty the stack at the bottom of Pattern */
	uhead = n ;

	for (k = kend ; k > kstart ; k--)
	{
	    /* Pattern [0..deg-1] is the pattern of row k of U */

	    /* -------------------------------------------------------------- */
	    /* make row k-1 of U in Pattern [0..deg-1] */
	    /* -------------------------------------------------------------- */

	    ASSERT (k >= 0 && k < npiv) ;
	    ulen = Uilen [k] ;
	    /* delete, and push on the stack */
	    for (j = 0 ; j < ulen ; j++)
	    {
		ASSERT (uhead >= deg) ;
		Pattern [--uhead] = Pattern [--deg] ;
	    }
	    DEBUG4 (("middle of chain for row of U "ID" deg "ID"\n", k, deg)) ;
	    ASSERT (deg >= 0) ;

	    pos = Upos [k] ;
	    if (pos != EMPTY)
	    {
		/* add the pivot column */
		DEBUG4 (("k "ID" add pivot entry at position "ID"\n", k, pos)) ;
		ASSERT (pos >= 0 && pos <= deg) ;
		Pattern [deg++] = Pattern [pos] ;
		Pattern [pos] = k ;
	    }
	}

	/* Pattern [0..deg-1] is now the pattern of the first row in Uchain */

	/* ------------------------------------------------------------------ */
	/* solve using this Uchain, in reverse order */
	/* ------------------------------------------------------------------ */

	DEBUG4 (("Unwinding Uchain\n")) ;
	for (k = kstart ; k <= kend ; k++)
	{

	    /* -------------------------------------------------------------- */
	    /* construct row k */
	    /* -------------------------------------------------------------- */

	    ASSERT (k >= 0 && k < npiv) ;
	    pos = Upos [k] ;
	    if (pos != EMPTY)
	    {
		/* remove the pivot column */
		DEBUG4 (("k "ID" add pivot entry at position "ID"\n", k, pos)) ;
		ASSERT (k > kstart) ;
		ASSERT (pos >= 0 && pos < deg) ;
		ASSERT (Pattern [pos] == k) ;
		Pattern [pos] = Pattern [--deg] ;
	    }

	    up = Uip [k] ;
	    ulen = Uilen [k] ;
	    if (k > kstart)
	    {
		/* concatenate the deleted pattern; pop from the stack */
		for (j = 0 ; j < ulen ; j++)
		{
		    ASSERT (deg <= uhead && uhead < n) ;
		    Pattern [deg++] = Pattern [uhead++] ;
		}
		DEBUG4 (("middle of chain, row of U "ID" deg "ID"\n", k, deg)) ;
		ASSERT (deg >= 0) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* use row k of U */
	    /* -------------------------------------------------------------- */

	    /* Go ahead and divide by zero if D [k] is zero. */
#ifdef CONJUGATE_SOLVE
	    /* xk = X [k] / conjugate (D [k]) ; */
	    DIV_CONJ (xk, X [k], D [k]) ;
#else
	    /* xk = X [k] / D [k] ; */
	    DIV (xk, X [k], D [k]) ;
#endif
	    X [k] = xk ;
	    if (IS_NONZERO (xk))
	    {
		if (k == kstart)
		{
		    up = -up ;
		    xp = (Entry *) (Numeric->Memory + up + UNITS (Int, ulen)) ;
		}
		else
		{
		    xp = (Entry *) (Numeric->Memory + up) ;
		}
		for (j = 0 ; j < deg ; j++)
		{
		    DEBUG4 (("  k "ID" col "ID" value", k, Pattern [j])) ;
		    EDEBUG4 (*xp) ;
		    DEBUG4 (("\n")) ;
#ifdef CONJUGATE_SOLVE
		    /* X [Pattern [j]] -= xk * conjugate (*xp) ; */
		    MULT_SUB_CONJ (X [Pattern [j]], xk, *xp) ;
#else
		    /* X [Pattern [j]] -= xk * (*xp) ; */
		    MULT_SUB (X [Pattern [j]], xk, *xp) ;
#endif
		    xp++ ;
		}
	    }
	}
	ASSERT (uhead == n) ;
    }

    for (k = npiv ; k < n ; k++)
    {
	/* This is an *** intentional *** divide-by-zero, to get Inf or Nan,
	 * as appropriate.  It is not a bug. */
	ASSERT (IS_ZERO (D [k])) ;
	/* For conjugate solve, D [k] == conjugate (D [k]), in this case */
	/* xk = X [k] / D [k] ; */
	DIV (xk, X [k], D [k]) ;
	X [k] = xk ;
    }

#ifndef NDEBUG
    for (j = 0 ; j < n ; j++)
    {
	DEBUG4 (("Utsolve done "ID": ", j)) ;
	EDEBUG4 (X [j]) ;
	DEBUG4 (("\n")) ;
    }
    DEBUG4 (("Utsolve done.\n")) ;
#endif

    return (DIV_FLOPS * ((double) n) + MULTSUB_FLOPS * ((double) Numeric->unz));
}
