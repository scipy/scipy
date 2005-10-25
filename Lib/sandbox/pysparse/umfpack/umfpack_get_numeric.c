/* ========================================================================== */
/* === UMFPACK_get_numeric ================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Gets the LU factors and the permutation vectors held in the
    Numeric object.  L is returned in sparse row form with sorted rows, U is
    returned in sparse column form with sorted columns, and P and Q are
    returned as permutation vectors.  See umfpack_get_numeric.h for a more
    detailed description.

    Returns TRUE if successful, FALSE if the Numeric object is invalid or
    if out of memory.

    Dynamic memory usage:  calls UMF_malloc twice, for a total space of
    2*n integers, and then frees all of it via UMF_free when done.

*/

#include "umf_internal.h"
#include "umf_valid_numeric.h"
#include "umf_malloc.h"
#include "umf_free.h"

#ifndef NDEBUG
PRIVATE Int init_count ;
#endif

PRIVATE void get_L
(
    Int Lp [ ],
    Int Lj [ ],
    double Lx [ ],
#ifdef COMPLEX
    double Lz [ ],
#endif
    NumericType *Numeric,
    Int Pattern [ ],
    Int Wi [ ]
) ;

PRIVATE void get_U
(
    Int Up [ ],
    Int Ui [ ],
    double Ux [ ],
#ifdef COMPLEX
    double Uz [ ],
#endif
    NumericType *Numeric,
    Int Pattern [ ],
    Int Wi [ ]
) ;

/* ========================================================================== */
/* === UMFPACK_get_numeric ================================================== */
/* ========================================================================== */

GLOBAL Int UMFPACK_get_numeric
(
    Int Lp [ ],
    Int Lj [ ],
    double Lx [ ],
#ifdef COMPLEX
    double Lz [ ],
#endif
    Int Up [ ],
    Int Ui [ ],
    double Ux [ ],
#ifdef COMPLEX
    double Uz [ ],
#endif
    Int P [ ],
    Int Q [ ],
    double Dx [ ],
#ifdef COMPLEX
    double Dz [ ],
#endif
    Int *p_do_recip,
    double Rs [ ],
    void *NumericHandle
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    NumericType *Numeric ;
    Int getL, getU, *Rperm, *Cperm, k, nn, n_row, n_col, *Wi, *Pattern,
	n_inner ;
    double *Rs1 ;
    Entry *D ;

#ifndef NDEBUG
    init_count = UMF_malloc_count ;
#endif

    Wi = (Int *) NULL ;
    Pattern = (Int *) NULL ;

    /* ---------------------------------------------------------------------- */
    /* check input parameters */
    /* ---------------------------------------------------------------------- */

    Numeric = (NumericType *) NumericHandle ;
    if (!UMF_valid_numeric (Numeric))
    {
	return (UMFPACK_ERROR_invalid_Numeric_object) ;
    }

    n_row = Numeric->n_row ;
    n_col = Numeric->n_col ;
    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

#ifdef COMPLEX
    getL = Lp && Lj && Lx && Lz ;
    getU = Up && Ui && Ux && Uz ;
#else
    getL = Lp && Lj && Lx ;
    getU = Up && Ui && Ux ;
#endif

    if (getL || getU)
    {
	Wi = (Int *) UMF_malloc (nn, sizeof (Int)) ;
	Pattern = (Int *) UMF_malloc (nn, sizeof (Int)) ;
	if (!Wi || !Pattern)
	{
	    (void) UMF_free ((void *) Wi) ;
	    (void) UMF_free ((void *) Pattern) ;
	    ASSERT (UMF_malloc_count == init_count) ;
	    DEBUGm4 (("out of memory: get numeric\n")) ;
	    return (UMFPACK_ERROR_out_of_memory) ;
	}
	ASSERT (UMF_malloc_count == init_count + 2) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get contents of Numeric */
    /* ---------------------------------------------------------------------- */

    if (P != (Int *) NULL)
    {
	Rperm = Numeric->Rperm ;
	for (k = 0 ; k < n_row ; k++)
	{
	    P [k] = Rperm [k] ;
	}
    }

    if (Q != (Int *) NULL)
    {
	Cperm = Numeric->Cperm ;
	for (k = 0 ; k < n_col ; k++)
	{
	    Q [k] = Cperm [k] ;
	}
    }

    if (getL)
    {
	get_L (Lp, Lj, Lx,
#ifdef COMPLEX
	    Lz,
#endif
	    Numeric, Pattern, Wi) ;
    }

    if (getU)
    {
	get_U (Up, Ui, Ux,
#ifdef COMPLEX
	    Uz,
#endif
	    Numeric, Pattern, Wi) ;
    }

    if (Dx != (double *) NULL)
    {
	D = Numeric->D ;
	for (k = 0 ; k < n_inner ; k++)
	{
	    Dx [k] = REAL_COMPONENT (D [k]) ;
	}
    }

#ifdef COMPLEX
    if (Dz != (double *) NULL)
    {
	D = Numeric->D ;
	for (k = 0 ; k < n_inner ; k++)
	{
	    Dz [k] = IMAG_COMPONENT (D [k]) ;
	}
    }
#endif

    /* return the flag stating whether the scale factors are to be multiplied,
     * or divided.   If do_recip is TRUE, multiply.  Otherwise, divided.
     * If NRECIPROCAL is defined at compile time, the scale factors are always
     * to be used by dividing.
     */
    if (p_do_recip != (Int *) NULL)
    {
#ifndef NRECIPROCAL
	*p_do_recip = Numeric->do_recip ;
#else
	*p_do_recip = FALSE ;
#endif
    }

    if (Rs != (double *) NULL)
    {
	Rs1 = Numeric->Rs ;
	if (Rs1 == (double *) NULL)
	{
	    /* R is the identity matrix.  */
	    for (k = 0 ; k < n_row ; k++)
	    {
		Rs [k] = 1.0 ;
	    }
	}
	else
	{
	    for (k = 0 ; k < n_row ; k++)
	    {
		Rs [k] = Rs1 [k] ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* free the workspace */
    /* ---------------------------------------------------------------------- */

    (void) UMF_free ((void *) Wi) ;
    (void) UMF_free ((void *) Pattern) ;
    ASSERT (UMF_malloc_count == init_count) ;

    return (UMFPACK_OK) ;
}


/* ========================================================================== */
/* === get_L ================================================================ */
/* ========================================================================== */

/*
    The matrix L is stored in the following arrays in the Numeric object:

	Int Lpos [0..npiv]
	Int Lip [0..npiv], index into Numeric->Memory
	Int Lilen [0..npiv]
	Unit *(Numeric->Memory), pointer to memory space holding row indices
		and numerical values

    where npiv is the number of pivot entries found.  If A is n_row-by-n_col,
    then npiv <= MIN (n_row,n_col).

    Let L_k denote the pattern of entries in column k of L (excluding the
    diagonal).

    An Lchain is a sequence of columns of L whose nonzero patterns are related.
    The start of an Lchain is denoted by a negative value of Lip [k].

    To obtain L_k:

    (1)	If column k starts an Lchain, then L_k is stored in its entirety.
	|Lip [k]| is an index into Numeric->Memory for the integer row indices
	in L_k.  The number of entries in the column is |L_k| = Lilen [k].
	This defines the pattern of the "leading" column of this chain.
	Lpos [k] is not used for the first column in the chain.  Column zero
	is always a leading column.

    (2) If column k does not start an Lchain, then L_k is represented as a
	superset of L_k-1.  Define Lnew_k such that (L_k-1 - {k} union Lnew_k)
	= L_k, where Lnew_k and (L_k-1)-{k} are disjoint.  Lnew_k are the
	entries in L_k that are not in L_k-1.  Lpos [k] holds the position of
	pivot row index k in the prior pattern L_k-1 (if it is present), so
	that the set subtraction (L_k-1)-{k} can be computed quickly, when
	computing the pattern of L_k from L_k-1.  The number of new entries in
	L_k is stored in Lilen [k] = |Lnew_k|.

	Note that this means we must have the pattern L_k-1 to compute L_k.

    In both cases (1) and (2), we obtain the pattern L_k.

    The numerical values are stored in Numeric->Memory, starting at the index
    |Lip [k]| + Lilen [k].  It is stored in the same order as the entries
    in L_k, after L_k is obtained from cases (1) or (2), above.

    The advantage of using this "packed" data structure is that it can
    dramatically reduce the amount of storage needed for the pattern of L.
    The disadvantage is that it can be difficult for the user to access,
    and it does not match the sparse matrix data structure used in MATLAB.
    Thus, this routine is provided to create a conventional sparse matrix
    data structure for L, in sparse-row form.  A row-form of L appears to
    MATLAB to be a column-oriented from of the transpose of L.  If you would
    like a column-form of L, then use UMFPACK_transpose (an example of this
    is in umfpackmex.c).

*/
/* ========================================================================== */

PRIVATE void get_L
(
    Int Lp [ ],		/* of size n_row+1 */
    Int Lj [ ],		/* of size lnz, where lnz = Lp [n_row] */
    double Lx [ ],	/* of size lnz */
#ifdef COMPLEX
    double Lz [ ],	/* of size lnz */
#endif
    NumericType *Numeric,
    Int Pattern [ ],	/* workspace of size n_row */
    Int Wi [ ]		/* workspace of size n_row */
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int deg, *ip, j, row, n_row, n_col, n_inner, *Lpos, *Lilen, *Lip, p, llen,
	lnz2, lp, newLchain, k, pos, npiv, *Li, n1 ;
    Entry *xp, value, *Lval ;

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    DEBUG4 (("get_L start:\n")) ;
    n_row = Numeric->n_row ;
    n_col = Numeric->n_col ;
    n_inner = MIN (n_row, n_col) ;
    npiv = Numeric->npiv ;
    n1 = Numeric->n1 ;
    Lpos = Numeric->Lpos ;
    Lilen = Numeric->Lilen ;
    Lip = Numeric->Lip ;
    deg = 0 ;

    /* ---------------------------------------------------------------------- */
    /* count the nonzeros in each row of L */
    /* ---------------------------------------------------------------------- */

#pragma ivdep
    for (row = 0 ; row < n_inner ; row++)
    {
	/* include the diagonal entry in the row counts */
	Wi [row] = 1 ;
    }
#pragma ivdep
    for (row = n_inner ; row < n_row ; row++)
    {
	Wi [row] = 0 ;
    }

    /* singletons */
    for (k = 0 ; k < n1 ; k++)
    {
	DEBUG4 (("Singleton k "ID"\n", k)) ;
	deg = Lilen [k] ;
	if (deg > 0)
	{
	    lp = Lip [k] ;
	    Li = (Int *) (Numeric->Memory + lp) ;
	    lp += UNITS (Int, deg) ;
	    Lval = (Entry *) (Numeric->Memory + lp) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		row = Li [j] ;
		value = Lval [j] ;
		DEBUG4 (("  row "ID"  k "ID" value", row, k)) ;
		EDEBUG4 (value) ;
		DEBUG4 (("\n")) ;
		if (IS_NONZERO (value))
		{
		    Wi [row]++ ;
		}
	    }
	}
    }

    /* non-singletons */
    for (k = n1 ; k < npiv ; k++)
    {

	/* ------------------------------------------------------------------ */
	/* make column of L in Pattern [0..deg-1] */
	/* ------------------------------------------------------------------ */

	lp = Lip [k] ;
	newLchain = (lp < 0) ;
	if (newLchain)
	{
	    lp = -lp ;
	    deg = 0 ;
	    DEBUG4 (("start of chain for column of L\n")) ;
	}

	/* remove pivot row */
	pos = Lpos [k] ;
	if (pos != EMPTY)
	{
	    DEBUG4 (("  k "ID" removing row "ID" at position "ID"\n",
	    k, Pattern [pos], pos)) ;
	    ASSERT (!newLchain) ;
	    ASSERT (deg > 0) ;
	    ASSERT (pos >= 0 && pos < deg) ;
	    ASSERT (Pattern [pos] == k) ;
	    Pattern [pos] = Pattern [--deg] ;
	}

	/* concatenate the pattern */
	ip = (Int *) (Numeric->Memory + lp) ;
	llen = Lilen [k] ;
	for (j = 0 ; j < llen ; j++)
	{
	    row = *ip++ ;
	    DEBUG4 (("  row "ID"  k "ID"\n", row, k)) ;
	    ASSERT (row > k && row < n_row) ;
	    Pattern [deg++] = row ;
	}

	xp = (Entry *) (Numeric->Memory + lp + UNITS (Int, llen)) ;

	for (j = 0 ; j < deg ; j++)
	{
	    DEBUG4 (("  row "ID"  k "ID" value", Pattern [j], k)) ;
	    row = Pattern [j] ;
	    value = *xp++ ;
	    EDEBUG4 (value) ;
	    DEBUG4 (("\n")) ;
	    if (IS_NONZERO (value))
	    {
		Wi [row]++ ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* construct the final row form of L */
    /* ---------------------------------------------------------------------- */

    /* create the row pointers */
    lnz2 = 0 ;
    for (row = 0 ; row < n_row ; row++)
    {
	Lp [row] = lnz2 ;
	lnz2 += Wi [row] ;
	Wi [row] = Lp [row] ;
    }
    Lp [n_row] = lnz2 ;
    ASSERT (Numeric->lnz + n_inner == lnz2) ;

    /* add entries from the rows of L (singletons) */
    for (k = 0 ; k < n1 ; k++)
    {
	DEBUG4 (("Singleton k "ID"\n", k)) ;
	deg = Lilen [k] ;
	if (deg > 0)
	{
	    lp = Lip [k] ;
	    Li = (Int *) (Numeric->Memory + lp) ;
	    lp += UNITS (Int, deg) ;
	    Lval = (Entry *) (Numeric->Memory + lp) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		row = Li [j] ;
		value = Lval [j] ;
		DEBUG4 (("  row "ID"  k "ID" value", row, k)) ;
		EDEBUG4 (value) ;
		DEBUG4 (("\n")) ;
		if (IS_NONZERO (value))
		{
		    p = Wi [row]++ ;
		    Lj [p] = k ;
		    Lx [p] = REAL_COMPONENT (value) ;
#ifdef COMPLEX
		    Lz [p] = IMAG_COMPONENT (value) ;
#endif
		}
	    }
	}
    }

    /* add entries from the rows of L (non-singletons) */
    for (k = n1 ; k < npiv ; k++)
    {

	/* ------------------------------------------------------------------ */
	/* make column of L in Pattern [0..deg-1] */
	/* ------------------------------------------------------------------ */

	lp = Lip [k] ;
	newLchain = (lp < 0) ;
	if (newLchain)
	{
	    lp = -lp ;
	    deg = 0 ;
	    DEBUG4 (("start of chain for column of L\n")) ;
	}

	/* remove pivot row */
	pos = Lpos [k] ;
	if (pos != EMPTY)
	{
	    DEBUG4 (("  k "ID" removing row "ID" at position "ID"\n",
	    k, Pattern [pos], pos)) ;
	    ASSERT (!newLchain) ;
	    ASSERT (deg > 0) ;
	    ASSERT (pos >= 0 && pos < deg) ;
	    ASSERT (Pattern [pos] == k) ;
	    Pattern [pos] = Pattern [--deg] ;
	}

	/* concatenate the pattern */
	ip = (Int *) (Numeric->Memory + lp) ;
	llen = Lilen [k] ;
	for (j = 0 ; j < llen ; j++)
	{
	    row = *ip++ ;
	    DEBUG4 (("  row "ID"  k "ID"\n", row, k)) ;
	    ASSERT (row > k) ;
	    Pattern [deg++] = row ;
	}

	xp = (Entry *) (Numeric->Memory + lp + UNITS (Int, llen)) ;

	for (j = 0 ; j < deg ; j++)
	{
	    DEBUG4 (("  row "ID"  k "ID" value", Pattern [j], k)) ;
	    row = Pattern [j] ;
	    value = *xp++ ;
	    EDEBUG4 (value) ;
	    DEBUG4 (("\n")) ;
	    if (IS_NONZERO (value))
	    {
		p = Wi [row]++ ;
		Lj [p] = k ;
		Lx [p] = REAL_COMPONENT (value) ;
#ifdef COMPLEX
		Lz [p] = IMAG_COMPONENT (value) ;
#endif
	    }
	}
    }

    /* add all of the diagonal entries (L is unit diagonal) */
    for (row = 0 ; row < n_inner ; row++)
    {
	p = Wi [row]++ ;
	Lj [p] = row ;
	Lx [p] = 1. ;
#ifdef COMPLEX
	Lz [p] = 0. ;
#endif
	ASSERT (Wi [row] == Lp [row+1]) ;
    }

#ifndef NDEBUG
    DEBUG6 (("L matrix (stored by rows):")) ;
    UMF_dump_col_matrix (Lx,
#ifdef COMPLEX
	Lz,
#endif
	Lj, Lp, n_inner, n_row, Numeric->lnz+n_inner) ;
#endif

    DEBUG4 (("get_L done:\n")) ;
}


/* ========================================================================== */
/* === get_U ================================================================ */
/* ========================================================================== */

/*
    The matrix U is stored in the following arrays in the Numeric object:

	Int Upos [0..npiv]
	Int Uip [0..npiv], index into Numeric->Memory
	Int Uilen [0..npiv]
	Unit *(Numeric->Memory), pointer to memory space holding column indices
		and numerical values

    where npiv is the number of pivot entries found.  If A is n_row-by-n_col,
    then npiv <= MIN (n_row,n_col).

    Let U_k denote the pattern of entries in row k of U (excluding the
    diagonal).

    A Uchain is a sequence of columns of U whose nonzero patterns are related.
    The start of a Uchain is denoted by a negative value of Uip [k].

    To obtain U_k-1:

    (1) If row k is the start of a Uchain then Uip [k] is negative and |Uip [k]|
	is an index into Numeric->Memory for the integer column indices in
	U_k-1.  The number of entries in the row is |U_k-1| = Uilen [k].  This
	defines the pattern of the "trailing" row of this chain that ends at
	row k-1.


    (2) If row k is not the start of a Uchain, then U_k-1 is a subset of U_k.
	The indices in U_k are arranged so that last Uilen [k] entries of
	U_k are those indices not in U_k-1.  Next, the pivot column index k is
	added if it appears in row U_k-1 (it never appears in U_k).  Upos [k]
	holds the position of pivot column index k in the pattern U_k-1 (if it
	is present), so that the set union (U_k-1)+{k} can be computed quickly,
	when computing the pattern of U_k-1 from U_k.

	Note that this means we must have the pattern U_k to compute L_k-1.

    In both cases (1) and (2), we obtain the pattern U_k.

    The numerical values are stored in Numeric->Memory.  If k is the start of a
    Uchain, then the offset is |Uip [k]| plus the size of the space needed to
    store the pattern U_k-1.  Otherwise, Uip [k] is the offset itself of the
    numerical values, since in this case no pattern is stored.
    The numerical values are stored in the same order as the entries in U_k,
    after U_k is obtained from cases (1) or (2), above.

    The advantage of using this "packed" data structure is that it can
    dramatically reduce the amount of storage needed for the pattern of U.
    The disadvantage is that it can be difficult for the user to access,
    and it does not match the sparse matrix data structure used in MATLAB.
    Thus, this routine is provided to create a conventional sparse matrix
    data structure for U, in sparse-column form.

*/
/* ========================================================================== */

PRIVATE void get_U
(
    Int Up [ ],		/* of size n_col+1 */
    Int Ui [ ],		/* of size unz, where unz = Up [n_col] */
    double Ux [ ],	/* of size unz */
#ifdef COMPLEX
    double Uz [ ],	/* of size unz */
#endif
    NumericType *Numeric,
    Int Pattern [ ],	/* workspace of size n_col */
    Int Wi [ ]		/* workspace of size n_col */
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int deg, j, *ip, col, *Upos, *Uilen, *Uip, n_col, ulen, *Usi,
	unz2, p, k, up, newUchain, pos, npiv, n1 ;
    Entry *xp, *D, value, *Uval ;

#ifndef NDEBUG
    Int nnzpiv = 0 ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    DEBUG4 (("get_U start:\n")) ;
    n_col = Numeric->n_col ;
    n1 = Numeric->n1 ;
    npiv = Numeric->npiv ;
    Upos = Numeric->Upos ;
    Uilen = Numeric->Uilen ;
    Uip = Numeric->Uip ;
    D = Numeric->D ;

    /* ---------------------------------------------------------------------- */
    /* count the nonzeros in each column of U */
    /* ---------------------------------------------------------------------- */

    for (col = 0 ; col < npiv ; col++)
    {
	/* include the diagonal entry in the column counts */
	DEBUG4 (("D ["ID"] = ", col)) ;
	EDEBUG4 (D [col]) ;
	Wi [col] = IS_NONZERO (D [col]) ;
	DEBUG4 ((" is nonzero: "ID"\n", Wi [col])) ;
#ifndef NDEBUG
	nnzpiv += IS_NONZERO (D [col]) ;
#endif
    }
    DEBUG4 (("nnzpiv "ID" "ID"\n", nnzpiv, Numeric->nnzpiv)) ;
    ASSERT (nnzpiv == Numeric->nnzpiv) ;
    for (col = npiv ; col < n_col ; col++)
    {
	/* diagonal entries are zero for structurally singular part */
	Wi [col] = 0 ;
    }

    deg = Numeric->ulen ;
    if (deg > 0)
    {
	/* make last pivot row of U (singular matrices only) */
	DEBUG0 (("Last pivot row of U: ulen "ID"\n", deg)) ;
	for (j = 0 ; j < deg ; j++)
	{
	    Pattern [j] = Numeric->Upattern [j] ;
	    DEBUG0 (("    column "ID"\n", Pattern [j])) ;
	}
    }

    /* non-singletons */
    for (k = npiv-1 ; k >= n1 ; k--)
    {

	/* ------------------------------------------------------------------ */
	/* use row k of U */
	/* ------------------------------------------------------------------ */

	up = Uip [k] ;
	ulen = Uilen [k] ;
	newUchain = (up < 0) ;
	if (newUchain)
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
	    DEBUG4 (("  k "ID" col "ID" value\n", k, Pattern [j])) ;
	    col = Pattern [j] ;
	    ASSERT (col >= 0 && col < n_col) ;
	    value = *xp++ ;
	    EDEBUG4 (value) ;
	    DEBUG4 (("\n")) ;
	    if (IS_NONZERO (value))
	    {
		Wi [col]++ ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* make row k-1 of U in Pattern [0..deg-1] */
	/* ------------------------------------------------------------------ */

	if (k == n1) break ;

	if (newUchain)
	{
	    /* next row is a new Uchain */
	    deg = ulen ;
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
	else
	{
	    deg -= ulen ;
	    DEBUG4 (("middle of chain for row of U "ID" deg "ID"\n", k-1, deg));
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
    }

    /* singletons */
    for (k = n1 - 1 ; k >= 0 ; k--)
    {
	deg = Uilen [k] ;
	DEBUG4 (("Singleton k "ID"\n", k)) ;
	if (deg > 0)
	{
	    up = Uip [k] ;
	    Usi = (Int *) (Numeric->Memory + up) ;
	    up += UNITS (Int, deg) ;
	    Uval = (Entry *) (Numeric->Memory + up) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		col = Usi [j] ;
		value = Uval [j] ;
		DEBUG4 (("  k "ID" col "ID" value", k, col)) ;
		EDEBUG4 (value) ;
		DEBUG4 (("\n")) ;
		if (IS_NONZERO (value))
		{
		    Wi [col]++ ;
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* construct the final column form of U */
    /* ---------------------------------------------------------------------- */

    /* create the column pointers */
    unz2 = 0 ;
    for (col = 0 ; col < n_col ; col++)
    {
	Up [col] = unz2 ;
	unz2 += Wi [col] ;
    }
    Up [n_col] = unz2 ;
    DEBUG1 (("Numeric->unz "ID"  npiv "ID" nnzpiv "ID" unz2 "ID"\n",
	Numeric->unz, npiv, Numeric->nnzpiv, unz2)) ;
    ASSERT (Numeric->unz + Numeric->nnzpiv == unz2) ;

    for (col = 0 ; col < n_col ; col++)
    {
	Wi [col] = Up [col+1] ;
    }

    /* add all of the diagonal entries */
    for (col = 0 ; col < npiv ; col++)
    {
	if (IS_NONZERO (D [col]))
	{
	    p = --(Wi [col]) ;
	    Ui [p] = col ;
	    Ux [p] = REAL_COMPONENT (D [col]) ;
#ifdef COMPLEX
	    Uz [p] = IMAG_COMPONENT (D [col]) ;
#endif
	}
    }

    /* add all the entries from the rows of U */

    deg = Numeric->ulen ;
    if (deg > 0)
    {
	/* make last pivot row of U (singular matrices only) */
	for (j = 0 ; j < deg ; j++)
	{
	    Pattern [j] = Numeric->Upattern [j] ;
	}
    }

    /* non-singletons */
    for (k = npiv-1 ; k >= n1 ; k--)
    {

	/* ------------------------------------------------------------------ */
	/* use row k of U */
	/* ------------------------------------------------------------------ */

	up = Uip [k] ;
	ulen = Uilen [k] ;
	newUchain = (up < 0) ;
	if (newUchain)
	{
	    up = -up ;
	    xp = (Entry *) (Numeric->Memory + up + UNITS (Int, ulen)) ;
	}
	else
	{
	    xp = (Entry *) (Numeric->Memory + up) ;
	}

	xp += deg ;
	for (j = deg-1 ; j >= 0 ; j--)
	{
	    DEBUG4 (("  k "ID" col "ID" value", k, Pattern [j])) ;
	    col = Pattern [j] ;
	    ASSERT (col >= 0 && col < n_col) ;
	    value = *(--xp) ;
	    EDEBUG4 (value) ;
	    DEBUG4 (("\n")) ;
	    if (IS_NONZERO (value))
	    {
		p = --(Wi [col]) ;
		Ui [p] = k ;
		Ux [p] = REAL_COMPONENT (value) ;
#ifdef COMPLEX
		Uz [p] = IMAG_COMPONENT (value) ;
#endif
	    }
	}

	/* ------------------------------------------------------------------ */
	/* make row k-1 of U in Pattern [0..deg-1] */
	/* ------------------------------------------------------------------ */

	if (newUchain)
	{
	    /* next row is a new Uchain */
	    deg = ulen ;
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
	else
	{
	    deg -= ulen ;
	    DEBUG4 (("middle of chain for row of U "ID" deg "ID"\n", k-1, deg));
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
    }

    /* singletons */
    for (k = n1 - 1 ; k >= 0 ; k--)
    {
	deg = Uilen [k] ;
	DEBUG4 (("Singleton k "ID"\n", k)) ;
	if (deg > 0)
	{
	    up = Uip [k] ;
	    Usi = (Int *) (Numeric->Memory + up) ;
	    up += UNITS (Int, deg) ;
	    Uval = (Entry *) (Numeric->Memory + up) ;
	    for (j = 0 ; j < deg ; j++)
	    {
		col = Usi [j] ;
		value = Uval [j] ;
		DEBUG4 (("  k "ID" col "ID" value", k, col)) ;
		EDEBUG4 (value) ;
		DEBUG4 (("\n")) ;
		if (IS_NONZERO (value))
		{
		    p = --(Wi [col]) ;
		    Ui [p] = k ;
		    Ux [p] = REAL_COMPONENT (value) ;
#ifdef COMPLEX
		    Uz [p] = IMAG_COMPONENT (value) ;
#endif
		}
	    }
	}
    }

#ifndef NDEBUG
    DEBUG6 (("U matrix:")) ;
    UMF_dump_col_matrix (Ux,
#ifdef COMPLEX
	Uz,
#endif
	Ui, Up, Numeric->n_row, n_col, Numeric->unz + Numeric->nnzpiv) ;
#endif

}
