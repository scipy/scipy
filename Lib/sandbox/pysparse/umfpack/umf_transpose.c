/* ========================================================================== */
/* === UMF_transpose ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*  Not user-callable.  Computes a permuted transpose, R = (A (P,Q(1:nq)))' in
	MATLAB notation, where R is in column-form.  A is n_row-by-n_col, the
	row-form matrix R is n_row-by-nq, where nq <= n_col.  A may be singular.
	The complex version can do transpose (') or array transpose (.').

	Uses Gustavson's method (Two Fast Algorithms for Sparse Matrices:
	Multiplication and Permuted Transposition, ACM Trans. on Math. Softw.,
	vol 4, no 3, pp. 250-269).
*/

#include "umf_internal.h"
#include "umf_is_permutation.h"

GLOBAL Int UMF_transpose
(
    Int n_row,			/* A is n_row-by-n_col */
    Int n_col,
    const Int Ap [ ],		/* size n_col+1 */
    const Int Ai [ ],		/* size nz = Ap [n_col] */
    const double Ax [ ],	/* size nz if present */

    const Int P [ ],	/* P [k] = i means original row i is kth row in A(P,Q)*/
			/* P is identity if not present */
			/* size n_row, if present */

    const Int Q [ ],	/* Q [k] = j means original col j is kth col in A(P,Q)*/
			/* Q is identity if not present */
			/* size nq, if present */
    Int nq,		/* size of Q, ignored if Q is (Int *) NULL */

			/* output matrix: Rp, Ri, Rx, and Rz: */
    Int Rp [ ],		/* size n_row+1 */
    Int Ri [ ],		/* size nz */
    double Rx [ ],	/* size nz, if present */

    Int W [ ],		/* size max (n_row,n_col) workspace */

    Int check		/* if true, then check inputs */
#ifdef COMPLEX
    , const double Az [ ]	/* size nz */
    , double Rz [ ]		/* size nz */
    , Int do_conjugate		/* if true, then do conjugate transpose */
				/* otherwise, do array transpose */
#endif
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int i, j, k, p, bp, newj, do_values ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    Int nz ;
    ASSERT (n_col >= 0) ;
    nz = (Ap != (Int *) NULL) ? Ap [n_col] : 0 ;
    DEBUG2 (("UMF_transpose: "ID"-by-"ID" nz "ID"\n", n_row, n_col, nz)) ;
#endif

    if (check)
    {
	/* UMFPACK_symbolic skips this check */
	/* UMFPACK_transpose always does this check */
	if (!Ai || !Ap || !Ri || !Rp || !W)
	{
	    return (UMFPACK_ERROR_argument_missing) ;
	}
	if (n_row <= 0 || n_col <= 0)		/* n_row,n_col must be > 0 */
	{
	    return (UMFPACK_ERROR_n_nonpositive) ;
	}
	if (!UMF_is_permutation (P, W, n_row, n_row) ||
	    !UMF_is_permutation (Q, W, nq, nq))
	{
	    return (UMFPACK_ERROR_invalid_permutation) ;
	}
	if (!AMD_valid (n_row, n_col, Ap, Ai))
	{
	    return (UMFPACK_ERROR_invalid_matrix) ;
	}
    }

#ifndef NDEBUG
    DEBUG2 (("UMF_transpose, input matrix:\n")) ;
    UMF_dump_col_matrix (Ax,
#ifdef COMPLEX
	Az,
#endif
	Ai, Ap, n_row, n_col, nz) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* count the entries in each row of A */
    /* ---------------------------------------------------------------------- */

    /* use W as workspace for RowCount */

    for (i = 0 ; i < n_row ; i++)
    {
	W [i] = 0 ;
	Rp [i] = 0 ;
    }

    if (Q != (Int *) NULL)
    {
	for (newj = 0 ; newj < nq ; newj++)
	{
	    j = Q [newj] ;
	    ASSERT (j >= 0 && j < n_col) ;
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		i = Ai [p] ;
		ASSERT (i >= 0 && i < n_row) ;
		W [i]++ ;
	    }
	}
    }
    else
    {
	for (j = 0 ; j < n_col ; j++)
	{
	    for (p = Ap [j] ; p < Ap [j+1] ; p++)
	    {
		i = Ai [p] ;
		ASSERT (i >= 0 && i < n_row) ;
		W [i]++ ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* compute the row pointers for R = A (P,Q) */
    /* ---------------------------------------------------------------------- */

    if (P != (Int *) NULL)
    {
	Rp [0] = 0 ;
	for (k = 0 ; k < n_row ; k++)
	{
	    i = P [k] ;
	    ASSERT (i >= 0 && i < n_row) ;
	    Rp [k+1] = Rp [k] + W [i] ;
	}
	for (k = 0 ; k < n_row ; k++)
	{
	    i = P [k] ;
	    ASSERT (i >= 0 && i < n_row) ;
	    W [i] = Rp [k] ;
	}
    }
    else
    {
	Rp [0] = 0 ;
	for (i = 0 ; i < n_row ; i++)
	{
	    Rp [i+1] = Rp [i] + W [i] ;
	}
	for (i = 0 ; i < n_row ; i++)
	{
	    W [i] = Rp [i] ;
	}
    }
    ASSERT (Rp [n_row] <= Ap [n_col]) ;

    /* at this point, W holds the permuted row pointers */

    /* ---------------------------------------------------------------------- */
    /* construct the row form of B */
    /* ---------------------------------------------------------------------- */

    do_values = Ax && Rx ;
#ifdef COMPLEX
    do_values = do_values && Az && Rz ;
#endif

#ifdef COMPLEX
    if (do_conjugate && do_values)
    {
	if (Q != (Int *) NULL)
	{

		/* R = A (P,Q)' */
		for (newj = 0 ; newj < nq ; newj++)
		{
		    j = Q [newj] ;
		    ASSERT (j >= 0 && j < n_col) ;
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			bp = W [Ai [p]]++ ;
			Ri [bp] = newj ;
			Rx [bp] = Ax [p] ;
			Rz [bp] = -Az [p] ;
		    }
		}

	}
	else
	{

		/* R = A (P,:)' */
		for (j = 0 ; j < n_col ; j++)
		{
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			bp = W [Ai [p]]++ ;
			Ri [bp] = j ;
			Rx [bp] = Ax [p] ;
			Rz [bp] = -Az [p] ;
		    }
		}

	}
    }
    else
#endif
    {
	if (Q != (Int *) NULL)
	{
	    if (do_values)
	    {

		/* R = A (P,Q).' */
		for (newj = 0 ; newj < nq ; newj++)
		{
		    j = Q [newj] ;
		    ASSERT (j >= 0 && j < n_col) ;
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			bp = W [Ai [p]]++ ;
			Ri [bp] = newj ;
			Rx [bp] = Ax [p] ;
#ifdef COMPLEX
			Rz [bp] = Az [p] ;
#endif
		    }
		}

	    }
	    else
	    {

		/* R = pattern of A (P,Q).' */
		for (newj = 0 ; newj < nq ; newj++)
		{
		    j = Q [newj] ;
		    ASSERT (j >= 0 && j < n_col) ;
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			Ri [W [Ai [p]]++] = newj ;
		    }
		}

	    }
	}
	else
	{
	    if (do_values)
	    {

		/* R = A (P,:).' */
		for (j = 0 ; j < n_col ; j++)
		{
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			bp = W [Ai [p]]++ ;
			Ri [bp] = j ;
			Rx [bp] = Ax [p] ;
#ifdef COMPLEX
			Rz [bp] = Az [p] ;
#endif
		    }
		}

	    }
	    else
	    {

		/* R = pattern of A (P,:).' */
		for (j = 0 ; j < n_col ; j++)
		{
		    for (p = Ap [j] ; p < Ap [j+1] ; p++)
		    {
			Ri [W [Ai [p]]++] = j ;
		    }
		}

	    }
	}

    }

#ifndef NDEBUG
    for (k = 0 ; k < n_row ; k++)
    {
	if (P != (Int *) NULL)
	{
	    i = P [k] ;
	}
	else
	{
	    i = k ;
	}
	DEBUG3 ((ID":  W[i] "ID" Rp[k+1] "ID"\n", i, W [i], Rp [k+1])) ;
	ASSERT (W [i] == Rp [k+1]) ;
    }
    DEBUG2 (("UMF_transpose, output matrix:\n")) ;
    UMF_dump_col_matrix (Rx,
#ifdef COMPLEX
	Rz,
#endif
	Ri, Rp, n_col, n_row, Rp [n_row]) ;
    ASSERT (AMD_valid (n_col, n_row, Rp, Ri)) ;
#endif

    return (UMFPACK_OK) ;
}
