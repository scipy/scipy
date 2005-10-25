/* ========================================================================== */
/* === UMF_solve ============================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Not user-callable.  Solves a linear system using the numerical factorization
    computed by UMFPACK_numeric.  No workspace is dynamically allocated.  Counts
    flops, but excludes floating-point comparisons (thus real abs (...) are
    zero flops, but complex abs (...) takes 6 flops).

    Returns UMFPACK_OK if successful, UMFPACK_ERROR_argument_missing if
    required arguments are missing, UMFPACK_ERROR_invalid_system if the sys
    string is not valid or if the matrix A is not square.

    Uses the sparse backward error method of Arioli, Demmel, and Duff
    (Solving sparse linear systems with sparse backward error, SIAM J. Matrix
    Analysis and Applic., vol 10, pp. 165-190).
*/

#include "umf_internal.h"
#include "umf_lsolve.h"
#include "umf_usolve.h"
#include "umf_ltsolve.h"
#include "umf_utsolve.h"
#include "umf_report_vector.h"

PRIVATE Int do_step
(
    double omega [3],
    Int step,
    const double B2 [ ],
    Entry X [ ],
    const Entry W [ ],
    const double Y [ ],
    const double Z2 [ ],
    Entry S [ ],
    Int n,
    double Info [UMFPACK_INFO]
) ;

/* ========================================================================== */
/* === UMF_solve ============================================================ */
/* ========================================================================== */

GLOBAL Int UMF_solve
(
    Int sys,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
    double Xx [ ],
    const double Bx [ ],
#ifdef COMPLEX
    const double Az [ ],
    double Xz [ ],
    const double Bz [ ],
#endif
    NumericType *Numeric,
    Int irstep,
    double Info [UMFPACK_INFO],
    Int Pattern [ ],		/* size n */
    double SolveWork [ ]	/* if irstep>0 real:  size 5*n.  complex:10*n */
				/* otherwise   real:  size   n.  complex: 4*n */
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int *Rperm, *Cperm, i, n, p, step, j, nz, status, p2, do_scale ;
    Entry axx, wi, xj, zi, xi, aij, *W, *Z, *S, *X, bi ;
    double omega [3], d, *Z2, *Y, z2i, yi, *B2, *Rs, flops ;

#ifndef NRECIPROCAL
    Int do_recip = Numeric->do_recip ;
#endif

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    UMF_dump_lu (Numeric) ;
    ASSERT (Numeric && Xx && Bx && Pattern && SolveWork && Info) ;
#ifdef COMPLEX
    ASSERT (Xz && Bz) ;
#endif
#endif

    nz = 0 ;
    omega [0] = 0. ;
    omega [1] = 0. ;
    omega [2] = 0. ;
    Rperm = Numeric->Rperm ;
    Cperm = Numeric->Cperm ;
    Rs = Numeric->Rs ;		/* row scale factors */
    do_scale = (Rs != (double *) NULL) ;
    flops = 0 ;
    Info [UMFPACK_SOLVE_FLOPS] = 0 ;
    Info [UMFPACK_IR_TAKEN] = 0 ;
    Info [UMFPACK_IR_ATTEMPTED] = 0 ;

    /* UMFPACK_solve does not call this routine if A is rectangular */
    ASSERT (Numeric->n_row == Numeric->n_col) ;
    n = Numeric->n_row ;
    if (Numeric->nnzpiv < n
	|| SCALAR_IS_ZERO (Numeric->rcond) || SCALAR_IS_NAN (Numeric->rcond))
    {
	/* Note that systems involving just L return UMFPACK_OK, even if */
	/* A is singular (L is always has a unit diagonal). */
	DEBUGm4 (("Note, matrix is singular in umf_solve\n")) ;
	status = UMFPACK_WARNING_singular_matrix ;
	irstep = 0 ;
    }
    else
    {
	status = UMFPACK_OK ;
    }
    irstep = MAX (0, irstep) ;			/* make sure irstep is >= 0 */

    W = (Entry *) SolveWork ;			/* Entry W [0..n-1] */

    Z = (Entry *) NULL ;	/* unused if no iterative refinement */
    S = (Entry *) NULL ;
    Y = (double *) NULL ;
    Z2 = (double *) NULL ;
    B2 = (double *) NULL ;

#ifdef COMPLEX
    X = (Entry *) (SolveWork + 2*n) ;		/* Entry X [0..n-1] */
    if (irstep > 0)
    {
	if (!Ap || !Ai || !Ax || !Az)
	{
	    return (UMFPACK_ERROR_argument_missing) ;
	}
	Z = (Entry *) (SolveWork + 4*n) ;	/* Entry Z [0..n-1] */
	S = (Entry *) (SolveWork + 6*n) ;	/* Entry S [0..n-1] */
	Y = (double *) (SolveWork + 8*n) ;	/* double Y [0..n-1] */
	B2 = (double *) (SolveWork + 9*n) ;	/* double B2 [0..n-1] */
	Z2 = (double *) Z ;		/* double Z2 [0..n-1], equiv. to Z */
    }
#else
    X = (Entry *) Xx ;				/* Entry X [0..n-1] */
    if (irstep > 0)
    {
	if (!Ap || !Ai || !Ax)
	{
	    return (UMFPACK_ERROR_argument_missing) ;
	}
	Z = (Entry *) (SolveWork + n) ;		/* Entry Z [0..n-1] */
	S = (Entry *) (SolveWork + 2*n) ;	/* Entry S [0..n-1] */
	Y = (double *) (SolveWork + 3*n) ;	/* double Y [0..n-1] */
	B2 = (double *) (SolveWork + 4*n) ;	/* double B2 [0..n-1] */
	Z2 = (double *) Z ;		/* double Z2 [0..n-1], equiv. to Z */
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* determine which system to solve */
    /* ---------------------------------------------------------------------- */

    if (sys == UMFPACK_A)
    {

	/* ------------------------------------------------------------------ */
	/* solve A x = b with optional iterative refinement */
	/* ------------------------------------------------------------------ */

	if (irstep > 0)
	{

	    /* -------------------------------------------------------------- */
	    /* using iterative refinement:  compute Y and B2 */
	    /* -------------------------------------------------------------- */

	    nz = Ap [n] ;
	    Info [UMFPACK_NZ] = nz ;

	    /* A is stored by column */
	    /* Y (i) = ||R A_i||, 1-norm of row i of R A */
	    for (i = 0 ; i < n ; i++)
	    {
		Y [i] = 0. ;
	    }
	    flops += (ABS_FLOPS + 1) * nz ;
	    p2 = Ap [n] ;
	    for (p = 0 ; p < p2 ; p++)
	    {
		/* Y [Ai [p]] += ABS (Ax [p]) ; */
		ASSIGN (aij, Ax [p], Az [p]) ;
		ABS (d, aij) ;
		Y [Ai [p]] += d ;
	    }

	    /* B2 = abs (B) */
	    flops += ABS_FLOPS * n ;
	    for (i = 0 ; i < n ; i++)
	    {
		/* B2 [i] = ABS (B [i]) ; */
		ASSIGN (bi, Bx [i], Bz [i]) ;
		ABS (B2 [i], bi) ;
	    }

	    /* scale Y and B2. */
	    if (do_scale)
	    {
		/* Y = R Y */
		/* B2 = R B2 */
#ifndef NRECIPROCAL
		if (do_recip)
		{
		    /* multiply by the scale factors */
		    for (i = 0 ; i < n ; i++)
		    {
			Y [i]  *= Rs [i] ;
			B2 [i] *= Rs [i] ;
		    }
		}
		else
#endif
		{
		    /* divide by the scale factors */
		    for (i = 0 ; i < n ; i++)
		    {
			Y [i]  /= Rs [i] ;
			B2 [i] /= Rs [i] ;
		    }
		}

		flops += 2 * n ;
	    }

	}

	for (step = 0 ; step <= irstep ; step++)
	{

	    /* -------------------------------------------------------------- */
	    /* Solve A x = b (step 0): */
	    /*  x = Q (U \ (L \ (P R b))) */
	    /* and then perform iterative refinement (step > 0): */
	    /*  x = x + Q (U \ (L \ (P R (b - A x)))) */
	    /* -------------------------------------------------------------- */

	    if (step == 0)
	    {
		if (do_scale)
		{
		    /* W = P R b, using X as workspace, since Z is not
		     * allocated if irstep = 0. */
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			/* multiply by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    ASSIGN (X [i], Bx [i], Bz [i]) ;
			    SCALE_RECIP (X [i], Rs [i]) ;
			}
		    }
		    else
#endif
		    {
			/* divide by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    ASSIGN (X [i], Bx [i], Bz [i]) ;
			    SCALE_DIV (X [i], Rs [i]) ;
			}
		    }
		    flops += SCALE_FLOPS * n ;
		    for (i = 0 ; i < n ; i++)
		    {
			W [i] = X [Rperm [i]] ;
		    }
		}
		else
		{
		    /* W = P b, since the row scaling R = I */
		    for (i = 0 ; i < n ; i++)
		    {
			/* W [i] = B [Rperm [i]] ; */
			ASSIGN (W [i], Bx [Rperm [i]], Bz [Rperm [i]]) ;
		    }
		}
	    }
	    else
	    {
		for (i = 0 ; i < n ; i++)
		{
		    /* Z [i] = B [i] ; */
		    ASSIGN (Z [i], Bx [i], Bz [i]) ;
		}
		flops += MULTSUB_FLOPS * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    xi = X [i] ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			/* Z [Ai [p]] -= Ax [p] * xi ; */
			ASSIGN (aij, Ax [p], Az [p]) ;
			MULT_SUB (Z [Ai [p]], aij, xi) ;
		    }
		}
		/* scale, Z = R Z */
		if (do_scale)
		{
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			/* multiply by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_RECIP (Z [i], Rs [i]) ;
			}
		    }
		    else
#endif
		    {
			/* divide by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_DIV (Z [i], Rs [i]) ;
			}
		    }
		    flops += SCALE_FLOPS * n ;
		}
		for (i = 0 ; i < n ; i++)
		{
		    W [i] = Z [Rperm [i]] ;
		}
	    }

	    flops += UMF_lsolve (Numeric, W, Pattern) ;
	    flops += UMF_usolve (Numeric, W, Pattern) ;

	    if (step == 0)
	    {
		for (i = 0 ; i < n ; i++)
		{
		    X [Cperm [i]] = W [i] ;
		}
	    }
	    else
	    {
		flops += ASSEMBLE_FLOPS * n ;
		for (i = 0 ; i < n ; i++)
		{
		    /* X [Cperm [i]] += W [i] ; */
		    ASSEMBLE (X [Cperm [i]], W [i]) ;
		}
	    }

	    /* -------------------------------------------------------------- */
	    /* sparse backward error estimate */
	    /* -------------------------------------------------------------- */

	    if (irstep > 0)
	    {

		/* ---------------------------------------------------------- */
		/* A is stored by column */
		/* W (i) = R (b - A x)_i, residual */
		/* Z2 (i) = R (|A||x|)_i */
		/* ---------------------------------------------------------- */

		for (i = 0 ; i < n ; i++)
		{
		    /* W [i] = B [i] ; */
		    ASSIGN (W [i], Bx [i], Bz [i]) ;
		    Z2 [i] = 0. ;
		}
		flops += (MULT_FLOPS + DECREMENT_FLOPS + ABS_FLOPS + 1) * nz ;
		for (j = 0 ; j < n ; j++)
		{
		    xj = X [j] ;
		    p2 = Ap [j+1] ;
		    for (p = Ap [j] ; p < p2 ; p++)
		    {
			i = Ai [p] ;

			/* axx = Ax [p] * xj ; */
			ASSIGN (aij, Ax [p], Az [p]) ;
			MULT (axx, aij, xj) ;

			/* W [i] -= axx ; */
			DECREMENT (W [i], axx) ;

			/* Z2 [i] += ABS (axx) ; */
			ABS (d, axx) ;
			Z2 [i] += d ;
		    }
		}

		/* scale W and Z2 */
		if (do_scale)
		{
		    /* Z2 = R Z2 */
		    /* W = R W */
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			/* multiply by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_RECIP (W [i], Rs [i]) ;
			    Z2 [i] *= Rs [i] ;
			}
		    }
		    else
#endif
		    {
			/* divide by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_DIV (W [i], Rs [i]) ;
			    Z2 [i] /= Rs [i] ;
			}
		    }
		    flops += (SCALE_FLOPS + 1) * n ;
		}

		flops += (2*ABS_FLOPS + 5) * n ;
		if (do_step (omega, step, B2, X, W, Y, Z2, S, n, Info))
		{
		    /* iterative refinement is done */
		    break ;
		}

	    }

	}

    }
    else if (sys == UMFPACK_At)
    {

	/* ------------------------------------------------------------------ */
	/* solve A' x = b with optional iterative refinement */
	/* ------------------------------------------------------------------ */

	/* A' is the complex conjugate transpose */

	if (irstep > 0)
	{

	    /* -------------------------------------------------------------- */
	    /* using iterative refinement:  compute Y */
	    /* -------------------------------------------------------------- */

	    nz = Ap [n] ;
	    Info [UMFPACK_NZ] = nz ;

	    /* A' is stored by row */
	    /* Y (i) = ||(A' R)_i||, 1-norm of row i of A' R */

	    if (do_scale)
	    {
		flops += (ABS_FLOPS + 2) * nz ;
#ifndef NRECIPROCAL
		if (do_recip)
		{
		    /* multiply by the scale factors */
		    for (i = 0 ; i < n ; i++)
		    {
			yi = 0. ;
			p2 = Ap [i+1] ;
			for (p = Ap [i] ; p < p2 ; p++)
			{
			    /* yi += ABS (Ax [p]) * Rs [Ai [p]] ; */
			    /* note that abs (aij) is the same as
			     * abs (conj (aij)) */
			    ASSIGN (aij, Ax [p], Az [p]) ;
			    ABS (d, aij) ;
			    yi += (d * Rs [Ai [p]]) ;
			}
			Y [i] = yi ;
		    }
		}
		else
#endif
		{
		    /* divide by the scale factors */
		    for (i = 0 ; i < n ; i++)
		    {
			yi = 0. ;
			p2 = Ap [i+1] ;
			for (p = Ap [i] ; p < p2 ; p++)
			{
			    /* yi += ABS (Ax [p]) / Rs [Ai [p]] ; */
			    /* note that abs (aij) is the same as
			     * abs (conj (aij)) */
			    ASSIGN (aij, Ax [p], Az [p]) ;
			    ABS (d, aij) ;
			    yi += (d / Rs [Ai [p]]) ;
			}
			Y [i] = yi ;
		    }
		}
	    }
	    else
	    {
		/* no scaling */
		flops += (ABS_FLOPS + 1) * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    yi = 0. ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			/* yi += ABS (Ax [p]) ; */
			/* note that abs (aij) is the same as
			 * abs (conj (aij)) */
			ASSIGN (aij, Ax [p], Az [p]) ;
			ABS (d, aij) ;
			yi += d ;
		    }
		    Y [i] = yi ;
		}
	    }

	    /* B2 = abs (B) */
	    for (i = 0 ; i < n ; i++)
	    {
		/* B2 [i] = ABS (B [i]) ; */
		ASSIGN (bi, Bx [i], Bz [i]) ;
		ABS (B2 [i], bi) ;
	    }

	}

	for (step = 0 ; step <= irstep ; step++)
	{

	    /* -------------------------------------------------------------- */
	    /* Solve A' x = b (step 0): */
	    /*	x = R P' (L' \ (U' \ (Q' b))) */
	    /* and then perform iterative refinement (step > 0): */
	    /*	x = x + R P' (L' \ (U' \ (Q' (b - A' x)))) */
	    /* -------------------------------------------------------------- */

	    if (step == 0)
	    {
		/* W = Q' b */
		for (i = 0 ; i < n ; i++)
		{
		    /* W [i] = B [Cperm [i]] ; */
		    ASSIGN (W [i], Bx [Cperm [i]], Bz [Cperm [i]]) ;
		}
	    }
	    else
	    {
		/* Z = b - A' x */
		for (i = 0 ; i < n ; i++)
		{
		    /* Z [i] = B [i] ; */
		    ASSIGN (Z [i], Bx [i], Bz [i]) ;
		}
		flops += MULTSUB_FLOPS * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    zi = Z [i] ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			/* zi -= conjugate (Ax [p]) * X [Ai [p]] ; */
			ASSIGN (aij, Ax [p], Az [p]) ;
			MULT_SUB_CONJ (zi, X [Ai [p]], aij) ;
		    }
		    Z [i] = zi ;
		}
		/* W = Q' Z */
		for (i = 0 ; i < n ; i++)
		{
		    W [i] = Z [Cperm [i]] ;
		}
	    }

	    flops += UMF_uhsolve (Numeric, W, Pattern) ;
	    flops += UMF_lhsolve (Numeric, W, Pattern) ;

	    if (step == 0)
	    {

		/* X = R P' W */
		/* do not use Z, since it isn't allocated if irstep = 0 */

		/* X = P' W */
		for (i = 0 ; i < n ; i++)
		{
		    X [Rperm [i]] = W [i] ;
		}
		if (do_scale)
		{
		    /* X = R X */
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			/* multiply by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_RECIP (X [i], Rs [i]) ;
			}
		    }
		    else
#endif
		    {
			/* divide by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_DIV (X [i], Rs [i]) ;
			}
		    }
		    flops += SCALE_FLOPS * n ;
		}

	    }
	    else
	    {

		/* Z = P' W */
		for (i = 0 ; i < n ; i++)
		{
		    Z [Rperm [i]] = W [i] ;
		}
		if (do_scale)
		{
		    /* Z = R Z */
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			/* multiply by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_RECIP (Z [i], Rs [i]) ;
			}
		    }
		    else
#endif
		    {
			/* divide by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_DIV (Z [i], Rs [i]) ;
			}
		    }
		    flops += SCALE_FLOPS * n ;
		}

		flops += ASSEMBLE_FLOPS * n ;
		/* X += Z */
		for (i = 0 ; i < n ; i++)
		{
		    /* X [i] += W [i] ; */
		    ASSEMBLE (X [i], W [i]) ;
		}
	    }

	    /* -------------------------------------------------------------- */
	    /* sparse backward error estimate */
	    /* -------------------------------------------------------------- */

	    if (irstep > 0)
	    {

		/* ---------------------------------------------------------- */
		/* A' is stored by row */
		/* W (i) = (b - A' x)_i, residual */
		/* Z2 (i) = (|A'||x|)_i */
		/* ---------------------------------------------------------- */

		flops += (MULT_FLOPS + DECREMENT_FLOPS + ABS_FLOPS + 1) * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    /* wi = B [i] ; */
		    ASSIGN (wi, Bx [i], Bz [i]) ;
		    z2i = 0. ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			/* axx = conjugate (Ax [p]) * X [Ai [p]] ; */
			ASSIGN (aij, Ax [p], Az [p]) ;
			MULT_CONJ (axx, X [Ai [p]], aij) ;

			/* wi -= axx ; */
			DECREMENT (wi, axx) ;

			/* z2i += ABS (axx) ; */
			ABS (d, axx) ;
			z2i += d ;
		    }
		    W [i] = wi ;
		    Z2 [i] = z2i ;
		}

		flops += (2*ABS_FLOPS + 5) * n ;
		if (do_step (omega, step, B2, X, W, Y, Z2, S, n, Info))
		{
		    /* iterative refinement is done */
		    break ;
		}

	    }

	}

    }
    else if (sys == UMFPACK_Aat)
    {

	/* ------------------------------------------------------------------ */
	/* solve A.' x = b with optional iterative refinement */
	/* ------------------------------------------------------------------ */

	/* A' is the array transpose */

	if (irstep > 0)
	{

	    /* -------------------------------------------------------------- */
	    /* using iterative refinement:  compute Y */
	    /* -------------------------------------------------------------- */

	    nz = Ap [n] ;
	    Info [UMFPACK_NZ] = nz ;

	    /* A.' is stored by row */
	    /* Y (i) = ||(A.' R)_i||, 1-norm of row i of A.' R */

	    if (do_scale)
	    {
		flops += (ABS_FLOPS + 2) * nz ;
#ifndef NRECIPROCAL
		if (do_recip)
		{
		    /* multiply by the scale factors */
		    for (i = 0 ; i < n ; i++)
		    {
			yi = 0. ;
			p2 = Ap [i+1] ;
			for (p = Ap [i] ; p < p2 ; p++)
			{
			    /* yi += ABS (Ax [p]) * Rs [Ai [p]] ; */
			    /* note that A.' is the array transpose,
			     * so no conjugate */
			    ASSIGN (aij, Ax [p], Az [p]) ;
			    ABS (d, aij) ;
			    yi += (d * Rs [Ai [p]]) ;
			}
			Y [i] = yi ;
		    }
		}
		else
#endif
		{
		    /* divide by the scale factors */
		    for (i = 0 ; i < n ; i++)
		    {
			yi = 0. ;
			p2 = Ap [i+1] ;
			for (p = Ap [i] ; p < p2 ; p++)
			{
			    /* yi += ABS (Ax [p]) / Rs [Ai [p]] ; */
			    /* note that A.' is the array transpose,
			     * so no conjugate */
			    ASSIGN (aij, Ax [p], Az [p]) ;
			    ABS (d, aij) ;
			    yi += (d / Rs [Ai [p]]) ;
			}
			Y [i] = yi ;
		    }
		}
	    }
	    else
	    {
		/* no scaling */
		flops += (ABS_FLOPS + 1) * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    yi = 0. ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			/* yi += ABS (Ax [p]) */
			/* note that A.' is the array transpose,
			 * so no conjugate */
			ASSIGN (aij, Ax [p], Az [p]) ;
			ABS (d, aij) ;
			yi += d ;
		    }
		    Y [i] = yi ;
		}
	    }

	    /* B2 = abs (B) */
	    for (i = 0 ; i < n ; i++)
	    {
		/* B2 [i] = ABS (B [i]) ; */
		ASSIGN (bi, Bx [i], Bz [i]) ;
		ABS (B2 [i], bi) ;
	    }

	}

	for (step = 0 ; step <= irstep ; step++)
	{

	    /* -------------------------------------------------------------- */
	    /* Solve A.' x = b (step 0): */
	    /*	x = R P' (L.' \ (U.' \ (Q' b))) */
	    /* and then perform iterative refinement (step > 0): */
	    /*	x = x + R P' (L.' \ (U.' \ (Q' (b - A.' x)))) */
	    /* -------------------------------------------------------------- */

	    if (step == 0)
	    {
		/* W = Q' b */
		for (i = 0 ; i < n ; i++)
		{
		    /* W [i] = B [Cperm [i]] ; */
		    ASSIGN (W [i], Bx [Cperm [i]], Bz [Cperm [i]]) ;
		}
	    }
	    else
	    {
		/* Z = b - A.' x */
		for (i = 0 ; i < n ; i++)
		{
		    /* Z [i] = B [i] ; */
		    ASSIGN (Z [i], Bx [i], Bz [i]) ;
		}
		flops += MULTSUB_FLOPS * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    zi = Z [i] ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			/* zi -= Ax [p] * X [Ai [p]] ; */
			ASSIGN (aij, Ax [p], Az [p]) ;
			MULT_SUB (zi, aij, X [Ai [p]]) ;
		    }
		    Z [i] = zi ;
		}
		/* W = Q' Z */
		for (i = 0 ; i < n ; i++)
		{
		    W [i] = Z [Cperm [i]] ;
		}
	    }

	    flops += UMF_utsolve (Numeric, W, Pattern) ;
	    flops += UMF_ltsolve (Numeric, W, Pattern) ;

	    if (step == 0)
	    {

		/* X = R P' W */
		/* do not use Z, since it isn't allocated if irstep = 0 */

		/* X = P' W */
		for (i = 0 ; i < n ; i++)
		{
		    X [Rperm [i]] = W [i] ;
		}
		if (do_scale)
		{
		    /* X = R X */
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			/* multiply by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_RECIP (X [i], Rs [i]) ;
			}
		    }
		    else
#endif
		    {
			/* divide by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_DIV (X [i], Rs [i]) ;
			}
		    }
		    flops += SCALE_FLOPS * n ;
		}

	    }
	    else
	    {

		/* Z = P' W */
		for (i = 0 ; i < n ; i++)
		{
		    Z [Rperm [i]] = W [i] ;
		}
		if (do_scale)
		{
		    /* Z = R Z */
#ifndef NRECIPROCAL
		    if (do_recip)
		    {
			/* multiply by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_RECIP (Z [i], Rs [i]) ;
			}
		    }
		    else
#endif
		    {
			/* divide by the scale factors */
			for (i = 0 ; i < n ; i++)
			{
			    SCALE_DIV (Z [i], Rs [i]) ;
			}
		    }
		    flops += SCALE_FLOPS * n ;
		}

		flops += ASSEMBLE_FLOPS * n ;
		/* X += Z */
		for (i = 0 ; i < n ; i++)
		{
		    /* X [i] += W [i] ; */
		    ASSEMBLE (X [i], W [i]) ;
		}
	    }

	    /* -------------------------------------------------------------- */
	    /* sparse backward error estimate */
	    /* -------------------------------------------------------------- */

	    if (irstep > 0)
	    {

		/* ---------------------------------------------------------- */
		/* A.' is stored by row */
		/* W (i) = (b - A.' x)_i, residual */
		/* Z (i) = (|A.'||x|)_i */
		/* ---------------------------------------------------------- */

		flops += (MULT_FLOPS + DECREMENT_FLOPS + ABS_FLOPS + 1) * nz ;
		for (i = 0 ; i < n ; i++)
		{
		    /* wi = B [i] ; */
		    ASSIGN (wi, Bx [i], Bz [i]) ;
		    z2i = 0. ;
		    p2 = Ap [i+1] ;
		    for (p = Ap [i] ; p < p2 ; p++)
		    {
			/* axx = Ax [p] * X [Ai [p]] ; */
			ASSIGN (aij, Ax [p], Az [p]) ;
			MULT (axx, aij, X [Ai [p]]) ;

			/* wi -= axx ; */
			DECREMENT (wi, axx) ;

			/* z2i += ABS (axx) ; */
			ABS (d, axx) ;
			z2i += d ;
		    }
		    W [i] = wi ;
		    Z2 [i] = z2i ;
		}

		flops += (2*ABS_FLOPS + 5) * n ;
		if (do_step (omega, step, B2, X, W, Y, Z2, S, n, Info))
		{
		    /* iterative refinement is done */
		    break ;
		}

	    }

	}

    }
    else if (sys == UMFPACK_Pt_L)
    {

	/* ------------------------------------------------------------------ */
	/* Solve P'Lx=b:  x = L \ Pb */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < n ; i++)
	{
	    /* X [i] = B [Rperm [i]] ; */
	    ASSIGN (X [i], Bx [Rperm [i]], Bz [Rperm [i]]) ;
	}
	flops = UMF_lsolve (Numeric, X, Pattern) ;
	status = UMFPACK_OK ;

    }
    else if (sys == UMFPACK_L)
    {

	/* ------------------------------------------------------------------ */
	/* Solve Lx=b:  x = L \ b */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < n ; i++)
	{
	    /* X [i] = B [i] ; */
	    ASSIGN (X [i], Bx [i], Bz [i]) ;
	}
	flops = UMF_lsolve (Numeric, X, Pattern) ;
	status = UMFPACK_OK ;

    }
    else if (sys == UMFPACK_Lt_P)
    {

	/* ------------------------------------------------------------------ */
	/* Solve L'Px=b:  x = P' (L' \ b) */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < n ; i++)
	{
	    /* W [i] = B [i] ; */
	    ASSIGN (W [i], Bx [i], Bz [i]) ;
	}
	flops = UMF_lhsolve (Numeric, W, Pattern) ;
	for (i = 0 ; i < n ; i++)
	{
	    X [Rperm [i]] = W [i] ;
	}
	status = UMFPACK_OK ;

    }
    else if (sys == UMFPACK_Lat_P)
    {

	/* ------------------------------------------------------------------ */
	/* Solve L.'Px=b:  x = P' (L.' \ b) */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < n ; i++)
	{
	    /* W [i] = B [i] ; */
	    ASSIGN (W [i], Bx [i], Bz [i]) ;
	}
	flops = UMF_ltsolve (Numeric, W, Pattern) ;
	for (i = 0 ; i < n ; i++)
	{
	    X [Rperm [i]] = W [i] ;
	}
	status = UMFPACK_OK ;

    }
    else if (sys == UMFPACK_Lt)
    {

	/* ------------------------------------------------------------------ */
	/* Solve L'x=b:  x = L' \ b */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < n ; i++)
	{
	    /* X [i] = B [i] ; */
	    ASSIGN (X [i], Bx [i], Bz [i]) ;
	}
	flops = UMF_lhsolve (Numeric, X, Pattern) ;
	status = UMFPACK_OK ;

    }
    else if (sys == UMFPACK_Lat)
    {

	/* ------------------------------------------------------------------ */
	/* Solve L.'x=b:  x = L.' \ b */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < n ; i++)
	{
	    /* X [i] = B [i] ; */
	    ASSIGN (X [i], Bx [i], Bz [i]) ;
	}
	flops = UMF_ltsolve (Numeric, X, Pattern) ;
	status = UMFPACK_OK ;

    }
    else if (sys == UMFPACK_U_Qt)
    {

	/* ------------------------------------------------------------------ */
	/* Solve UQ'x=b:  x = Q (U \ b) */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < n ; i++)
	{
	    /* W [i] = B [i] ; */
	    ASSIGN (W [i], Bx [i], Bz [i]) ;
	}
	flops = UMF_usolve (Numeric, W, Pattern) ;
	for (i = 0 ; i < n ; i++)
	{
	    X [Cperm [i]] = W [i] ;
	}

    }
    else if (sys == UMFPACK_U)
    {

	/* ------------------------------------------------------------------ */
	/* Solve Ux=b:  x = U \ b */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < n ; i++)
	{
	    /* X [i] = B [i] ; */
	    ASSIGN (X [i], Bx [i], Bz [i]) ;
	}
	flops = UMF_usolve (Numeric, X, Pattern) ;

    }
    else if (sys == UMFPACK_Q_Ut)
    {

	/* ------------------------------------------------------------------ */
	/* Solve QU'x=b:  x = U' \ Q'b */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < n ; i++)
	{
	    /* X [i] = B [Cperm [i]] ; */
	    ASSIGN (X [i], Bx [Cperm [i]], Bz [Cperm [i]]) ;
	}
	flops = UMF_uhsolve (Numeric, X, Pattern) ;

    }
    else if (sys == UMFPACK_Q_Uat)
    {

	/* ------------------------------------------------------------------ */
	/* Solve QU.'x=b:  x = U.' \ Q'b */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < n ; i++)
	{
	    /* X [i] = B [Cperm [i]] ; */
	    ASSIGN (X [i], Bx [Cperm [i]], Bz [Cperm [i]]) ;
	}
	flops = UMF_utsolve (Numeric, X, Pattern) ;

    }
    else if (sys == UMFPACK_Ut)
    {

	/* ------------------------------------------------------------------ */
	/* Solve U'x=b:  x = U' \ b */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < n ; i++)
	{
	    /* X [i] = B [i] ; */
	    ASSIGN (X [i], Bx [i], Bz [i]) ;
	}
	flops = UMF_uhsolve (Numeric, X, Pattern) ;

    }
    else if (sys == UMFPACK_Uat)
    {

	/* ------------------------------------------------------------------ */
	/* Solve U'x=b:  x = U' \ b */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < n ; i++)
	{
	    /* X [i] = B [i] ; */
	    ASSIGN (X [i], Bx [i], Bz [i]) ;
	}
	flops = UMF_utsolve (Numeric, X, Pattern) ;

    }
    else
    {
	return (UMFPACK_ERROR_invalid_system) ;
    }

#ifdef COMPLEX
    /* copy the solution back, from Entry X [ ] to double Xx [ ] and Xz [ ] */
    for (i = 0 ; i < n ; i++)
    {
	Xx [i] = REAL_COMPONENT (X [i]) ;
	Xz [i] = IMAG_COMPONENT (X [i]) ;
    }
#endif

    /* return UMFPACK_OK, or UMFPACK_WARNING_singular_matrix */
    /* Note that systems involving just L will return UMFPACK_OK */
    Info [UMFPACK_SOLVE_FLOPS] = flops ;
    return (status) ;
}


/* ========================================================================== */
/* === do_step ============================================================== */
/* ========================================================================== */

/* Perform one step of iterative refinement, for A x = b or A' x = b */

PRIVATE Int do_step		/* return TRUE if iterative refinement done */
(
    double omega [3],
    Int step,			/* which step of iterative refinement to do */
    const double B2 [ ],	/* abs (B) */
    Entry X [ ],
    const Entry W [ ],
    const double Y [ ],
    const double Z2 [ ],
    Entry S [ ],
    Int n,
    double Info [UMFPACK_INFO]
)
{
    double last_omega [3], tau, nctau, d1, wd1, d2, wd2, xi, yix, wi, xnorm ;
    Int i ;

    /* DBL_EPSILON is a standard ANSI C term defined in <float.h> */
    /* It is the smallest positive x such that 1.0+x != 1.0 */

    nctau = 1000 * n * DBL_EPSILON ;
    DEBUG0 (("do_step start: nctau = %30.20e\n", nctau)) ;
    ASSERT (UMF_report_vector (n, (double *) X, (double *) NULL, UMF_debug,
	FALSE, FALSE) == UMFPACK_OK) ;

    /* for approximate flop count, assume d1 > tau is always true */
    /* flops += (2*ABS_FLOPS + 5) * n ; (done in UMF_solve, above) */

    /* ---------------------------------------------------------------------- */
    /* save the last iteration in case we need to reinstate it */
    /* ---------------------------------------------------------------------- */

    last_omega [0] = omega [0] ;
    last_omega [1] = omega [1] ;
    last_omega [2] = omega [2] ;

    /* ---------------------------------------------------------------------- */
    /* compute sparse backward errors: omega [1] and omega [2] */
    /* ---------------------------------------------------------------------- */

    /* xnorm = ||x|| maxnorm */
    xnorm = 0.0 ;
    for (i = 0 ; i < n ; i++)
    {
	/* xi = ABS (X [i]) ; */
	ABS (xi, X [i]) ;
	if (SCALAR_IS_NAN (xi))
	{
	    xnorm = xi ;
	    break ;
	}
	/* no NaN's to consider here: */
	xnorm = MAX (xnorm, xi) ;
    }

    omega [1] = 0. ;
    omega [2] = 0. ;
    for (i = 0 ; i < n ; i++)
    {
	yix = Y [i] * xnorm ;
	tau = (yix + B2 [i]) * nctau ;
	d1 = Z2 [i] + B2 [i] ;
	/* wi = ABS (W [i]) ; */
	ABS (wi, W [i]) ;
	if (SCALAR_IS_NAN (d1))
	{
	    omega [1] = d1 ;
	    omega [2] = d1 ;
	    break ;
	}
	if (SCALAR_IS_NAN (tau))
	{
	    omega [1] = tau ;
	    omega [2] = tau ;
	    break ;
	}
	if (d1 > tau)		/* a double relop, but no NaN's here */
	{
	    wd1 = wi / d1 ;
	    omega [1] = MAX (omega [1], wd1) ;
	}
	else if (tau > 0.0)	/* a double relop, but no NaN's here */
	{
	    d2 = Z2 [i] + yix ;
	    wd2 = wi / d2 ;
	    omega [2] = MAX (omega [2], wd2) ;
	}
    }

    omega [0] = omega [1] + omega [2] ;
    Info [UMFPACK_OMEGA1] = omega [1] ;
    Info [UMFPACK_OMEGA2] = omega [2] ;

    /* ---------------------------------------------------------------------- */
    /* stop the iterations if the backward error is small, or NaN */
    /* ---------------------------------------------------------------------- */

    Info [UMFPACK_IR_TAKEN] = step ;
    Info [UMFPACK_IR_ATTEMPTED] = step ;

    if (SCALAR_IS_NAN (omega [0]))
    {
	DEBUG0 (("omega[0] is NaN - done.\n")) ;
	ASSERT (UMF_report_vector (n, (double *) X, (double *) NULL, UMF_debug,
	    FALSE, FALSE) == UMFPACK_OK) ;
	return (TRUE) ;
    }

    if (omega [0] < DBL_EPSILON)    /* double relop, but no NaN case here */
    {
	DEBUG0 (("omega[0] too small - done.\n")) ;
	ASSERT (UMF_report_vector (n, (double *) X, (double *) NULL, UMF_debug,
	    FALSE, FALSE) == UMFPACK_OK) ;
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* stop if insufficient decrease in omega */
    /* ---------------------------------------------------------------------- */

    /* double relop, but no NaN case here: */
    if (step > 0 && omega [0] > last_omega [0] / 2)
    {
	DEBUG0 (("stop refinement\n")) ;
	if (omega [0] > last_omega [0])
	{
	    /* last iteration better than this one, reinstate it */
	    DEBUG0 (("last iteration better\n")) ;
	    for (i = 0 ; i < n ; i++)
	    {
		X [i] = S [i] ;
	    }
	    Info [UMFPACK_OMEGA1] = last_omega [1] ;
	    Info [UMFPACK_OMEGA2] = last_omega [2] ;
	}
	Info [UMFPACK_IR_TAKEN] = step - 1 ;
	ASSERT (UMF_report_vector (n, (double *) X, (double *) NULL, UMF_debug,
	    FALSE, FALSE) == UMFPACK_OK) ;
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* save current solution in case we need to reinstate */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < n ; i++)
    {
	S [i] = X [i] ;
    }

    /* ---------------------------------------------------------------------- */
    /* iterative refinement continues */
    /* ---------------------------------------------------------------------- */

    ASSERT (UMF_report_vector (n, (double *) X, (double *) NULL, UMF_debug,
	FALSE, FALSE) == UMFPACK_OK) ;
    return (FALSE) ;
}
