/* ========================================================================== */
/* === UMFPACK_scale ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Applies the scale factors computed during numerical
    factorization to a vector. See umfpack_scale.h for more details.

    The LU factorization is L*U = P*R*A*Q, where P and Q are permutation
    matrices, and R is diagonal.  This routine computes X = R * B using the
    matrix R stored in the Numeric object.

    Returns FALSE if any argument is invalid, TRUE otherwise.

    If R not present in the Numeric object, then R = I and no floating-point
    work is done.  B is simply copied into X.
*/

#include "umf_internal.h"
#include "umf_valid_numeric.h"

GLOBAL Int UMFPACK_scale
(
    double Xx [ ],
#ifdef COMPLEX
    double Xz [ ],
#endif
    const double Bx [ ],
#ifdef COMPLEX
    const double Bz [ ],
#endif
    void *NumericHandle
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    NumericType *Numeric ;
    Int n, i ;
    double *Rs ;

    Numeric = (NumericType *) NumericHandle ;
    if (!UMF_valid_numeric (Numeric))
    {
	return (UMFPACK_ERROR_invalid_Numeric_object) ;
    }

    n = Numeric->n_row ;
    Rs = Numeric->Rs ;

    if (!Xx || !Bx
#ifdef COMPLEX
	|| !Xz || !Bz
#endif
    )
    {
	return (UMFPACK_ERROR_argument_missing) ;
    }

    /* ---------------------------------------------------------------------- */
    /* X = R*B or R\B */
    /* ---------------------------------------------------------------------- */

    if (Rs != (double *) NULL)
    {
#ifndef NRECIPROCAL
	if (Numeric->do_recip)
	{
	    /* multiply by the scale factors */
	    for (i = 0 ; i < n ; i++)
	    {
		Xx [i] = Bx [i] * Rs [i] ;
#ifdef COMPLEX
		Xz [i] = Bz [i] * Rs [i] ;
#endif
	    }
	}
	else
#endif
	{
	    /* divide by the scale factors */
	    for (i = 0 ; i < n ; i++)
	    {
		Xx [i] = Bx [i] / Rs [i] ;
#ifdef COMPLEX
		Xz [i] = Bz [i] / Rs [i] ;
#endif
	    }
	}
    }
    else
    {
	/* no scale factors, just copy B into X */
	for (i = 0 ; i < n ; i++)
	{
	    Xx [i] = Bx [i] ;
#ifdef COMPLEX
	    Xz [i] = Bz [i] ;
#endif
	}
    }

    return (UMFPACK_OK) ;
}
