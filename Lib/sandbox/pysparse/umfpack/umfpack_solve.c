/* ========================================================================== */
/* === UMFPACK_solve ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Solves a linear system using the numerical factorization
    computed by UMFPACK_numeric.  See umfpack_solve.h for more details.

    For umfpack_*_solve:
	Dynamic memory usage:  UMFPACK_solve calls UMF_malloc twice, for
	workspace of size c*n*sizeof(double) + n*sizeof(Int), where c is
	defined below.  On return, all of this workspace is free'd via UMF_free.

    For umfpack_*_wsolve:
	No dynamic memory usage.  Input arrays are used for workspace instead.
	Pattern is a workspace of size n Integers.  The double array W must be
	at least of size c*n, where c is defined below.

    If iterative refinement is requested, and Ax=b, A'x=b or A.'x=b is being
    solved, and the matrix A is not singular, then c is 5 for the real version
    and 10 for the complex version.  Otherwise, c is 1 for the real version and
    4 for the complex version.
*/

#include "umf_internal.h"
#include "umf_valid_numeric.h"
#include "umf_solve.h"

#ifndef WSOLVE
#include "umf_malloc.h"
#include "umf_free.h"
#ifndef NDEBUG
PRIVATE Int init_count ;
#endif
#endif

GLOBAL Int
#ifdef WSOLVE
UMFPACK_wsolve
#else
UMFPACK_solve
#endif
(
    Int sys,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    double Xx [ ],
#ifdef COMPLEX
    double Xz [ ],
#endif
    const double Bx [ ],
#ifdef COMPLEX
    const double Bz [ ],
#endif
    void *NumericHandle,
    const double Control [UMFPACK_CONTROL],
    double User_Info [UMFPACK_INFO]
#ifdef WSOLVE
    , Int Pattern [ ],
    double W [ ]
#endif
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    NumericType *Numeric ;
    Int n, i, irstep, status ;
    double Info2 [UMFPACK_INFO], *Info, /* tstart, tend */ stats [2] ;
#ifndef WSOLVE
    Int *Pattern, wsize ;
    double *W ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get the amount of time used by the process so far */
    /* ---------------------------------------------------------------------- */

    umfpack_tic (stats) ;

#ifndef WSOLVE
#ifndef NDEBUG
    init_count = UMF_malloc_count ;
#endif
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    irstep = GET_CONTROL (UMFPACK_IRSTEP, UMFPACK_DEFAULT_IRSTEP) ;

    if (User_Info != (double *) NULL)
    {
	/* return Info in user's array */
	Info = User_Info ;
	/* clear the parts of Info that are set by UMFPACK_solve */
	for (i = UMFPACK_IR_TAKEN ; i <= UMFPACK_SOLVE_TIME ; i++)
	{
	    Info [i] = EMPTY ;
	}
    }
    else
    {
	/* no Info array passed - use local one instead */
	Info = Info2 ;
	for (i = 0 ; i < UMFPACK_INFO ; i++)
	{
	    Info [i] = EMPTY ;
	}
    }

    Info [UMFPACK_STATUS] = UMFPACK_OK ;
    Info [UMFPACK_SOLVE_FLOPS] = 0 ;

    Numeric = (NumericType *) NumericHandle ;
    if (!UMF_valid_numeric (Numeric))
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_invalid_Numeric_object ;
	return (UMFPACK_ERROR_invalid_Numeric_object) ;
    }

    Info [UMFPACK_NROW] = Numeric->n_row ;
    Info [UMFPACK_NCOL] = Numeric->n_col ;

    if (Numeric->n_row != Numeric->n_col)
    {
	/* only square systems can be handled */
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_invalid_system ;
	return (UMFPACK_ERROR_invalid_system) ;
    }
    n = Numeric->n_row ;
    if (Numeric->nnzpiv < n
	|| SCALAR_IS_ZERO (Numeric->rcond) || SCALAR_IS_NAN (Numeric->rcond))
    {
	/* turn off iterative refinement if A is singular */
	/* or if U has NaN's on the diagonal. */
	irstep = 0 ;
    }

    if (!Xx || !Bx
#ifdef COMPLEX
	|| !Xz || !Bz
#endif
    )
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_argument_missing ;
	return (UMFPACK_ERROR_argument_missing) ;
    }

    if (sys >= UMFPACK_Pt_L)
    {
	/* no iterative refinement except for nonsingular Ax=b, A'x=b, A.'x=b */
	irstep = 0 ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate or check the workspace */
    /* ---------------------------------------------------------------------- */

#ifdef WSOLVE

    if (!W || !Pattern)
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_argument_missing ;
	return (UMFPACK_ERROR_argument_missing) ;
    }

#else

#ifdef COMPLEX
    if (irstep > 0)
    {
	wsize = 10*n ;		/* W, X, Z, S, Y, B2 */
    }
    else
    {
	wsize = 4*n ;		/* W, X */
    }
#else
    if (irstep > 0)
    {
	wsize = 5*n ;		/* W, Z, S, Y, B2 */
    }
    else
    {
	wsize = n ;		/* W */
    }
#endif

    Pattern = (Int *) UMF_malloc (n, sizeof (Int)) ;
    W = (double *) UMF_malloc (wsize, sizeof (double)) ;
    if (!W || !Pattern)
    {
	DEBUGm4 (("out of memory: solve work\n")) ;
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	(void) UMF_free ((void *) W) ;
	(void) UMF_free ((void *) Pattern) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

#endif	/* WSOLVE */

    /* ---------------------------------------------------------------------- */
    /* solve the system */
    /* ---------------------------------------------------------------------- */

    status = UMF_solve (sys, Ap, Ai, Ax, Xx, Bx,
#ifdef COMPLEX
	Az, Xz, Bz,
#endif
	Numeric, irstep, Info, Pattern, W) ;

    /* ---------------------------------------------------------------------- */
    /* free the workspace (if allocated) */
    /* ---------------------------------------------------------------------- */

#ifndef WSOLVE
    (void) UMF_free ((void *) W) ;
    (void) UMF_free ((void *) Pattern) ;
    ASSERT (UMF_malloc_count == init_count) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get the time used by UMFPACK_*solve */
    /* ---------------------------------------------------------------------- */

    Info [UMFPACK_STATUS] = status ;
    if (status >= 0)
    {
	umfpack_toc (stats) ;
	Info [UMFPACK_SOLVE_WALLTIME] = stats [0] ;
	Info [UMFPACK_SOLVE_TIME] = stats [1] ;
    }

    return (status) ;
}
