/* ========================================================================== */
/* === UMF_scale ============================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Divide a vector of stride 1 by the pivot value. */

#include "umf_internal.h"

GLOBAL void UMF_scale
(
    Int n,
    Entry pivot,
    Entry X [ ]
)
{
    Int i ;
    Entry x ;
    double s ;

    /* ---------------------------------------------------------------------- */
    /* compute the approximate absolute value of the pivot, and select method */
    /* ---------------------------------------------------------------------- */

    APPROX_ABS (s, pivot) ;

    if (s < RECIPROCAL_TOLERANCE || IS_NAN (pivot))
    {
	/* ------------------------------------------------------------------ */
	/* tiny, or zero, pivot case */
	/* ------------------------------------------------------------------ */

	/* The pivot is tiny, or NaN.  Do not divide zero by the pivot value,
	 * and do not multiply by 1/pivot, either. */

	for (i = 0 ; i < n ; i++)
	{
	    /* X [i] /= pivot ; */
	    x = X [i] ;
	    if (IS_NONZERO (x))
	    {
		DIV (X [i], x, pivot) ;
	    }
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* normal case.  select the x/pivot or x * (1/pivot) method */
	/* ------------------------------------------------------------------ */

	/* The pivot is not tiny, and is not NaN.   Don't bother to check for
	 * zeros in the pivot column, X. */

#if !defined (NRECIPROCAL) && !(defined (__GNUC__) && defined (COMPLEX))

	    /* -------------------------------------------------------------- */
	    /* multiply x by (1/pivot) */
	    /* -------------------------------------------------------------- */

	    /* Slightly less accurate, but faster.  It allows the use of
	     * the level-1 BLAS dscal or zscal routine.  This not used when
	     * UMFPACK is used in MATLAB (either as a built-in routine, or as
	     * a mexFunction).
	     *
	     * Using gcc version 3.2 can cause the following code to fail for
	     * some complex matrices (not all), with or without the BLAS.  This
	     * was found in Red Hat Linux 7.3 on a Dell Latitude C840 with a
	     * Pentium 4M.  Thus, this code is not used when gcc is used, for
	     * the complex case.
	     *
	     * It works just fine with Intel's icc compiler, version 7.0.
	     */

	    /* pivot = 1 / pivot */
	    RECIPROCAL (pivot) ;

#if defined (USE_NO_BLAS)
	    for (i = 0 ; i < n ; i++)
	    {
		/* X [i] *= pivot ; */
		x = X [i] ;
		MULT (X [i], x, pivot) ;
	    }
#else
	    BLAS_SCAL (n, pivot, X) ;
#endif

#else

	    /* -------------------------------------------------------------- */
	    /* divide x by the pivot */
	    /* -------------------------------------------------------------- */

	    /* This is slightly more accurate, particularly if the pivot column
	     * consists of only IEEE subnormals.  Always do this if UMFPACK is
	     * being compiled as a built-in routine or mexFunction in MATLAB,
	     * or if gcc is being used with complex matrices. */

	    for (i = 0 ; i < n ; i++)
	    {
		/* X [i] /= pivot ; */
		x = X [i] ;
		DIV (X [i], x, pivot) ;
	    }

#endif

    }
}
