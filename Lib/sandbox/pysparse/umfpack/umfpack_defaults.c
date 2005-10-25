/* ========================================================================== */
/* === UMFPACK_defaults ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Sets default control parameters.  See umfpack_defaults.h
    for details.
*/

#include "umf_internal.h"

GLOBAL void UMFPACK_defaults
(
    double Control [UMFPACK_CONTROL]
)
{
    Int i ;

    if (!Control)
    {
	/* silently return if no Control array */
	return ;
    }

    for (i = 0 ; i < UMFPACK_CONTROL ; i++)
    {
	Control [i] = 0 ;
    }

    /* ---------------------------------------------------------------------- */
    /* default control settings: can be modified at run-time */
    /* ---------------------------------------------------------------------- */

    /* used in UMFPACK_report_* routines: */
    Control [UMFPACK_PRL] = UMFPACK_DEFAULT_PRL ;

    /* used in UMFPACK_*symbolic: */
    Control [UMFPACK_DENSE_ROW] = UMFPACK_DEFAULT_DENSE_ROW ;
    Control [UMFPACK_DENSE_COL] = UMFPACK_DEFAULT_DENSE_COL ;
    Control [UMFPACK_AMD_DENSE] = UMFPACK_DEFAULT_AMD_DENSE ;
    Control [UMFPACK_STRATEGY] = UMFPACK_DEFAULT_STRATEGY ;
    Control [UMFPACK_2BY2_TOLERANCE] = UMFPACK_DEFAULT_2BY2_TOLERANCE ;
    Control [UMFPACK_AGGRESSIVE] = UMFPACK_DEFAULT_AGGRESSIVE ;

    /* used in UMFPACK_numeric: */
    Control [UMFPACK_PIVOT_TOLERANCE] = UMFPACK_DEFAULT_PIVOT_TOLERANCE ;
    Control [UMFPACK_SYM_PIVOT_TOLERANCE] = UMFPACK_DEFAULT_SYM_PIVOT_TOLERANCE;
    Control [UMFPACK_BLOCK_SIZE] = UMFPACK_DEFAULT_BLOCK_SIZE ;
    Control [UMFPACK_ALLOC_INIT] = UMFPACK_DEFAULT_ALLOC_INIT ;
    Control [UMFPACK_FRONT_ALLOC_INIT] = UMFPACK_DEFAULT_FRONT_ALLOC_INIT ;
    Control [UMFPACK_SCALE] = UMFPACK_DEFAULT_SCALE ;

    /* used in UMFPACK_*solve: */
    Control [UMFPACK_IRSTEP] = UMFPACK_DEFAULT_IRSTEP ;

    /* ---------------------------------------------------------------------- */
    /* compile-time settings: cannot be modified at run-time */
    /* ---------------------------------------------------------------------- */

#ifdef USE_NO_BLAS
    /* do not use the BLAS - use in-line C code instead */
    Control [UMFPACK_COMPILED_WITH_BLAS] = 0 ;
#else
    /* use externally-provided BLAS (dgemm, dger, dgemv, zgemm, zgeru, zgemv) */
    Control [UMFPACK_COMPILED_WITH_BLAS] = 1 ;
#endif

#ifdef MATLAB_MEX_FILE
    /* use mxMalloc, mxFree, mxRealloc, and mexPrintf */
    /* use mxAssert if debugging is enabled */
    Control [UMFPACK_COMPILED_FOR_MATLAB] = 1 ;
#else
#ifdef MATHWORKS
    /* use internal utMalloc, utFree, utRealloc, and utPrintf routines. */
    /* use utDivideComplex and utFdlibm_hypot for complex version. */
    /* use utAssert if debugging is enabled. */
    Control [UMFPACK_COMPILED_FOR_MATLAB] = 2 ;
#else
    /* use ANSI C malloc, free, realloc, and print */
    /* use ANSI C assert if debugging is enabled */
    Control [UMFPACK_COMPILED_FOR_MATLAB] = 0 ;
#endif
#endif

#ifndef NPOSIX
    /* uses the POSIX sysconf ( ) and times ( ) routines in UMFPACK_tic, toc */
    Control [UMFPACK_COMPILED_WITH_GETRUSAGE] = 2 ;
#else
#ifdef GETRUSAGE
    /* uses the non-standard getrusage to get CPU time (Solaris) */
    Control [UMFPACK_COMPILED_WITH_GETRUSAGE] = 1 ;
#else
    /* uses the ANSI standard clock routine to get CPU time */
    /* this may wrap around */
    Control [UMFPACK_COMPILED_WITH_GETRUSAGE] = 0 ;
#endif
#endif

#ifndef NDEBUG
    /* UMFPACK is compiled in debug mode. */
    /* This is exceedingly slow. */
    DEBUG0 (("UMFPACK is running in debug mode.  This is very slow!\n")) ;
    Control [UMFPACK_COMPILED_IN_DEBUG_MODE] = 1 ;
#else
    /* UMFPACK is compiled in normal (non-debug) mode */
    Control [UMFPACK_COMPILED_IN_DEBUG_MODE] = 0 ;
#endif
}
