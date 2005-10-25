/* ========================================================================== */
/* === umfpack_defaults ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

void umfpack_di_defaults
(
    double Control [UMFPACK_CONTROL]
) ;

void umfpack_dl_defaults
(
    double Control [UMFPACK_CONTROL]
) ;

void umfpack_zi_defaults
(
    double Control [UMFPACK_CONTROL]
) ;

void umfpack_zl_defaults
(
    double Control [UMFPACK_CONTROL]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    double Control [UMFPACK_CONTROL] ;
    umfpack_di_defaults (Control) ;

double long Syntax:

    #include "umfpack.h"
    double Control [UMFPACK_CONTROL] ;
    umfpack_dl_defaults (Control) ;

complex int Syntax:

    #include "umfpack.h"
    double Control [UMFPACK_CONTROL] ;
    umfpack_zi_defaults (Control) ;

complex long Syntax:

    #include "umfpack.h"
    double Control [UMFPACK_CONTROL] ;
    umfpack_zl_defaults (Control) ;

Purpose:

    Sets the default control parameter settings.

    NOTE:  new control parameters have been added to the Control array for
    Version 4.1.  These entries were unused in Version 4.0.  The default block
    size for the BLAS has increased from 24 to 32.  Some rarely used control
    parameters have been removed (those that controlled relaxed amalgamation).

Arguments:

    double Control [UMFPACK_CONTROL] ;	Output argument.

	Control is set to the default control parameter settings.  You can
	then modify individual settings by changing specific entries in the
	Control array.  If Control is a (double *) NULL pointer, then
	umfpack_*_defaults returns silently (no error is generated, since
	passing a NULL pointer for Control to any UMFPACK routine is valid).
*/
