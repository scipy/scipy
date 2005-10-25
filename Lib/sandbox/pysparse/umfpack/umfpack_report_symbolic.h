/* ========================================================================== */
/* === umfpack_report_symbolic ============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_report_symbolic
(
    void *Symbolic,
    const double Control [UMFPACK_CONTROL]
) ;

long umfpack_dl_report_symbolic
(
    void *Symbolic,
    const double Control [UMFPACK_CONTROL]
) ;

int umfpack_zi_report_symbolic
(
    void *Symbolic,
    const double Control [UMFPACK_CONTROL]
) ;

long umfpack_zl_report_symbolic
(
    void *Symbolic,
    const double Control [UMFPACK_CONTROL]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    double Control [UMFPACK_CONTROL] ;
    int status ;
    status = umfpack_di_report_symbolic (Symbolic, Control) ;

double long Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    double Control [UMFPACK_CONTROL] ;
    long status ;
    status = umfpack_dl_report_symbolic (Symbolic, Control) ;

complex int Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    double Control [UMFPACK_CONTROL] ;
    int status ;
    status = umfpack_zi_report_symbolic (Symbolic, Control) ;

complex long Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    double Control [UMFPACK_CONTROL] ;
    long status ;
    status = umfpack_zl_report_symbolic (Symbolic, Control) ;

Purpose:

    Verifies and prints a Symbolic object.  This routine checks the object more
    carefully than the computational routines.  Normally, this check is not
    required, since umfpack_*_*symbolic either returns (void *) NULL, or a valid
    Symbolic object.  However, if you suspect that your own code has corrupted
    the Symbolic object (by overruning memory bounds, for example), then this
    routine might be able to detect a corrupted Symbolic object.  Since this is
    a complex object, not all such user-generated errors are guaranteed to be
    caught by this routine.

Returns:

    UMFPACK_OK if Control [UMFPACK_PRL] is <= 2 (no inputs are checked).

    Otherwise:

    UMFPACK_OK if the Symbolic object is valid.
    UMFPACK_ERROR_invalid_Symbolic_object if the Symbolic object is invalid.
    UMFPACK_ERROR_out_of_memory if out of memory.

Arguments:

    void *Symbolic ;			Input argument, not modified.

	The Symbolic object, which holds the symbolic factorization computed by
	umfpack_*_*symbolic.

    double Control [UMFPACK_CONTROL] ;	Input argument, not modified.

	If a (double *) NULL pointer is passed, then the default control
	settings are used.  Otherwise, the settings are determined from the
	Control array.  See umfpack_*_defaults on how to fill the Control
	array with the default settings.  If Control contains NaN's, the
	defaults are used.  The following Control parameters are used:

	Control [UMFPACK_PRL]:  printing level.

	    2 or less: no output.  returns silently without checking anything.
	    3: fully check input, and print a short summary of its status
	    4: as 3, but print first few entries of the input
	    5: as 3, but print all of the input
	    Default: 1
*/
