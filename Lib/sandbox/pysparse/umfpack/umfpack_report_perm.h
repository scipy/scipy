/* ========================================================================== */
/* === umfpack_report_perm ================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_report_perm
(
    int np,
    const int Perm [ ],
    const double Control [UMFPACK_CONTROL]
) ;

long umfpack_dl_report_perm
(
    long np,
    const long Perm [ ],
    const double Control [UMFPACK_CONTROL]
) ;

int umfpack_zi_report_perm
(
    int np,
    const int Perm [ ],
    const double Control [UMFPACK_CONTROL]
) ;

long umfpack_zl_report_perm
(
    long np,
    const long Perm [ ],
    const double Control [UMFPACK_CONTROL]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    int np, *Perm, status ;
    double Control [UMFPACK_CONTROL] ;
    status = umfpack_di_report_perm (np, Perm, Control) ;

double long Syntax:

    #include "umfpack.h"
    long np, *Perm, status ;
    double Control [UMFPACK_CONTROL] ;
    status = umfpack_dl_report_perm (np, Perm, Control) ;

complex int Syntax:

    #include "umfpack.h"
    int np, *Perm, status ;
    double Control [UMFPACK_CONTROL] ;
    status = umfpack_zi_report_perm (np, Perm, Control) ;

complex long Syntax:

    #include "umfpack.h"
    long np, *Perm, status ;
    double Control [UMFPACK_CONTROL] ;
    status = umfpack_zl_report_perm (np, Perm, Control) ;

Purpose:

    Verifies and prints a permutation vector.

Returns:

    UMFPACK_OK if Control [UMFPACK_PRL] <= 2 (the input is not checked).

    Otherwise:
    UMFPACK_OK if the permutation vector is valid (this includes that case
	when Perm is (Int *) NULL, which is not an error condition).
    UMFPACK_ERROR_n_nonpositive if np <= 0.
    UMFPACK_ERROR_out_of_memory if out of memory.
    UMFPACK_ERROR_invalid_permutation if Perm is not a valid permutation vector.

Arguments:

    Int np ;		Input argument, not modified.

	Perm is an integer vector of size np.  Restriction: np > 0.

    Int Perm [np] ;	Input argument, not modified.

	A permutation vector of size np.  If Perm is not present (an (Int *)
	NULL pointer), then it is assumed to be the identity permutation.  This
	is consistent with its use as an input argument to umfpack_*_qsymbolic,
	and is not an error condition.  If Perm is present, the entries in Perm
	must range between 0 and np-1, and no duplicates may exist.

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
