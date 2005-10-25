/* ========================================================================== */
/* === umfpack_report_vector ================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_report_vector
(
    int n,
    const double X [ ],
    const double Control [UMFPACK_CONTROL]
) ;

long umfpack_dl_report_vector
(
    long n,
    const double X [ ],
    const double Control [UMFPACK_CONTROL]
) ;

int umfpack_zi_report_vector
(
    int n,
    const double Xx [ ], const double Xz [ ],
    const double Control [UMFPACK_CONTROL]
) ;

long umfpack_zl_report_vector
(
    long n,
    const double Xx [ ], const double Xz [ ],
    const double Control [UMFPACK_CONTROL]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    int n, status ;
    double *X, Control [UMFPACK_CONTROL] ;
    status = umfpack_di_report_vector (n, X, Control) ;

double long Syntax:

    #include "umfpack.h"
    long n, status ;
    double *X, Control [UMFPACK_CONTROL] ;
    status = umfpack_dl_report_vector (n, X, Control) ;

complex int Syntax:

    #include "umfpack.h"
    int n, status ;
    double *Xx, *Xz, Control [UMFPACK_CONTROL] ;
    status = umfpack_zi_report_vector (n, Xx, Xz, Control) ;

complex long Syntax:

    #include "umfpack.h"
    long n, status ;
    double *Xx, *Xz, Control [UMFPACK_CONTROL] ;
    status = umfpack_zl_report_vector (n, Xx, Xz, Control) ;

Purpose:

    Verifies and prints a dense vector.

Returns:

    UMFPACK_OK if Control [UMFPACK_PRL] <= 2 (the input is not checked).

    Otherwise:
    
    UMFPACK_OK if the vector is valid.
    UMFPACK_ERROR_argument_missing if X or Xx is missing.
    UMFPACK_ERROR_n_nonpositive if n <= 0.

Arguments:

    Int n ;		Input argument, not modified.

	X is a real or complex vector of size n.  Restriction: n > 0.

    double X [n] ;	Input argument, not modified.  For real versions.

	A real vector of size n.  X must not be (double *) NULL.

    double Xx [n or 2*n] ; Input argument, not modified.  For complex versions.
    double Xz [n or 0] ;   Input argument, not modified.  For complex versions.

	A complex vector of size n, in one of two storage formats.
	Xx must not be (double *) NULL.

	If Xz is not (double *) NULL, then Xx [i] is the real part of X (i) and
	Xz [i] is the imaginary part of X (i).  Both vectors are of length n.
	This is the "split" form of the complex vector X.

	If Xz is (double *) NULL, then Xx holds both real and imaginary parts,
	where Xx [2*i] is the real part of X (i) and Xx [2*i+1] is the imaginary
	part of X (i).  Xx is of length 2*n doubles.  If you have an ANSI C99
	compiler with the intrinsic double _Complex type, then Xx can be of
	type double _Complex in the calling routine and typecast to (double *)
	when passed to umfpack_*_report_vector (this is untested, however).
	This is the "merged" form of the complex vector X.

	Future work:  all complex routines in UMFPACK could use this same
	strategy for their complex arguments.  The split format is useful for
	MATLAB, which holds its real and imaginary parts in seperate arrays.
	The merged format is compatible with the intrinsic double _Complex
	type in ANSI C99, and is also compatible with SuperLU's method of
	storing complex matrices.  In the current version, only 
	umfpack_*_report_vector supports both formats.

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
