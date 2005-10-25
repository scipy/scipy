/* ========================================================================== */
/* === umfpack_solve ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_solve
(
    int sys,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    double X [ ],
    const double B [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

long umfpack_dl_solve
(
    long sys,
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ],
    double X [ ],
    const double B [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

int umfpack_zi_solve
(
    int sys,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ], const double Az [ ],
    double Xx [ ],	 double Xz [ ],
    const double Bx [ ], const double Bz [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

long umfpack_zl_solve
(
    long sys,
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ], const double Az [ ],
    double Xx [ ],	 double Xz [ ],
    const double Bx [ ], const double Bz [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    int status, *Ap, *Ai, sys ;
    double *B, *X, *Ax, Info [UMFPACK_INFO], Control [UMFPACK_CONTROL] ;
    status = umfpack_di_solve (sys, Ap, Ai, Ax, X, B, Numeric, Control, Info) ;

double long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    long status, *Ap, *Ai, sys ;
    double *B, *X, *Ax, Info [UMFPACK_INFO], Control [UMFPACK_CONTROL] ;
    status = umfpack_dl_solve (sys, Ap, Ai, Ax, X, B, Numeric, Control, Info) ;

complex int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    int status, *Ap, *Ai, sys ;
    double *Bx, *Bz, *Xx, *Xz, *Ax, *Az, Info [UMFPACK_INFO],
	Control [UMFPACK_CONTROL] ;
    status = umfpack_zi_solve (sys, Ap, Ai, Ax, Az, Xx, Xz, Bx, Bz, Numeric,
	Control, Info) ;

complex long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    long status, *Ap, *Ai, sys ;
    double *Bx, *Bz, *Xx, *Xz, *Ax, *Az, Info [UMFPACK_INFO],
	Control [UMFPACK_CONTROL] ;
    status = umfpack_zl_solve (sys, Ap, Ai, Ax, Az, Xx, Xz, Bx, Bz, Numeric,
	Control, Info) ;

Purpose:

    Given LU factors computed by umfpack_*_numeric (PAQ=LU, PRAQ=LU, or
    P(R\A)Q=LU) and the right-hand-side, B, solve a linear system for the
    solution X.  Iterative refinement is optionally performed.  Only square
    systems are handled.  Singular matrices result in a divide-by-zero for all
    systems except those involving just the matrix L.  Iterative refinement is
    not performed for singular matrices.  In the discussion below, n is equal
    to n_row and n_col, because only square systems are handled.

Returns:

    The status code is returned.  See Info [UMFPACK_STATUS], below.

Arguments:

    Int sys ;		Input argument, not modified.

	Defines which system to solve.  (') is the linear algebraic transpose
	(complex conjugate if A is complex), and (.') is the array transpose.

	    sys value	    system solved
	    UMFPACK_A       Ax=b
	    UMFPACK_At      A'x=b
	    UMFPACK_Aat     A.'x=b
	    UMFPACK_Pt_L    P'Lx=b
	    UMFPACK_L       Lx=b
	    UMFPACK_Lt_P    L'Px=b
	    UMFPACK_Lat_P   L.'Px=b
	    UMFPACK_Lt      L'x=b
	    UMFPACK_U_Qt    UQ'x=b
	    UMFPACK_U       Ux=b
	    UMFPACK_Q_Ut    QU'x=b
	    UMFPACK_Q_Uat   QU.'x=b
	    UMFPACK_Ut      U'x=b
	    UMFPACK_Uat     U.'x=b

	Iterative refinement can be optionally performed when sys is any of
	the following:

	    UMFPACK_A       Ax=b
	    UMFPACK_At      A'x=b
	    UMFPACK_Aat     A.'x=b

	For the other values of the sys argument, iterative refinement is not
	performed (Control [UMFPACK_IRSTEP], Ap, Ai, Ax, and Az are ignored).

	Earlier versions used a string argument for sys.  It was changed to an
	integer to make it easier for a Fortran code to call UMFPACK.

    Int Ap [n+1] ;	Input argument, not modified.
    Int Ai [nz] ;	Input argument, not modified.
    double Ax [nz] ;	Input argument, not modified.
    double Az [nz] ;	Input argument, not modified, for complex versions.

	If iterative refinement is requested (Control [UMFPACK_IRSTEP] >= 1,
	Ax=b, A'x=b, or A.'x=b is being solved, and A is nonsingular), then
	these arrays must be identical to the same ones passed to
	umfpack_*_numeric.  The umfpack_*_solve routine does not check the
	contents of these arguments, so the results are undefined if Ap, Ai, Ax,
	and/or Az are modified between the calls the umfpack_*_numeric and
	umfpack_*_solve.  These three arrays do not need to be present (NULL
	pointers can be passed) if Control [UMFPACK_IRSTEP] is zero, or if a
	system other than Ax=b, A'x=b, or A.'x=b is being solved, or if A is
	singular, since in each of these cases A is not accessed.

	Future complex version:  if Ax is present and Az is NULL, then both real
	and imaginary parts will be contained in Ax[0..2*nz-1], with Ax[2*k]
	and Ax[2*k+1] being the real and imaginary part of the kth entry.

    double X [n] ;	Output argument.
    or:
    double Xx [n] ;	Output argument, real part.
    double Xz [n] ;	Output argument, imaginary part.

	The solution to the linear system, where n = n_row = n_col is the
	dimension of the matrices A, L, and U.

	Future complex version:  if Xx is present and Xz is NULL, then both real
	and imaginary parts will be returned in Xx[0..2*n-1], with Xx[2*k] and
	Xx[2*k+1] being the real and imaginary part of the kth entry.

    double B [n] ;	Input argument, not modified.
    or:
    double Bx [n] ;	Input argument, not modified, real part.
    double Bz [n] ;	Input argument, not modified, imaginary part.

	The right-hand side vector, b, stored as a conventional array of size n
	(or two arrays of size n for complex versions).  This routine does not
	solve for multiple right-hand-sides, nor does it allow b to be stored in
	a sparse-column form.

	Future complex version:  if Bx is present and Bz is NULL, then both real
	and imaginary parts will be contained in Bx[0..2*n-1], with Bx[2*k]
	and Bx[2*k+1] being the real and imaginary part of the kth entry.

    void *Numeric ;		Input argument, not modified.

	Numeric must point to a valid Numeric object, computed by
	umfpack_*_numeric.

    double Control [UMFPACK_CONTROL] ;	Input argument, not modified.

	If a (double *) NULL pointer is passed, then the default control
	settings are used.  Otherwise, the settings are determined from the
	Control array.  See umfpack_*_defaults on how to fill the Control
	array with the default settings.  If Control contains NaN's, the
	defaults are used.  The following Control parameters are used:

	Control [UMFPACK_IRSTEP]:  The maximum number of iterative refinement
	    steps to attempt.  A value less than zero is treated as zero.  If
	    less than 1, or if Ax=b, A'x=b, or A.'x=b is not being solved, or
	    if A is singular, then the Ap, Ai, Ax, and Az arguments are not
	    accessed.  Default: 2.

    double Info [UMFPACK_INFO] ;	Output argument.

	Contains statistics about the solution factorization.  If a
	(double *) NULL pointer is passed, then no statistics are returned in
	Info (this is not an error condition).  The following statistics are
	computed in umfpack_*_solve:

	Info [UMFPACK_STATUS]: status code.  This is also the return value,
	    whether or not Info is present.

	    UMFPACK_OK

		The linear system was successfully solved.

	    UMFPACK_WARNING_singular_matrix

		A divide-by-zero occured.  Your solution will contain Inf's
		and/or NaN's.  Some parts of the solution may be valid.  For
		example, solving Ax=b with

		A = [2 0]  b = [ 1 ]  returns x = [ 0.5 ]
		    [0 0]      [ 0 ]              [ Inf ]

	    UMFPACK_ERROR_out_of_memory

		Insufficient memory to solve the linear system.

	    UMFPACK_ERROR_argument_missing

		One or more required arguments are missing.  The B, X, (or
		Bx, Bz, Xx and Xz for the complex versions) arguments
		are always required.  Info and Control are not required.  Ap,
		Ai, Ax (and Az for complex versions) are required if Ax=b,
		A'x=b, A.'x=b is to be solved, the (default) iterative
		refinement is requested, and the matrix A is nonsingular.

	    UMFPACK_ERROR_invalid_system

		The sys argument is not valid, or the matrix A is not square.

	    UMFPACK_ERROR_invalid_Numeric_object

		The Numeric object is not valid.

	Info [UMFPACK_NROW], Info [UMFPACK_NCOL]:
		The dimensions of the matrix A (L is n_row-by-n_inner and
		U is n_inner-by-n_col, with n_inner = min(n_row,n_col)).

	Info [UMFPACK_NZ]:  the number of entries in the input matrix, Ap [n],
	    if iterative refinement is requested (Ax=b, A'x=b, or A.'x=b is
	    being solved, Control [UMFPACK_IRSTEP] >= 1, and A is nonsingular).

	Info [UMFPACK_IR_TAKEN]:  The number of iterative refinement steps
	    effectively taken.  The number of steps attempted may be one more
	    than this; the refinement algorithm backtracks if the last
	    refinement step worsens the solution.

	Info [UMFPACK_IR_ATTEMPTED]:   The number of iterative refinement steps
	    attempted.  The number of times a linear system was solved is one
	    more than this (once for the initial Ax=b, and once for each Ay=r
	    solved for each iterative refinement step attempted).

	Info [UMFPACK_OMEGA1]:  sparse backward error estimate, omega1, if
	    iterative refinement was performed, or -1 if iterative refinement
	    not performed.

	Info [UMFPACK_OMEGA2]:  sparse backward error estimate, omega2, if
	    iterative refinement was performed, or -1 if iterative refinement
	    not performed.

	Info [UMFPACK_SOLVE_FLOPS]:  the number of floating point operations
	    performed to solve the linear system.  This includes the work
	    taken for all iterative refinement steps, including the backtrack
	    (if any).

	Info [UMFPACK_SOLVE_TIME]:  The time taken, in seconds.

	------------------------------------------------------------------------
	The following statistic was added to Version 4.1:
	------------------------------------------------------------------------

        Info [UMFPACK_SOLVE_WALLTIME]:  The wallclock time taken, in seconds.

	Only the above listed Info [...] entries are accessed.  The remaining
	entries of Info are not accessed or modified by umfpack_*_solve.
	Future versions might modify different parts of Info.
*/
