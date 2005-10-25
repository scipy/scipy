/* ========================================================================== */
/* === umfpack_wsolve ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_wsolve
(
    int sys,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    double X [ ],
    const double B [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO],
    int Wi [ ],
    double W [ ]
) ;

long umfpack_dl_wsolve
(
    long sys,
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ],
    double X [ ],
    const double B [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO],
    long Wi [ ],
    double W [ ]
) ;

int umfpack_zi_wsolve
(
    int sys,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ], const double Az [ ],
    double Xx [ ],	 double Xz [ ],
    const double Bx [ ], const double Bz [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO],
    int Wi [ ],
    double W [ ]
) ;

long umfpack_zl_wsolve
(
    long sys,
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ], const double Az [ ],
    double Xx [ ],	 double Xz [ ],
    const double Bx [ ], const double Bz [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO],
    long Wi [ ],
    double W [ ]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    int status, *Ap, *Ai, *Wi, sys ;
    double *B, *X, *Ax, *W, Info [UMFPACK_INFO], Control [UMFPACK_CONTROL] ;
    status = umfpack_di_wsolve (sys, Ap, Ai, Ax, X, B, Numeric,
	Control, Info, Wi, W) ;

double long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    long status, *Ap, *Ai, *Wi, sys ;
    double *B, *X, *Ax, *W, Info [UMFPACK_INFO], Control [UMFPACK_CONTROL] ;
    status = umfpack_dl_wsolve (sys, Ap, Ai, Ax, X, B, Numeric,
	Control, Info, Wi, W) ;

complex int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    int status, *Ap, *Ai, *Wi, sys ;
    double *Bx, *Bz, *Xx, *Xz, *Ax, *Az, *W,
	Info [UMFPACK_INFO], Control [UMFPACK_CONTROL] ;
    status = umfpack_zi_wsolve (sys, Ap, Ai, Ax, Az, Xx, Xz, Bx, Bz, Numeric,
	Control, Info, Wi, W) ;

complex long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    long status, *Ap, *Ai, *Wi, sys ;
    double *Bx, *Bz, *Xx, *Xz, *Ax, *Az, *W,
	Info [UMFPACK_INFO], Control [UMFPACK_CONTROL] ;
    status = umfpack_zl_wsolve (sys, Ap, Ai, Ax, Az, Xx, Xz, Bx, Bz, Numeric,
	Control, Info, Wi, W) ;

Purpose:

    Given LU factors computed by umfpack_*_numeric (PAQ=LU) and the
    right-hand-side, B, solve a linear system for the solution X.  Iterative
    refinement is optionally performed.  This routine is identical to
    umfpack_*_solve, except that it does not dynamically allocate any workspace.
    When you have many linear systems to solve, this routine is faster than
    umfpack_*_solve, since the workspace (Wi, W) needs to be allocated only
    once, prior to calling umfpack_*_wsolve.

Returns:

    The status code is returned.  See Info [UMFPACK_STATUS], below.

Arguments:

    Int sys ;		Input argument, not modified.
    Int Ap [n+1] ;	Input argument, not modified.
    Int Ai [nz] ;	Input argument, not modified.
    double Ax [nz] ;	Input argument, not modified.
    double X [n] ;	Output argument.
    double B [n] ;	Input argument, not modified.
    void *Numeric ;	Input argument, not modified.
    double Control [UMFPACK_CONTROL] ;	Input argument, not modified.
    double Info [UMFPACK_INFO] ;	Output argument.

    for complex versions:
    double Az [nz] ;	Input argument, not modified, imaginary part
    double Xx [n] ;	Output argument, real part.
    double Xz [n] ;	Output argument, imaginary part
    double Bx [n] ;	Input argument, not modified, real part
    double Bz [n] ;	Input argument, not modified, imaginary part

	The above arguments are identical to umfpack_*_solve, except that the
	error code UMFPACK_ERROR_out_of_memory will not be returned in
	Info [UMFPACK_STATUS], since umfpack_*_wsolve does not allocate any
	memory.

    Int Wi [n] ;		Workspace.
    double W [c*n] ;		Workspace, where c is defined below.

	The Wi and W arguments are workspace used by umfpack_*_wsolve.  They
	need not be initialized on input, and their contents are undefined on
	output.  The size of W depends on whether or not iterative refinement is
	used, and which version (real or complex) is called.  Iterative
	refinement is performed if Ax=b, A'x=b, or A.'x=b is being solved,
	Control [UMFPACK_IRSTEP] > 0, and A is nonsingular.  The size of W is
	given below:

				no iter.	with iter.
				refinement	refinement
	umfpack_di_wsolve	n		5*n
	umfpack_dl_wsolve	n		5*n
	umfpack_zi_wsolve	4*n		10*n
	umfpack_zl_wsolve	4*n		10*n
*/
