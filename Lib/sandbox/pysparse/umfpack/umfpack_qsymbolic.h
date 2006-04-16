/* ========================================================================== */
/* === umfpack_qsymbolic ==================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_qsymbolic
(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    const int Qinit [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

long umfpack_dl_qsymbolic
(
    long n_row,
    long n_col,
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ],
    const long Qinit [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

int umfpack_zi_qsymbolic
(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ], const double Az [ ],
    const int Qinit [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

long umfpack_zl_qsymbolic
(
    long n_row,
    long n_col,
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ], const double Az [ ],
    const long Qinit [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    int n_row, n_col, *Ap, *Ai, *Qinit, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO], *Ax ;
    status = umfpack_di_qsymbolic (n_row, n_col, Ap, Ai, Ax, Qinit,
	&Symbolic, Control, Info) ;

double long Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    long n_row, n_col, *Ap, *Ai, *Qinit, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO], *Ax ;
    status = umfpack_dl_qsymbolic (n_row, n_col, Ap, Ai, Ax, Qinit,
	&Symbolic, Control, Info) ;

complex int Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    int n_row, n_col, *Ap, *Ai, *Qinit, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO], *Ax, *Az ;
    status = umfpack_zi_qsymbolic (n_row, n_col, Ap, Ai, Ax, Az, Qinit,
	&Symbolic, Control, Info) ;

complex long Syntax:

    #include "umfpack.h"
    void *Symbolic ;
    long n_row, n_col, *Ap, *Ai, *Qinit, status ;
    double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO], *Ax, *Az ;
    status = umfpack_zl_qsymbolic (n_row, n_col, Ap, Ai, Ax, Az, Qinit,
	&Symbolic, Control, Info) ;

Purpose:

    Given the nonzero pattern of a sparse matrix A in column-oriented form, and
    a sparsity preserving column pre-ordering Qinit, umfpack_*_qsymbolic
    performs the symbolic factorization of A*Qinit (or A (:,Qinit) in MATLAB
    notation).  This is identical to umfpack_*_symbolic, except that neither
    COLAMD nor AMD are called and the user input column order Qinit is used
    instead.  Note that in general, the Qinit passed to umfpack_*_qsymbolic
    can differ from the final Q found in umfpack_*_numeric.  The unsymmetric
    strategy will perform a column etree postordering done in
    umfpack_*_qsymbolic and sparsity-preserving modifications are made within
    each frontal matrix during umfpack_*_numeric.  The symmetric and 2-by-2
    strategies will preserve Qinit, unless the matrix is structurally singular.

    See umfpack_*_symbolic for more information.

    *** WARNING ***  A poor choice of Qinit can easily cause umfpack_*_numeric
    to use a huge amount of memory and do a lot of work.  The "default" symbolic
    analysis method is umfpack_*_symbolic, not this routine.  If you use this
    routine, the performance of UMFPACK is your responsibility;  UMFPACK will
    not try to second-guess a poor choice of Qinit.

Returns:

    The value of Info [UMFPACK_STATUS]; see umfpack_*_symbolic.
    Also returns UMFPACK_ERROR_invalid_permuation if Qinit is not a valid
    permutation vector.

Arguments:

    All arguments are the same as umfpack_*_symbolic, except for the following:

    Int Qinit [n_col] ;		Input argument, not modified.

	The user's fill-reducing initial column pre-ordering.  This must be a
	permutation of 0..n_col-1.  If Qinit [k] = j, then column j is the kth
	column of the matrix A (:,Qinit) to be factorized.  If Qinit is an
	(Int *) NULL pointer, then COLAMD or AMD are called instead.

*/
