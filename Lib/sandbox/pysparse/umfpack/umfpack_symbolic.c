/* ========================================================================== */
/* === UMFPACK_symbolic ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Performs a symbolic factorization.
    See umfpack_symbolic.h for details.
*/

#include "umf_internal.h"

GLOBAL Int UMFPACK_symbolic
(
    Int n_row,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    void **SymbolicHandle,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
)
{
    Int *Qinit = (Int *) NULL ;
    return (UMFPACK_qsymbolic (n_row, n_col, Ap, Ai, Ax,
#ifdef COMPLEX
	Az,
#endif
	Qinit, SymbolicHandle, Control, Info)) ;
}
