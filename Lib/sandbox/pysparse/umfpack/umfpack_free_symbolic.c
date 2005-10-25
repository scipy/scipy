/* ========================================================================== */
/* === UMFPACK_free_symbolic ================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  See umfpack_free_symbolic.h for details.
    All 10 objects comprising the Symbolic object are free'd via UMF_free.
*/

#include "umf_internal.h"
#include "umf_free.h"

GLOBAL void UMFPACK_free_symbolic
(
    void **SymbolicHandle
)
{

    SymbolicType *Symbolic ;
    if (!SymbolicHandle)
    {
	return ;
    }
    Symbolic = *((SymbolicType **) SymbolicHandle) ;
    if (!Symbolic)
    {
	return ;
    }

    (void) UMF_free ((void *) Symbolic->Cperm_init) ;
    (void) UMF_free ((void *) Symbolic->Rperm_init) ;
    (void) UMF_free ((void *) Symbolic->Front_npivcol) ;
    (void) UMF_free ((void *) Symbolic->Front_parent) ;
    (void) UMF_free ((void *) Symbolic->Front_1strow) ;
    (void) UMF_free ((void *) Symbolic->Front_leftmostdesc) ;
    (void) UMF_free ((void *) Symbolic->Chain_start) ;
    (void) UMF_free ((void *) Symbolic->Chain_maxrows) ;
    (void) UMF_free ((void *) Symbolic->Chain_maxcols) ;
    (void) UMF_free ((void *) Symbolic->Cdeg) ;
    (void) UMF_free ((void *) Symbolic->Rdeg) ;

    /* only when dense rows are present */
    (void) UMF_free ((void *) Symbolic->Esize) ;

    /* only when diagonal pivoting is prefered */
    (void) UMF_free ((void *) Symbolic->Diagonal_map) ;

    (void) UMF_free ((void *) Symbolic) ;
    *SymbolicHandle = (void *) NULL ;
}
