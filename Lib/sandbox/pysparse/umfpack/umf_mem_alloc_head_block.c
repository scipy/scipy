/* ========================================================================== */
/* === UMF_mem_alloc_head_block ============================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* The UMF_mem_* routines manage the Numeric->Memory memory space. */

/* allocate nunits from head of Numeric->Memory.  No header allocated. */
/* Returns the index into Numeric->Memory if successful, or 0 on failure. */

#include "umf_internal.h"

GLOBAL Int UMF_mem_alloc_head_block
(
    NumericType *Numeric,
    Int nunits
)
{
    Int p, usage ;
    DEBUG2 (("GET  BLOCK: from head, size "ID" ", nunits)) ;

    ASSERT (Numeric != (NumericType *) NULL) ;
    ASSERT (Numeric->Memory != (Unit *) NULL) ;

#ifndef NDEBUG
    if (UMF_allocfail)
    {
	/* pretend to fail, to test garbage_collection */
	DEBUGm2 (("UMF_mem_alloc_head_block: pretend to fail\n")) ;
	UMF_allocfail = FALSE ;	/* don't fail the next time */
	return (0) ;
    }
#endif

    if (nunits > (Numeric->itail - Numeric->ihead))
    {
	DEBUG2 ((" failed\n")) ;
	return (0) ;
    }

    /* return p as an offset from Numeric->Memory */
    p = Numeric->ihead ;
    Numeric->ihead += nunits ;

    DEBUG2 (("p: "ID"\n", p)) ;
    usage = Numeric->ihead + Numeric->tail_usage ;
    Numeric->max_usage = MAX (Numeric->max_usage, usage) ;
    return (p) ;
}
