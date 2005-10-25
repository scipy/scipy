/* ========================================================================== */
/* === UMF_mem_free_tail_block ============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* The UMF_mem_* routines manage the Numeric->Memory memory space. */

/* free a block from the tail of Numeric->memory */

#include "umf_internal.h"

GLOBAL void UMF_mem_free_tail_block
(
    NumericType *Numeric,
    Int i
)
{
    Unit *pprev, *pnext, *p, *pbig ;
    Int sprev ;

    ASSERT (Numeric != (NumericType *) NULL) ;
    ASSERT (Numeric->Memory != (Unit *) NULL) ;
    if (i == EMPTY || i == 0) return ;	/* already deallocated */

    /* ---------------------------------------------------------------------- */
    /* get the block */
    /* ---------------------------------------------------------------------- */

    p = Numeric->Memory + i ;

    p-- ;	/* get the corresponding header */
    DEBUG2 (("free block: p: "ID, (Int) (p-Numeric->Memory))) ;
    ASSERT (p >= Numeric->Memory + Numeric->itail) ;
    ASSERT (p < Numeric->Memory + Numeric->size) ;
    ASSERT (p->header.size > 0) ;		/* block not already free */
    ASSERT (p->header.prevsize >= 0) ;

    Numeric->tail_usage -= p->header.size + 1 ;

    /* ---------------------------------------------------------------------- */
    /* merge with next free block, if any */
    /* ---------------------------------------------------------------------- */

    pnext = p + 1 + p->header.size ;
    DEBUG2 (("size: "ID" next: "ID" ", p->header.size,
	(Int) (pnext-Numeric->Memory))) ;
    ASSERT (pnext < Numeric->Memory + Numeric->size) ;
    ASSERT (pnext->header.prevsize == p->header.size) ;
    ASSERT (pnext->header.size != 0) ;

    if (pnext->header.size < 0)
    {
	/* next block is also free - merge with current block */
	p->header.size += (-(pnext->header.size)) + 1 ;
	DEBUG2 ((" NEXT FREE ")) ;
    }

    /* ---------------------------------------------------------------------- */
    /* merge with previous free block, if any */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    if (p == Numeric->Memory + Numeric->itail)
    {
	DEBUG2 ((" at top of tail ")) ;
	ASSERT (p->header.prevsize == 0) ;
    }
#endif

    if (p > Numeric->Memory + Numeric->itail)
    {
	ASSERT (p->header.prevsize > 0) ;
	pprev = p - 1 - p->header.prevsize ;
	DEBUG2 ((" prev: "ID" ", (Int) (pprev-Numeric->Memory))) ;
	ASSERT (pprev >= Numeric->Memory + Numeric->itail) ;
	sprev = pprev->header.size ;
	if (sprev < 0)
	{
	    /* previous block is also free - merge it with current block */
	    ASSERT (p->header.prevsize == -sprev) ;
	    pprev->header.size = p->header.size + (-sprev) + 1 ;
	    p = pprev ;
	    DEBUG2 ((" PREV FREE ")) ;
	    /* note that p may now point to Numeric->itail */
	}
#ifndef NDEBUG
	else
	{
	    ASSERT (p->header.prevsize == sprev) ;
	}
#endif
    }

    /* ---------------------------------------------------------------------- */
    /* free the block, p */
    /* ---------------------------------------------------------------------- */

    pnext = p + 1 + p->header.size ;
    ASSERT (pnext < Numeric->Memory + Numeric->size) ;

    if (p == Numeric->Memory + Numeric->itail)
    {
	/* top block in list is freed */
	Numeric->itail = pnext - Numeric->Memory ;
	pnext->header.prevsize = 0 ;
	DEBUG2 ((" NEW TAIL : "ID" ", Numeric->itail)) ;
	ASSERT (pnext->header.size > 0) ;
	if (Numeric->ibig != EMPTY && Numeric->ibig <= Numeric->itail)
	{
	    /* the big free block is now above the tail */
	    Numeric->ibig = EMPTY ;
	}
    }
    else
    {
	/* keep track of the biggest free block seen */
	if (Numeric->ibig == EMPTY)
	{
	    Numeric->ibig = p - Numeric->Memory ;
	}
	else
	{
	    pbig = Numeric->Memory + Numeric->ibig ;
	    if (-(pbig->header.size) < p->header.size)
	    {
		Numeric->ibig = p - Numeric->Memory ;
	    }
	}
	/* flag the block as free, somewhere in the middle of the tail */
	pnext->header.prevsize = p->header.size ;
	p->header.size = -(p->header.size) ;
    }

    DEBUG2 (("new p: "ID" freesize: "ID"\n", (Int) (p-Numeric->Memory),
	-(p->header.size))) ;

}
