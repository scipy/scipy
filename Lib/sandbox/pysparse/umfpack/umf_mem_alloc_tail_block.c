/* ========================================================================== */
/* === UMF_mem_alloc_tail_block ============================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* The UMF_mem_* routines manage the Numeric->Memory memory space. */

#include "umf_internal.h"

/* allocate nunits from tail of Numeric->Memory */
/* (requires nunits+1, for header). */
/* Returns the index into Numeric->Memory if successful, or 0 on failure. */

GLOBAL Int UMF_mem_alloc_tail_block
(
    NumericType *Numeric,
    Int nunits
)
{
    Int bigsize, usage ;
    Unit *p, *pnext, *pbig ;

    ASSERT (Numeric != (NumericType *) NULL) ;
    ASSERT (Numeric->Memory != (Unit *) NULL) ;

#ifndef NDEBUG
    if (UMF_allocfail)
    {
	/* pretend to fail, to test garbage_collection */
	DEBUGm2 (("UMF_mem_alloc_tail_block: pretend to fail\n")) ;
	UMF_allocfail = FALSE ;	/* don't fail the next time */
	return (0) ;
    }
    DEBUG2 (("UMF_mem_alloc_tail_block, size: "ID" + 1 = "ID":  ",
	nunits, nunits+1)) ;
#endif

    bigsize = 0 ;
    pbig = (Unit *) NULL ;

    ASSERT (nunits > 0) ;	/* size must be positive */
    if (Numeric->ibig != EMPTY)
    {
	ASSERT (Numeric->ibig > Numeric->itail) ;
	ASSERT (Numeric->ibig < Numeric->size) ;
	pbig = Numeric->Memory + Numeric->ibig ;
	bigsize = -pbig->header.size ;
	ASSERT (bigsize > 0) ;	/* Numeric->ibig is free */
	ASSERT (pbig->header.prevsize >= 0) ;	/* prev. is not free */
    }

    if (pbig && bigsize >= nunits)
    {

	/* use the biggest block, somewhere in middle of memory */
	p = pbig ;
	pnext = p + 1 + bigsize ;
	/* next is in range */
	ASSERT (pnext < Numeric->Memory + Numeric->size) ;
	/* prevsize of next = this size */
	ASSERT (pnext->header.prevsize == bigsize) ;
	/* next is not free */
	ASSERT (pnext->header.size > 0) ;
	bigsize -= nunits + 1 ;

	if (bigsize < 4)
	{
	    /* internal fragmentation would be too small */
	    /* allocate the entire free block */
	    p->header.size = -p->header.size ;
	    DEBUG2 (("GET  BLOCK: p: "ID" size: "ID", all of big: "ID" size: "
		ID"\n", (Int) (p-Numeric->Memory), nunits, Numeric->ibig,
		p->header.size)) ;
	    /* no more biggest block */
	    Numeric->ibig = EMPTY ;

	}
	else
	{

	    /* allocate just the first nunits Units of the free block */
	    p->header.size = nunits ;
	    /* make a new free block */
	    Numeric->ibig += nunits + 1 ;
	    pbig = Numeric->Memory + Numeric->ibig ;
	    pbig->header.size = -bigsize ;
	    pbig->header.prevsize = nunits ;
	    pnext->header.prevsize = bigsize ;
	    DEBUG2 (("GET  BLOCK: p: "ID" size: "ID", some of big: "ID" left: "
		ID"\n", (Int) (p-Numeric->Memory), nunits, Numeric->ibig,
		bigsize)) ;
	}

    }
    else
    {

	/* allocate from the top of tail */
	pnext = Numeric->Memory + Numeric->itail ;
	DEBUG2 (("GET  BLOCK: from tail ")) ;
	if ((nunits + 1) > (Numeric->itail - Numeric->ihead))
	{
	    DEBUG2 (("\n")) ;
	    return (0) ;
	}
	Numeric->itail -= (nunits + 1) ;
	p = Numeric->Memory + Numeric->itail ;
	p->header.size = nunits ;
	p->header.prevsize = 0 ;
	pnext->header.prevsize = nunits ;
	DEBUG2 (("p: "ID" size: "ID", new tail "ID"\n",
	    (Int) (p-Numeric->Memory), nunits, Numeric->itail)) ;
    }

    Numeric->tail_usage += p->header.size + 1 ;
    usage = Numeric->ihead + Numeric->tail_usage ;
    Numeric->max_usage = MAX (Numeric->max_usage, usage) ;

#ifndef NDEBUG
    UMF_debug -= 10 ;
    UMF_dump_memory (Numeric) ;
    UMF_debug += 10 ;
#endif

    /* p points to the header.  Add one to point to the usable block itself. */
    /* return the offset into Numeric->Memory */
    return ((p - Numeric->Memory) + 1) ;
}
