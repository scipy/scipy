/* ========================================================================== */
/* === UMF_grow_front ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Current frontal matrix is too small.  Make it bigger. */

#include "umf_internal.h"
#include "umf_realloc.h"
#include "umf_mem_free_tail_block.h"
#include "umf_mem_alloc_tail_block.h"
#include "umf_get_memory.h"

GLOBAL Int UMF_grow_front
(
    NumericType *Numeric,
    Int fnr2,		/* desired size is fnr2-by-fnc2 */
    Int fnc2,
    WorkType *Work,
    Int do_what		/* -1: UMF_start_front
			 * 0:  UMF_init_front, do not recompute Fcpos
			 * 1:  UMF_extend_front
			 * 2:  UMF_init_front, recompute Fcpos */
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Entry *Fcold, *Fcnew ;
    Int j, i, col, *Fcpos, *Fcols, fnrows_max, fncols_max, fnr_curr, nb,
	fnrows_new, fncols_new, fnr_min, fnc_min, minsize,
	newsize, fnrows, fncols, *E, eloc ;
    double s ;

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    if (do_what != -1) UMF_debug++ ;
    DEBUG0 (("\n\n====================GROW FRONT: do_what: "ID"\n", do_what)) ;
    if (do_what != -1) UMF_debug-- ;
    ASSERT (Work->do_grow) ;
    ASSERT (Work->fnpiv == 0) ;
#endif

    Fcols = Work->Fcols ;
    Fcpos = Work->Fcpos ;
    E = Work->E ;

    /* ---------------------------------------------------------------------- */
    /* The current front is too small, find the new size */
    /* ---------------------------------------------------------------------- */

    /* maximum size of frontal matrix for this chain */
    nb = Work->nb ;
    fnrows_max = Work->fnrows_max + nb ;
    fncols_max = Work->fncols_max + nb ;
    ASSERT (fnrows_max >= 0 && (fnrows_max % 2) == 1) ;
    DEBUG0 (("Max     size: "ID"-by-"ID" (incl. "ID" pivot block\n",
	fnrows_max, fncols_max, nb)) ;

    /* current dimensions of frontal matrix: fnr-by-fnc */
    DEBUG0 (("Current : "ID"-by-"ID" (excl "ID" pivot blocks)\n",
		Work->fnr_curr, Work->fnc_curr, nb)) ;
    ASSERT (Work->fnr_curr >= 0) ;
    ASSERT ((Work->fnr_curr % 2 == 1) || Work->fnr_curr == 0) ;

    /* required dimensions of frontal matrix: fnr_min-by-fnc_min */
    fnrows_new = Work->fnrows_new + 1 ;
    fncols_new = Work->fncols_new + 1 ;
    ASSERT (fnrows_new >= 0) ;
    if (fnrows_new % 2 == 0) fnrows_new++ ;
    fnrows_new += nb ;
    fncols_new += nb ;
    fnr_min = MIN (fnrows_new, fnrows_max) ;
    fnc_min = MIN (fncols_new, fncols_max) ;
    minsize = fnr_min * fnc_min ;
    if (INT_OVERFLOW ((double) fnr_min * (double) fnc_min * sizeof (Entry)))
    {
	/* :: the minimum front size is bigger than the integer maximum :: */
	return (FALSE) ;
    }
    ASSERT (fnr_min >= 0) ;
    ASSERT (fnr_min % 2 == 1) ;

    DEBUG0 (("Min     : "ID"-by-"ID"\n", fnr_min, fnc_min)) ;

    /* grow the front to fnr2-by-fnc2, but no bigger than the maximum,
     * and no smaller than the minumum. */
    DEBUG0 (("Desired : ("ID"+"ID")-by-("ID"+"ID")\n", fnr2, nb, fnc2, nb)) ;
    fnr2 += nb ;
    fnc2 += nb ;
    ASSERT (fnr2 >= 0) ;
    if (fnr2 % 2 == 0) fnr2++ ;
    fnr2 = MAX (fnr2, fnr_min) ;
    fnc2 = MAX (fnc2, fnc_min) ;
    fnr2 = MIN (fnr2, fnrows_max) ;
    fnc2 = MIN (fnc2, fncols_max) ;
    DEBUG0 (("Try     : "ID"-by-"ID"\n", fnr2, fnc2)) ;
    ASSERT (fnr2 >= 0) ;
    ASSERT (fnr2 % 2 == 1) ;

    s = ((double) fnr2) * ((double) fnc2) ;
    if (INT_OVERFLOW (s * sizeof (Entry)))
    {
	/* :: frontal matrix size int overflow :: */
	/* the desired front size is bigger than the integer maximum */
	/* compute a such that a*a*s < Int_MAX / sizeof (Entry) */
	double a = 0.9 * sqrt ((Int_MAX / sizeof (Entry)) / s) ;
	fnr2 = MAX (fnr_min, a * fnr2) ;
	fnc2 = MAX (fnc_min, a * fnc2) ;
	/* the new frontal size is a*r*a*c = a*a*s */
	newsize = fnr2 * fnc2 ;
	ASSERT (fnr2 >= 0) ;
	if (fnr2 % 2 == 0) fnr2++ ;
	fnc2 = newsize / fnr2 ;
    }

    fnr2 = MAX (fnr2, fnr_min) ;
    fnc2 = MAX (fnc2, fnc_min) ;
    newsize = fnr2 * fnc2 ;

    ASSERT (fnr2 >= 0) ;
    ASSERT (fnr2 % 2 == 1) ;
    ASSERT (fnr2 >= fnr_min) ;
    ASSERT (fnc2 >= fnc_min) ;
    ASSERT (newsize >= minsize) ;

    /* ---------------------------------------------------------------------- */
    /* free the current front if it is empty of any numerical values */
    /* ---------------------------------------------------------------------- */

    if (E [0] && do_what != 1)
    {
	/* free the current front, if it exists and has nothing in it */
	DEBUG0 (("Freeing empty front\n")) ;
	UMF_mem_free_tail_block (Numeric, E [0]) ;
	E [0] = 0 ;
	Work->Flublock = (Entry *) NULL ;
	Work->Flblock  = (Entry *) NULL ;
	Work->Fublock  = (Entry *) NULL ;
	Work->Fcblock  = (Entry *) NULL ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate the new front, doing garbage collection if necessary */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    UMF_allocfail = FALSE ;
    if (UMF_gprob > 0)  /* a double relop, but ignore NaN case */
    {
	double rrr = ((double) (rand ( ))) / (((double) RAND_MAX) + 1) ;
	DEBUG1 (("Check random %e %e\n", rrr, UMF_gprob)) ;
	UMF_allocfail = rrr < UMF_gprob ;
	if (UMF_allocfail) DEBUGm2 (("Random garbage collection (grow)\n")) ;
    }
#endif

    DEBUG0 (("Attempt size: "ID"-by-"ID"\n", fnr2, fnc2)) ;
    eloc = UMF_mem_alloc_tail_block (Numeric, UNITS (Entry, newsize)) ;

    if (!eloc)
    {
	/* Do garbage collection, realloc, and try again. Compact the current
	 * contribution block in the front to fnrows-by-fncols.  Note that
	 * there are no pivot rows/columns in current front.  Do not recompute
	 * Fcpos in UMF_garbage_collection. */
	DEBUGm3 (("get_memory from umf_grow_front\n")) ;
	if (!UMF_get_memory (Numeric, Work, 1 + UNITS (Entry, newsize),
	    Work->fnrows, Work->fncols, FALSE))
	{
	    /* :: out of memory in umf_grow_front :: */
	    return (FALSE) ;	/* out of memory */
	}
	DEBUG0 (("Attempt size: "ID"-by-"ID" again\n", fnr2, fnc2)) ;
	eloc = UMF_mem_alloc_tail_block (Numeric, UNITS (Entry, newsize)) ;
    }

    /* try again with something smaller */
    while ((fnr2 != fnr_min || fnc2 != fnc_min) && !eloc)
    {
	fnr2 = MIN (fnr2 - 2, fnr2 * UMF_REALLOC_REDUCTION) ;
	fnc2 = MIN (fnc2 - 2, fnc2 * UMF_REALLOC_REDUCTION) ;
	ASSERT (fnr_min >= 0) ;
	ASSERT (fnr_min % 2 == 1) ;
	fnr2 = MAX (fnr_min, fnr2) ;
	fnc2 = MAX (fnc_min, fnc2) ;
	ASSERT (fnr2 >= 0) ;
	if (fnr2 % 2 == 0) fnr2++ ;
	newsize = fnr2 * fnc2 ;
	DEBUGm3 (("Attempt smaller size: "ID"-by-"ID" minsize "ID"-by-"ID"\n",
	    fnr2, fnc2, fnr_min, fnc_min)) ;
	eloc = UMF_mem_alloc_tail_block (Numeric, UNITS (Entry, newsize)) ;
    }

    /* try again with the smallest possible size */
    if (!eloc)
    {
	fnr2 = fnr_min ;
	fnc2 = fnc_min ;
	newsize = minsize ;
	DEBUG0 (("Attempt minsize: "ID"-by-"ID"\n", fnr2, fnc2)) ;
	eloc = UMF_mem_alloc_tail_block (Numeric, UNITS (Entry, newsize)) ;
    }

    if (!eloc)
    {
	/* out of memory */
	return (FALSE) ;
    }

    ASSERT (fnr2 >= 0) ;
    ASSERT (fnr2 % 2 == 1) ;
    ASSERT (fnr2 >= fnr_min && fnc2 >= fnc_min) ;

    /* ---------------------------------------------------------------------- */
    /* copy the old frontal matrix into the new one */
    /* ---------------------------------------------------------------------- */

    /* old contribution block (if any) */
    fnr_curr = Work->fnr_curr ;	    /* garbage collection can change fn*_curr */
    ASSERT (fnr_curr >= 0) ;
    ASSERT ((fnr_curr % 2 == 1) || fnr_curr == 0) ;
    fnrows = Work->fnrows ;
    fncols = Work->fncols ;
    Fcold = Work->Fcblock ;

    /* remove nb from the sizes */
    fnr2 -= nb ;
    fnc2 -= nb ;

    /* new frontal matrix */
    Work->Flublock = (Entry *) (Numeric->Memory + eloc) ;
    Work->Flblock  = Work->Flublock + nb * nb ;
    Work->Fublock  = Work->Flblock  + nb * fnr2 ;
    Work->Fcblock  = Work->Fublock  + nb * fnc2 ;
    Fcnew = Work->Fcblock ;

    if (E [0])
    {
	/* copy the old contribution block into the new one */
	for (j = 0 ; j < fncols ; j++)
	{
	    col = Fcols [j] ;
	    DEBUG1 (("copy col "ID" \n",col)) ;
	    ASSERT (col >= 0 && col < Work->n_col) ;
	    for (i = 0 ; i < fnrows ; i++)
	    {
		Fcnew [i] = Fcold [i] ;
	    }
	    Fcnew += fnr2 ;
	    Fcold += fnr_curr ;
	    DEBUG1 (("new offset col "ID" "ID"\n",col, j * fnr2)) ;
	    Fcpos [col] = j * fnr2 ;
	}
    }
    else if (do_what == 2)
    {
	/* just find the new column offsets */
	for (j = 0 ; j < fncols ; j++)
	{
	    col = Fcols [j] ;
	    DEBUG1 (("new offset col "ID" "ID"\n",col, j * fnr2)) ;
	    Fcpos [col] = j * fnr2 ;
	}
    }

    /* free the old frontal matrix */
    UMF_mem_free_tail_block (Numeric, E [0]) ;

    /* ---------------------------------------------------------------------- */
    /* new frontal matrix size */
    /* ---------------------------------------------------------------------- */

    E [0] = eloc ;
    Work->fnr_curr = fnr2 ;	    /* C block is fnr2-by-fnc2 */
    Work->fnc_curr = fnc2 ;
    Work->fcurr_size = newsize ;    /* including LU, L, U, and C blocks */
    Work->do_grow = FALSE ;	    /* the front has just been grown */

    ASSERT (Work->fnr_curr >= 0) ;
    ASSERT (Work->fnr_curr % 2 == 1) ;
    DEBUG0 (("Newly grown front: "ID"+"ID" by "ID"+"ID"\n", Work->fnr_curr,
	nb, Work->fnc_curr, nb)) ;
    return (TRUE) ;
}
