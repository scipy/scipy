/* ========================================================================== */
/* === UMF_start_front ====================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Allocate the initial frontal matrix working array for a single chain.  The
 * front does not have to be big enough, since if it's too small it will get
 * reallocated.  The size computed here is just an estimate. */

#include "umf_internal.h"
#include "umf_grow_front.h"

GLOBAL Int UMF_start_front    /* returns TRUE if successful, FALSE otherwise */
(
    Int chain,
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
)
{
    Int fnrows_max, fncols_max, fnr2, fnc2, fsize, fcurr_size, maxfrsize,
	overflow, nb, f, cdeg ;
    double maxbytes ;

    nb = Symbolic->nb ;
    fnrows_max = Symbolic->Chain_maxrows [chain] ;
    fncols_max = Symbolic->Chain_maxcols [chain] ;

    DEBUGm2 (("Start Front for chain "ID".  fnrows_max "ID" fncols_max "ID"\n",
	chain, fnrows_max, fncols_max)) ;

    Work->fnrows_max = fnrows_max ;
    Work->fncols_max = fncols_max ;
    Work->any_skip = FALSE ;

    maxbytes = sizeof (Entry) *
	(double) (fnrows_max + nb) * (double) (fncols_max + nb) ;
    fcurr_size = Work->fcurr_size ;

    if (Symbolic->prefer_diagonal)
    {
	/* Get a rough upper bound on the degree of the first pivot column in
	 * this front.  Note that Col_degree is not maintained if diagonal
	 * pivoting is preferred.  For most matrices, the first pivot column
	 * of the first frontal matrix of a new chain has only one tuple in
	 * it anyway, so this bound is exact in that case. */
	Int col, tpi, e, *E, *Col_tuples, *Col_tlen, *Cols ;
	Tuple *tp, *tpend ;
	Unit *Memory, *p ;
	Element *ep ;
	E = Work->E ;
	Memory = Numeric->Memory ;
	Col_tuples = Numeric->Lip ;
	Col_tlen = Numeric->Lilen ;
	col = Work->nextcand ;
	tpi = Col_tuples [col] ;
	tp = (Tuple *) Memory + tpi ;
	tpend = tp + Col_tlen [col] ;
	cdeg = 0 ;
	DEBUGm3 (("\n=============== start front: col "ID" tlen "ID"\n",
		col, Col_tlen [col])) ;
	for ( ; tp < tpend ; tp++)
	{
	    DEBUG1 (("Tuple ("ID","ID")\n", tp->e, tp->f)) ;
	    e = tp->e ;
	    if (!E [e]) continue ;
	    f = tp->f ;
	    p = Memory + E [e] ;
	    ep = (Element *) p ;
	    p += UNITS (Element, 1) ;
	    Cols = (Int *) p ;
	    if (Cols [f] == EMPTY) continue ;
	    DEBUG1 (("  nrowsleft "ID"\n", ep->nrowsleft)) ;
	    cdeg += ep->nrowsleft ;
	}
#ifndef NDEBUG
	DEBUGm3 (("start front cdeg: "ID" col "ID"\n", cdeg, col)) ;
	UMF_dump_rowcol (1, Numeric, Work, col, FALSE) ;
#endif

	/* cdeg is now the rough upper bound on the degree of the next pivot
	 * column. */

	/* If AMD was called, we know the maximum number of nonzeros in any
	 * column of L.  Use this as an upper bound for cdeg, but add 2 to
	 * account for a small amount of off-diagonal pivoting. */
	if (Symbolic->amd_dmax > 0)
	{
	    cdeg = MIN (cdeg, Symbolic->amd_dmax) ;
	}

	/* Increase it to account for larger columns later on.
	 * Also ensure that it's larger than zero. */
	cdeg += 2 ;

	/* cdeg cannot be larger than fnrows_max */
	cdeg = MIN (cdeg, fnrows_max) ;

    }
    else
    {
	/* don't do the above cdeg computation */
	cdeg = 0 ;
    }

    DEBUGm2 (("fnrows max "ID" fncols_max "ID"\n", fnrows_max, fncols_max)) ;

    /* the current frontal matrix is empty */
    ASSERT (Work->fnrows == 0 && Work->fncols == 0 && Work->fnpiv == 0) ;

    /* maximum row dimension is always odd, to avoid bad cache effects */
    ASSERT (fnrows_max >= 0) ;
    ASSERT (fnrows_max % 2 == 1) ;

    /* ----------------------------------------------------------------------
     * allocate working array for current frontal matrix:
     * minimum size: 1-by-1
     * maximum size: fnrows_max-by-fncols_max
     * desired size:
     *
     *   if Numeric->front_alloc_init >= 0:
     *
     *	    for unsymmetric matrices:
     *	    Numeric->front_alloc_init * (fnrows_max-by-fncols_max)
     *
     *	    for symmetric matrices (diagonal pivoting preference, actually):
     *	    Numeric->front_alloc_init * (fnrows_max-by-fncols_max), or
     *	    cdeg*cdeg, whichever is smaller.
     *
     *   if Numeric->front_alloc_init < 0:
     *	    allocate a front of size -Numeric->front_alloc_init.
     *
     * Allocate the whole thing if it's small (less than 2*nb^2).  Make sure the
     * leading dimension of the frontal matrix is odd.
     *
     * Also allocate the nb-by-nb LU block, the dr-by-nb L block, and the
     * nb-by-dc U block.
     * ---------------------------------------------------------------------- */

    /* get the maximum front size, avoiding integer overflow */
    overflow = INT_OVERFLOW (maxbytes) ;
    if (overflow)
    {
	/* :: int overflow, max front size :: */
	maxfrsize = Int_MAX / sizeof (Entry) ;
    }
    else
    {
	maxfrsize = (fnrows_max + nb) * (fncols_max + nb) ;
    }
    ASSERT (!INT_OVERFLOW ((double) maxfrsize * sizeof (Entry))) ;

    if (Numeric->front_alloc_init < 0)
    {
	/* allocate a front of -Numeric->front_alloc_init entries */
	fsize = -Numeric->front_alloc_init ;
	fsize = MAX (1, fsize) ;
    }
    else
    {
	if (INT_OVERFLOW (Numeric->front_alloc_init * maxbytes))
	{
	    /* :: int overflow, requested front size :: */
	    fsize = Int_MAX / sizeof (Entry) ;
	}
	else
	{
	    fsize = Numeric->front_alloc_init * maxfrsize ;
	}

	if (cdeg > 0)
	{
	    /* diagonal pivoting is in use.  cdeg was computed above */
	    Int fsize2 ;

	    /* add the L and U blocks */
	    cdeg += nb ;

	    if (INT_OVERFLOW (((double) cdeg * (double) cdeg) * sizeof (Entry)))
	    {
		/* :: int overflow, symmetric front size :: */
		fsize2 = Int_MAX / sizeof (Entry) ;
	    }
	    else
	    {
		fsize2 = MAX (cdeg * cdeg, fcurr_size) ;
	    }
	    fsize = MIN (fsize, fsize2) ;
	}
    }

    fsize = MAX (fsize, 2*nb*nb) ;

    /* fsize and maxfrsize are now safe from integer overflow.  They both
     * include the size of the pivot blocks. */
    ASSERT (!INT_OVERFLOW ((double) fsize * sizeof (Entry))) ;

    Work->fnrows_new = 0 ;
    Work->fncols_new = 0 ;

    /* desired size is fnr2-by-fnc2 (includes L and U blocks): */
    DEBUGm2 (("    fsize "ID"  fcurr_size "ID"\n", fsize, fcurr_size)) ;
    DEBUGm2 (("    maxfrsize "ID"  fnr_curr "ID" fnc_curr "ID"\n", maxfrsize,
	Work->fnr_curr, Work->fnc_curr)) ;

    if (fsize >= maxfrsize && !overflow)
    {
	/* max working array is small, allocate all of it */
	fnr2 = fnrows_max + nb ;
	fnc2 = fncols_max + nb ;
	fsize = maxfrsize ;
	DEBUGm1 (("   sufficient for ("ID"+"ID")-by-("ID"+"ID")\n",
	    fnrows_max, nb, fncols_max, nb)) ;
    }
    else
    {
	/* allocate a smaller working array */
	if (fnrows_max <= fncols_max)
	{
	    fnr2 = (Int) sqrt ((double) fsize) ;
	    /* make sure fnr2 is odd */
	    fnr2 = MAX (fnr2, 1) ;
	    if (fnr2 % 2 == 0) fnr2++ ;
	    fnr2 = MIN (fnr2, fnrows_max + nb) ;
	    fnc2 = fsize / fnr2 ;
	}
	else
	{
	    fnc2 = (Int) sqrt ((double) fsize) ;
	    fnc2 = MIN (fnc2, fncols_max + nb) ;
	    fnr2 = fsize / fnc2 ;
	    /* make sure fnr2 is odd */
	    fnr2 = MAX (fnr2, 1) ;
	    if (fnr2 % 2 == 0)
	    {
		fnr2++ ;
		fnc2 = fsize / fnr2 ;
	    }
	}
	DEBUGm1 (("   smaller "ID"-by-"ID"\n", fnr2, fnc2)) ;
    }
    fnr2 = MIN (fnr2, fnrows_max + nb) ;
    fnc2 = MIN (fnc2, fncols_max + nb) ;
    ASSERT (fnr2 % 2 == 1) ;
    ASSERT (fnr2 * fnc2 <= fsize) ;

    fnr2 -= nb ;
    fnc2 -= nb ;
    ASSERT (fnr2 >= 0) ;
    ASSERT (fnc2 >= 0) ;

    if (fsize > fcurr_size)
    {
	DEBUGm1 (("   Grow front \n")) ;
	Work->do_grow = TRUE ;
	if (!UMF_grow_front (Numeric, fnr2, fnc2, Work, -1))
	{
	    /* since the minimum front size is 1-by-1, it would be nearly
	     * impossible to run out of memory here. */
	    DEBUGm4 (("out of memory: start front\n")) ;
	    return (FALSE) ;
	}
    }
    else
    {
	/* use the existing front */
	DEBUGm1 (("   existing front ok\n")) ;
	Work->fnr_curr = fnr2 ;
	Work->fnc_curr = fnc2 ;
	Work->Flblock  = Work->Flublock + nb * nb ;
	Work->Fublock  = Work->Flblock  + nb * fnr2 ;
	Work->Fcblock  = Work->Fublock  + nb * fnc2 ;
    }
    ASSERT (Work->Flblock  == Work->Flublock + Work->nb*Work->nb) ;
    ASSERT (Work->Fublock  == Work->Flblock  + Work->fnr_curr*Work->nb) ;
    ASSERT (Work->Fcblock  == Work->Fublock  + Work->nb*Work->fnc_curr) ;
    return (TRUE) ;
}
