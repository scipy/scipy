
/*! @file zsnode_bmod.c
 * \brief Performs numeric block updates within the relaxed snode.
 *
 * <pre>
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 * Copyright (c) 1994 by Xerox Corporation.  All rights reserved.
 *
 * THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
 * EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
 * 
 * Permission is hereby granted to use or copy this program for any
 * purpose, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is
 * granted, provided the above notices are retained, and a notice that
 * the code was modified is included with the above copyright notice.
 * </pre>
 */


#include "slu_zdefs.h"


/*! \brief Performs numeric block updates within the relaxed snode. 
 */
int
zsnode_bmod (
	    const int  jcol,	  /* in */
	    const int  jsupno,    /* in */
	    const int  fsupc,     /* in */
	    doublecomplex     *dense,    /* in */
	    doublecomplex     *tempv,    /* working array */
	    GlobalLU_t *Glu,      /* modified */
	    SuperLUStat_t *stat   /* output */
	    )
{
#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
    _fcd ftcs1 = _cptofcd("L", strlen("L")),
	 ftcs2 = _cptofcd("N", strlen("N")),
	 ftcs3 = _cptofcd("U", strlen("U"));
#endif
    int            incx = 1, incy = 1;
    doublecomplex         alpha = {-1.0, 0.0},  beta = {1.0, 0.0};
#endif

    doublecomplex   comp_zero = {0.0, 0.0};
    int            luptr, nsupc, nsupr, nrow;
    int            isub, irow, i, iptr; 
    register int   ufirst, nextlu;
    int            *lsub, *xlsub;
    doublecomplex         *lusup;
    int            *xlusup;
    flops_t *ops = stat->ops;

    lsub    = Glu->lsub;
    xlsub   = Glu->xlsub;
    lusup   = (doublecomplex *) Glu->lusup;
    xlusup  = Glu->xlusup;

    nextlu = xlusup[jcol];
    
    /*
     *	Process the supernodal portion of L\U[*,j]
     */
    for (isub = xlsub[fsupc]; isub < xlsub[fsupc+1]; isub++) {
  	irow = lsub[isub];
	lusup[nextlu] = dense[irow];
        dense[irow] = comp_zero;
	++nextlu;
    }

    xlusup[jcol + 1] = nextlu;	/* Initialize xlusup for next column */
    
    if ( fsupc < jcol ) {

	luptr = xlusup[fsupc];
	nsupr = xlsub[fsupc+1] - xlsub[fsupc];
	nsupc = jcol - fsupc;	/* Excluding jcol */
	ufirst = xlusup[jcol];	/* Points to the beginning of column
				   jcol in supernode L\U(jsupno). */
	nrow = nsupr - nsupc;

	ops[TRSV] += 4 * nsupc * (nsupc - 1);
	ops[GEMV] += 8 * nrow * nsupc;

#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
	CTRSV( ftcs1, ftcs2, ftcs3, &nsupc, &lusup[luptr], &nsupr, 
	      &lusup[ufirst], &incx );
	CGEMV( ftcs2, &nrow, &nsupc, &alpha, &lusup[luptr+nsupc], &nsupr, 
		&lusup[ufirst], &incx, &beta, &lusup[ufirst+nsupc], &incy );
#else
	ztrsv_( "L", "N", "U", &nsupc, &lusup[luptr], &nsupr, 
	      &lusup[ufirst], &incx );
	zgemv_( "N", &nrow, &nsupc, &alpha, &lusup[luptr+nsupc], &nsupr, 
		&lusup[ufirst], &incx, &beta, &lusup[ufirst+nsupc], &incy );
#endif
#else
	zlsolve ( nsupr, nsupc, &lusup[luptr], &lusup[ufirst] );
	zmatvec ( nsupr, nrow, nsupc, &lusup[luptr+nsupc], 
			&lusup[ufirst], &tempv[0] );

        /* Scatter tempv[*] into lusup[*] */
	iptr = ufirst + nsupc;
	for (i = 0; i < nrow; i++) {
	    z_sub(&lusup[iptr], &lusup[iptr], &tempv[i]);
            ++iptr;
	    tempv[i] = comp_zero;
	}
#endif

    }

    return 0;
}
