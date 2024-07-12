/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file dsnode_bmod.c
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


#include "slu_ddefs.h"


/*! \brief Performs numeric block updates within the relaxed snode. 
 */
int
dsnode_bmod (
	    const int  jcol,	  /* in */
	    const int  jsupno,    /* in */
	    const int  fsupc,     /* in */
	    double     *dense,    /* in */
	    double     *tempv,    /* working array */
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
    double         alpha = -1.0, beta = 1.0;
#endif

    int     nsupc, nsupr, nrow;
    int_t   isub, irow;
    int_t   ufirst, nextlu;
    int_t   *lsub, *xlsub;
    double *lusup;
    int_t   *xlusup, luptr;
    flops_t *ops = stat->ops;

    lsub    = Glu->lsub;
    xlsub   = Glu->xlsub;
    lusup   = (double *) Glu->lusup;
    xlusup  = Glu->xlusup;

    nextlu = xlusup[jcol];
    
    /*
     *	Process the supernodal portion of L\U[*,j]
     */
    for (isub = xlsub[fsupc]; isub < xlsub[fsupc+1]; isub++) {
  	irow = lsub[isub];
	lusup[nextlu] = dense[irow];
	dense[irow] = 0;
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

	ops[TRSV] += nsupc * (nsupc - 1);
	ops[GEMV] += 2 * nrow * nsupc;

#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
	STRSV( ftcs1, ftcs2, ftcs3, &nsupc, &lusup[luptr], &nsupr, 
	      &lusup[ufirst], &incx );
	SGEMV( ftcs2, &nrow, &nsupc, &alpha, &lusup[luptr+nsupc], &nsupr, 
		&lusup[ufirst], &incx, &beta, &lusup[ufirst+nsupc], &incy );
#else
#if SCIPY_FIX
       if (nsupr < nsupc) {
           /* Fail early rather than passing in invalid parameters to TRSV. */
           ABORT("failed to factorize matrix");
       }
#endif
	dtrsv_( "L", "N", "U", &nsupc, &lusup[luptr], &nsupr, 
	      &lusup[ufirst], &incx );
	dgemv_( "N", &nrow, &nsupc, &alpha, &lusup[luptr+nsupc], &nsupr, 
		&lusup[ufirst], &incx, &beta, &lusup[ufirst+nsupc], &incy );
#endif
#else
	dlsolve ( nsupr, nsupc, &lusup[luptr], &lusup[ufirst] );
	dmatvec ( nsupr, nrow, nsupc, &lusup[luptr+nsupc], 
			&lusup[ufirst], &tempv[0] );

	int_t i, iptr; 
        /* Scatter tempv[*] into lusup[*] */
	iptr = ufirst + nsupc;
	for (i = 0; i < nrow; i++) {
	    lusup[iptr++] -= tempv[i];
	    tempv[i] = 0.0;
	}
#endif

    }

    return 0;
}
