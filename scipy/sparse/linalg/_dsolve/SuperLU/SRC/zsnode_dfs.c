/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file zsnode_dfs.c
 * \brief Determines the union of row structures of columns within the relaxed node
 *
 * <pre>
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
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

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *    zsnode_dfs() - Determine the union of the row structures of those 
 *    columns within the relaxed snode.
 *    Note: The relaxed snodes are leaves of the supernodal etree, therefore, 
 *    the portion outside the rectangular supernode must be zero.
 *
 * Return value
 * ============
 *     0   success;
 *    >0   number of bytes allocated when run out of memory.
 * </pre>
 */

int_t
zsnode_dfs (
	   const int   jcol,	    /* in - start of the supernode */
	   const int   kcol, 	    /* in - end of the supernode */
	   const int_t *asub,        /* in */
	   const int_t *xa_begin,    /* in */
	   const int_t *xa_end,      /* in */
	   int_t      *xprune,      /* out */
	   int        *marker,      /* modified */
	   GlobalLU_t *Glu          /* modified */
	   )
{

    int_t i, k, ifrom, ito, nextl, new_next, nzlmax;
    int   nsuper, krow, kmark;
    int_t mem_error;
    int   *xsup, *supno;
    int_t *lsub, *xlsub;
    
    xsup    = Glu->xsup;
    supno   = Glu->supno;
    lsub    = Glu->lsub;
    xlsub   = Glu->xlsub;
    nzlmax  = Glu->nzlmax;

    nsuper = ++supno[jcol];	/* Next available supernode number */
    nextl = xlsub[jcol];

    for (i = jcol; i <= kcol; i++) {
	/* For each nonzero in A[*,i] */
	for (k = xa_begin[i]; k < xa_end[i]; k++) {	
	    krow = asub[k];
	    kmark = marker[krow];
	    if ( kmark != kcol ) { /* First time visit krow */
		marker[krow] = kcol;
		lsub[nextl++] = krow;
		if ( nextl >= nzlmax ) {
		    mem_error = zLUMemXpand(jcol, nextl, LSUB, &nzlmax, Glu);
		    if ( mem_error ) return (mem_error);
		    lsub = Glu->lsub;
		}
	    }
    	}
	supno[i] = nsuper;
    }

    /* Supernode > 1, then make a copy of the subscripts for pruning */
    if ( jcol < kcol ) {
	new_next = nextl + (nextl - xlsub[jcol]);
	while ( new_next > nzlmax ) {
	    mem_error = zLUMemXpand(jcol, nextl, LSUB, &nzlmax, Glu);
	    if ( mem_error ) return (mem_error);
	    lsub = Glu->lsub;
	}
	ito = nextl;
	for (ifrom = xlsub[jcol]; ifrom < nextl; )
	    lsub[ito++] = lsub[ifrom++];	
        for (i = jcol+1; i <= kcol; i++) xlsub[i] = nextl;
	nextl = ito;
    }

    xsup[nsuper+1] = kcol + 1;
    supno[kcol+1]  = nsuper;
    xprune[kcol]   = nextl;
    xlsub[kcol+1]  = nextl;

    return 0;
}

