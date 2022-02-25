/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file ilu_zsnode_dfs.c
 * \brief Determines the union of row structures of columns within the relaxed node
 *
 * <pre>
 * -- SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Laboratory
 * June 30, 2009
 * </pre>
 */

#include "slu_zdefs.h"

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *    ilu_zsnode_dfs() - Determine the union of the row structures of those
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

int
ilu_zsnode_dfs(
	   const int  jcol,	    /* in - start of the supernode */
	   const int  kcol,	    /* in - end of the supernode */
	   const int  *asub,	    /* in */
	   const int  *xa_begin,    /* in */
	   const int  *xa_end,	    /* in */
	   int	      *marker,	    /* modified */
	   GlobalLU_t *Glu	    /* modified */
	   )
{

    register int i, k, nextl;
    int 	 nsuper, krow, kmark, mem_error;
    int 	 *xsup, *supno;
    int 	 *lsub, *xlsub;
    int 	 nzlmax;

    xsup    = Glu->xsup;
    supno   = Glu->supno;
    lsub    = Glu->lsub;
    xlsub   = Glu->xlsub;
    nzlmax  = Glu->nzlmax;

    nsuper = ++supno[jcol];	/* Next available supernode number */
    nextl = xlsub[jcol];

    for (i = jcol; i <= kcol; i++)
    {
	/* For each nonzero in A[*,i] */
	for (k = xa_begin[i]; k < xa_end[i]; k++)
	{
	    krow = asub[k];
	    kmark = marker[krow];
	    if ( kmark != kcol )
	    { /* First time visit krow */
		marker[krow] = kcol;
		lsub[nextl++] = krow;
		if ( nextl >= nzlmax )
		{
		    if ( (mem_error = zLUMemXpand(jcol, nextl, LSUB, &nzlmax,
			    Glu)) != 0)
			return (mem_error);
		    lsub = Glu->lsub;
		}
	    }
	}
	supno[i] = nsuper;
    }

    /* Supernode > 1 */
    if ( jcol < kcol )
	for (i = jcol+1; i <= kcol; i++) xlsub[i] = nextl;

    xsup[nsuper+1] = kcol + 1;
    supno[kcol+1]  = nsuper;
    xlsub[kcol+1]  = nextl;

    return 0;
}
