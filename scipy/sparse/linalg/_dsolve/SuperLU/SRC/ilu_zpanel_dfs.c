/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file ilu_zpanel_dfs.c
 * \brief Peforms a symbolic factorization on a panel of symbols and
 * record the entries with maximum absolute value in each column
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
 *
 *   Performs a symbolic factorization on a panel of columns [jcol, jcol+w).
 *
 *   A supernode representative is the last column of a supernode.
 *   The nonzeros in U[*,j] are segments that end at supernodal
 *   representatives.
 *
 *   The routine returns one list of the supernodal representatives
 *   in topological order of the dfs that generates them. This list is
 *   a superset of the topological order of each individual column within
 *   the panel.
 *   The location of the first nonzero in each supernodal segment
 *   (supernodal entry location) is also returned. Each column has a
 *   separate list for this purpose.
 *
 *   Two marker arrays are used for dfs:
 *     marker[i] == jj, if i was visited during dfs of current column jj;
 *     marker1[i] >= jcol, if i was visited by earlier columns in this panel;
 *
 *   marker: A-row --> A-row/col (0/1)
 *   repfnz: SuperA-col --> PA-row
 *   parent: SuperA-col --> SuperA-col
 *   xplore: SuperA-col --> index to L-structure
 * </pre>
 */
void
ilu_zpanel_dfs(
   const int  m,	   /* in - number of rows in the matrix */
   const int  w,	   /* in */
   const int  jcol,	   /* in */
   SuperMatrix *A,	   /* in - original matrix */
   int	      *perm_r,	   /* in */
   int	      *nseg,	   /* out */
   doublecomplex     *dense,	   /* out */
   double     *amax,	   /* out - max. abs. value of each column in panel */
   int	      *panel_lsub, /* out */
   int	      *segrep,	   /* out */
   int	      *repfnz,	   /* out */
   int	      *marker,	   /* out */
   int	      *parent,	   /* working array */
   int_t      *xplore,	   /* working array */
   GlobalLU_t *Glu	   /* modified */
)
{

    NCPformat *Astore;
    doublecomplex    *a;
    int_t     *asub;
    int_t     *xa_begin, *xa_end;
    int       krep, chperm, chmark, chrep, oldrep, kchild, myfnz;
    int       krow, kmark, kperm, kpar;
    int_t     xdfs, maxdfs, k;
    int       jj;	   /* index through each column in the panel */
    int       *marker1;    /* marker1[jj] >= jcol if vertex jj was visited
			      by a previous column within this panel. */
    int       *repfnz_col; /* start of each column in the panel */
    doublecomplex    *dense_col;  /* start of each column in the panel */
    int_t     nextl_col;   /* next available position in panel_lsub[*,jj] */
    int       *xsup, *supno;
    int_t     *lsub, *xlsub;
    double    *amax_col;
    register double tmp;

    /* Initialize pointers */
    Astore     = A->Store;
    a	       = Astore->nzval;
    asub       = Astore->rowind;
    xa_begin   = Astore->colbeg;
    xa_end     = Astore->colend;
    marker1    = marker + m;
    repfnz_col = repfnz;
    dense_col  = dense;
    amax_col   = amax;
    *nseg      = 0;
    xsup       = Glu->xsup;
    supno      = Glu->supno;
    lsub       = Glu->lsub;
    xlsub      = Glu->xlsub;

    /* For each column in the panel */
    for (jj = jcol; jj < jcol + w; jj++) {
	nextl_col = (jj - jcol) * m;

#ifdef CHK_DFS
	printf("\npanel col %d: ", jj);
#endif

	*amax_col = 0.0;
	/* For each nonz in A[*,jj] do dfs */
	for (k = xa_begin[jj]; k < xa_end[jj]; k++) {
	    krow = asub[k];
	    tmp = z_abs1(&a[k]);
	    if (tmp > *amax_col) *amax_col = tmp;
	    dense_col[krow] = a[k];
	    kmark = marker[krow];
	    if ( kmark == jj )
		continue;     /* krow visited before, go to the next nonzero */

	    /* For each unmarked nbr krow of jj
	     * krow is in L: place it in structure of L[*,jj]
	     */
	    marker[krow] = jj;
	    kperm = perm_r[krow];

	    if ( kperm == EMPTY ) {
		panel_lsub[nextl_col++] = krow; /* krow is indexed into A */
	    }
	    /*
	     * krow is in U: if its supernode-rep krep
	     * has been explored, update repfnz[*]
	     */
	    else {
		
		krep = xsup[supno[kperm]+1] - 1;
		myfnz = repfnz_col[krep];
		
#ifdef CHK_DFS
		printf("krep %d, myfnz %d, perm_r[%d] %d\n", krep, myfnz, krow, kperm);
#endif
		if ( myfnz != EMPTY ) { /* Representative visited before */
		    if ( myfnz > kperm ) repfnz_col[krep] = kperm;
		    /* continue; */
		}
		else {
		    /* Otherwise, perform dfs starting at krep */
		    oldrep = EMPTY;
		    parent[krep] = oldrep;
		    repfnz_col[krep] = kperm;
		    xdfs = xlsub[xsup[supno[krep]]];
		    maxdfs = xlsub[krep + 1];

#ifdef CHK_DFS
		    printf("  xdfs %d, maxdfs %d: ", xdfs, maxdfs);
		    for (i = xdfs; i < maxdfs; i++) printf(" %d", lsub[i]);
		    printf("\n");
#endif
		    do {
			/*
			 * For each unmarked kchild of krep
			 */
			while ( xdfs < maxdfs ) {

			    kchild = lsub[xdfs];
			    xdfs++;
			    chmark = marker[kchild];

			    if ( chmark != jj ) { /* Not reached yet */
				marker[kchild] = jj;
				chperm = perm_r[kchild];

				/* Case kchild is in L: place it in L[*,j] */
				if ( chperm == EMPTY ) {
				    panel_lsub[nextl_col++] = kchild;
				}
				/* Case kchild is in U:
				 *   chrep = its supernode-rep. If its rep has
				 *   been explored, update its repfnz[*]
				 */
				else {
				    
				    chrep = xsup[supno[chperm]+1] - 1;
				    myfnz = repfnz_col[chrep];
#ifdef CHK_DFS
				    printf("chrep %d,myfnz %d,perm_r[%d] %d\n",chrep,myfnz,kchild,chperm);
#endif
				    if ( myfnz != EMPTY ) { /* Visited before */
					if ( myfnz > chperm )
					    repfnz_col[chrep] = chperm;
				    }
				    else {
					/* Cont. dfs at snode-rep of kchild */
					xplore[krep] = xdfs;
					oldrep = krep;
					krep = chrep; /* Go deeper down G(L) */
					parent[krep] = oldrep;
					repfnz_col[krep] = chperm;
					xdfs = xlsub[xsup[supno[krep]]];
					maxdfs = xlsub[krep + 1];
#ifdef CHK_DFS
					printf("  xdfs %d, maxdfs %d: ", xdfs, maxdfs);
					for (i = xdfs; i < maxdfs; i++) printf(" %d", lsub[i]);
					printf("\n");
#endif
				    } /* else */

				} /* else */

			    } /* if... */

			} /* while xdfs < maxdfs */

			/* krow has no more unexplored nbrs:
			 *    Place snode-rep krep in postorder DFS, if this
			 *    segment is seen for the first time. (Note that
			 *    "repfnz[krep]" may change later.)
			 *    Backtrack dfs to its parent.
			 */
			if ( marker1[krep] < jcol ) {
			    segrep[*nseg] = krep;
			    ++(*nseg);
			    marker1[krep] = jj;
			}

			kpar = parent[krep]; /* Pop stack, mimic recursion */
			if ( kpar == EMPTY ) break; /* dfs done */
			krep = kpar;
			xdfs = xplore[krep];
			maxdfs = xlsub[krep + 1];

#ifdef CHK_DFS
			printf("  pop stack: krep %d,xdfs %d,maxdfs %d: ", krep,xdfs,maxdfs);
			for (i = xdfs; i < maxdfs; i++) printf(" %d", lsub[i]);
			printf("\n");
#endif
		    } while ( kpar != EMPTY ); /* do-while - until empty stack */

		} /* else */
		
	    } /* else */

	} /* for each nonz in A[*,jj] */

	repfnz_col += m;    /* Move to next column */
	dense_col += m;
	amax_col++;

    } /* for jj ... */

}
