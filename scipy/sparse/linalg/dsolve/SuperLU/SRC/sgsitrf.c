/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file sgsitrf.c
 * \brief Computes an ILU factorization of a general sparse matrix
 *
 * <pre>
 * -- SuperLU routine (version 4.1) --
 * Lawrence Berkeley National Laboratory.
 * June 30, 2009
 *
 * </pre>
 */

#include "slu_sdefs.h"

#ifdef DEBUG
int num_drop_L;
#endif

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * SGSITRF computes an ILU factorization of a general sparse m-by-n
 * matrix A using partial pivoting with row interchanges.
 * The factorization has the form
 *     Pr * A = L * U
 * where Pr is a row permutation matrix, L is lower triangular with unit
 * diagonal elements (lower trapezoidal if A->nrow > A->ncol), and U is upper
 * triangular (upper trapezoidal if A->nrow < A->ncol).
 *
 * See supermatrix.h for the definition of 'SuperMatrix' structure.
 *
 * Arguments
 * =========
 *
 * options (input) superlu_options_t*
 *	   The structure defines the input parameters to control
 *	   how the ILU decomposition will be performed.
 *
 * A	    (input) SuperMatrix*
 *	    Original matrix A, permuted by columns, of dimension
 *	    (A->nrow, A->ncol). The type of A can be:
 *	    Stype = SLU_NCP; Dtype = SLU_S; Mtype = SLU_GE.
 *
 * relax    (input) int
 *	    To control degree of relaxing supernodes. If the number
 *	    of nodes (columns) in a subtree of the elimination tree is less
 *	    than relax, this subtree is considered as one supernode,
 *	    regardless of the row structures of those columns.
 *
 * panel_size (input) int
 *	    A panel consists of at most panel_size consecutive columns.
 *
 * etree    (input) int*, dimension (A->ncol)
 *	    Elimination tree of A'*A.
 *	    Note: etree is a vector of parent pointers for a forest whose
 *	    vertices are the integers 0 to A->ncol-1; etree[root]==A->ncol.
 *	    On input, the columns of A should be permuted so that the
 *	    etree is in a certain postorder.
 *
 * work     (input/output) void*, size (lwork) (in bytes)
 *	    User-supplied work space and space for the output data structures.
 *	    Not referenced if lwork = 0;
 *
 * lwork   (input) int
 *	   Specifies the size of work array in bytes.
 *	   = 0:  allocate space internally by system malloc;
 *	   > 0:  use user-supplied work array of length lwork in bytes,
 *		 returns error if space runs out.
 *	   = -1: the routine guesses the amount of space needed without
 *		 performing the factorization, and returns it in
 *		 *info; no other side effects.
 *
 * perm_c   (input) int*, dimension (A->ncol)
 *	    Column permutation vector, which defines the
 *	    permutation matrix Pc; perm_c[i] = j means column i of A is
 *	    in position j in A*Pc.
 *	    When searching for diagonal, perm_c[*] is applied to the
 *	    row subscripts of A, so that diagonal threshold pivoting
 *	    can find the diagonal of A, rather than that of A*Pc.
 *
 * perm_r   (input/output) int*, dimension (A->nrow)
 *	    Row permutation vector which defines the permutation matrix Pr,
 *	    perm_r[i] = j means row i of A is in position j in Pr*A.
 *	    If options->Fact = SamePattern_SameRowPerm, the pivoting routine
 *	       will try to use the input perm_r, unless a certain threshold
 *	       criterion is violated. In that case, perm_r is overwritten by
 *	       a new permutation determined by partial pivoting or diagonal
 *	       threshold pivoting.
 *	    Otherwise, perm_r is output argument;
 *
 * L	    (output) SuperMatrix*
 *	    The factor L from the factorization Pr*A=L*U; use compressed row
 *	    subscripts storage for supernodes, i.e., L has type:
 *	    Stype = SLU_SC, Dtype = SLU_S, Mtype = SLU_TRLU.
 *
 * U	    (output) SuperMatrix*
 *	    The factor U from the factorization Pr*A*Pc=L*U. Use column-wise
 *	    storage scheme, i.e., U has types: Stype = SLU_NC,
 *	    Dtype = SLU_S, Mtype = SLU_TRU.
 *
 * Glu      (input/output) GlobalLU_t *
 *          If options->Fact == SamePattern_SameRowPerm, it is an input;
 *              The matrix A will be factorized assuming that a 
 *              factorization of a matrix with the same sparsity pattern
 *              and similar numerical values was performed prior to this one.
 *              Therefore, this factorization will reuse both row and column
 *		scaling factors R and C, both row and column permutation
 *		vectors perm_r and perm_c, and the L & U data structures
 *		set up from the previous factorization.
 *          Otherwise, it is an output.
 *
 * stat     (output) SuperLUStat_t*
 *	    Record the statistics on runtime and floating-point operation count.
 *	    See slu_util.h for the definition of 'SuperLUStat_t'.
 *
 * info     (output) int*
 *	    = 0: successful exit
 *	    < 0: if info = -i, the i-th argument had an illegal value
 *	    > 0: if info = i, and i is
 *	       <= A->ncol: number of zero pivots. They are replaced by small
 *		  entries according to options->ILU_FillTol.
 *	       > A->ncol: number of bytes allocated when memory allocation
 *		  failure occurred, plus A->ncol. If lwork = -1, it is
 *		  the estimated amount of space needed, plus A->ncol.
 *
 * ======================================================================
 *
 * Local Working Arrays:
 * ======================
 *   m = number of rows in the matrix
 *   n = number of columns in the matrix
 *
 *   marker[0:3*m-1]: marker[i] = j means that node i has been
 *	reached when working on column j.
 *	Storage: relative to original row subscripts
 *	NOTE: There are 4 of them:
 *	      marker/marker1 are used for panel dfs, see (ilu_)dpanel_dfs.c;
 *	      marker2 is used for inner-factorization, see (ilu)_dcolumn_dfs.c;
 *	      marker_relax(has its own space) is used for relaxed supernodes.
 *
 *   parent[0:m-1]: parent vector used during dfs
 *	Storage: relative to new row subscripts
 *
 *   xplore[0:m-1]: xplore[i] gives the location of the next (dfs)
 *	unexplored neighbor of i in lsub[*]
 *
 *   segrep[0:nseg-1]: contains the list of supernodal representatives
 *	in topological order of the dfs. A supernode representative is the
 *	last column of a supernode.
 *	The maximum size of segrep[] is n.
 *
 *   repfnz[0:W*m-1]: for a nonzero segment U[*,j] that ends at a
 *	supernodal representative r, repfnz[r] is the location of the first
 *	nonzero in this segment.  It is also used during the dfs: repfnz[r]>0
 *	indicates the supernode r has been explored.
 *	NOTE: There are W of them, each used for one column of a panel.
 *
 *   panel_lsub[0:W*m-1]: temporary for the nonzeros row indices below
 *	the panel diagonal. These are filled in during dpanel_dfs(), and are
 *	used later in the inner LU factorization within the panel.
 *	panel_lsub[]/dense[] pair forms the SPA data structure.
 *	NOTE: There are W of them.
 *
 *   dense[0:W*m-1]: sparse accumulating (SPA) vector for intermediate values;
 *		   NOTE: there are W of them.
 *
 *   tempv[0:*]: real temporary used for dense numeric kernels;
 *	The size of this array is defined by NUM_TEMPV() in slu_util.h.
 *	It is also used by the dropping routine ilu_ddrop_row().
 * </pre>
 */

void
sgsitrf(superlu_options_t *options, SuperMatrix *A, int relax, int panel_size,
	int *etree, void *work, int lwork, int *perm_c, int *perm_r,
	SuperMatrix *L, SuperMatrix *U, 
    	GlobalLU_t *Glu, /* persistent to facilitate multiple factorizations */
	SuperLUStat_t *stat, int *info)
{
    /* Local working arrays */
    NCPformat *Astore;
    int       *iperm_r = NULL; /* inverse of perm_r; used when
				  options->Fact == SamePattern_SameRowPerm */
    int       *iperm_c; /* inverse of perm_c */
    int       *swap, *iswap; /* swap is used to store the row permutation
				during the factorization. Initially, it is set
				to iperm_c (row indeces of Pc*A*Pc').
				iswap is the inverse of swap. After the
				factorization, it is equal to perm_r. */
    int       *iwork;
    float   *swork;
    int       *segrep, *repfnz, *parent, *xplore;
    int       *panel_lsub; /* dense[]/panel_lsub[] pair forms a w-wide SPA */
    int       *marker, *marker_relax;
    float    *dense, *tempv;
    int       *relax_end, *relax_fsupc;
    float    *a;
    int       *asub;
    int       *xa_begin, *xa_end;
    int       *xsup, *supno;
    int       *xlsub, *xlusup, *xusub;
    int       nzlumax;
    float    *amax; 
    float    drop_sum;
    float alpha, omega;  /* used in MILU, mimicing DRIC */
    float    *swork2;	   /* used by the second dropping rule */

    /* Local scalars */
    fact_t    fact = options->Fact;
    double    diag_pivot_thresh = options->DiagPivotThresh;
    double    drop_tol = options->ILU_DropTol; /* tau */
    double    fill_ini = options->ILU_FillTol; /* tau^hat */
    double    gamma = options->ILU_FillFactor;
    int       drop_rule = options->ILU_DropRule;
    milu_t    milu = options->ILU_MILU;
    double    fill_tol;
    int       pivrow;	/* pivotal row number in the original matrix A */
    int       nseg1;	/* no of segments in U-column above panel row jcol */
    int       nseg;	/* no of segments in each U-column */
    register int jcol;
    register int kcol;	/* end column of a relaxed snode */
    register int icol;
    register int i, k, jj, new_next, iinfo;
    int       m, n, min_mn, jsupno, fsupc, nextlu, nextu;
    int       w_def;	/* upper bound on panel width */
    int       usepr, iperm_r_allocated = 0;
    int       nnzL, nnzU;
    int       *panel_histo = stat->panel_histo;
    flops_t   *ops = stat->ops;

    int       last_drop;/* the last column which the dropping rules applied */
    int       quota;
    int       nnzAj;	/* number of nonzeros in A(:,1:j) */
    int       nnzLj, nnzUj;
    double    tol_L = drop_tol, tol_U = drop_tol;
    float zero = 0.0;
    float one = 1.0;

    /* Executable */	   
    iinfo    = 0;
    m	     = A->nrow;
    n	     = A->ncol;
    min_mn   = SUPERLU_MIN(m, n);
    Astore   = A->Store;
    a	     = Astore->nzval;
    asub     = Astore->rowind;
    xa_begin = Astore->colbeg;
    xa_end   = Astore->colend;

    /* Allocate storage common to the factor routines */
    *info = sLUMemInit(fact, work, lwork, m, n, Astore->nnz, panel_size,
		       gamma, L, U, Glu, &iwork, &swork);
    if ( *info ) return;

    xsup    = Glu->xsup;
    supno   = Glu->supno;
    xlsub   = Glu->xlsub;
    xlusup  = Glu->xlusup;
    xusub   = Glu->xusub;

    SetIWork(m, n, panel_size, iwork, &segrep, &parent, &xplore,
	     &repfnz, &panel_lsub, &marker_relax, &marker);
    sSetRWork(m, panel_size, swork, &dense, &tempv);

    usepr = (fact == SamePattern_SameRowPerm);
    if ( usepr ) {
	/* Compute the inverse of perm_r */
	iperm_r = (int *) intMalloc(m);
	for (k = 0; k < m; ++k) iperm_r[perm_r[k]] = k;
	iperm_r_allocated = 1;
    }

    iperm_c = (int *) intMalloc(n);
    for (k = 0; k < n; ++k) iperm_c[perm_c[k]] = k;
    swap = (int *)intMalloc(n);
    for (k = 0; k < n; k++) swap[k] = iperm_c[k];
    iswap = (int *)intMalloc(n);
    for (k = 0; k < n; k++) iswap[k] = perm_c[k];
    amax = (float *) floatMalloc(panel_size);
    if (drop_rule & DROP_SECONDARY)
	swork2 = (float *)floatMalloc(n);
    else
	swork2 = NULL;

    nnzAj = 0;
    nnzLj = 0;
    nnzUj = 0;
    last_drop = SUPERLU_MAX(min_mn - 2 * sp_ienv(7), (int)(min_mn * 0.95));
    alpha = pow((double)n, -1.0 / options->ILU_MILU_Dim);

    /* Identify relaxed snodes */
    relax_end = (int *) intMalloc(n);
    relax_fsupc = (int *) intMalloc(n);
    if ( options->SymmetricMode == YES )
	ilu_heap_relax_snode(n, etree, relax, marker, relax_end, relax_fsupc);
    else
	ilu_relax_snode(n, etree, relax, marker, relax_end, relax_fsupc);

    ifill (perm_r, m, EMPTY);
    ifill (marker, m * NO_MARKER, EMPTY);
    supno[0] = -1;
    xsup[0]  = xlsub[0] = xusub[0] = xlusup[0] = 0;
    w_def    = panel_size;

    /* Mark the rows used by relaxed supernodes */
    ifill (marker_relax, m, EMPTY);
    i = mark_relax(m, relax_end, relax_fsupc, xa_begin, xa_end,
	         asub, marker_relax);
#if ( PRNTlevel >= 1)
    printf("%d relaxed supernodes.\n", i);
#endif

    /*
     * Work on one "panel" at a time. A panel is one of the following:
     *	   (a) a relaxed supernode at the bottom of the etree, or
     *	   (b) panel_size contiguous columns, defined by the user
     */
    for (jcol = 0; jcol < min_mn; ) {

	if ( relax_end[jcol] != EMPTY ) { /* start of a relaxed snode */
	    kcol = relax_end[jcol];	  /* end of the relaxed snode */
	    panel_histo[kcol-jcol+1]++;

	    /* Drop small rows in the previous supernode. */
	    if (jcol > 0 && jcol < last_drop) {
		int first = xsup[supno[jcol - 1]];
		int last = jcol - 1;
		int quota;

		/* Compute the quota */
		if (drop_rule & DROP_PROWS)
		    quota = gamma * Astore->nnz / m * (m - first) / m
			    * (last - first + 1);
		else if (drop_rule & DROP_COLUMN) {
		    int i;
		    quota = 0;
		    for (i = first; i <= last; i++)
			quota += xa_end[i] - xa_begin[i];
		    quota = gamma * quota * (m - first) / m;
		} else if (drop_rule & DROP_AREA)
		    quota = gamma * nnzAj * (1.0 - 0.5 * (last + 1.0) / m)
			    - nnzLj;
		else
		    quota = m * n;
		fill_tol = pow(fill_ini, 1.0 - 0.5 * (first + last) / min_mn);

		/* Drop small rows */
		i = ilu_sdrop_row(options, first, last, tol_L, quota, &nnzLj,
				  &fill_tol, Glu, tempv, swork2, 0);
		/* Reset the parameters */
		if (drop_rule & DROP_DYNAMIC) {
		    if (gamma * nnzAj * (1.0 - 0.5 * (last + 1.0) / m)
			     < nnzLj)
			tol_L = SUPERLU_MIN(1.0, tol_L * 2.0);
		    else
			tol_L = SUPERLU_MAX(drop_tol, tol_L * 0.5);
		}
		if (fill_tol < 0) iinfo -= (int)fill_tol;
#ifdef DEBUG
		num_drop_L += i * (last - first + 1);
#endif
	    }

	    /* --------------------------------------
	     * Factorize the relaxed supernode(jcol:kcol)
	     * -------------------------------------- */
	    /* Determine the union of the row structure of the snode */
	    if ( (*info = ilu_ssnode_dfs(jcol, kcol, asub, xa_begin, xa_end,
					 marker, Glu)) != 0 )
		return;

	    nextu    = xusub[jcol];
	    nextlu   = xlusup[jcol];
	    jsupno   = supno[jcol];
	    fsupc    = xsup[jsupno];
	    new_next = nextlu + (xlsub[fsupc+1]-xlsub[fsupc])*(kcol-jcol+1);
	    nzlumax = Glu->nzlumax;
	    while ( new_next > nzlumax ) {
		if ((*info = sLUMemXpand(jcol, nextlu, LUSUP, &nzlumax, Glu)))
		    return;
	    }

	    for (icol = jcol; icol <= kcol; icol++) {
		xusub[icol+1] = nextu;

		amax[0] = 0.0;
		/* Scatter into SPA dense[*] */
		for (k = xa_begin[icol]; k < xa_end[icol]; k++) {
		    register float tmp = fabs(a[k]);
		    if (tmp > amax[0]) amax[0] = tmp;
		    dense[asub[k]] = a[k];
		}
		nnzAj += xa_end[icol] - xa_begin[icol];
		if (amax[0] == 0.0) {
		    amax[0] = fill_ini;
#if ( PRNTlevel >= 1)
		    printf("Column %d is entirely zero!\n", icol);
		    fflush(stdout);
#endif
		}

		/* Numeric update within the snode */
		ssnode_bmod(icol, jsupno, fsupc, dense, tempv, Glu, stat);

		if (usepr) pivrow = iperm_r[icol];
		fill_tol = pow(fill_ini, 1.0 - (double)icol / (double)min_mn);
		if ( (*info = ilu_spivotL(icol, diag_pivot_thresh, &usepr,
					  perm_r, iperm_c[icol], swap, iswap,
					  marker_relax, &pivrow,
                                          amax[0] * fill_tol, milu, zero,
                                          Glu, stat)) ) {
		    iinfo++;
		    marker[pivrow] = kcol;
		}

	    }

	    jcol = kcol + 1;

	} else { /* Work on one panel of panel_size columns */

	    /* Adjust panel_size so that a panel won't overlap with the next
	     * relaxed snode.
	     */
	    panel_size = w_def;
	    for (k = jcol + 1; k < SUPERLU_MIN(jcol+panel_size, min_mn); k++)
		if ( relax_end[k] != EMPTY ) {
		    panel_size = k - jcol;
		    break;
		}
	    if ( k == min_mn ) panel_size = min_mn - jcol;
	    panel_histo[panel_size]++;

	    /* symbolic factor on a panel of columns */
	    ilu_spanel_dfs(m, panel_size, jcol, A, perm_r, &nseg1,
                          dense, amax, panel_lsub, segrep, repfnz,
                          marker, parent, xplore, Glu);

	    /* numeric sup-panel updates in topological order */
	    spanel_bmod(m, panel_size, jcol, nseg1, dense,
			tempv, segrep, repfnz, Glu, stat);

	    /* Sparse LU within the panel, and below panel diagonal */
	    for (jj = jcol; jj < jcol + panel_size; jj++) {

		k = (jj - jcol) * m; /* column index for w-wide arrays */

		nseg = nseg1;	/* Begin after all the panel segments */

		nnzAj += xa_end[jj] - xa_begin[jj];

		if ((*info = ilu_scolumn_dfs(m, jj, perm_r, &nseg,
					     &panel_lsub[k], segrep, &repfnz[k],
					     marker, parent, xplore, Glu)))
		    return;

		/* Numeric updates */
		if ((*info = scolumn_bmod(jj, (nseg - nseg1), &dense[k],
					  tempv, &segrep[nseg1], &repfnz[k],
					  jcol, Glu, stat)) != 0) return;

		/* Make a fill-in position if the column is entirely zero */
		if (xlsub[jj + 1] == xlsub[jj]) {
		    register int i, row;
		    int nextl;
		    int nzlmax = Glu->nzlmax;
		    int *lsub = Glu->lsub;
		    int *marker2 = marker + 2 * m;

		    /* Allocate memory */
		    nextl = xlsub[jj] + 1;
		    if (nextl >= nzlmax) {
			int error = sLUMemXpand(jj, nextl, LSUB, &nzlmax, Glu);
			if (error) { *info = error; return; }
			lsub = Glu->lsub;
		    }
		    xlsub[jj + 1]++;
		    assert(xlusup[jj]==xlusup[jj+1]);
		    xlusup[jj + 1]++;
		    ((float *) Glu->lusup)[xlusup[jj]] = zero;

		    /* Choose a row index (pivrow) for fill-in */
		    for (i = jj; i < n; i++)
			if (marker_relax[swap[i]] <= jj) break;
		    row = swap[i];
		    marker2[row] = jj;
		    lsub[xlsub[jj]] = row;
#ifdef DEBUG
		    printf("Fill col %d.\n", jj);
		    fflush(stdout);
#endif
		}

		/* Computer the quota */
		if (drop_rule & DROP_PROWS)
		    quota = gamma * Astore->nnz / m * jj / m;
		else if (drop_rule & DROP_COLUMN)
		    quota = gamma * (xa_end[jj] - xa_begin[jj]) *
			    (jj + 1) / m;
		else if (drop_rule & DROP_AREA)
		    quota = gamma * 0.9 * nnzAj * 0.5 - nnzUj;
		else
		    quota = m;

		/* Copy the U-segments to ucol[*] and drop small entries */
		if ((*info = ilu_scopy_to_ucol(jj, nseg, segrep, &repfnz[k],
					       perm_r, &dense[k], drop_rule,
					       milu, amax[jj - jcol] * tol_U,
					       quota, &drop_sum, &nnzUj, Glu,
					       swork2)) != 0)
		    return;

		/* Reset the dropping threshold if required */
		if (drop_rule & DROP_DYNAMIC) {
		    if (gamma * 0.9 * nnzAj * 0.5 < nnzLj)
			tol_U = SUPERLU_MIN(1.0, tol_U * 2.0);
		    else
			tol_U = SUPERLU_MAX(drop_tol, tol_U * 0.5);
		}

		if (drop_sum != zero)
		{
		    if (drop_sum > zero)
			omega = SUPERLU_MIN(2.0 * (1.0 - alpha)
				* amax[jj - jcol] / drop_sum, one);
		    else
			omega = SUPERLU_MAX(2.0 * (1.0 - alpha)
				* amax[jj - jcol] / drop_sum, -one);
		    drop_sum *= omega;
                }
		if (usepr) pivrow = iperm_r[jj];
		fill_tol = pow(fill_ini, 1.0 - (double)jj / (double)min_mn);
		if ( (*info = ilu_spivotL(jj, diag_pivot_thresh, &usepr, perm_r,
					  iperm_c[jj], swap, iswap,
					  marker_relax, &pivrow,
					  amax[jj - jcol] * fill_tol, milu,
					  drop_sum, Glu, stat)) ) {
		    iinfo++;
		    marker[m + pivrow] = jj;
		    marker[2 * m + pivrow] = jj;
		}

		/* Reset repfnz[] for this column */
		resetrep_col (nseg, segrep, &repfnz[k]);

		/* Start a new supernode, drop the previous one */
		if (jj > 0 && supno[jj] > supno[jj - 1] && jj < last_drop) {
		    int first = xsup[supno[jj - 1]];
		    int last = jj - 1;
		    int quota;

		    /* Compute the quota */
		    if (drop_rule & DROP_PROWS)
			quota = gamma * Astore->nnz / m * (m - first) / m
				* (last - first + 1);
		    else if (drop_rule & DROP_COLUMN) {
			int i;
			quota = 0;
			for (i = first; i <= last; i++)
			    quota += xa_end[i] - xa_begin[i];
			quota = gamma * quota * (m - first) / m;
		    } else if (drop_rule & DROP_AREA)
			quota = gamma * nnzAj * (1.0 - 0.5 * (last + 1.0)
				/ m) - nnzLj;
		    else
			quota = m * n;
		    fill_tol = pow(fill_ini, 1.0 - 0.5 * (first + last) /
			    (double)min_mn);

		    /* Drop small rows */
		    i = ilu_sdrop_row(options, first, last, tol_L, quota,
				      &nnzLj, &fill_tol, Glu, tempv, swork2,
				      1);

		    /* Reset the parameters */
		    if (drop_rule & DROP_DYNAMIC) {
			if (gamma * nnzAj * (1.0 - 0.5 * (last + 1.0) / m)
				< nnzLj)
			    tol_L = SUPERLU_MIN(1.0, tol_L * 2.0);
			else
			    tol_L = SUPERLU_MAX(drop_tol, tol_L * 0.5);
		    }
		    if (fill_tol < 0) iinfo -= (int)fill_tol;
#ifdef DEBUG
		    num_drop_L += i * (last - first + 1);
#endif
		} /* if start a new supernode */

	    } /* for */

	    jcol += panel_size; /* Move to the next panel */

	} /* else */

    } /* for */

    *info = iinfo;

    if ( m > n ) {
	k = 0;
	for (i = 0; i < m; ++i)
	    if ( perm_r[i] == EMPTY ) {
		perm_r[i] = n + k;
		++k;
	    }
    }

    ilu_countnz(min_mn, &nnzL, &nnzU, Glu);
    fixupL(min_mn, perm_r, Glu);

    sLUWorkFree(iwork, swork, Glu); /* Free work space and compress storage */

    if ( fact == SamePattern_SameRowPerm ) {
	/* L and U structures may have changed due to possibly different
	   pivoting, even though the storage is available.
	   There could also be memory expansions, so the array locations
	   may have changed, */
	((SCformat *)L->Store)->nnz = nnzL;
	((SCformat *)L->Store)->nsuper = Glu->supno[n];
	((SCformat *)L->Store)->nzval = (float *) Glu->lusup;
	((SCformat *)L->Store)->nzval_colptr = Glu->xlusup;
	((SCformat *)L->Store)->rowind = Glu->lsub;
	((SCformat *)L->Store)->rowind_colptr = Glu->xlsub;
	((NCformat *)U->Store)->nnz = nnzU;
	((NCformat *)U->Store)->nzval = (float *) Glu->ucol;
	((NCformat *)U->Store)->rowind = Glu->usub;
	((NCformat *)U->Store)->colptr = Glu->xusub;
    } else {
	sCreate_SuperNode_Matrix(L, A->nrow, min_mn, nnzL,
              (float *) Glu->lusup, Glu->xlusup,
              Glu->lsub, Glu->xlsub, Glu->supno, Glu->xsup,
	      SLU_SC, SLU_S, SLU_TRLU);
	sCreate_CompCol_Matrix(U, min_mn, min_mn, nnzU,
	      (float *) Glu->ucol, Glu->usub, Glu->xusub,
	      SLU_NC, SLU_S, SLU_TRU);
    }

    ops[FACT] += ops[TRSV] + ops[GEMV];
    stat->expansions = --(Glu->num_expansions);

    if ( iperm_r_allocated ) SUPERLU_FREE (iperm_r);
    SUPERLU_FREE (iperm_c);
    SUPERLU_FREE (relax_end);
    SUPERLU_FREE (swap);
    SUPERLU_FREE (iswap);
    SUPERLU_FREE (relax_fsupc);
    SUPERLU_FREE (amax);
    if ( swork2 ) SUPERLU_FREE (swork2);

}
