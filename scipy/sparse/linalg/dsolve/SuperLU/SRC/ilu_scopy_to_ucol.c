/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file ilu_scopy_to_ucol.c
 * \brief Copy a computed column of U to the compressed data structure
 * and drop some small entries
 *
 * <pre>
 * -- SuperLU routine (version 4.1) --
 * Lawrence Berkeley National Laboratory
 * November, 2010
 * </pre>
 */

#include "slu_sdefs.h"

#ifdef DEBUG
int num_drop_U;
#endif

extern void scopy_(int *, float [], int *, float [], int *);

#if 0
static float *A;  /* used in _compare_ only */
static int _compare_(const void *a, const void *b)
{
    register int *x = (int *)a, *y = (int *)b;
    register double xx = fabs(A[*x]), yy = fabs(A[*y]);
    if (xx > yy) return -1;
    else if (xx < yy) return 1;
    else return 0;
}
#endif

int
ilu_scopy_to_ucol(
	      int	 jcol,	   /* in */
	      int	 nseg,	   /* in */
	      int	 *segrep,  /* in */
	      int	 *repfnz,  /* in */
	      int	 *perm_r,  /* in */
	      float	 *dense,   /* modified - reset to zero on return */
	      int  	 drop_rule,/* in */
	      milu_t	 milu,	   /* in */
	      double	 drop_tol, /* in */
	      int	 quota,    /* maximum nonzero entries allowed */
	      float	 *sum,	   /* out - the sum of dropped entries */
	      int	 *nnzUj,   /* in - out */
	      GlobalLU_t *Glu,	   /* modified */
	      float	 *work	   /* working space with minimum size n,
				    * used by the second dropping rule */
	      )
{
/*
 * Gather from SPA dense[*] to global ucol[*].
 */
    int       ksub, krep, ksupno;
    int       i, k, kfnz, segsze;
    int       fsupc, isub, irow;
    int       jsupno, nextu;
    int       new_next, mem_error;
    int       *xsup, *supno;
    int       *lsub, *xlsub;
    float    *ucol;
    int       *usub, *xusub;
    int       nzumax;
    int       m; /* number of entries in the nonzero U-segments */
    register float d_max = 0.0, d_min = 1.0 / smach("Safe minimum");
    register double tmp;
    float zero = 0.0;
    int i_1 = 1;

    xsup    = Glu->xsup;
    supno   = Glu->supno;
    lsub    = Glu->lsub;
    xlsub   = Glu->xlsub;
    ucol    = (float *) Glu->ucol;
    usub    = Glu->usub;
    xusub   = Glu->xusub;
    nzumax  = Glu->nzumax;

    *sum = zero;
    if (drop_rule == NODROP) {
	drop_tol = -1.0, quota = Glu->n;
    }

    jsupno = supno[jcol];
    nextu  = xusub[jcol];
    k = nseg - 1;
    for (ksub = 0; ksub < nseg; ksub++) {
	krep = segrep[k--];
	ksupno = supno[krep];

	if ( ksupno != jsupno ) { /* Should go into ucol[] */
	    kfnz = repfnz[krep];
	    if ( kfnz != EMPTY ) {	/* Nonzero U-segment */

		fsupc = xsup[ksupno];
		isub = xlsub[fsupc] + kfnz - fsupc;
		segsze = krep - kfnz + 1;

		new_next = nextu + segsze;
		while ( new_next > nzumax ) {
		    if ((mem_error = sLUMemXpand(jcol, nextu, UCOL, &nzumax,
			    Glu)) != 0)
			return (mem_error);
		    ucol = Glu->ucol;
		    if ((mem_error = sLUMemXpand(jcol, nextu, USUB, &nzumax,
			    Glu)) != 0)
			return (mem_error);
		    usub = Glu->usub;
		    lsub = Glu->lsub;
		}

		for (i = 0; i < segsze; i++) {
		    irow = lsub[isub++];
		    tmp = fabs(dense[irow]);

		    /* first dropping rule */
		    if (quota > 0 && tmp >= drop_tol) {
			if (tmp > d_max) d_max = tmp;
			if (tmp < d_min) d_min = tmp;
			usub[nextu] = perm_r[irow];
			ucol[nextu] = dense[irow];
			nextu++;
		    } else {
			switch (milu) {
			    case SMILU_1:
			    case SMILU_2:
				*sum += dense[irow];
				break;
			    case SMILU_3:
				/* *sum += fabs(dense[irow]);*/
				*sum += tmp;
				break;
			    case SILU:
			    default:
				break;
			}
#ifdef DEBUG
			num_drop_U++;
#endif
		    }
		    dense[irow] = zero;
		}

	    }

	}

    } /* for each segment... */

    xusub[jcol + 1] = nextu;	  /* Close U[*,jcol] */
    m = xusub[jcol + 1] - xusub[jcol];

    /* second dropping rule */
    if (drop_rule & DROP_SECONDARY && m > quota) {
	register double tol = d_max;
	register int m0 = xusub[jcol] + m - 1;

	if (quota > 0) {
	    if (drop_rule & DROP_INTERP) {
		d_max = 1.0 / d_max; d_min = 1.0 / d_min;
		tol = 1.0 / (d_max + (d_min - d_max) * quota / m);
	    } else {
		scopy_(&m, &ucol[xusub[jcol]], &i_1, work, &i_1);
		tol = sqselect(m, work, quota);
#if 0
		A = &ucol[xusub[jcol]];
		for (i = 0; i < m; i++) work[i] = i;
		qsort(work, m, sizeof(int), _compare_);
		tol = fabs(usub[xusub[jcol] + work[quota]]);
#endif
	    }
	}
	for (i = xusub[jcol]; i <= m0; ) {
	    if (fabs(ucol[i]) <= tol) {
		switch (milu) {
		    case SMILU_1:
		    case SMILU_2:
			*sum += ucol[i];
			break;
		    case SMILU_3:
			*sum += fabs(ucol[i]);
			break;
		    case SILU:
		    default:
			break;
		}
		ucol[i] = ucol[m0];
		usub[i] = usub[m0];
		m0--;
		m--;
#ifdef DEBUG
		num_drop_U++;
#endif
		xusub[jcol + 1]--;
		continue;
	    }
	    i++;
	}
    }

    if (milu == SMILU_2) *sum = fabs(*sum);

    *nnzUj += m;

    return 0;
}
