
/*! @file ilu_ccopy_to_ucol.c
 * \brief Copy a computed column of U to the compressed data structure
 * and drop some small entries
 *
 * <pre>
 * -- SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Laboratory
 * June 30, 2009
 * </pre>
 */

#include "slu_cdefs.h"

#ifdef DEBUG
int num_drop_U;
#endif

static complex *A;  /* used in _compare_ only */
static int _compare_(const void *a, const void *b)
{
    register int *x = (int *)a, *y = (int *)b;
    register float xx = slu_c_abs1(&A[*x]), yy = slu_c_abs1(&A[*y]);
    if (xx > yy) return -1;
    else if (xx < yy) return 1;
    else return 0;
}


int
ilu_ccopy_to_ucol(
	      int	 jcol,	   /* in */
	      int	 nseg,	   /* in */
	      int	 *segrep,  /* in */
	      int	 *repfnz,  /* in */
	      int	 *perm_r,  /* in */
	      complex	 *dense,   /* modified - reset to zero on return */
	      int  	 drop_rule,/* in */
	      milu_t	 milu,	   /* in */
	      double	 drop_tol, /* in */
	      int	 quota,    /* maximum nonzero entries allowed */
	      complex	 *sum,	   /* out - the sum of dropped entries */
	      int	 *nnzUj,   /* in - out */
	      GlobalLU_t *Glu,	   /* modified */
	      int	 *work	   /* working space with minimum size n,
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
    complex    *ucol;
    int       *usub, *xusub;
    int       nzumax;
    int       m; /* number of entries in the nonzero U-segments */
    register float d_max = 0.0, d_min = 1.0 / dlamch_("Safe minimum");
    register double tmp;
    complex zero = {0.0, 0.0};

    xsup    = Glu->xsup;
    supno   = Glu->supno;
    lsub    = Glu->lsub;
    xlsub   = Glu->xlsub;
    ucol    = Glu->ucol;
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
		    if ((mem_error = cLUMemXpand(jcol, nextu, UCOL, &nzumax,
			    Glu)) != 0)
			return (mem_error);
		    ucol = Glu->ucol;
		    if ((mem_error = cLUMemXpand(jcol, nextu, USUB, &nzumax,
			    Glu)) != 0)
			return (mem_error);
		    usub = Glu->usub;
		    lsub = Glu->lsub;
		}

		for (i = 0; i < segsze; i++) {
		    irow = lsub[isub++];
         	    tmp = slu_c_abs1(&dense[irow]);

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
                                c_add(sum, sum, &dense[irow]);
				break;
			    case SMILU_3:
				/* *sum += fabs(dense[irow]);*/
				sum->r += tmp;
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
		A = &ucol[xusub[jcol]];
		for (i = 0; i < m; i++) work[i] = i;
		qsort(work, m, sizeof(int), _compare_);
		tol = fabs(usub[xusub[jcol] + work[quota]]);
	    }
	}
	for (i = xusub[jcol]; i <= m0; ) {
	    if (slu_c_abs1(&ucol[i]) <= tol) {
		switch (milu) {
		    case SMILU_1:
		    case SMILU_2:
			c_add(sum, sum, &ucol[i]);
			break;
		    case SMILU_3:
			sum->r += tmp;
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

    if (milu == SMILU_2) {
        sum->r = slu_c_abs1(sum); sum->i = 0.0;
    }
    if (milu == SMILU_3) sum->i = 0.0;

    *nnzUj += m;

    return 0;
}
