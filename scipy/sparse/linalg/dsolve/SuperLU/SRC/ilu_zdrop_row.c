
/*! @file ilu_zdrop_row.c
 * \brief Drop small rows from L
 *
 * <pre>
 * -- SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Laboratory.
 * June 30, 2009
 * <\pre>
 */

#include <math.h>
#include <stdlib.h>
#include "slu_zdefs.h"

extern void zswap_(int *, doublecomplex [], int *, doublecomplex [], int *);
extern void zaxpy_(int *, doublecomplex *, doublecomplex [], int *, doublecomplex [], int *);

static double *A;  /* used in _compare_ only */
static int _compare_(const void *a, const void *b)
{
    register int *x = (int *)a, *y = (int *)b;
    if (A[*x] - A[*y] > 0.0) return -1;
    else if (A[*x] - A[*y] < 0.0) return 1;
    else return 0;
}

/*! \brief
 * <pre>
 * Purpose
 * =======
 *    ilu_zdrop_row() - Drop some small rows from the previous 
 *    supernode (L-part only).
 * </pre>
 */
int ilu_zdrop_row(
	superlu_options_t *options, /* options */
	int    first,	    /* index of the first column in the supernode */
	int    last,	    /* index of the last column in the supernode */
	double drop_tol,    /* dropping parameter */
	int    quota,	    /* maximum nonzero entries allowed */
	int    *nnzLj,	    /* in/out number of nonzeros in L(:, 1:last) */
	double *fill_tol,   /* in/out - on exit, fill_tol=-num_zero_pivots,
			     * does not change if options->ILU_MILU != SMILU1 */
	GlobalLU_t *Glu,    /* modified */
	double dwork[],     /* working space with minimum size last-first+1 */
	int    iwork[],     /* working space with minimum size m - n,
			     * used by the second dropping rule */
	int    lastc	    /* if lastc == 0, there is nothing after the
			     * working supernode [first:last];
			     * if lastc == 1, there is one more column after
			     * the working supernode. */ )
{
    register int i, j, k, m1;
    register int nzlc; /* number of nonzeros in column last+1 */
    register int xlusup_first, xlsub_first;
    int m, n; /* m x n is the size of the supernode */
    int r = 0; /* number of dropped rows */
    register double *temp;
    register doublecomplex *lusup = Glu->lusup;
    register int *lsub = Glu->lsub;
    register int *xlsub = Glu->xlsub;
    register int *xlusup = Glu->xlusup;
    register double d_max = 0.0, d_min = 1.0;
    int    drop_rule = options->ILU_DropRule;
    milu_t milu = options->ILU_MILU;
    norm_t nrm = options->ILU_Norm;
    doublecomplex zero = {0.0, 0.0};
    doublecomplex one = {1.0, 0.0};
    doublecomplex none = {-1.0, 0.0};
    int inc_diag; /* inc_diag = m + 1 */
    int nzp = 0;  /* number of zero pivots */

    xlusup_first = xlusup[first];
    xlsub_first = xlsub[first];
    m = xlusup[first + 1] - xlusup_first;
    n = last - first + 1;
    m1 = m - 1;
    inc_diag = m + 1;
    nzlc = lastc ? (xlusup[last + 2] - xlusup[last + 1]) : 0;
    temp = dwork - n;

    /* Quick return if nothing to do. */
    if (m == 0 || m == n || drop_rule == NODROP)
    {
	*nnzLj += m * n;
	return 0;
    }

    /* basic dropping: ILU(tau) */
    for (i = n; i <= m1; )
    {
	/* the average abs value of ith row */
	switch (nrm)
	{
	    case ONE_NORM:
		temp[i] = dzasum_(&n, &lusup[xlusup_first + i], &m) / (double)n;
		break;
	    case TWO_NORM:
		temp[i] = dznrm2_(&n, &lusup[xlusup_first + i], &m)
		    / sqrt((double)n);
		break;
	    case INF_NORM:
	    default:
		k = izamax_(&n, &lusup[xlusup_first + i], &m) - 1;
		temp[i] = z_abs1(&lusup[xlusup_first + i + m * k]);
		break;
	}

	/* drop small entries due to drop_tol */
	if (drop_rule & DROP_BASIC && temp[i] < drop_tol)
	{
	    r++;
	    /* drop the current row and move the last undropped row here */
	    if (r > 1) /* add to last row */
	    {
		/* accumulate the sum (for MILU) */
		switch (milu)
		{
		    case SMILU_1:
		    case SMILU_2:
			zaxpy_(&n, &one, &lusup[xlusup_first + i], &m,
				&lusup[xlusup_first + m - 1], &m);
			break;
		    case SMILU_3:
			for (j = 0; j < n; j++)
			    lusup[xlusup_first + (m - 1) + j * m].r +=
				    z_abs1(&lusup[xlusup_first + i + j * m]);
			break;
		    case SILU:
		    default:
			break;
		}
		zcopy_(&n, &lusup[xlusup_first + m1], &m,
                       &lusup[xlusup_first + i], &m);
	    } /* if (r > 1) */
	    else /* move to last row */
	    {
		zswap_(&n, &lusup[xlusup_first + m1], &m,
			&lusup[xlusup_first + i], &m);
		if (milu == SMILU_3)
		    for (j = 0; j < n; j++) {
			lusup[xlusup_first + m1 + j * m].r =
				z_abs1(&lusup[xlusup_first + m1 + j * m]);
			lusup[xlusup_first + m1 + j * m].i = 0.0;
                    }
	    }
	    lsub[xlsub_first + i] = lsub[xlsub_first + m1];
	    m1--;
	    continue;
	} /* if dropping */
	else
	{
	    if (temp[i] > d_max) d_max = temp[i];
	    if (temp[i] < d_min) d_min = temp[i];
	}
	i++;
    } /* for */

    /* Secondary dropping: drop more rows according to the quota. */
    quota = ceil((double)quota / (double)n);
    if (drop_rule & DROP_SECONDARY && m - r > quota)
    {
	register double tol = d_max;

	/* Calculate the second dropping tolerance */
	if (quota > n)
	{
	    if (drop_rule & DROP_INTERP) /* by interpolation */
	    {
		d_max = 1.0 / d_max; d_min = 1.0 / d_min;
		tol = 1.0 / (d_max + (d_min - d_max) * quota / (m - n - r));
	    }
	    else /* by quick sort */
	    {
		register int *itemp = iwork - n;
		A = temp;
		for (i = n; i <= m1; i++) itemp[i] = i;
		qsort(iwork, m1 - n + 1, sizeof(int), _compare_);
		tol = temp[iwork[quota]];
	    }
	}

	for (i = n; i <= m1; )
	{
	    if (temp[i] <= tol)
	    {
		register int j;
		r++;
		/* drop the current row and move the last undropped row here */
		if (r > 1) /* add to last row */
		{
		    /* accumulate the sum (for MILU) */
		    switch (milu)
		    {
			case SMILU_1:
			case SMILU_2:
			    zaxpy_(&n, &one, &lusup[xlusup_first + i], &m,
				    &lusup[xlusup_first + m - 1], &m);
			    break;
			case SMILU_3:
			    for (j = 0; j < n; j++)
				lusup[xlusup_first + (m - 1) + j * m].r +=
   				  z_abs1(&lusup[xlusup_first + i + j * m]);
			    break;
			case SILU:
			default:
			    break;
		    }
		    zcopy_(&n, &lusup[xlusup_first + m1], &m,
			    &lusup[xlusup_first + i], &m);
		} /* if (r > 1) */
		else /* move to last row */
		{
		    zswap_(&n, &lusup[xlusup_first + m1], &m,
			    &lusup[xlusup_first + i], &m);
		    if (milu == SMILU_3)
			for (j = 0; j < n; j++) {
			    lusup[xlusup_first + m1 + j * m].r =
				    z_abs1(&lusup[xlusup_first + m1 + j * m]);
			    lusup[xlusup_first + m1 + j * m].i = 0.0;
                        }
		}
		lsub[xlsub_first + i] = lsub[xlsub_first + m1];
		m1--;
		temp[i] = temp[m1];

		continue;
	    }
	    i++;

	} /* for */

    } /* if secondary dropping */

    for (i = n; i < m; i++) temp[i] = 0.0;

    if (r == 0)
    {
	*nnzLj += m * n;
	return 0;
    }

    /* add dropped entries to the diagnal */
    if (milu != SILU)
    {
	register int j;
	doublecomplex t;
	for (j = 0; j < n; j++)
	{
	    zd_mult(&t, &lusup[xlusup_first + (m - 1) + j * m],
                               MILU_ALPHA);
 	    switch (milu)
	    {
		case SMILU_1:
		    if ( !(z_eq(&t, &none)) ) {
                        z_add(&t, &t, &one);
                        zz_mult(&lusup[xlusup_first + j * inc_diag],
			                  &lusup[xlusup_first + j * inc_diag],
                                          &t);
                    }
		    else
		    {
                        zd_mult(
                                &lusup[xlusup_first + j * inc_diag],
			        &lusup[xlusup_first + j * inc_diag],
                                *fill_tol);
#ifdef DEBUG
			printf("[1] ZERO PIVOT: FILL col %d.\n", first + j);
			fflush(stdout);
#endif
			nzp++;
		    }
		    break;
		case SMILU_2:
                    zd_mult(&lusup[xlusup_first + j * inc_diag],
                                          &lusup[xlusup_first + j * inc_diag],
                                          1.0 + z_abs1(&t));
		    break;
		case SMILU_3:
                    z_add(&t, &t, &one);
                    zz_mult(&lusup[xlusup_first + j * inc_diag],
	                              &lusup[xlusup_first + j * inc_diag],
                                      &t);
		    break;
		case SILU:
		default:
		    break;
	    }
	}
	if (nzp > 0) *fill_tol = -nzp;
    }

    /* Remove dropped entries from the memory and fix the pointers. */
    m1 = m - r;
    for (j = 1; j < n; j++)
    {
	register int tmp1, tmp2;
	tmp1 = xlusup_first + j * m1;
	tmp2 = xlusup_first + j * m;
	for (i = 0; i < m1; i++)
	    lusup[i + tmp1] = lusup[i + tmp2];
    }
    for (i = 0; i < nzlc; i++)
	lusup[xlusup_first + i + n * m1] = lusup[xlusup_first + i + n * m];
    for (i = 0; i < nzlc; i++)
	lsub[xlsub[last + 1] - r + i] = lsub[xlsub[last + 1] + i];
    for (i = first + 1; i <= last + 1; i++)
    {
	xlusup[i] -= r * (i - first);
	xlsub[i] -= r;
    }
    if (lastc)
    {
	xlusup[last + 2] -= r * n;
	xlsub[last + 2] -= r;
    }

    *nnzLj += (m - r) * n;
    return r;
}
