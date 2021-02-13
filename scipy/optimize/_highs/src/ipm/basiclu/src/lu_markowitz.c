/*
 * lu_markowitz.c
 *
 * Copyright (C) 2016-2019  ERGO-Code
 *
 * Search for pivot element with small Markowitz cost. An eligible pivot
 * must be nonzero and satisfy
 *
 *  (1) abs(piv) >= abstol,
 *  (2) abs(piv) >= reltol * max[pivot column].
 *
 *  From all eligible pivots search for one that minimizes
 *
 *   mc := (nnz[pivot row] - 1) * (nnz[pivot column] - 1).
 *
 * The search is terminated when maxsearch rows or columns with eligible pivots
 * have been searched (if not before). The row and column of the cheapest one
 * found is stored in pivot_row and pivot_col.
 *
 * When the active submatrix contains columns with column count = 0, then such a
 * column is chosen immediately and pivot_row = -1 is returned. Otherwise, when
 * the Markowitz search does not find a pivot that is nonzero and >= abstol,
 * then pivot_col = pivot_row = -1 is returned. (The latter cannot happen in the
 * current version of the code because lu_pivot() erases columns of the active
 * submatrix whose maximum absolute value drops below abstol.)
 *
 * The Markowitz search is implemented as described in [1].
 *
 * [1] U. Suhl, L. Suhl, "Computing Sparse LU Factorizations for Large-Scale
 *     Linear Programming Bases", ORSA Journal on Computing (1990)
 *
 */

#include "lu_internal.h"
#include "lu_list.h"

lu_int lu_markowitz(struct lu *this)
{
    const lu_int m                  = this->m;
    const lu_int *Wbegin            = this->Wbegin;
    const lu_int *Wend              = this->Wend;
    const lu_int *Windex            = this->Windex;
    const double *Wvalue            = this->Wvalue;
    const lu_int *colcount_flink    = this->colcount_flink;
    lu_int *rowcount_flink          = this->rowcount_flink;
    lu_int *rowcount_blink          = this->rowcount_blink;
    const double *colmax            = this->col_pivot;
    const double abstol             = this->abstol;
    const double reltol             = this->reltol;
    const lu_int maxsearch          = this->maxsearch;
    const lu_int search_rows        = this->search_rows;
    const lu_int nz_start           = search_rows ?
        MIN(this->min_colnz, this->min_rownz) : this->min_colnz;

    lu_int i, j, pos, where, inext, nz, pivot_row, pivot_col;
    lu_int nsearch, cheap, found, min_colnz, min_rownz;
    double cmx, x, tol, tic[2];

    /* integers for Markowitz cost must be 64 bit to prevent overflow */
    const int_least64_t M = m;
    int_least64_t nz1, nz2, mc, MC;

    pivot_row = -1;             /* row of best pivot so far */
    pivot_col = -1;             /* col of best pivot so far */
    MC = M*M;                   /* Markowitz cost of best pivot so far */
    nsearch = 0;                /* count rows/columns searched */
    min_colnz = -1;             /* minimum col count in active submatrix */
    min_rownz = -1;             /* minimum row count in active submatrix */
    assert(nz_start >= 1);

    /* If the active submatrix contains empty columns, choose one and return
       with pivot_row = -1. */
    if (colcount_flink[m] != m)
    {
        pivot_col = colcount_flink[m];
        assert(pivot_col >= 0 && pivot_col < m);
        assert(Wend[pivot_col] == Wbegin[pivot_col]);
        goto done;
    }

    for (nz = nz_start; nz <= m; nz++)
    {
        /* Search columns with nz nonzeros. */
        for (j = colcount_flink[m+nz]; j < m; j = colcount_flink[j])
        {
            if (min_colnz == -1)
                min_colnz = nz;
            assert(Wend[j] - Wbegin[j] == nz);
            cmx = colmax[j];
            assert(cmx >= 0);
            if (!cmx || cmx < abstol)
                continue;
            tol = fmax(abstol, reltol*cmx);
            for (pos = Wbegin[j]; pos < Wend[j]; pos++)
            {
                x = fabs(Wvalue[pos]);
                if (!x || x < tol)
                    continue;
                i = Windex[pos];
                assert(i >= 0 && i < m);
                nz1 = nz;
                nz2 = Wend[m+i] - Wbegin[m+i];
                assert(nz2 >= 1);
                mc = (nz1-1) * (nz2-1);
                if (mc < MC)
                {
                    MC = mc;
                    pivot_row = i;
                    pivot_col = j;
                    if (search_rows && MC <= (nz1-1)*(nz1-1))
                        goto done;
                }
            }
            /* We have seen at least one eligible pivot in column j. */
            assert(MC < M*M);
            if (++nsearch >= maxsearch)
                goto done;
        }
        assert(j == m+nz);

        if (!search_rows)
            continue;

        /* Search rows with nz nonzeros. */
        for (i = rowcount_flink[m+nz]; i < m; i = inext)
        {
            if (min_rownz == -1)
                min_rownz = nz;
            /* rowcount_flink[i] might be changed below, so keep a copy */
            inext = rowcount_flink[i];
            assert(Wend[m+i] - Wbegin[m+i] == nz);
            cheap = 0;          /* row has entries with Markowitz cost < MC? */
            found = 0;          /* eligible pivot found? */
            for (pos = Wbegin[m+i]; pos < Wend[m+i]; pos++)
            {
                j = Windex[pos];
                assert(j >= 0 && j < m);
                nz1 = nz;
                nz2 = Wend[j] - Wbegin[j];
                assert(nz2 >= 1);
                mc = (nz1-1) * (nz2-1);
                if (mc >= MC)
                    continue;
                cheap = 1;
                cmx = colmax[j];
                assert(cmx >= 0);
                if (!cmx || cmx < abstol)
                    continue;
                /* find position of pivot in column file */
                for (where = Wbegin[j]; Windex[where] != i; where++)
                    assert(where < Wend[j] - 1);
                x = fabs(Wvalue[where]);
                if (x >= abstol && x >= reltol*cmx)
                {
                    found = 1;
                    MC = mc;
                    pivot_row = i;
                    pivot_col = j;
                    if (MC <= nz1*(nz1-1))
                        goto done;
                }
            }
            /* If row i has cheap entries but none of them is numerically
             * acceptable, then don't search the row again until updated. */
            if (cheap && !found)
            {
                lu_list_move(i, m+1, rowcount_flink, rowcount_blink, m, NULL);
            }
            else
            {
                assert(MC < M*M);
                if (++nsearch >= maxsearch)
                    goto done;
            }
        }
        assert(i == m+nz);
    }

done:
    this->pivot_row = pivot_row;
    this->pivot_col = pivot_col;
    this->nsearch_pivot += nsearch;
    if (min_colnz >= 0)
        this->min_colnz = min_colnz;
    if (min_rownz >= 0)
        this->min_rownz = min_rownz;
    return BASICLU_OK;
}
