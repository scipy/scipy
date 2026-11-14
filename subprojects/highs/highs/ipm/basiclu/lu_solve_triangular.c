/*
 * lu_solve_triangular.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "ipm/basiclu/lu_internal.h"


/**
 * lu_solve_triangular() - substitution with triangular matrix
 *
 * The symbolic nonzero pattern of the solution must be given in topological
 * order in pattern_symb[0..nz_symb-1]. On return pattern[0..nz-1] holds the
 * nonzero pattern of the solution after dropping numerical zeros; nz is
 * returned. pattern and pattern_symb can point to the same array.
 *
 * Entries in the solution that are less than or equal to droptol are set to
 * zero. When droptol is zero or negative, then no entries will be set to zero.
 *
 * Note: The nonzero pattern of the solution never includes zeros. That means,
 *       even if droptol is negative, the output pattern is not identical to
 *       the symbolic pattern when exact cancellation happens.
 *
 * The pivot elements must be stored separately to the matrix. When pivot is
 * NULL, then the pivot elements are assumed to be 1. The matrix is given in
 * parallel arrays index, value. When end is not NULL, column j has elements
 *
 *   index[begin[j]..end[j]-1], value[begin[j]..end[j]-1].
 *
 * When end is NULL, then each column must be terminated by a negative index.
 *
 */
lu_int lu_solve_triangular
(
    const lu_int nz_symb,
    const lu_int *pattern_symb,
    const lu_int *begin,
    const lu_int *end,
    const lu_int *index,
    const double *value,
    const double *pivot,
    const double droptol,
    double *lhs,                /* solution overwrites RHS */
    lu_int *pattern,
    lu_int *flops               /* add flop count */
)
{
    lu_int i, ipivot, pos, n, nz = 0, flop_count = 0;
    double x;

    if (pivot && end)
    {
        for (n = 0; n < nz_symb; n++)
        {
            ipivot = pattern_symb[n];
            if (lhs[ipivot])
            {
                x = lhs[ipivot] /= pivot[ipivot];
                flop_count++;
                for (pos = begin[ipivot]; pos < end[ipivot]; pos++)
                {
                    i = index[pos];
                    lhs[i] -= x * value[pos];
                    flop_count++;
                }
                if (fabs(x) > droptol)
                    pattern[nz++] = ipivot;
                else
                    lhs[ipivot] = 0.0;
            }
        }
    }
    else if (pivot)
    {
        for (n = 0; n < nz_symb; n++)
        {
            ipivot = pattern_symb[n];
            if (lhs[ipivot])
            {
                x = lhs[ipivot] /= pivot[ipivot];
                flop_count++;
                for (pos = begin[ipivot]; (i = index[pos]) >= 0; pos++)
                {
                    lhs[i] -= x * value[pos];
                    flop_count++;
                }
                if (fabs(x) > droptol)
                    pattern[nz++] = ipivot;
                else
                    lhs[ipivot] = 0.0;
            }
        }
    }
    else if (end)
    {
        for (n = 0; n < nz_symb; n++)
        {
            ipivot = pattern_symb[n];
            if (lhs[ipivot])
            {
                x = lhs[ipivot];
                for (pos = begin[ipivot]; pos < end[ipivot]; pos++)
                {
                    i = index[pos];
                    lhs[i] -= x * value[pos];
                    flop_count++;
                }
                if (fabs(x) > droptol)
                    pattern[nz++] = ipivot;
                else
                    lhs[ipivot] = 0.0;
            }
        }
    }
    else
    {
        for (n = 0; n < nz_symb; n++)
        {
            ipivot = pattern_symb[n];
            if (lhs[ipivot])
            {
                x = lhs[ipivot];
                for (pos = begin[ipivot]; (i = index[pos]) >= 0; pos++)
                {
                    lhs[i] -= x * value[pos];
                    flop_count++;
                }
                if (fabs(x) > droptol)
                    pattern[nz++] = ipivot;
                else
                    lhs[ipivot] = 0.0;
            }
        }
    }

    *flops += flop_count;
    return nz;
}
