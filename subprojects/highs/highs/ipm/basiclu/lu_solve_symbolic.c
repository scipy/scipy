/*
 * lu_solve_symbolic.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * Symbolic solve with triangular matrix
 *
 */

#include "ipm/basiclu/lu_internal.h"

/*
 * lu_solve_symbolic()
 *
 * The pattern of the right-hand side is given in irhs[0..nrhs-1]. The
 * pattern of the solution is returned in ilhs[top..m-1] in topological
 * order, top is the function return value.
 *
 * When end is not NULL, then the pattern of column j of the matrix must be
 * given in index[begin[j]..end[j]-1]. When end is NULL, then each column must
 * be terminated by a negative index.
 *
 * The method is due to J. Gilbert and T. Peierls, "Sparse partial pivoting
 * in time proportional to arithmetic operations", (1988).
 */

lu_int lu_solve_symbolic
(
    const lu_int m,
    const lu_int *begin,
    const lu_int *end,
    const lu_int *index,
    const lu_int nrhs,
    const lu_int *irhs,
    lu_int *ilhs,
    lu_int *pstack,             /* size m workspace */
    lu_int *marked,             /* marked[i] != M on entry */
    const lu_int M
)
{
    lu_int i, n, top = m;

    for (n = 0; n < nrhs; n++)
        if (marked[i = irhs[n]] != M)
            top = lu_dfs(i, begin, end, index, top, ilhs, pstack, marked, M);

    return top;
}
