/*
 * lu_garbage_perm.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "ipm/basiclu/lu_internal.h"

/*
 * lu_garbage_perm()
 *
 * The sequence of pivot columns and pivot rows is stored in
 *
 *  pivotcol[0..pivotlen-1], pivotrow[0..pivotlen-1],
 *
 * where pivotlen >= m. When pivotlen > m, then the arrays contain duplicates.
 * For each index its last occurrence in the arrays is its position in the pivot
 * sequence and occurrences before mark unused slots.
 *
 * This routine removes duplicates and compresses the indices such that
 * pivotlen == m.
 */

void lu_garbage_perm(struct lu *this)
{
    const lu_int m      = this->m;
    lu_int pivotlen     = this->pivotlen;
    lu_int *pivotcol    = this->pivotcol;
    lu_int *pivotrow    = this->pivotrow;
    lu_int *marked      = this->marked;

    lu_int j, get, put, marker;

    if (pivotlen > m)
    {
        marker = ++this->marker;
        put = pivotlen;
        for (get = pivotlen-1; get >= 0; get--)
        {
            if (marked[j = pivotcol[get]] != marker)
            {
                marked[j] = marker;
                pivotcol[--put] = j;
                pivotrow[put] = pivotrow[get];
            }
        }
        assert(put+m == pivotlen);
        memmove(pivotcol, pivotcol+put, m*sizeof(lu_int));
        memmove(pivotrow, pivotrow+put, m*sizeof(lu_int));
        this->pivotlen = m;
    }
}
