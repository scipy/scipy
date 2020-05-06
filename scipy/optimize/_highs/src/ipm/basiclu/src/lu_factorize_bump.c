/* 
 * lu_factorize_bump.c
 *
 * Copyright (C) 2016-2019  ERGO-Code
 *
 * Bump factorization driver routine
 *
 */

#include "lu_internal.h"
#include "lu_list.h"

lu_int lu_factorize_bump(struct lu *this)
{
    const lu_int m          = this->m;
    lu_int *colcount_flink  = this->colcount_flink;
    lu_int *colcount_blink  = this->colcount_blink;
    lu_int *pinv            = this->pinv;
    lu_int *qinv            = this->qinv;
    lu_int status = BASICLU_OK;

    while (this->rank + this->rankdef < m)
    {
        /*
         * Find pivot element. Markowitz search need not be called if the
         * previous call to lu_pivot() returned for reallocation. In this case
         * this->pivot_col is valid.
         */
        if (this->pivot_col < 0)
            lu_markowitz(this);
        assert(this->pivot_col >= 0);

        if (this->pivot_row < 0)
        {
            /* Eliminate empty column without choosing a pivot. */
            lu_list_remove(colcount_flink, colcount_blink, this->pivot_col);
            this->pivot_col = -1;
            this->rankdef++;
        }
        else
        {
            /* Eliminate pivot. This may require reallocation. */
            assert(pinv[this->pivot_row] == -1);
            assert(qinv[this->pivot_col] == -1);
            status = lu_pivot(this);
            if (status != BASICLU_OK)
                break;
            pinv[this->pivot_row] = this->rank;
            qinv[this->pivot_col] = this->rank;
            this->pivot_col = -1;
            this->pivot_row = -1;
            this->rank++;
        }
    }
    return status;
}
