/*
 * lu_solve_dense.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "ipm/basiclu/lu_internal.h"

void lu_solve_dense(struct lu *this, const double *rhs, double *lhs, char trans)
{
    const lu_int m                  = this->m;
    const lu_int nforrest           = this->nforrest;
    const lu_int *p                 = this->p;
    const lu_int *eta_row           = this->eta_row;
    const lu_int *pivotcol          = this->pivotcol;
    const lu_int *pivotrow          = this->pivotrow;
    const lu_int *Lbegin_p          = this->Lbegin_p;
    const lu_int *Ltbegin_p         = this->Ltbegin_p;
    const lu_int *Ubegin            = this->Ubegin;
    const lu_int *Rbegin            = this->Rbegin;
    const lu_int *Wbegin            = this->Wbegin;
    const lu_int *Wend              = this->Wend;
    const double *col_pivot         = this->col_pivot;
    const double *row_pivot         = this->row_pivot;
    const lu_int *Lindex            = this->Lindex;
    const double *Lvalue            = this->Lvalue;
    const lu_int *Uindex            = this->Uindex;
    const double *Uvalue            = this->Uvalue;
    const lu_int *Windex            = this->Windex;
    const double *Wvalue            = this->Wvalue;
    double *work1                   = this->work1;

    lu_int i, k, t, ipivot, jpivot, pos;
    double x;

    lu_garbage_perm(this);
    assert(this->pivotlen == m);

    if (trans == 't' || trans == 'T')
    {
        /* ----------------------- */
        /* Solve transposed system */
        /* ----------------------- */

        memcpy(work1, rhs, m*sizeof(double));

        /* Solve with U'. */
        for (k = 0; k < m; k++)
        {
            jpivot = pivotcol[k];
            ipivot = pivotrow[k];
            x = work1[jpivot] / col_pivot[jpivot];
            for (pos = Wbegin[jpivot]; pos < Wend[jpivot]; pos++)
            {
                work1[Windex[pos]] -= x * Wvalue[pos];
            }
            lhs[ipivot] = x;
        }

        /* Solve with update ETAs backwards. */
        for (t = nforrest-1; t >= 0; t--)
        {
            ipivot = eta_row[t];
            x = lhs[ipivot];
            for (pos = Rbegin[t]; pos < Rbegin[t+1]; pos++)
            {
                i = Lindex[pos];
                lhs[i] -= x * Lvalue[pos];
            }
        }

        /* Solve with L'. */
        for (k = m-1; k >= 0; k--)
        {
            x = 0.0;
            for (pos = Lbegin_p[k]; (i = Lindex[pos]) >= 0; pos++)
            {
                x += lhs[i] * Lvalue[pos];
            }
            lhs[p[k]] -= x;
        }
    }
    else
    {
        /* -------------------- */
        /* Solve forward system */
        /* -------------------- */

        memcpy(work1, rhs, m*sizeof(double));

        /* Solve with L. */
        for (k = 0; k < m; k++)
        {
            x = 0.0;
            for (pos = Ltbegin_p[k]; (i = Lindex[pos]) >= 0; pos++)
            {
                x += work1[i] * Lvalue[pos];
            }
            work1[p[k]] -= x;
        }

        /* Solve with update ETAs. */
        pos = Rbegin[0];
        for (t = 0; t < nforrest; t++)
        {
            ipivot = eta_row[t];
            x = 0.0;
            for ( ; pos < Rbegin[t+1]; pos++)
            {
                x += work1[Lindex[pos]] * Lvalue[pos];
            }
            work1[ipivot] -= x;
        }

        /* Solve with U. */
        for (k = m-1; k >= 0; k--)
        {
            jpivot = pivotcol[k];
            ipivot = pivotrow[k];
            x = work1[ipivot] / row_pivot[ipivot];
            for (pos = Ubegin[ipivot]; (i = Uindex[pos]) >= 0; pos++)
            {
                work1[i] -= x * Uvalue[pos];
            }
            lhs[jpivot] = x;
        }
    }
}
