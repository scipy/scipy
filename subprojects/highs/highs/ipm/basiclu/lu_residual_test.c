/*
 * lu_residual_test.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * Stability test of fresh LU factorization based on relative residual.
 *
 */

#include "ipm/basiclu/lu_internal.h"

static double lu_onenorm(lu_int m, const double *x)
{
    lu_int i;
    double d = 0;
    for (i = 0; i < m; i++)
    {
        d += fabs(x[i]);
    }
    return d;
}

void lu_residual_test(
    struct lu *this, const lu_int *Bbegin, const lu_int *Bend, const lu_int *Bi,
    const double *Bx)
{
    const lu_int m                  = this->m;
    const lu_int rank               = this->rank;
    const lu_int *p                 = this->p;
    const lu_int *pivotcol          = this->pivotcol;
    const lu_int *pivotrow          = this->pivotrow;
    const lu_int *Lbegin_p          = this->Lbegin_p;
    const lu_int *Ltbegin_p         = this->Ltbegin_p;
    const lu_int *Ubegin            = this->Ubegin;
    const double *row_pivot         = this->row_pivot;
    const lu_int *Lindex            = this->Lindex;
    const double *Lvalue            = this->Lvalue;
    const lu_int *Uindex            = this->Uindex;
    const double *Uvalue            = this->Uvalue;
    double *rhs                     = this->work0;
    double *lhs                     = this->work1;

    lu_int i, k, ipivot, jpivot, pos;
    double norm_ftran, norm_ftran_res, norm_btran, norm_btran_res, d;

    assert(this->nupdate == 0);

    /* --------------------------------- */
    /* Residual Test with Forward System */
    /* --------------------------------- */

    /* Compute lhs = L\rhs and build rhs on-the-fly. */
    for (k = 0; k < m; k++)
    {
        d = 0.0;
        for (pos = Ltbegin_p[k]; (i = Lindex[pos]) >= 0; pos++)
        {
            d += lhs[i] * Lvalue[pos];
        }
        ipivot = p[k];
        rhs[ipivot] = d <= 0.0 ? 1.0 : -1.0;
        lhs[ipivot] = rhs[ipivot] - d;
    }

    /* Overwrite lhs by U\lhs. */
    for (k = m-1; k >= 0; k--)
    {
        ipivot = pivotrow[k];
        d = lhs[ipivot] /= row_pivot[ipivot];
        for (pos = Ubegin[ipivot]; (i = Uindex[pos]) >= 0; pos++)
        {
            lhs[i] -= d * Uvalue[pos];
        }
    }

    /* Overwrite rhs by the residual rhs-B*lhs. */
    for (k = 0; k < rank; k++)
    {
        ipivot = pivotrow[k];
        jpivot = pivotcol[k];
        d = lhs[ipivot];
        for (pos = Bbegin[jpivot]; pos < Bend[jpivot]; pos++)
        {
            rhs[Bi[pos]] -= d * Bx[pos];
        }
    }
    for (k = rank; k < m; k++)
    {
        ipivot = pivotrow[k];
        rhs[ipivot] -= lhs[ipivot];
    }
    norm_ftran = lu_onenorm(m, lhs);
    norm_ftran_res = lu_onenorm(m, rhs);

    /* ---------------------------------- */
    /* Residual Test with Backward System */
    /* ---------------------------------- */

    /* Compute lhs = U'\rhs and build rhs on-the-fly. */
    for (k = 0; k < m; k++)
    {
        ipivot = pivotrow[k];
        d = 0.0;
        for (pos = Ubegin[ipivot]; (i = Uindex[pos]) >= 0; pos++)
        {
            d += lhs[i] * Uvalue[pos];
        }
        rhs[ipivot] = d <= 0.0 ? 1.0 : -1.0;
        lhs[ipivot] = (rhs[ipivot] - d) / row_pivot[ipivot];
    }

    /* Overwrite lhs by L'\lhs. */
    for (k = m-1; k >= 0; k--)
    {
        d = 0.0;
        for (pos = Lbegin_p[k]; (i = Lindex[pos]) >= 0; pos++)
        {
            d += lhs[i] * Lvalue[pos];
        }
        lhs[p[k]] -= d;
    }

    /* Overwrite rhs by the residual rhs-B'*lhs. */
    for (k = 0; k < rank; k++)
    {
        ipivot = pivotrow[k];
        jpivot = pivotcol[k];
        d = 0.0;
        for (pos = Bbegin[jpivot]; pos < Bend[jpivot]; pos++)
        {
            d += lhs[Bi[pos]] * Bx[pos];
        }
        rhs[ipivot] -= d;
    }
    for (k = rank; k < m; k++)
    {
        ipivot = pivotrow[k];
        rhs[ipivot] -= lhs[ipivot];
    }
    norm_btran = lu_onenorm(m, lhs);
    norm_btran_res = lu_onenorm(m, rhs);

    /* -------- */
    /* Finalize */
    /* -------- */

    lu_matrix_norm(this, Bbegin, Bend, Bi, Bx);
    assert(this->onenorm > 0.0);
    assert(this->infnorm > 0.0);
    this->residual_test = fmax(norm_ftran_res / (m + this->onenorm*norm_ftran),
                               norm_btran_res / (m + this->infnorm*norm_btran));

    /* reset workspace */
    for (i = 0; i < m; i++) rhs[i] = 0;
}
