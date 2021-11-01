/*
 * lu_matrix_norm.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * Computes the 1-norm and infinity-norm of the matrix that was freshly
 * factorized. Unit cols inserted by the factorization are handled implicitly.
 *
 */

#include "lu_internal.h"

void lu_matrix_norm(
    struct lu *this, const lu_int *Bbegin, const lu_int *Bend, const lu_int *Bi,
    const double *Bx)
{
    const lu_int m          = this->m;
    const lu_int rank       = this->rank;
    const lu_int *pivotcol  = this->pivotcol;
    const lu_int *pivotrow  = this->pivotrow;
    double *rowsum          = this->work1;
    lu_int ipivot, jpivot, i, k, pos;
    double onenorm, infnorm, colsum;

    assert(this->nupdate == 0);

    for (i = 0; i < m; i++) rowsum[i] = 0;
    onenorm = 0;
    infnorm = 0;
    for (k = 0; k < rank; k++)
    {
        jpivot = pivotcol[k];
        colsum = 0;
        for (pos = Bbegin[jpivot]; pos < Bend[jpivot]; pos++)
        {
            colsum += fabs(Bx[pos]);
            rowsum[Bi[pos]] += fabs(Bx[pos]);
        }
        onenorm = fmax(onenorm, colsum);
    }
    for (k = rank; k < m; k++)
    {
        ipivot = pivotrow[k];
        rowsum[ipivot] += 1;
        onenorm = fmax(onenorm, 1);
    }
    for (i = 0; i < m; i++) infnorm = fmax(infnorm, rowsum[i]);

    this->onenorm = onenorm;
    this->infnorm = infnorm;
}
