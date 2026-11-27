/*
 * basiclu_get_factors.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "ipm/basiclu/lu_internal.h"

lu_int basiclu_get_factors
(
    lu_int istore[],
    double xstore[],
    lu_int Li[],
    double Lx[],
    lu_int Ui[],
    double Ux[],
    lu_int Wi[],
    double Wx[],
    lu_int rowperm[],
    lu_int colperm[],
    lu_int Lcolptr[],
    lu_int Lrowidx[],
    double Lvalue_[],
    lu_int Ucolptr[],
    lu_int Urowidx[],
    double Uvalue_[]
)
{
    struct lu this;
    lu_int m, status;

    status = lu_load(&this, istore, xstore, Li, Lx, Ui, Ux, Wi, Wx);
    if (status != BASICLU_OK)
        return status;
    if (this.nupdate != 0)
    {
        status = BASICLU_ERROR_invalid_call;
        return lu_save(&this, istore, xstore, status);
    }
    m = this.m;

    if (rowperm)
        memcpy(rowperm, this.pivotrow, m*sizeof(lu_int));
    if (colperm)
        memcpy(colperm, this.pivotcol, m*sizeof(lu_int));

    if (Lcolptr && Lrowidx && Lvalue_)
    {
        const lu_int *Lbegin_p  = this.Lbegin_p;
        const lu_int *Ltbegin_p = this.Ltbegin_p;
        const lu_int *Lindex    = this.Lindex;
        const double *Lvalue    = this.Lvalue;
        const lu_int *p         = this.p;
        lu_int *colptr          = this.iwork1; /* size m workspace */
        lu_int i, k, put, pos;

        /*
         * L[:,k] will hold the elimination factors from the k-th pivot step.
         * First set the column pointers and store the unit diagonal elements
         * at the front of each column. Then scatter each row of L' into the
         * columnwise L so that the row indices become sorted.
         */
        put = 0;
        for (k = 0; k < m; k++)
        {
            Lcolptr[k] = put;
            Lrowidx[put] = k;
            Lvalue_[put++] = 1.0;
            colptr[p[k]] = put; /* next free position in column */
            put += Lbegin_p[k+1] - Lbegin_p[k] - 1;
            /* subtract 1 because internal storage uses (-1) terminators */
        }
        Lcolptr[m] = put;
        assert(put == this.Lnz+m);

        for (k = 0; k < m; k++)
        {
            for (pos = Ltbegin_p[k]; (i = Lindex[pos]) >= 0; pos++)
            {
                put = colptr[i]++;
                Lrowidx[put] = k;
                Lvalue_[put] = Lvalue[pos];
            }
        }

        #ifndef NDEBUG
        for (k = 0; k < m; k++)
        {
            assert(colptr[p[k]] == Lcolptr[k+1]);
        }
        #endif
    }

    if (Ucolptr && Urowidx && Uvalue_)
    {
        const lu_int *Wbegin    = this.Wbegin;
        const lu_int *Wend      = this.Wend;
        const lu_int *Windex    = this.Windex;
        const double *Wvalue    = this.Wvalue;
        const double *col_pivot = this.col_pivot;
        const lu_int *pivotcol  = this.pivotcol;
        lu_int *colptr          = this.iwork1; /* size m workspace */
        lu_int j, k, put, pos;

        /*
         * U[:,k] will hold the column of B from the k-th pivot step.
         * First set the column pointers and store the pivot element at the end
         * of each column. Then scatter each row of U' into the columnwise U so
         * that the row indices become sorted.
         */
        memset(colptr, 0, m*sizeof(lu_int)); /* column counts */
        for (j = 0; j < m; j++)
        {
            for (pos = Wbegin[j]; pos < Wend[j]; pos++)
                colptr[Windex[pos]]++;
        }
        put = 0;
        for (k = 0; k < m; k++) /* set column pointers */
        {
            j = pivotcol[k];
            Ucolptr[k] = put;
            put += colptr[j];
            colptr[j] = Ucolptr[k]; /* next free position in column */
            Urowidx[put] = k;
            Uvalue_[put++] = col_pivot[j];
        }
        Ucolptr[m] = put;
        assert(put == this.Unz+m);
        for (k = 0; k < m; k++) /* scatter row k */
        {
            j = pivotcol[k];
            for (pos = Wbegin[j]; pos < Wend[j]; pos++)
            {
                put = colptr[Windex[pos]]++;
                Urowidx[put] = k;
                Uvalue_[put] = Wvalue[pos];
            }
        }

        #ifndef NDEBUG
        for (k = 0; k < m; k++) assert(colptr[pivotcol[k]] == Ucolptr[k+1]-1);
        for (k = 0; k < m; k++) assert(Urowidx[Ucolptr[k+1]-1] == k);
        #endif
    }

    return BASICLU_OK;
}
