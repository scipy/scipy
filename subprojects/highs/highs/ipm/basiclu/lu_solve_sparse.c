/*
 * lu_solve_sparse.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "ipm/basiclu/lu_internal.h"

void lu_solve_sparse(
    struct lu *this, const lu_int nrhs, const lu_int *irhs, const double *xrhs,
    lu_int *p_nlhs, lu_int *ilhs, double *xlhs, char trans)
{
    const lu_int m                  = this->m;
    const lu_int nforrest           = this->nforrest;
    const lu_int pivotlen           = this->pivotlen;
    const lu_int nz_sparse          = this->sparse_thres * m;
    const double droptol            = this->droptol;
    const lu_int *p                 = this->p;
    const lu_int *pmap              = this->pmap;
    const lu_int *qmap              = this->qmap;
    const lu_int *eta_row           = this->eta_row;
    const lu_int *pivotcol          = this->pivotcol;
    const lu_int *pivotrow          = this->pivotrow;
    const lu_int *Lbegin            = this->Lbegin;
    const lu_int *Ltbegin           = this->Ltbegin;
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
    lu_int *marked                  = this->marked;

    lu_int i, j, k, n, t, top, pos, ipivot, jpivot, nz, nz_symb, M;
    double x;

    lu_int Lflops = 0, Uflops = 0, Rflops = 0;

    if (trans == 't' || trans == 'T')
    {
        /* ----------------------- */
        /* Solve transposed system */
        /* ----------------------- */

        lu_int *pattern_symb = this->iwork1;
        lu_int *pattern      = this->iwork1 + m;
        double *work         = this->work0;
        lu_int *pstack       = (void *) this->work1;
        assert(sizeof(lu_int) <= sizeof(double));

        /*
         * Sparse triangular solve with U'.
         * Solution scattered into work, indices in pattern[0..nz-1].
         */
        M = ++this->marker;
        top = lu_solve_symbolic(m, Wbegin, Wend, Windex, nrhs, irhs,
                                pattern_symb, pstack, marked, M);
        nz_symb = m-top;

        for (n = 0; n < nrhs; n++)
            work[irhs[n]] = xrhs[n];
        nz = lu_solve_triangular(nz_symb, pattern_symb+top, Wbegin, Wend,
                                 Windex, Wvalue, col_pivot, droptol, work,
                                 pattern, &Uflops);

        /*
         * Permute solution into xlhs.
         * Map pattern from column indices to row indices.
         */
        M = ++this->marker;
        for (n = 0; n < nz; n++)
        {
            j = pattern[n];
            i = pmap[j];
            pattern[n] = i;
            xlhs[i] = work[j];
            work[j] = 0;
            marked[i] = M;
        }

        /*
         * Solve with update etas.
         * Append fill-in to pattern.
         */
        for (t = nforrest-1; t >= 0; t--)
        {
            ipivot = eta_row[t];
            if (xlhs[ipivot])
            {
                x = xlhs[ipivot];
                for (pos = Rbegin[t]; pos < Rbegin[t+1]; pos++)
                {
                    i = Lindex[pos];
                    if (marked[i] != M)
                    {
                        marked[i] = M;
                        pattern[nz++] = i;
                    }
                    xlhs[i] -= x * Lvalue[pos];
                    Rflops++;
                }
            }
        }

        if (nz <= nz_sparse)
        {
            /*
             * Sparse triangular solve with L'.
             * Solution scattered into xlhs, indices in ilhs[0..nz-1].
             */
            M = ++this->marker;
            top = lu_solve_symbolic(m, Ltbegin, NULL, Lindex, nz, pattern,
                                    pattern_symb, pstack, marked, M);
            nz_symb = m-top;

            nz = lu_solve_triangular(nz_symb, pattern_symb+top, Ltbegin, NULL,
                                      Lindex, Lvalue, NULL, droptol, xlhs, ilhs,
                                      &Lflops);
            *p_nlhs = nz;
        }
        else
        {
            /*
             * Sequential triangular solve with L'.
             * Solution scattered into xlhs, indices in ilhs[0..nz-1].
             */
            nz = 0;
            for (k = m-1; k >= 0; k--)
            {
                ipivot = p[k];
                if (xlhs[ipivot])
                {
                    x = xlhs[ipivot];
                    for (pos = Ltbegin_p[k]; (i = Lindex[pos]) >= 0; pos++)
                    {
                        xlhs[i] -= x * Lvalue[pos];
                        Lflops++;
                    }
                    if (fabs(x) > droptol)
                        ilhs[nz++] = ipivot;
                    else
                        xlhs[ipivot] = 0.0;
                }
            }
            *p_nlhs = nz;
        }
    }
    else
    {
        /* -------------------- */
        /* Solve forward system */
        /* -------------------- */

        lu_int *pattern_symb = this->iwork1;
        lu_int *pattern      = this->iwork1 + m;
        double *work         = this->work0;
        lu_int *pstack       = (void *) this->work1;
        assert(sizeof(lu_int) <= sizeof(double));

        /*
         * Sparse triangular solve with L.
         * Solution scattered into work, indices in pattern[0..nz-1].
         */
        M = ++this->marker;
        top = lu_solve_symbolic(m, Lbegin, NULL, Lindex, nrhs, irhs,
                                pattern_symb, pstack, marked, M);
        nz_symb = m-top;

        for (n = 0; n < nrhs; n++)
            work[irhs[n]] = xrhs[n];
        nz = lu_solve_triangular(nz_symb, pattern_symb+top, Lbegin, NULL,
                                 Lindex, Lvalue, NULL, droptol, work, pattern,
                                 &Lflops);

        /* unmark cancellation */
        if (nz < nz_symb)
        {
            for (t = top, n = 0; n < nz; t++)
            {
                i = pattern_symb[t];
                if (i == pattern[n])
                    n++;
                else
                    marked[i]--;
            }
            for ( ; t < m; t++)
                marked[pattern_symb[t]]--;
        }

        /*
         * Solve with update etas.
         * Append fill-in to pattern.
         */
        pos = Rbegin[0];
        for (t = 0; t < nforrest; t++)
        {
            ipivot = eta_row[t];
            x = 0.0;
            for ( ; pos < Rbegin[t+1]; pos++)
            {
                x += work[Lindex[pos]] * Lvalue[pos];
            }
            work[ipivot] -= x;
            if (x && marked[ipivot] != M)
            {
                marked[ipivot] = M;
                pattern[nz++] = ipivot;
            }
        }
        Rflops += Rbegin[nforrest] - Rbegin[0];

        if (nz <= nz_sparse)
        {
            /*
             * Sparse triangular solve with U.
             * Solution scattered into work, indices in ilhs[0..nz-1].
             */
            M = ++this->marker;
            top = lu_solve_symbolic(m, Ubegin, NULL, Uindex, nz, pattern,
                                    pattern_symb, pstack, marked, M);
            nz_symb = m-top;

            nz = lu_solve_triangular(nz_symb, pattern_symb+top, Ubegin, NULL,
                                     Uindex, Uvalue, row_pivot, droptol, work,
                                     ilhs, &Uflops);

            /*
             * Permute solution into xlhs.
             * Map pattern from row indices to column indices.
             */
            for (n = 0; n < nz; n++)
            {
                i = ilhs[n];
                j = qmap[i];
                ilhs[n] = j;
                xlhs[j] = work[i];
                work[i] = 0;
            }
        }
        else
        {
            /*
             * Sequential triangular solve with U.
             * Solution computed in work and permuted into xlhs.
             * Pattern (in column indices) stored in ilhs[0..nz-1].
             */
            nz = 0;
            for (k = pivotlen-1; k >= 0; k--)
            {
                ipivot = pivotrow[k];
                jpivot = pivotcol[k];
                if (work[ipivot])
                {
                    x = work[ipivot] / row_pivot[ipivot];
                    work[ipivot] = 0.0;
                    for (pos = Ubegin[ipivot]; (i = Uindex[pos]) >= 0; pos++)
                    {
                        work[i] -= x * Uvalue[pos];
                        Uflops++;
                    }
                    if (fabs(x) > droptol)
                    {
                        ilhs[nz++] = jpivot;
                        xlhs[jpivot] = x;
                    }
                }
            }
        }
        *p_nlhs = nz;
    }

    this->Lflops += Lflops;
    this->Uflops += Uflops;
    this->Rflops += Rflops;
    this->update_cost_numer += Rflops;
}
