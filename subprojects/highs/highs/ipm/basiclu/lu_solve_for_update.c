/*
 * lu_solve_for_update.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "ipm/basiclu/lu_internal.h"

lu_int lu_solve_for_update(
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
    lu_int *eta_row                 = this->eta_row;
    const lu_int *pivotcol          = this->pivotcol;
    const lu_int *pivotrow          = this->pivotrow;
    const lu_int *Lbegin            = this->Lbegin;
    const lu_int *Ltbegin           = this->Ltbegin;
    const lu_int *Ltbegin_p         = this->Ltbegin_p;
    lu_int *Ubegin                  = this->Ubegin;
    lu_int *Rbegin                  = this->Rbegin;
    const lu_int *Wbegin            = this->Wbegin;
    const lu_int *Wend              = this->Wend;
    const double *col_pivot         = this->col_pivot;
    const double *row_pivot         = this->row_pivot;
    lu_int *Lindex                  = this->Lindex;
    double *Lvalue                  = this->Lvalue;
    lu_int *Uindex                  = this->Uindex;
    double *Uvalue                  = this->Uvalue;
    const lu_int *Windex            = this->Windex;
    const double *Wvalue            = this->Wvalue;
    lu_int *marked                  = this->marked;

    lu_int i, j, k, n, t, top, pos, put, ipivot, jpivot, nz, nz_symb, M,
        room, need, jbegin, jend;
    double x, xdrop, pivot;

    const lu_int want_solution = p_nlhs && ilhs && xlhs;
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

        jpivot = irhs[0];
        ipivot = pmap[jpivot];
        jbegin = Wbegin[jpivot];
        jend   = Wend[jpivot];

        /*
         * Compute row eta vector.
         * Symbolic pattern in pattern_symb[top..m-1], indices of (actual)
         * nonzeros in pattern[0..nz-1], values scattered into work.
         * We do not drop small elements to zero, but the symbolic and the
         * numeric pattern will still be different when we have exact
         * cancellation.
         */
        M = ++this->marker;
        top = lu_solve_symbolic(m, Wbegin, Wend, Windex, jend-jbegin,
                                Windex+jbegin, pattern_symb, pstack, marked, M);
        nz_symb = m-top;

        /* reallocate if not enough memory in Li, Lx (where we store R) */
        room = this->Lmem - Rbegin[nforrest];
        if (room < nz_symb)
        {
            this->addmemL = nz_symb - room;
            return BASICLU_REALLOCATE;
        }

        for (pos = jbegin; pos < jend; pos++)
        {
            work[Windex[pos]] = Wvalue[pos];
        }
        lu_solve_triangular(nz_symb, pattern_symb+top, Wbegin, Wend, Windex,
                            Wvalue, col_pivot, 0.0, work, pattern, &Uflops);

        /*
         * Compress row eta into L, pattern mapped from column to row indices.
         * The triangularity test in lu_update requires the symbolic pattern.
         */
        put = Rbegin[nforrest];
        for (t = top; t < m; t++)
        {
            j = pattern_symb[t];
            i = pmap[j];
            Lindex[put] = i;
            Lvalue[put++] = work[j];
            work[j] = 0;
        }
        Rbegin[nforrest+1] = put;
        eta_row[nforrest] = ipivot;
        this->btran_for_update = jpivot;

        if (!want_solution)
            goto done;

        /*
         * Scatter the row eta into xlhs and scale it to become the solution
         * to U^{-1}*[unit vector]. Now we can drop small entries to zero and
         * recompute the numerical pattern.
         */
        M = ++this->marker;
        pattern[0] = ipivot;
        marked[ipivot] = M;
        pivot = col_pivot[jpivot];
        xlhs[ipivot] = 1.0 / pivot;

        xdrop = droptol * fabs(pivot);
        nz = 1;
        for (pos = Rbegin[nforrest]; pos < Rbegin[nforrest+1]; pos++)
        {
            if (fabs(Lvalue[pos]) > xdrop)
            {
                pattern[nz++] = i = Lindex[pos];
                marked[i] = M;
                xlhs[i] = -Lvalue[pos] / pivot;
            }
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

        /* reallocate if not enough memory in U */
        room = this->Umem - Ubegin[m];
        need = nz+1;
        if (room < need)
        {
            for (n = 0; n < nz; n++)
                work[pattern[n]] = 0.0;
            this->addmemU = need - room;
            return BASICLU_REALLOCATE;
        }

        /* Compress spike into U. */
        put = Ubegin[m];
        for (n = 0; n < nz; n++)
        {
            i = pattern[n];
            Uindex[put] = i;
            Uvalue[put++] = work[i];
            if (!want_solution)
                work[i] = 0.0;
        }
        Uindex[put++] = -1;     /* terminate column */
        this->ftran_for_update = 0;

        if (!want_solution)
            goto done;

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

done:
    this->Lflops += Lflops;
    this->Uflops += Uflops;
    this->Rflops += Rflops;
    this->update_cost_numer += Rflops;
    return BASICLU_OK;
}
