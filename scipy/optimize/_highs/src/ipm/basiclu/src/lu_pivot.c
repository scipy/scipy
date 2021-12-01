/* 
 * lu_pivot.c
 *
 * Copyright (C) 2016-2019  ERGO-Code
 *
 * Pivot elimination from active submatrix
 *
 * lu_pivot() is the only routine callable from extern. It branches out the
 * implementation of the pivot operation. The pivot operation removes row
 * this->pivot_row and column this->pivot_col from the active submatrix and
 * applies a rank-1 update to the remaining active submatrix. It updates the
 * row and column counts and the column maxima.
 *
 * Each pivot elimination adds one column to L and one row to U. On entry
 * Lbegin_p[rank] and Ubegin[rank] must point to the next position in Lindex,
 * Lvalue respectively Uindex, Uvalue. On return the column in L is terminated
 * by index -1 and Lbegin_p[rank+1] and Ubegin[rank+1] point to the next free
 * position.
 *
 * The routines check if memory is sufficient and request reallocation before
 * manipulating data structures otherwise.
 *
 * Updating the columns of the active submatrix is implemented like in the
 * Coin-OR code Clp (J. Forrest) by compressing unmodified entries and then
 * appending the entries updated or filled-in by the pivot column. Compared
 * to the technique from the Suhl/Suhl paper, this method often produces
 * sparser factors. I assume the reason is tie breaking in the Markowitz
 * search and that in Forrest's method updated elements are moved to the end
 * of the column (and likewise for rows).
 */

#include "lu_internal.h"
#include "lu_list.h"
#include "lu_file.h"

/*
 * MAXROW_SMALL is the maximum number of off-diagonal elements in the pivot
 * column handled by lu_pivot_small(). lu_pivot_small() uses int64_t integers
 * for bit masking. Since each row to be updated requires one bit, the routine
 * can handle pivot operations for up to 64 rows (excluding pivot row).
 *
 * Since int64_t is optional in the C99 standard, using it limits portability
 * of the code. However, using a fixed threshold to switch between
 * lu_pivot_small() and lu_pivot_any() guarantees identical pivot operations
 * on all architectures. If int64_t does not exist, then the user can adapt it
 * by hand and is aware of it.
 */

#define MAXROW_SMALL 64

static lu_int lu_pivot_any(struct lu *this);
static lu_int lu_pivot_small(struct lu *this);
static lu_int lu_pivot_singleton_row(struct lu *this);
static lu_int lu_pivot_singleton_col(struct lu *this);
static lu_int lu_pivot_doubleton_col(struct lu *this);
static void lu_remove_col(struct lu *this, lu_int j);


/* ==========================================================================
   lu_pivot
   ========================================================================== */

lu_int lu_pivot(struct lu *this)
{
    const lu_int m          = this->m;
    const lu_int rank       = this->rank;
    const lu_int Lmem       = this->Lmem;
    const lu_int Umem       = this->Umem;
    const lu_int pivot_col  = this->pivot_col;
    const lu_int pivot_row  = this->pivot_row;
    const double *colmax    = this->col_pivot;
    const lu_int *Lbegin_p  = this->Lbegin_p;
    const lu_int *Ubegin    = this->Ubegin;
    const lu_int *Wbegin    = this->Wbegin;
    const lu_int *Wend      = this->Wend;
    const lu_int *Uindex    = this->Uindex;
    const lu_int nz_col     = Wend[pivot_col] - Wbegin[pivot_col];
    const lu_int nz_row     = Wend[m+pivot_row] - Wbegin[m+pivot_row];

    lu_int room, need, pos, j;
    lu_int status = BASICLU_OK;
    double tic[2];

    assert(nz_row >= 1);
    assert(nz_col >= 1);

    /* Check if room is available in L and U. */
    room = Lmem - Lbegin_p[rank];
    need = nz_col;          /* # off-diagonals in pivot col + end marker (-1) */
    if (room < need)
    {
        this->addmemL = need-room;
        status = BASICLU_REALLOCATE;
    }
    room = Umem - Ubegin[rank];
    need = nz_row-1;            /* # off-diagonals in pivot row */
    if (room < need)
    {
        this->addmemU = need-room;
        status = BASICLU_REALLOCATE;
    }
    if (status != BASICLU_OK)
        return status;

    /* Branch out implementation of pivot operation. */
    if (nz_row == 1)
    {
        status = lu_pivot_singleton_row(this);
    }
    else if (nz_col == 1)
    {
        status = lu_pivot_singleton_col(this);
    }
    else if (nz_col == 2)
    {
        status = lu_pivot_doubleton_col(this);
    }
    else if (nz_col-1 <= MAXROW_SMALL)
    {
        status = lu_pivot_small(this);
    }
    else
    {
        status = lu_pivot_any(this);
    }

    /* Remove all entries in columns whose maximum entry has dropped below
       absolute pivot tolerance. */
    if (status == BASICLU_OK)
    {
        for (pos = Ubegin[rank]; pos < Ubegin[rank+1]; pos++)
        {
            j = Uindex[pos];
            assert(j != pivot_col);
            if (colmax[j] == 0.0 || colmax[j] < this->abstol)
                lu_remove_col(this, j);
        }
    }

    this->factor_flops += (nz_col-1) * (nz_row-1);
    return status;
}


/* ==========================================================================
   lu_pivot_any
   ========================================================================== */

static lu_int lu_pivot_any(struct lu *this)
{
    const lu_int m          = this->m;
    const lu_int rank       = this->rank;
    const double droptol    = this->droptol;
    const lu_int pad        = this->pad;
    const double stretch    = this->stretch;
    const lu_int pivot_col  = this->pivot_col;
    const lu_int pivot_row  = this->pivot_row;
    lu_int *colcount_flink  = this->colcount_flink;
    lu_int *colcount_blink  = this->colcount_blink;
    lu_int *rowcount_flink  = this->rowcount_flink;
    lu_int *rowcount_blink  = this->rowcount_blink;
    double *colmax          = this->col_pivot;
    lu_int *Lbegin_p        = this->Lbegin_p;
    lu_int *Ubegin          = this->Ubegin;
    lu_int *Wbegin          = this->Wbegin;
    lu_int *Wend            = this->Wend;
    lu_int *Wflink          = this->Wflink;
    lu_int *Wblink          = this->Wblink;
    lu_int *Lindex          = this->Lindex;
    double *Lvalue          = this->Lvalue;
    lu_int *Uindex          = this->Uindex;
    double *Uvalue          = this->Uvalue;
    lu_int *Windex          = this->Windex;
    double *Wvalue          = this->Wvalue;
    lu_int *marked          = this->iwork0;
    double *work            = this->work0;
    
    lu_int cbeg = Wbegin[pivot_col]; /* changed by file compression */
    lu_int cend = Wend[pivot_col];
    lu_int rbeg = Wbegin[m+pivot_row];
    lu_int rend = Wend[m+pivot_row];
    const lu_int cnz1 = cend-cbeg-1;   /* nz in pivot column except pivot */
    const lu_int rnz1 = rend-rbeg-1;   /* nz in pivot row except pivot */

    lu_int i, j, pos, pos1, rpos, put, Uput, where, nz, *wi;
    lu_int grow, room, found, position;
    double pivot, a, x, cmx, xrj, *wx;

    /*
     * Check if room is available in W. At most each updated row and each
     * updated column will be reappended and filled-in with rnz1 respectively
     * cnz1 elements. Move pivot to the front of pivot row and pivot column.
     */
    grow = 0;
    where = -1;
    for (pos = cbeg; pos < cend; pos++)
    {
        if ((i = Windex[pos]) == pivot_row)
            where = pos;
        else
        {
            nz = Wend[m+i] - Wbegin[m+i];
            grow += nz+rnz1 + stretch*(nz+rnz1) + pad;
        }
    }
    assert(where >= 0);
    lu_iswap(Windex, cbeg, where);
    lu_fswap(Wvalue, cbeg, where);
    pivot = Wvalue[cbeg];
    assert(pivot);
    where = -1;
    for (rpos = rbeg; rpos < rend; rpos++)
    {
        if ((j = Windex[rpos]) == pivot_col)
            where = rpos;
        else
        {
            nz = Wend[j] - Wbegin[j];
            grow += nz+cnz1 + stretch*(nz+cnz1) + pad;
        }
    }
    assert(where >= 0);
    lu_iswap(Windex, rbeg, where);
    room = Wend[2*m] - Wbegin[2*m];
    if (grow > room)
    {
        lu_file_compress(2*m,Wbegin,Wend,Wflink,Windex,Wvalue,stretch,pad);
        cbeg = Wbegin[pivot_col];
        cend = Wend[pivot_col];
        rbeg = Wbegin[m+pivot_row];
        rend = Wend[m+pivot_row];
        room = Wend[2*m] - Wbegin[2*m];
        this->ngarbage++;
    }
    if (grow > room)
    {
        this->addmemW = grow-room;
        return BASICLU_REALLOCATE;
    }

    /* get pointer to U */
    Uput = Ubegin[rank];
    assert(Uput >= 0);
    assert(Uput < this->Umem);

    /* ---------------------------------------------------------------------- */
    /* Column file update */
    /* ---------------------------------------------------------------------- */

    /* For each row i to be updated set marked[i] > 0 to its position
       in the (packed) pivot column. */
    position = 1;
    for (pos = cbeg+1; pos < cend; pos++)
    {
        i = Windex[pos];
        marked[i] = position++;
    }

    wi = Windex + cbeg;
    wx = Wvalue + cbeg;
    for (rpos = rbeg+1; rpos < rend; rpos++)
    {
        j = Windex[rpos];
        assert(j != pivot_col);
        cmx = 0.0;              /* column maximum */

        /* Compress unmodified column entries. Store entries to be updated
           in workspace. Move pivot row entry to the front of column. */
        where = -1;
        put = pos1 = Wbegin[j];
        for (pos = pos1; pos < Wend[j]; pos++)
        {
            i = Windex[pos];
            if ((position = marked[i]) > 0)
            {
                assert(i != pivot_row);
                work[position] = Wvalue[pos];
            }
            else
            {
                assert(position == 0);
                if (i == pivot_row)
                    where = put;
                else if ((x = fabs(Wvalue[pos])) > cmx)
                    cmx = x;
                Windex[put] = Windex[pos];
                Wvalue[put++] = Wvalue[pos];
            }
        }
        assert(where >= 0);
        Wend[j] = put;
        lu_iswap(Windex, pos1, where);
        lu_fswap(Wvalue, pos1, where);
        xrj = Wvalue[pos1];     /* pivot row entry */

        /* Reappend column if no room for update. */
        room = Wbegin[Wflink[j]] - put;
        if (room < cnz1)
        {
            nz = Wend[j] - Wbegin[j];
            room = cnz1 + stretch*(nz+cnz1) + pad;
            lu_file_reappend(j, 2*m, Wbegin, Wend, Wflink, Wblink, Windex,
                             Wvalue, room);
            put = Wend[j];
            assert(Wbegin[Wflink[j]] - put == room);
            this->nexpand++;
        }

        /* Compute update in workspace and append to column. */
        a = xrj/pivot;
        for (pos = 1; pos <= cnz1; pos++)
            work[pos] -= a * wx[pos];
        for (pos = 1; pos <= cnz1; pos++)
        {
            Windex[put] = wi[pos];
            Wvalue[put++] = work[pos];
            if ((x = fabs(work[pos])) > cmx)
                cmx = x;
            work[pos] = 0.0;
        }
        Wend[j] = put;

        /* Write pivot row entry to U and remove from file. */
        if (fabs(xrj) > droptol)
        {
            assert(Uput < this->Umem);
            Uindex[Uput] = j;
            Uvalue[Uput++] = xrj;
        }
        assert(Windex[Wbegin[j]] == pivot_row);
        Wbegin[j]++;

        /* Move column to new list and update min_colnz. */
        nz = Wend[j] - Wbegin[j];
        lu_list_move(j, nz, colcount_flink, colcount_blink, m,
                     &this->min_colnz);

        colmax[j] = cmx;
    }
    for (pos = cbeg+1; pos < cend; pos++)
        marked[Windex[pos]] = 0;

    /* ---------------------------------------------------------------------- */
    /* Row file update */
    /* ---------------------------------------------------------------------- */

    for (rpos = rbeg; rpos < rend; rpos++)
        marked[Windex[rpos]] = 1;
    assert(marked[pivot_col] == 1);

    for (pos = cbeg+1; pos < cend; pos++)
    {
        i = Windex[pos];
        assert(i != pivot_row);

        /* Compress unmodified row entries (not marked). Remove
           overlap with pivot row, including pivot column entry. */
        found = 0;
        put = Wbegin[m+i];
        for (rpos = Wbegin[m+i]; rpos < Wend[m+i]; rpos++)
        {
            if ((j = Windex[rpos]) == pivot_col)
                found = 1;
            if (!marked[j])
                Windex[put++] = j;
        }
        assert(found);
        Wend[m+i] = put;

        /* Reappend row if no room for update. Append pattern of pivot row. */
        room = Wbegin[Wflink[m+i]] - put;
        if (room < rnz1)
        {
            nz = Wend[m+i] - Wbegin[m+i];
            room = rnz1 + stretch*(nz+rnz1) + pad;
            lu_file_reappend(m+i, 2*m, Wbegin, Wend, Wflink, Wblink, Windex,
                             Wvalue, room);
            put = Wend[m+i];
            assert(Wbegin[Wflink[m+i]] - put == room);
            this->nexpand++;
        }
        for (rpos = rbeg+1; rpos < rend; rpos++)
            Windex[put++] = Windex[rpos];
        Wend[m+i] = put;

        /* Move to new list. The row must be reinserted even if nz are
           unchanged since it might have been taken out in Markowitz search. */
        nz = Wend[m+i] - Wbegin[m+i];
        lu_list_move(i, nz, rowcount_flink, rowcount_blink, m,
                     &this->min_rownz);
    }
    for (rpos = rbeg; rpos < rend; rpos++)
        marked[Windex[rpos]] = 0;


    /* Store column in L. */
    put = Lbegin_p[rank];
    for (pos = cbeg+1; pos < cend; pos++)
    {
        x = Wvalue[pos] / pivot;
        if (fabs(x) > droptol)
        {
            Lindex[put] = Windex[pos];
            Lvalue[put++] = x;
        }
    }
    Lindex[put++] = -1;         /* terminate column */
    Lbegin_p[rank+1] = put;
    Ubegin[rank+1] = Uput;

    /*
     * Cleanup:
     * store pivot elemnt;
     * remove pivot colum from column file, pivot row from row file;
     * remove pivot column/row from column/row counts
     */
    colmax[pivot_col] = pivot;
    Wend[pivot_col] = cbeg;
    Wend[m+pivot_row] = rbeg;
    lu_list_remove(colcount_flink, colcount_blink, pivot_col);
    lu_list_remove(rowcount_flink, rowcount_blink, pivot_row);

    /* 
     * Check that row file and column file are consistent. Only use when
     * DEBUG_EXTRA since this check is really expensive.
     */
    #ifdef DEBUG_EXTRA
    assert(lu_file_diff(m, Wbegin+m, Wend+m, Wbegin, Wend, Windex, NULL) == 0);
    assert(lu_file_diff(m, Wbegin, Wend, Wbegin+m, Wend+m, Windex, NULL) == 0);
    #endif

    return BASICLU_OK;
}


/* ==========================================================================
   lu_pivot_small
   ========================================================================== */

static lu_int lu_pivot_small(struct lu *this)
{
    const lu_int m          = this->m;
    const lu_int rank       = this->rank;
    const double droptol    = this->droptol;
    const lu_int pad        = this->pad;
    const double stretch    = this->stretch;
    const lu_int pivot_col  = this->pivot_col;
    const lu_int pivot_row  = this->pivot_row;
    lu_int *colcount_flink  = this->colcount_flink;
    lu_int *colcount_blink  = this->colcount_blink;
    lu_int *rowcount_flink  = this->rowcount_flink;
    lu_int *rowcount_blink  = this->rowcount_blink;
    double *colmax          = this->col_pivot;
    lu_int *Lbegin_p        = this->Lbegin_p;
    lu_int *Ubegin          = this->Ubegin;
    lu_int *Wbegin          = this->Wbegin;
    lu_int *Wend            = this->Wend;
    lu_int *Wflink          = this->Wflink;
    lu_int *Wblink          = this->Wblink;
    lu_int *Lindex          = this->Lindex;
    double *Lvalue          = this->Lvalue;
    lu_int *Uindex          = this->Uindex;
    double *Uvalue          = this->Uvalue;
    lu_int *Windex          = this->Windex;
    double *Wvalue          = this->Wvalue;
    lu_int *marked          = this->iwork0;
    double *work            = this->work0;
    int64_t *cancelled      = (void *) this->row_pivot;
    
    lu_int cbeg = Wbegin[pivot_col]; /* changed by file compression */
    lu_int cend = Wend[pivot_col];
    lu_int rbeg = Wbegin[m+pivot_row];
    lu_int rend = Wend[m+pivot_row];
    const lu_int cnz1 = cend-cbeg-1;   /* nz in pivot column except pivot */
    const lu_int rnz1 = rend-rbeg-1;   /* nz in pivot row except pivot */

    lu_int i, j, pos, pos1, rpos, put, Uput, where, nz, *wi;
    lu_int grow, room, found, position, col_number;
    double pivot, a, x, cmx, xrj, *wx;
    int64_t mask;

    assert(cnz1 <= MAXROW_SMALL);

    /*
     * Check if room is available in W. At most each updated row and each
     * updated column will be reappended and filled-in with rnz1 respectively
     * cnz1 elements. Move pivot to the front of pivot row and pivot column.
      */
    grow = 0;
    where = -1;
    for (pos = cbeg; pos < cend; pos++)
    {
        if ((i = Windex[pos]) == pivot_row)
            where = pos;
        else
        {
            nz = Wend[m+i] - Wbegin[m+i];
            grow += nz+rnz1 + stretch*(nz+rnz1) + pad;
        }
    }
    assert(where >= 0);
    lu_iswap(Windex, cbeg, where);
    lu_fswap(Wvalue, cbeg, where);
    pivot = Wvalue[cbeg];
    assert(pivot);
    where = -1;
    for (rpos = rbeg; rpos < rend; rpos++)
    {
        if ((j = Windex[rpos]) == pivot_col)
            where = rpos;
        else
        {
            nz = Wend[j] - Wbegin[j];
            grow += nz+cnz1 + stretch*(nz+cnz1) + pad;
        }
    }
    assert(where >= 0);
    lu_iswap(Windex, rbeg, where);
    room = Wend[2*m] - Wbegin[2*m];
    if (grow > room)
    {
        lu_file_compress(2*m,Wbegin,Wend,Wflink,Windex,Wvalue,stretch,pad);
        cbeg = Wbegin[pivot_col];
        cend = Wend[pivot_col];
        rbeg = Wbegin[m+pivot_row];
        rend = Wend[m+pivot_row];
        room = Wend[2*m] - Wbegin[2*m];
        this->ngarbage++;
    }
    if (grow > room)
    {
        this->addmemW = grow-room;
        return BASICLU_REALLOCATE;
    }

    /* get pointer to U */
    Uput = Ubegin[rank];
    assert(Uput >= 0);
    assert(Uput < this->Umem);

    /* ---------------------------------------------------------------------- */
    /* Column file update */
    /* ---------------------------------------------------------------------- */

    /* For each row i to be updated set marked[i] > 0 to its position
       in the (packed) pivot column. */
    position = 1;
    for (pos = cbeg+1; pos < cend; pos++)
    {
        i = Windex[pos];
        marked[i] = position++;
    }

    wi = Windex + cbeg;
    wx = Wvalue + cbeg;
    col_number = 0;             /* mask cancelled[col_number] */
    for (rpos = rbeg+1; rpos < rend; rpos++, col_number++)
    {
        j = Windex[rpos];
        assert(j != pivot_col);
        cmx = 0.0;              /* column maximum */

        /* Compress unmodified column entries. Store entries to be updated
           in workspace. Move pivot row entry to the front of column. */
        where = -1;
        put = pos1 = Wbegin[j];
        for (pos = pos1; pos < Wend[j]; pos++)
        {
            i = Windex[pos];
            if ((position = marked[i]) > 0)
            {
                assert(i != pivot_row);
                work[position] = Wvalue[pos];
            }
            else
            {
                assert(position == 0);
                if (i == pivot_row)
                    where = put;
                else if ((x = fabs(Wvalue[pos])) > cmx)
                    cmx = x;
                Windex[put] = Windex[pos];
                Wvalue[put++] = Wvalue[pos];
            }
        }
        assert(where >= 0);
        Wend[j] = put;
        lu_iswap(Windex, pos1, where);
        lu_fswap(Wvalue, pos1, where);
        xrj = Wvalue[pos1];     /* pivot row entry */

        /* Reappend column if no room for update. */
        room = Wbegin[Wflink[j]] - put;
        if (room < cnz1)
        {
            nz = Wend[j] - Wbegin[j];
            room = cnz1 + stretch*(nz+cnz1) + pad;
            lu_file_reappend(j, 2*m, Wbegin, Wend, Wflink, Wblink, Windex,
                             Wvalue, room);
            put = Wend[j];
            assert(Wbegin[Wflink[j]] - put == room);
            this->nexpand++;
        }

        /* Compute update in workspace and append to column. */
        a = xrj/pivot;
        for (pos = 1; pos <= cnz1; pos++)
            work[pos] -= a * wx[pos];
        mask = 0;
        for (pos = 1; pos <= cnz1; pos++)
        {
            x = fabs(work[pos]);
            if (x > droptol)
            {
                Windex[put] = wi[pos];
                Wvalue[put++] = work[pos];
                if (x > cmx)
                    cmx = x;
            }
            else
            {
                /* cancellation in row wi[pos] */
                mask |= (int64_t) 1 << (pos-1);
            }
            work[pos] = 0.0;
        }
        Wend[j] = put;
        cancelled[col_number] = mask;

        /* Write pivot row entry to U and remove from file. */
        if (fabs(xrj) > droptol)
        {
            assert(Uput < this->Umem);
            Uindex[Uput] = j;
            Uvalue[Uput++] = xrj;
        }
        assert(Windex[Wbegin[j]] == pivot_row);
        Wbegin[j]++;

        /* Move column to new list and update min_colnz. */
        nz = Wend[j] - Wbegin[j];
        lu_list_move(j, nz, colcount_flink, colcount_blink, m,
                     &this->min_colnz);

        colmax[j] = cmx;
    }
    for (pos = cbeg+1; pos < cend; pos++)
        marked[Windex[pos]] = 0;

    /* ---------------------------------------------------------------------- */
    /* Row file update */
    /* ---------------------------------------------------------------------- */

    for (rpos = rbeg; rpos < rend; rpos++)
        marked[Windex[rpos]] = 1;
    assert(marked[pivot_col] == 1);

    mask = 1;
    for (pos = cbeg+1; pos < cend; pos++, mask <<= 1)
    {
        assert(mask);
        i = Windex[pos];
        assert(i != pivot_row);

        /* Compress unmodified row entries (not marked). Remove
           overlap with pivot row, including pivot column entry. */
        found = 0;
        put = Wbegin[m+i];
        for (rpos = Wbegin[m+i]; rpos < Wend[m+i]; rpos++)
        {
            if ((j = Windex[rpos]) == pivot_col)
                found = 1;
            if (!marked[j])
                Windex[put++] = j;
        }
        assert(found);
        Wend[m+i] = put;

        /* Reappend row if no room for update. Append pattern of pivot row. */
        room = Wbegin[Wflink[m+i]] - put;
        if (room < rnz1)
        {
            nz = Wend[m+i] - Wbegin[m+i];
            room = rnz1 + stretch*(nz+rnz1) + pad;
            lu_file_reappend(m+i, 2*m, Wbegin, Wend, Wflink, Wblink, Windex,
                             Wvalue, room);
            put = Wend[m+i];
            assert(Wbegin[Wflink[m+i]] - put == room);
            this->nexpand++;
        }

        col_number = 0;
        for (rpos = rbeg+1; rpos < rend; rpos++, col_number++)
        {
            if (! (cancelled[col_number] & mask))
            {
                Windex[put++] = Windex[rpos];
            }
        }
        Wend[m+i] = put;

        /* Move to new list. The row must be reinserted even if nz are
           unchanged since it might have been taken out in Markowitz search. */
        nz = Wend[m+i] - Wbegin[m+i];
        lu_list_move(i, nz, rowcount_flink, rowcount_blink, m,
                     &this->min_rownz);
    }
    for (rpos = rbeg; rpos < rend; rpos++)
        marked[Windex[rpos]] = 0;


    /* Store column in L. */
    put = Lbegin_p[rank];
    for (pos = cbeg+1; pos < cend; pos++)
    {
        x = Wvalue[pos] / pivot;
        if (fabs(x) > droptol)
        {
            Lindex[put] = Windex[pos];
            Lvalue[put++] = x;
        }
    }
    Lindex[put++] = -1;         /* terminate column */
    Lbegin_p[rank+1] = put;
    Ubegin[rank+1] = Uput;

    /*
     * Cleanup:
     * store pivot elemnt;
     * remove pivot colum from column file, pivot row from row file;
     * remove pivot column/row from column/row counts
     */
    colmax[pivot_col] = pivot;
    Wend[pivot_col] = cbeg;
    Wend[m+pivot_row] = rbeg;
    lu_list_remove(colcount_flink, colcount_blink, pivot_col);
    lu_list_remove(rowcount_flink, rowcount_blink, pivot_row);

    /*
     * Check that row file and column file are consistent. Only use when
     * DEBUG_EXTRA since this check is really expensive.
     */
    #ifdef DEBUG_EXTRA
    assert(lu_file_diff(m, Wbegin+m, Wend+m, Wbegin, Wend, Windex, NULL) == 0);
    assert(lu_file_diff(m, Wbegin, Wend, Wbegin+m, Wend+m, Windex, NULL) == 0);
    #endif

    return BASICLU_OK;
}


/* ==========================================================================
   lu_pivot_singleton_row
   ========================================================================== */

static lu_int lu_pivot_singleton_row(struct lu *this)
{
    const lu_int m          = this->m;
    const lu_int rank       = this->rank;
    const double droptol    = this->droptol;
    const lu_int pivot_col  = this->pivot_col;
    const lu_int pivot_row  = this->pivot_row;
    lu_int *colcount_flink  = this->colcount_flink;
    lu_int *colcount_blink  = this->colcount_blink;
    lu_int *rowcount_flink  = this->rowcount_flink;
    lu_int *rowcount_blink  = this->rowcount_blink;
    double *colmax          = this->col_pivot;
    lu_int *Lbegin_p        = this->Lbegin_p;
    lu_int *Ubegin          = this->Ubegin;
    lu_int *Wbegin          = this->Wbegin;
    lu_int *Wend            = this->Wend;
    lu_int *Lindex          = this->Lindex;
    double *Lvalue          = this->Lvalue;
    lu_int *Windex          = this->Windex;
    double *Wvalue          = this->Wvalue;
    
    const lu_int cbeg = Wbegin[pivot_col];
    const lu_int cend = Wend[pivot_col];
    const lu_int rbeg = Wbegin[m+pivot_row];
    const lu_int rend = Wend[m+pivot_row];
    const lu_int rnz1 = rend-rbeg-1;   /* nz in pivot row except pivot */

    lu_int i, pos, put, nz, where;
    double pivot, x;

    assert(rnz1 == 0);

    /* Find pivot. */
    for (where = cbeg; Windex[where] != pivot_row; where++)
        assert(where < cend-1);
    pivot = Wvalue[where];
    assert(pivot);

    /* Store column in L. */
    put = Lbegin_p[rank];
    for (pos = cbeg; pos < cend; pos++)
    {
        x = Wvalue[pos] / pivot;
        if (pos != where && fabs(x) > droptol)
        {
            Lindex[put] = Windex[pos];
            Lvalue[put++] = x;
        }
    }
    Lindex[put++] = -1;         /* terminate column */
    Lbegin_p[rank+1] = put;
    Ubegin[rank+1] = Ubegin[rank];

    /* Remove pivot column from row file. Update row lists. */
    for (pos = cbeg; pos < cend; pos++)
    {
        i = Windex[pos];
        if (i == pivot_row)
            continue;
        for (where = Wbegin[m+i]; Windex[where] != pivot_col; where++)
            assert(where < Wend[m+i]-1);
        Windex[where] = Windex[--Wend[m+i]];
        nz = Wend[m+i] - Wbegin[m+i];
        lu_list_move(i, nz, rowcount_flink, rowcount_blink, m,
                     &this->min_rownz);
    }

    /*
     * Cleanup:
     * store pivot elemnt;
     * remove pivot colum from column file, pivot row from row file;
     * remove pivot column/row from column/row counts
     */
    colmax[pivot_col] = pivot;
    Wend[pivot_col] = cbeg;
    Wend[m+pivot_row] = rbeg;
    lu_list_remove(colcount_flink, colcount_blink, pivot_col);
    lu_list_remove(rowcount_flink, rowcount_blink, pivot_row);

    return BASICLU_OK;
}


/* ==========================================================================
   lu_pivot_singleton_col
   ========================================================================== */

static lu_int lu_pivot_singleton_col(struct lu *this)
{
    const lu_int m          = this->m;
    const lu_int rank       = this->rank;
    const double droptol    = this->droptol;
    const lu_int pivot_col  = this->pivot_col;
    const lu_int pivot_row  = this->pivot_row;
    lu_int *colcount_flink  = this->colcount_flink;
    lu_int *colcount_blink  = this->colcount_blink;
    lu_int *rowcount_flink  = this->rowcount_flink;
    lu_int *rowcount_blink  = this->rowcount_blink;
    double *colmax          = this->col_pivot;
    lu_int *Lbegin_p        = this->Lbegin_p;
    lu_int *Ubegin          = this->Ubegin;
    lu_int *Wbegin          = this->Wbegin;
    lu_int *Wend            = this->Wend;
    lu_int *Lindex          = this->Lindex;
    lu_int *Uindex          = this->Uindex;
    double *Uvalue          = this->Uvalue;
    lu_int *Windex          = this->Windex;
    double *Wvalue          = this->Wvalue;
    
    const lu_int cbeg = Wbegin[pivot_col];
    const lu_int cend = Wend[pivot_col];
    const lu_int rbeg = Wbegin[m+pivot_row];
    const lu_int rend = Wend[m+pivot_row];
    const lu_int cnz1 = cend-cbeg-1;   /* nz in pivot column except pivot */

    lu_int j, pos, rpos, put, nz, where, found;
    double pivot, cmx, x, xrj;

    assert(cnz1 == 0);

    /* Remove pivot row from column file and store in U. Update column lists. */
    put = Ubegin[rank];
    pivot = Wvalue[cbeg];
    assert(pivot);
    found = 0;
    xrj = 0.0;                  /* initialize to make gcc happy */
    for (rpos = rbeg; rpos < rend; rpos++)
    {
        j = Windex[rpos];
        if (j == pivot_col)
        {
            found = 1;
            continue;
        }
        where = -1;
        cmx = 0.0;              /* column maximum */
        for (pos = Wbegin[j]; pos < Wend[j]; pos++)
        {
            if (Windex[pos] == pivot_row)
            {
                where = pos;
                xrj = Wvalue[pos];
            }
            else if ((x = fabs(Wvalue[pos])) > cmx)
                cmx = x;
        }
        assert(where >= 0);
        if (fabs(xrj) > droptol)
        {
            Uindex[put] = j;
            Uvalue[put++] = xrj;
        }
        Windex[where] = Windex[--Wend [j]];
        Wvalue[where] = Wvalue[Wend [j]];
        nz = Wend[j] - Wbegin[j];
        lu_list_move(j, nz, colcount_flink, colcount_blink, m,
                     &this->min_colnz);
        colmax[j] = cmx;
    }
    assert(found);
    Ubegin[rank+1] = put;

    /* Store empty column in L. */
    put = Lbegin_p[rank];
    Lindex[put++] = -1;         /* terminate column */
    Lbegin_p[rank+1] = put;

    /*
     * Cleanup:
     * store pivot elemnt;
     * remove pivot colum from column file, pivot row from row file;
     * remove pivot column/row from column/row counts
     */
    colmax[pivot_col] = pivot;
    Wend[pivot_col] = cbeg;
    Wend[m+pivot_row] = rbeg;
    lu_list_remove(colcount_flink, colcount_blink, pivot_col);
    lu_list_remove(rowcount_flink, rowcount_blink, pivot_row);

    return BASICLU_OK;
}


/* ==========================================================================
   lu_pivot_doubleton_col
   ========================================================================== */

static lu_int lu_pivot_doubleton_col(struct lu *this)
{
    const lu_int m          = this->m;
    const lu_int rank       = this->rank;
    const double droptol    = this->droptol;
    const lu_int pad        = this->pad;
    const double stretch    = this->stretch;
    const lu_int pivot_col  = this->pivot_col;
    const lu_int pivot_row  = this->pivot_row;
    lu_int *colcount_flink  = this->colcount_flink;
    lu_int *colcount_blink  = this->colcount_blink;
    lu_int *rowcount_flink  = this->rowcount_flink;
    lu_int *rowcount_blink  = this->rowcount_blink;
    double *colmax          = this->col_pivot;
    lu_int *Lbegin_p        = this->Lbegin_p;
    lu_int *Ubegin          = this->Ubegin;
    lu_int *Wbegin          = this->Wbegin;
    lu_int *Wend            = this->Wend;
    lu_int *Wflink          = this->Wflink;
    lu_int *Wblink          = this->Wblink;
    lu_int *Lindex          = this->Lindex;
    double *Lvalue          = this->Lvalue;
    lu_int *Uindex          = this->Uindex;
    double *Uvalue          = this->Uvalue;
    lu_int *Windex          = this->Windex;
    double *Wvalue          = this->Wvalue;
    lu_int *marked          = this->iwork0;
    
    lu_int cbeg = Wbegin[pivot_col]; /* changed by file compression */
    lu_int cend = Wend[pivot_col];
    lu_int rbeg = Wbegin[m+pivot_row];
    lu_int rend = Wend[m+pivot_row];
    const lu_int cnz1 = cend-cbeg-1;   /* nz in pivot column except pivot */
    const lu_int rnz1 = rend-rbeg-1;   /* nz in pivot row except pivot */

    lu_int j, pos, rpos, put, Uput, nz, nfill, where, where_pivot, where_other;
    lu_int other_row, grow, room, space, end, ncancelled;
    double pivot, other_value, xrj, cmx, x, xabs;

    assert(cnz1 == 1);

    /* Move pivot element to front of pivot column and pivot row. */
    if (Windex[cbeg] != pivot_row)
    {
        lu_iswap(Windex, cbeg, cbeg+1);
        lu_fswap(Wvalue, cbeg, cbeg+1);
    }
    assert(Windex[cbeg] == pivot_row);
    pivot = Wvalue[cbeg];
    assert(pivot);
    other_row = Windex[cbeg+1];
    other_value = Wvalue[cbeg+1];
    for (where = rbeg; Windex[where] != pivot_col; where++)
        assert(where < rend-1);
    lu_iswap(Windex, rbeg, where);

    /*
     * Check if room is available in W.
     * Columns can be updated in place but the updated row may need to be
     * expanded.
     */
    nz = Wend[m+other_row] - Wbegin[m+other_row];
    grow = nz+rnz1 + stretch*(nz+rnz1) + pad;
    room = Wend[2*m] - Wbegin[2*m];
    if (grow > room)
    {
        lu_file_compress(2*m,Wbegin,Wend,Wflink,Windex,Wvalue,stretch,pad);
        cbeg = Wbegin[pivot_col];
        cend = Wend[pivot_col];
        rbeg = Wbegin[m+pivot_row];
        rend = Wend[m+pivot_row];
        room = Wend[2*m] - Wbegin[2*m];
        this->ngarbage++;
    }
    if (grow > room)
    {
        this->addmemW = grow-room;
        return BASICLU_REALLOCATE;
    }

    /* ---------------------------------------------------------------------- */
    /* Column file update */
    /* ---------------------------------------------------------------------- */

    Uput = Ubegin[rank];
    put = rbeg+1;
    ncancelled = 0;
    for (rpos = rbeg+1; rpos < rend; rpos++)
    {
        j = Windex[rpos];
        assert(j != pivot_col);
        cmx = 0.0;              /* column maximum */

        /* Find position of pivot row entry and possibly other row entry in
           column j. */
        where_pivot = -1;
        where_other = -1;
        end = Wend[j];
        for (pos = Wbegin[j]; pos < end; pos++)
        {
            if (Windex[pos] == pivot_row)
                where_pivot = pos;
            else if (Windex[pos] == other_row)
                where_other = pos;
            else if ((x = fabs(Wvalue[pos])) > cmx)
                cmx = x;
        }
        assert(where_pivot >= 0);
        xrj = Wvalue[where_pivot];

        /* Store pivot row entry in U. */
        if (fabs(Wvalue[where_pivot]) > droptol)
        {
            Uindex[Uput] = j;
            Uvalue[Uput++] = Wvalue[where_pivot];
        }

        if (where_other == -1)
        {
            /* Compute fill-in element. */
            x = -xrj * (other_value / pivot);
            xabs = fabs(x);
            if (xabs > droptol)
            {
                /* Store fill-in where pivot row entry was. */
                Windex[where_pivot] = other_row;
                Wvalue[where_pivot] = x;
                Windex[put++] = j;
                if (xabs > cmx)
                    cmx = xabs;
            }
            else
            {
                /* Remove pivot row entry. */
                end = --Wend[j];
                Windex[where_pivot] = Windex[end];
                Wvalue[where_pivot] = Wvalue[end];

                /* Decrease column count. */
                nz = end - Wbegin[j];
                lu_list_move(j, nz, colcount_flink, colcount_blink, m,
                             &this->min_colnz);
            }
        }
        else
        {
            /* Remove pivot row entry and update other row entry. */
            end = --Wend[j];
            Windex[where_pivot] = Windex[end];
            Wvalue[where_pivot] = Wvalue[end];
            if (where_other == end)
                where_other = where_pivot;
            Wvalue[where_other] -= xrj * (other_value/pivot);

            /* If we have numerical cancellation, then remove the entry and mark
               the column. */
            x = fabs(Wvalue[where_other]);
            if (x <= droptol)
            {
                end = --Wend[j];
                Windex[where_other] = Windex[end];
                Wvalue[where_other] = Wvalue[end];
                marked[j] = 1;
                ncancelled++;
            }
            else if (x > cmx)
                cmx = x;

            /* Decrease column count. */
            nz = Wend[j] - Wbegin[j];
            lu_list_move(j, nz, colcount_flink, colcount_blink, m,
                         &this->min_colnz);
        }
        colmax[j] = cmx;
    }
    rend = put;
    Ubegin[rank+1] = Uput;
    

    /* ---------------------------------------------------------------------- */
    /* Row file update */
    /* ---------------------------------------------------------------------- */

    /* If we have numerical cancellation, then we have to remove these entries
       (marked) from the row pattern. In any case remove pivot column entry. */
    if (ncancelled)
    {
        assert(marked[pivot_col] == 0);
        marked[pivot_col] = 1;     /* treat as cancelled */
        put = Wbegin[m+other_row]; /* compress remaining entries */
        end = Wend[m+other_row];
        for (pos = put; pos < end; pos++)
        {
            j = Windex[pos];
            if (marked[j])
                marked[j] = 0;
            else
                Windex[put++] = j;
        }
        assert(end-put == ncancelled+1);
        Wend[m+other_row] = put;
    }
    else
    {
        for (where = Wbegin[m+other_row]; Windex[where] != pivot_col; where++)
            assert (where < Wend[m+other_row]-1);
        end = --Wend[m+other_row];
        Windex[where] = Windex[end];
    }

    /* Reappend row if no room for update. */
    nfill = rend - (rbeg+1);
    room = Wbegin[Wflink [m+other_row]] - Wend[m+other_row];
    if (nfill > room)
    {
        nz = Wend[m+other_row] - Wbegin[m+other_row];
        space = nfill + stretch*(nz+nfill) + pad;
        lu_file_reappend(m+other_row, 2*m, Wbegin, Wend, Wflink, Wblink,
                         Windex, Wvalue, space);
        this->nexpand++;
    }

    /* Append fill-in to row pattern. */
    put = Wend[m+other_row];
    for (pos = rbeg+1; pos < rend; pos++)
        Windex[put++] = Windex[pos];
    Wend[m+other_row] = put;

    /* Reinsert other row into row counts. */
    nz = Wend[m+other_row] - Wbegin[m+other_row];
    lu_list_move(other_row, nz, rowcount_flink, rowcount_blink, m,
                 &this->min_rownz);

    /* Store column in L. */
    put = Lbegin_p[rank];
    x = other_value / pivot;
    if (fabs(x) > droptol)
    {
        Lindex[put] = other_row;
        Lvalue[put++] = x;
    }
    Lindex[put++] = -1;         /* terminate column */
    Lbegin_p[rank+1] = put;

    /*
     * Cleanup:
     * store pivot elemnt;
     * remove pivot colum from column file, pivot row from row file;
     * remove pivot column/row from column/row counts
     */
    colmax[pivot_col] = pivot;
    Wend[pivot_col] = cbeg;
    Wend[m+pivot_row] = rbeg;
    lu_list_remove(colcount_flink, colcount_blink, pivot_col);
    lu_list_remove(rowcount_flink, rowcount_blink, pivot_row);

    /*
     * Check that row file and column file are consistent. Only use when
     * DEBUG_EXTRA since this check is really expensive.
     */
    #ifdef DEBUG_EXTRA
    assert(lu_file_diff(m, Wbegin+m, Wend+m, Wbegin, Wend, Windex, NULL) == 0);
    assert(lu_file_diff(m, Wbegin, Wend, Wbegin+m, Wend+m, Windex, NULL) == 0);
    #endif

    return BASICLU_OK;
}


/* ==========================================================================
   lu_remove_col
   ========================================================================== */

static void lu_remove_col(struct lu *this, lu_int j)
{
    const lu_int m          = this->m;
    lu_int *colcount_flink  = this->colcount_flink;
    lu_int *colcount_blink  = this->colcount_blink;
    lu_int *rowcount_flink  = this->rowcount_flink;
    lu_int *rowcount_blink  = this->rowcount_blink;
    double *colmax          = this->col_pivot;
    lu_int *Wbegin          = this->Wbegin;
    lu_int *Wend            = this->Wend;
    lu_int *Windex          = this->Windex;
    const lu_int cbeg       = Wbegin[j];
    const lu_int cend       = Wend[j];

    lu_int i, pos, nz, where;

    /* Remove column j from row file. */
    for (pos = cbeg; pos < cend; pos++)
    {
        i = Windex[pos];
        for (where = Wbegin[m+i]; Windex[where] != j; where++)
            assert(where < Wend[m+i]-1);
        Windex[where] = Windex[--Wend[m+i]];
        nz = Wend[m+i] - Wbegin[m+i];
        lu_list_move(i, nz, rowcount_flink, rowcount_blink, m,
                     &this->min_rownz);
    }

    /* Remove column j from column file. */
    colmax[j] = 0.0;
    Wend[j] = cbeg;
    lu_list_move(j, 0, colcount_flink, colcount_blink, m, &this->min_colnz);
}
