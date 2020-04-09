/*
 * lu_setup_bump.c
 *
 * Copyright (C) 2016-2019  ERGO-Code
 *
 * Setup data structures for bump factorization
 *
 */

#include "lu_internal.h"
#include "lu_list.h"
#include "lu_file.h"


/*
 * lu_setup_bump()
 *
 * The bump is composed of rows i and columns j for which pinv[i] < 0 and
 * qinv[j] < 0. For the factorization, the bump is stored in Windex, Wvalue
 * columnwise and additionally the nonzero pattern rowwise:
 *
 *  Wbegin[j]   points to the first element in column j.
 *  Wend[j]     points to one past the last element in colum j.
 *  Wbegin[m+i] points to the first element in row i.
 *  Wend[m+i]   points to one past the last element in row i.
 *
 *  Wflink, Wblink hold the 2*m lines in a double linked list in memory order.
 *
 * When a row or column is empty, then Wbegin == Wend. In the rowwise storage
 * the entries in Wvalue are undefined.
 *
 * The Markowitz search requires double linked lists of columns with equal
 * column counts and rows with equal row counts:
 *
 *  colcount_flink, colcount_blink
 *  rowcount_flink, rowcount_blink
 *
 * They organize m elements (cols/rows) in m+2 lists. Column j is in list
 * 0 <= nz <= m when it has nz nonzeros in the active submatrix. Row i can
 * alternatively be in list m+1 to exclude it temporarily from the search.
 * A column/row not in the active submatrix is not in any list.
 *
 * The Markowitz search also requires the maximum in each column of the
 * active submatrix. For column j the maximum is stored in col_pivot[j].
 * When j becomes pivot column, the maximum is replaced by the pivot element.
 *
 * Return:
 *
 *  BASICLU_REALLOCATE  require more memory in W
 *  BASICLU_OK
 */

lu_int lu_setup_bump(
    struct lu *this, const lu_int *Bbegin, const lu_int *Bend, const lu_int *Bi,
    const double *Bx)
{
    const lu_int m          = this->m;
    const lu_int rank       = this->rank;
    const lu_int Wmem       = this->Wmem;
    const lu_int Bnz        = this->matrix_nz;
    const lu_int Lnz        = this->Lbegin_p[rank] - rank;
    const lu_int Unz        = this->Ubegin[rank];
    const double abstol     = this->abstol;
    const lu_int pad        = this->pad;
    const double stretch    = this->stretch;
    lu_int *colcount_flink  = this->colcount_flink;
    lu_int *colcount_blink  = this->colcount_blink;
    lu_int *rowcount_flink  = this->rowcount_flink;
    lu_int *rowcount_blink  = this->rowcount_blink;
    const lu_int *pinv      = this->pinv;
    const lu_int *qinv      = this->qinv;
    lu_int *Wbegin          = this->Wbegin;
    lu_int *Wend            = this->Wend;
    lu_int *Wbegin2         = Wbegin + m; /* alias for row file */
    lu_int *Wend2           = Wend + m;
    lu_int *Wflink          = this->Wflink;
    lu_int *Wblink          = this->Wblink;
    lu_int *Windex          = this->Windex;
    double *Wvalue          = this->Wvalue;
    double *colmax          = this->col_pivot;
    lu_int *iwork0          = this->iwork0;

    lu_int bump_nz = Bnz-Lnz-Unz-rank; /* will change if columns are dropped */
    lu_int i, j, pos, put, cnz, rnz, need, min_rownz, min_colnz;
    double cmx;

    assert(Lnz >= 0);
    assert(Unz >= 0);
    assert(bump_nz >= 0);
    #ifndef NDEBUG
    for (i = 0; i < m; i++)
        assert(iwork0[i] == 0);
    #endif

    /*
     * Calculate memory and reallocate. For each row/column with nz nonzeros
     * add stretch*nz+pad elements extra space for fill-in.
     */
    need = bump_nz + stretch*bump_nz + (m-rank)*pad;
    need = 2*need;              /* rowwise + columnwise */
    if (need > Wmem)
    {
        this->addmemW = need-Wmem;
        return BASICLU_REALLOCATE;
    }

    lu_file_empty(2*m, Wbegin, Wend, Wflink, Wblink, Wmem);

    /*
     * Build columnwise storage. Build row counts in iwork0.
     */
    lu_list_init(colcount_flink, colcount_blink, m, m+2, &min_colnz);
    put = 0;
    for (j = 0; j < m; j++)
    {
        if (qinv[j] >= 0)
            continue;
        cnz = 0;                /* count nz per column */
        cmx = 0.0;              /* find column maximum */
        for (pos = Bbegin[j]; pos < Bend[j]; pos++)
        {
            i = Bi[pos];
            if (pinv[i] >= 0)
                continue;
            cmx = fmax(cmx, fabs(Bx[pos]));
            cnz++;
        }
        if (cmx == 0.0 || cmx < abstol)
        {
            /* Leave column of active submatrix empty. */
            colmax[j] = 0.0;
            lu_list_add(j, 0, colcount_flink, colcount_blink, m, &min_colnz);
            bump_nz -= cnz;
        }
        else
        {
            /* Copy column into active submatrix. */
            colmax[j] = cmx;
            lu_list_add(j, cnz, colcount_flink, colcount_blink, m, &min_colnz);
            Wbegin[j] = put;
            for (pos = Bbegin[j]; pos < Bend[j]; pos++)
            {
                i = Bi[pos];
                if (pinv[i] >= 0)
                    continue;
                Windex[put] = i;
                Wvalue[put++] = Bx[pos];
                iwork0[i]++;
            }
            Wend[j] = put;
            put += stretch*cnz + pad;
            /* reappend line to list end */
            lu_list_move(j, 0, Wflink, Wblink, 2*m, NULL);
        }
    }

    /*
     * Build rowwise storage (pattern only).
     */
    lu_list_init(rowcount_flink, rowcount_blink, m, m+2, &min_rownz);
    for (i = 0; i < m; i++)     /* set row pointers */
    {
        if (pinv[i] >= 0)
            continue;
        rnz = iwork0[i];
        iwork0[i] = 0;
        lu_list_add(i, rnz, rowcount_flink, rowcount_blink, m, &min_rownz);
        Wbegin2[i] = Wend2[i] = put;
        put += rnz;
        /* reappend line to list end */
        lu_list_move(m+i, 0, Wflink, Wblink, 2*m, NULL);
        put += stretch*rnz + pad;
    }
    for (j = 0; j < m; j++)     /* fill rows */
    {
        for (pos = Wbegin[j]; pos < Wend[j]; pos++)
        {
            i = Windex[pos];
            Windex[Wend2[i]++] = j;
        }
    }
    Wbegin[2*m] = put;          /* set beginning of free space */
    assert(Wbegin[2*m] <= Wend[2*m]);

    assert(lu_file_diff(m, Wbegin, Wend, Wbegin2, Wend2, Windex, NULL) == 0);
    assert(lu_file_diff(m, Wbegin2, Wend2, Wbegin, Wend, Windex, NULL) == 0);

    #ifndef NDEBUG
    for (i = 0; i < m; i++)
        assert(iwork0[i] == 0);
    #endif

    this->bump_nz = bump_nz;
    this->bump_size = m-rank;
    this->min_colnz = min_colnz;
    this->min_rownz = min_rownz;
    return BASICLU_OK;
}
