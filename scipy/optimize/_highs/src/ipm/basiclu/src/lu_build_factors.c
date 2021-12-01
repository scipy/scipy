/*
 * lu_build_factors.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * Build rowwise and columnwise form of L and U
 *
 */

#include "lu_internal.h"
#include "lu_list.h"
#include "lu_file.h"

/*
   BASICLU maintains the factorization in the form

     B = L * R^1 * R^2 * ... * R^{nforrest} * U,

   where L[p,p] is unit lower triangular and U[pivotrow,pivotcol] is upper
   triangular. After refactorization nforrest = 0 and p and pivotrow hold the
   same permutation. pivotrow and pivotcol are modified by updates, p is not.

   The permutations are stored as follows:
   ---------------------------------------

     p[0..m-1] is a vector.

     pivotrow[0..pivotlen-1],
     pivotcol[0..pivotlen-1]

       are vectors of length m <= pivotlen < 2*m which may contain duplicate
       indices. For each index its last occurance is its position in the pivot
       sequence, see lu_garbage_perm().

     pmap[0..m-1],
     qmap[0..m-1]

       are vectors such that i = pmap[j] and j = qmap[i] when element (i,j) of
       U is pivot element.

   The matrix L is stored as follows:
   ----------------------------------

     Lindex[0..Lnz+m-1],
     Lvalue[0..Lnz+m-1]

       hold L columnwise without the unit diagonal. Each column is terminated
       by index -1. Row indices are row indices of B.

     Lbegin[i]       points to the first element in column i.
     Lbegin_p[k]     points to the first element in column p[k].

     Lindex[Lnz+m..2*(Lnz+m)-1],
     Lvalue[Lnz+m..2*(Lnz+m)-1]

       hold L rowwise without the unit diagonal. Each row is terminated
       by index -1. Column indices are such that column i holds the elimination
       factors from the pivot step in which row i was pivot row.

     Ltbegin[i]      points to the first element in row i.
     Ltbegin_p[k]    points to the first element in row p[k].

   The matrices R^k are stored as follows:
   ---------------------------------------

     Lindex[Rbegin[k]..Rbegin[k+1]-1],
     Lvalue[Rbegin[k]..Rbegin[k+1]-1]

       hold the nontrivial column of R^k without the unit diagonal.
       Row indices are row indices of B. Rbegin[0] is one past the last
       element of the L storage.

     eta_row[k]

       holds the row index of the diagonal element of the nontrivial column
       of R^k. These are row indices of B.

   The matrix U is stored as follows:
   ----------------------------------

     Uindex[1..Unz+m],
     Uvalue[1..Unz+m]

       hold U columnwise without the pivot elements. Each column is terminated
       by index -1. Row indices are row indices of B. Updates will introduce
       gaps into the data structure.

       Uindex[0] stores index -1. All empty columns (ie. column in U with no
       off-diagonal elements) are stored in Uindex[0]. For each empty column
       the length of the data structure decreases by 1.

     Ubegin[i]       points to the first element in column qmap[i].

     Windex[..],
     Wvalue[..]

       hold U rowwise without the pivot elements. Column indices are column
       indices of B. The rows are stored in a dynamic file structure with gaps
       between them and out of order.

     Wbegin[j]       points to the first element in row pmap[j].
     Wend[j]         points to one past the last element in row pmap[j].
     Wflink, Wblink  double linked list of rows in memory order.

     col_pivot[0..m-1],
     row_pivot[0..m-1]

       hold the pivot elements by column and by row index.
*/

/*
 * lu_build_factors() - build data structures for L, R, U and permutations
 *
 * Return:
 *
 *  BASICLU_REALLOCATE  require more memory in L, U, and/or W
 *  BASICLU_OK
 */

lu_int lu_build_factors(struct lu *this)
{
    const lu_int m          = this->m;
    const lu_int rank       = this->rank;
    const lu_int Lmem       = this->Lmem;
    const lu_int Umem       = this->Umem;
    const lu_int Wmem       = this->Wmem;
    const lu_int pad        = this->pad;
    const double stretch    = this->stretch;
    lu_int *pinv            = this->pinv;
    lu_int *qinv            = this->qinv;
    lu_int *pmap            = this->pmap; /* shares memory with pinv */
    lu_int *qmap            = this->qmap; /* shares memory with qinv */
    lu_int *pivotcol        = this->pivotcol;
    lu_int *pivotrow        = this->pivotrow;
    lu_int *Lbegin          = this->Lbegin;
    lu_int *Lbegin_p        = this->Lbegin_p;
    lu_int *Ltbegin         = this->Ltbegin;
    lu_int *Ltbegin_p       = this->Ltbegin_p;
    lu_int *Ubegin          = this->Ubegin;
    lu_int *Rbegin          = this->Rbegin;
    lu_int *Wbegin          = this->Wbegin;
    lu_int *Wend            = this->Wend;
    lu_int *Wflink          = this->Wflink;
    lu_int *Wblink          = this->Wblink;
    double *col_pivot       = this->col_pivot;
    double *row_pivot       = this->row_pivot;
    lu_int *Lindex          = this->Lindex;
    double *Lvalue          = this->Lvalue;
    lu_int *Uindex          = this->Uindex;
    double *Uvalue          = this->Uvalue;
    lu_int *Windex          = this->Windex;
    double *Wvalue          = this->Wvalue;
    lu_int *iwork1          = this->iwork1;

    lu_int i, j, ipivot, jpivot, k, lrank, nz, Lnz, Unz, need, get, put, pos;
    double pivot, min_pivot, max_pivot;
    lu_int status = BASICLU_OK;

    /*
     * So far L is stored columnwise in Lindex, Lvalue and U stored rowwise
     * in Uindex, Uvalue. The factorization has computed rank columns of L
     * and rank rows of U. If rank < m, then the columns which have not been
     * pivotal will be removed from U.
     */
    Lnz = Lbegin_p[rank];
    Lnz -= rank;                /* because each column is terminated by -1 */
    Unz = Ubegin[rank];         /* might be decreased when rank < m */

    /*
     * Calculate memory and reallocate. The rowwise and columnwise storage of
     * L both need space for Lnz nonzeros + m terminators. The same for the
     * columnwise storage of U except that Uindex[0] = -1 is reserved to
     * accomodate pointers to empty columns. In the rowwise storage of U each
     * row with nz nonzeros is padded by stretch*nz + pad elements.
     */
    need = 2*(Lnz+m);
    if (Lmem < need)
    {
        this->addmemL = need-Lmem;
        status = BASICLU_REALLOCATE;
    }
    need = Unz+m+1;
    if (Umem < need)
    {
        this->addmemU = need-Umem;
        status = BASICLU_REALLOCATE;
    }
    need = Unz + stretch*Unz + m*pad;
    if (Wmem < need)
    {
        this->addmemW = need-Wmem;
        status = BASICLU_REALLOCATE;
    }
    if (status != BASICLU_OK)
        return status;

    /* ------------------ */
    /* Build permutations */
    /* ------------------ */

    /*
     * Append columns/rows which have not been pivotal to the end of the
     * pivot sequence. Build pivotrow, pivotcol as inverse of pinv, qinv.
     */
    #ifndef NDEBUG
    for (k = 0; k < m; k++)
        pivotrow[k] = -1;
    for (k = 0; k < m; k++)
        pivotcol[k] = -1;
    #endif

    lrank = rank;
    for (i = 0; i < m; i++)
    {
        if (pinv[i] < 0)
            pinv[i] = lrank++;
        pivotrow[pinv[i]] = i;
    }
    assert(lrank == m);
    lrank = rank;
    for (j = 0; j < m; j++)
    {
        if (qinv[j] < 0)
            qinv[j] = lrank++;
        pivotcol[qinv[j]] = j;
    }
    assert(lrank == m);

    #ifndef NDEBUG
    for (k = 0; k < m; k++)
        assert(pivotrow[k] >= 0);
    for (k = 0; k < m; k++)
        assert(pivotcol[k] >= 0);
    #endif

    /* Dependent columns get unit pivot elements. */
    for (k = rank; k < m; k++)
        col_pivot[pivotcol[k]] = 1.0;

    /* ----------------------- */
    /* Lower triangular factor */
    /* ----------------------- */

    /*
     * L columnwise. If rank < m, then complete with unit columns (no
     * off-diagonals, so nothing to store here).
     */
    put = Lbegin_p[rank];
    for (k = rank; k < m; k++)
    {
        Lindex[put++] = -1;
        Lbegin_p[k+1] = put;
    }
    assert(Lbegin_p[m] == Lnz+m);
    for (i = 0; i < m; i++)
        Lbegin[i] = Lbegin_p[pinv[i]];

    /*
     * L rowwise.
     */
    memset(iwork1, 0, m*sizeof(lu_int)); /* row counts */
    for (get = 0; get < Lnz+m; get++)
    {
        if ((i = Lindex[get]) >= 0)
            iwork1[i]++;
    }
    put = Lnz+m;                /* L rowwise starts here */
    for (k = 0; k < m; k++)
    {
        i = pivotrow[k];
        Ltbegin_p[k] = put;
        Ltbegin[i] = put;
        put += iwork1[i];
        Lindex[put++] = -1;     /* terminate row */
        iwork1[i]= Ltbegin_p[k];
    }
    assert(put == 2*(Lnz+m));
    for (k = 0; k < m; k++)     /* fill rows */
    {
        ipivot = pivotrow[k];
        for (get = Lbegin_p[k]; (i = Lindex[get]) >= 0; get++)
        {
            put = iwork1[i]++;  /* put into row i */
            Lindex[put] = ipivot;
            Lvalue[put] = Lvalue[get];
        }
    }

    #ifndef NDEBUG
    for (i = 0; i < m; i++)
        assert(Lindex[iwork1[i]] == -1);
    #endif
    Rbegin[0] = 2*(Lnz+m);      /* beginning of update etas */

    /* ----------------------- */
    /* Upper triangular factor */
    /* ----------------------- */

    /*
     * U rowwise.
     */
    lu_file_empty(m, Wbegin, Wend, Wflink, Wblink, Wmem);
    memset(iwork1, 0, m*sizeof(lu_int)); /* column counts */
    put = 0;

    /*
     * Use separate loops for full rank and rank deficient factorizations. In
     * the first case no elements are removed from U, so skip the test.
     */
    if (rank == m)
    {
        for (k = 0; k < m; k++)
        {
            jpivot = pivotcol[k];
            Wbegin[jpivot] = put;
            nz = 0;
            for (pos = Ubegin[k]; pos < Ubegin[k+1]; pos++)
            {
                j = Uindex[pos];
                Windex[put] = j;
                Wvalue[put++] = Uvalue[pos];
                iwork1[j]++;
                nz++;
            }
            Wend[jpivot] = put;
            put += stretch*nz + pad;
            lu_list_move(jpivot, 0, Wflink, Wblink, m, NULL);
        }
    }
    else
    {
        Unz = 0;                /* actual number of nonzeros */
        for (k = 0; k < rank; k++)
        {
            jpivot = pivotcol[k];
            Wbegin[jpivot] = put;
            nz = 0;
            for (pos = Ubegin[k]; pos < Ubegin[k+1]; pos++)
            {
                j = Uindex[pos];
                if (qinv[j] < rank)
                {
                    Windex[put] = j;
                    Wvalue[put++] = Uvalue[pos];
                    iwork1[j]++;
                    nz++;
                }
            }
            Wend[jpivot] = put;
            put += stretch*nz + pad;
            lu_list_move(jpivot, 0, Wflink, Wblink, m, NULL);
            Unz += nz;
        }
        for (k = rank; k < m; k++)
        {
            jpivot = pivotcol[k];
            Wbegin[jpivot] = put;
            Wend[jpivot] = put;
            put += pad;
            lu_list_move(jpivot, 0, Wflink, Wblink, m, NULL);
        }
    }
    assert(put <= Wend[m]);
    Wbegin[m] = put;            /* beginning of free space */

    /*
     * U columnwise.
     */
    Uindex[0] = -1;
    put = 1;
    for (k = 0; k < m; k++)     /* set column pointers */
    {
        j = pivotcol[k];
        i = pivotrow[k];
        nz = iwork1[j];
        if (nz == 0)
        {
            Ubegin[i] = 0;      /* empty columns all in position 0 */
        }
        else
        {
            Ubegin[i] = put;
            put += nz;
            Uindex[put++] = -1; /* terminate column */
        }
        iwork1[j] = Ubegin[i];
    }
    Ubegin[m] = put;
    for (k = 0; k < m; k++)     /* fill columns */
    {
        jpivot = pivotcol[k];
        i = pivotrow[k];
        for (pos = Wbegin[jpivot]; pos < Wend[jpivot]; pos++)
        {
            j = Windex[pos];
            put = iwork1[j]++;
            assert(put >= 1);
            Uindex[put] = i;
            Uvalue[put] = Wvalue[pos];
        }
    }

    #ifndef NDEBUG
    for (j = 0; j < m; j++)
        assert(Uindex[iwork1[j]] == -1);
    #endif

    /* -------------------- */
    /* Build pivot sequence */
    /* -------------------- */

    /* Build row-column mappings, overwriting pinv, qinv. */
    for (k = 0; k < m; k++)
    {
        i = pivotrow[k];
        j = pivotcol[k];
        pmap[j] = i;
        qmap[i] = j;
    }

    /* Build pivots by row index. */
    max_pivot = 0.0;
    min_pivot = INFINITY;
    for (i = 0; i < m; i++)
    {
        row_pivot[i] = col_pivot[qmap[i]];
        pivot = fabs(row_pivot[i]);
        max_pivot = fmax(pivot, max_pivot);
        min_pivot = fmin(pivot, min_pivot);
    }

    memcpy(this->p, pivotrow, m*sizeof(lu_int));

    this->min_pivot = min_pivot;
    this->max_pivot = max_pivot;
    this->pivotlen = m;
    this->Lnz = Lnz;
    this->Unz = Unz;
    this->Rnz = 0;
    return status;
}
