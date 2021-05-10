/*
 * lu_singletons.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * Build initial triangular factors
 *
 */

#include "lu_internal.h"

static lu_int singleton_cols
(
    const lu_int m,
    const lu_int *Bbegin,       /* B columnwise */
    const lu_int *Bend,
    const lu_int *Bi,
    const double *Bx,
    const lu_int *Btp,          /* B rowwise */
    const lu_int *Bti,
    const double *Btx,
    lu_int *Up,
    lu_int *Ui,
    double *Ux,
    lu_int *Lp,
    lu_int *Li,
    double *Lx,
    double *col_pivot,
    lu_int *pinv,
    lu_int *qinv,
    lu_int *iset,               /* size m workspace */
    lu_int *queue,              /* size m workspace */
    lu_int rank,
    double abstol
);

static lu_int singleton_rows
(
    const lu_int m,
    const lu_int *Bbegin,       /* B columnwise */
    const lu_int *Bend,
    const lu_int *Bi,
    const double *Bx,
    const lu_int *Btp,          /* B rowwise */
    const lu_int *Bti,
    const double *Btx,
    lu_int *Up,
    lu_int *Ui,
    double *Ux,
    lu_int *Lp,
    lu_int *Li,
    double *Lx,
    double *col_pivot,
    lu_int *pinv,
    lu_int *qinv,
    lu_int *iset,               /* size m workspace */
    lu_int *queue,              /* size m workspace */
    lu_int rank,
    double abstol
);


/*
 * lu_singletons()
 *
 * Initialize the data structures which store the LU factors during
 * factorization and eliminate pivots with Markowitz cost zero.
 *
 * During factorization the inverse pivot sequence is recorded in pinv, qinv:
 *
 *   pinv[i] >=  0   if row i was pivot row in stage pinv[i]
 *   pinv[i] == -1   if row i has not been pivot row yet
 *   qinv[j] >=  0   if col j was pivot col in stage qinv[j]
 *   qinv[j] == -1   if col j has not been pivot col yet
 *
 * The lower triangular factor is composed columnwise in Lindex, Lvalue.
 * The upper triangular factor is composed rowwise in Uindex, Uvalue.
 * After rank steps of factorization:
 *
 *   Lbegin_p[rank] is the next unused position in Lindex, Lvalue.
 *
 *   Lindex[Lbegin_p[k]..], Lvalue[Lbegin_p[k]..] for 0 <= k < rank
 *   stores the column of L computed in stage k without the unit diagonal.
 *   The column is terminated by a negative index.
 *
 *   Ubegin[rank] is the next unused position in Uindex, Uvalue.
 *
 *   Uindex[Ubegin[k]..Ubegin[k+1]-1], Uvalue[Ubegin[k]..Ubegin[k+1]-1]
 *   stores the row of U computed in stage k without the pivot element.
 *
 * lu_singletons() does rank >= 0 steps of factorization until no singletons are
 * left. We can either eliminate singleton columns before singleton rows or vice
 * versa. When nzbias >= 0, then eliminate singleton columns first to keep L
 * sparse. Otherwise eliminate singleton rows first. The resulting permutations
 * P, Q (stored in inverse form) make PBQ' of the form
 *
 *          \uuuuuuuuuuuuuuuuuuuuuuu
 *           \u                    u
 *            \u                   u
 *             \u                  u
 *              \u                 u
 *  PBQ' =       \uuuuuuu__________u               singleton columns before
 *                \     |          |               singleton rows
 *                l\    |          |
 *                ll\   |          |
 *                l l\  |   BUMP   |
 *                l  l\ |          |
 *                lllll\|__________|
 *
 *          \
 *          l\
 *          ll\
 *          l l\
 *          l  l\
 *          l   l\       __________
 *  PBQ' =  l    l\uuuuu|          |               singleton rows before
 *          l    l \u  u|          |               singleton columns
 *          l    l  \u u|          |
 *          l    l   \uu|   BUMP   |
 *          l    l    \u|          |
 *          llllll     \|__________|
 *
 * Off-diagonals from singleton columns (u) are stored in U, off-diagonals from
 * singleton rows (l) are stored in L and divided by the diagonal. Diagonals (\)
 * are stored in col_pivot.
 *
 * Do not pivot on elements which are zero or less than abstol in magnitude.
 * When such pivots occur, the row/column remains in the active submatrix and
 * the bump factorization will detect the singularity.
 *
 * Return:
 *
 *  BASICLU_REALLOCATE              less than nnz(B) memory in L, U or W
 *  BASICLU_ERROR_invalid_argument  matrix B is invalid (negative number of
 *                                  entries in column, index out of range,
 *                                  duplicates)
 *  BASICLU_OK
 */

lu_int lu_singletons(
    struct lu *this, const lu_int *Bbegin, const lu_int *Bend, const lu_int *Bi,
    const double *Bx)
{
    const lu_int m      = this->m;
    const lu_int Lmem   = this->Lmem;
    const lu_int Umem   = this->Umem;
    const lu_int Wmem   = this->Wmem;
    const double abstol = this->abstol;
    const lu_int nzbias = this->nzbias;
    lu_int *pinv        = this->pinv;
    lu_int *qinv        = this->qinv;
    lu_int *Lbegin_p    = this->Lbegin_p;
    lu_int *Ubegin      = this->Ubegin;
    double *col_pivot   = this->col_pivot;
    lu_int *Lindex      = this->Lindex;
    double *Lvalue      = this->Lvalue;
    lu_int *Uindex      = this->Uindex;
    double *Uvalue      = this->Uvalue;
    lu_int *iwork1      = this->iwork1;
    lu_int *iwork2      = iwork1 + m;

    lu_int *Btp         = this->Wbegin; /* build B rowwise in W */
    lu_int *Bti         = this->Windex;
    double *Btx         = this->Wvalue;

    lu_int i, j, pos, put, rank, Bnz, ok;
    double tic[2];

    /* -------------------------------- */
    /* Check matrix and build transpose */
    /* -------------------------------- */

    /* Check pointers and count nnz(B). */
    Bnz = 0;
    ok = 1;
    for (j = 0; j < m && ok; j++)
    {
        if (Bend[j] < Bbegin[j])
            ok = 0;
        else
            Bnz += Bend[j] - Bbegin[j];
    }
    if (!ok)
        return BASICLU_ERROR_invalid_argument;

    /* Check if sufficient memory in L, U, W. */
    ok = 1;
    if (Lmem < Bnz) { this->addmemL = Bnz-Lmem; ok = 0; }
    if (Umem < Bnz) { this->addmemU = Bnz-Umem; ok = 0; }
    if (Wmem < Bnz) { this->addmemW = Bnz-Wmem; ok = 0; }
    if (!ok)
        return BASICLU_REALLOCATE;

    /* Count nz per row, check indices. */
    memset(iwork1, 0, m*sizeof(lu_int)); /* row counts */
    ok = 1;
    for (j = 0; j < m && ok; j++)
    {
        for (pos = Bbegin[j]; pos < Bend[j] && ok; pos++)
        {
            i = Bi[pos];
            if (i < 0 || i >= m)
                ok = 0;
            else
                iwork1[i]++;
        }
    }
    if (!ok)
        return BASICLU_ERROR_invalid_argument;

    /* Pack matrix rowwise, check for duplicates. */
    put = 0;
    for (i = 0; i < m; i++)     /* set row pointers */
    {
        Btp[i] = put;
        put += iwork1[i];
        iwork1[i] = Btp[i];
    }
    Btp[m] = put;
    assert(put == Bnz);
    ok = 1;
    for (j = 0; j < m; j++)     /* fill rows */
    {
        for (pos = Bbegin[j]; pos < Bend[j]; pos++)
        {
            i = Bi[pos];
            put = iwork1[i]++;
            Bti[put] = j;
            Btx[put] = Bx [pos];
            if (put > Btp[i] && Bti[put-1] == j)
                ok = 0;
        }
    }
    if (!ok)
        return BASICLU_ERROR_invalid_argument;

    /* ---------------- */
    /* Pivot singletons */
    /* ---------------- */

    /* No pivot rows or pivot columns so far. */
    for (i = 0; i < m; i++)
        pinv[i] = -1;
    for (j = 0; j < m; j++)
        qinv[j] = -1;

    if (nzbias >= 0)            /* put more in U */
    {
        Lbegin_p[0] = Ubegin[0] = rank = 0;

        rank = singleton_cols(m, Bbegin, Bend, Bi, Bx, Btp, Bti, Btx,
                              Ubegin, Uindex, Uvalue, Lbegin_p, Lindex, Lvalue,
                              col_pivot, pinv, qinv, iwork1, iwork2, rank,
                              abstol);

        rank = singleton_rows(m, Bbegin, Bend, Bi, Bx, Btp, Bti, Btx,
                              Ubegin, Uindex, Uvalue, Lbegin_p, Lindex, Lvalue,
                              col_pivot, pinv, qinv, iwork1, iwork2, rank,
                              abstol);
    }
    else                        /* put more in L */
    {
        Lbegin_p[0] = Ubegin[0] = rank = 0;

        rank = singleton_rows(m, Bbegin, Bend, Bi, Bx, Btp, Bti, Btx,
                              Ubegin, Uindex, Uvalue, Lbegin_p, Lindex, Lvalue,
                              col_pivot, pinv, qinv, iwork1, iwork2, rank,
                              abstol);

        rank = singleton_cols(m, Bbegin, Bend, Bi, Bx, Btp, Bti, Btx,
                              Ubegin, Uindex, Uvalue, Lbegin_p, Lindex, Lvalue,
                              col_pivot, pinv, qinv, iwork1, iwork2, rank,
                              abstol);
    }

    /* pinv, qinv were used as nonzero counters. Reset to -1 if not pivoted. */
    for (i = 0; i < m; i++)
        if (pinv[i] < 0)
            pinv[i] = -1;
    for (j = 0; j < m; j++)
        if (qinv[j] < 0)
            qinv[j] = -1;

    this->matrix_nz = Bnz;
    this->rank = rank;
    return BASICLU_OK;
}


/*
 * singleton_cols()
 *
 * The method successively removes singleton cols from an active submatrix.
 * The active submatrix is composed of columns j for which qinv[j] < 0 and
 * rows i for which pinv[i] < 0. When removing a singleton column and its
 * associated row generates new singleton columns, these are appended to a
 * queue. The method stops when the active submatrix has no more singleton
 * columns.
 *
 * For each active column j iset[j] is the XOR of row indices in the column
 * in the active submatrix. For a singleton column, this is its single row
 * index. The technique is due to J. Gilbert and described in [1], ex 3.7.
 *
 * For each eliminated column its associated row is stored in U without the
 * pivot element. The pivot elements are stored in col_pivot. For each
 * eliminated pivot an empty column is appended to L.
 *
 * Pivot elements which are zero or less than abstol, and empty columns in
 * the active submatrix are not eliminated. In these cases the matrix is
 * numerically or structurally singular and the bump factorization handles
 * it. (We want singularities at the end of the pivot sequence.)
 *
 * [1] T. Davis, "Direct methods for sparse linear systems"
 */

static lu_int singleton_cols
(
    const lu_int m,
    const lu_int *Bbegin,       /* B columnwise */
    const lu_int *Bend,
    const lu_int *Bi,
    const double *Bx,
    const lu_int *Btp,          /* B rowwise */
    const lu_int *Bti,
    const double *Btx,
    lu_int *Up,
    lu_int *Ui,
    double *Ux,
    lu_int *Lp,
    lu_int *Li,
    double *Lx,
    double *col_pivot,
    lu_int *pinv,
    lu_int *qinv,
    lu_int *iset,               /* size m workspace */
    lu_int *queue,              /* size m workspace */
    lu_int rank,
    double abstol
)
{
    lu_int i, j, j2, nz, pos, put, end, front, tail, rk = rank;
    double piv;

    /* Build index sets and initialize queue. */
    tail = 0;
    for (j = 0; j < m; j++)
    {
        if (qinv[j] < 0)
        {
            nz = Bend[j] - Bbegin[j];
            i = 0;
            for (pos = Bbegin[j]; pos < Bend[j]; pos++)
                i ^= Bi[pos];   /* put row into set j */
            iset[j] = i;
            qinv[j] = -nz-1;    /* use as nonzero counter */
            if (nz == 1)
                queue[tail++] = j;
        }
    }

    /* Eliminate singleton columns. */
    put = Up [rank];
    for (front = 0; front < tail; front++)
    {
        j = queue[front];
        assert(qinv[j] == -2 || qinv[j] == -1);
        if (qinv[j] == -1)
            continue;           /* empty column in active submatrix */
        i = iset[j];
        assert(i >= 0 && i < m);
        assert(pinv[i] < 0);
        end = Btp[i+1];
        for (pos = Btp[i]; Bti[pos] != j; pos++) /* find pivot */
            assert(pos < end-1);
        piv = Btx[pos];
        if (!piv || fabs(piv) < abstol)
            continue;           /* skip singularity */

        /* Eliminate pivot. */
        qinv[j] = rank;
        pinv[i] = rank;
        for (pos = Btp[i]; pos < end; pos++)
        {
            j2 = Bti[pos];
            if (qinv[j2] < 0)
                /* test is mandatory because the initial active submatrix may
                   not be the entire matrix (rows eliminated before) */
            {
                Ui[put] = j2;
                Ux[put++] = Btx[pos];
                iset[j2] ^= i;  /* remove i from set j2 */
                if (++qinv[j2] == -2)
                    queue[tail++] = j2; /* new singleton */
            }
        }
        Up[rank+1] = put;
        col_pivot[j] = piv;
        rank++;
    }

    /* Put empty columns into L. */
    pos = Lp[rk];
    for ( ; rk < rank; rk++)
    {
        Li[pos++] = -1;
        Lp[rk+1] = pos;
    }
    return rank;
}


/*
 * singleton_rows()
 *
 * Analogeous singleton_cols except that for each singleton row the
 * associated column is stored in L and divided by the pivot element. The
 * pivot element is stored in col_pivot.
 */

static lu_int singleton_rows
(
    const lu_int m,
    const lu_int *Bbegin,       /* B columnwise */
    const lu_int *Bend,
    const lu_int *Bi,
    const double *Bx,
    const lu_int *Btp,          /* B rowwise */
    const lu_int *Bti,
    const double *Btx,
    lu_int *Up,
    lu_int *Ui,
    double *Ux,
    lu_int *Lp,
    lu_int *Li,
    double *Lx,
    double *col_pivot,
    lu_int *pinv,
    lu_int *qinv,
    lu_int *iset,               /* size m workspace */
    lu_int *queue,              /* size m workspace */
    lu_int rank,
    double abstol
)
{
    lu_int i, j, i2, nz, pos, put, end, front, tail, rk = rank;
    double piv;

    /* Build index sets and initialize queue. */
    tail = 0;
    for (i = 0; i < m; i++)
    {
        if (pinv[i] < 0)
        {
            nz = Btp[i+1] - Btp[i];
            j = 0;
            for (pos = Btp[i]; pos < Btp[i+1]; pos++)
                j ^= Bti[pos];  /* put column into set i */
            iset[i] = j;
            pinv[i] = -nz-1;    /* use as nonzero counter */
            if (nz == 1)
                queue[tail++] = i;
        }
    }

    /* Eliminate singleton rows. */
    put = Lp[rank];
    for (front = 0; front < tail; front++)
    {
        i = queue[front];
        assert(pinv[i] == -2 || pinv[i] == -1);
        if (pinv[i] == -1)
            continue;           /* empty column in active submatrix */
        j = iset [i];
        assert(j >= 0 && j < m);
        assert(qinv[j] < 0);
        end = Bend[j];
        for (pos = Bbegin[j]; Bi[pos] != i; pos++) /* find pivot */
            assert(pos < end-1);
        piv = Bx[pos];
        if (!piv || fabs(piv) < abstol)
            continue;           /* skip singularity */

        /* Eliminate pivot. */
        qinv[j] = rank;
        pinv[i] = rank;
        for (pos = Bbegin[j]; pos < end; pos++)
        {
            i2 = Bi[pos];
            if (pinv[i2] < 0)
                /* test is mandatory because the initial active submatrix may
                   not be the entire matrix (columns eliminated before) */
            {
                Li[put] = i2;
                Lx[put++] = Bx[pos] / piv;
                iset[i2] ^= j;  /* remove j from set i2 */
                if (++pinv[i2] == -2)
                    queue[tail++] = i2; /* new singleton */
            }
        }
        Li[put++] = -1;         /* terminate column */
        Lp[rank+1] = put;
        col_pivot[j] = piv;
        rank++;
    }

    /* Put empty rows into U. */
    pos = Up[rk];
    for ( ; rk < rank; rk++)
        Up[rk+1] = pos;

    return rank;
}
