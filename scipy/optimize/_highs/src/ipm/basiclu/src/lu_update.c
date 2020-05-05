/*
 * lu_update.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * Forrest-Tomlin update with reordering
 *
 */

#include "lu_internal.h"
#include "lu_list.h"
#include "lu_file.h"

#define GAP (-1)
#define FLIP(i) (-(i)-1)


/**
 * find()
 *
 * Find position of index j in index[start..end-1].
 * If end < 0, then the search stops at the first nonnegative index.
 * Return end if not found.
 *
 */
static lu_int find(
    const lu_int j, const lu_int *index, lu_int start, const lu_int end)
{
    if (end >= 0)
    {
        while (start < end && index[start] != j)
            start++;
        return start;
    }
    else
    {
        while (index[start] != j && index[start] >= 0)
            start++;
        return index[start] == j ? start : end;
    }
}


/**
 * bfs_path()
 *
 * Find a path from j0 to j0 by a breadth first search. When top < m is
 * returned, then the indices in the path (excluding the final j0) are
 * jlist[top..m-1]. When top == m is returned, then no such path exists.
 *
 * The neighbours of node j are index[begin[j]..end[j]-1]. On entry
 * marked[j] >= 0 for all nodes j. On return some elements of marked
 * are set to zero.
 *
 */
static lu_int bfs_path
(
    const lu_int m,             /* graph has m nodes */
    const lu_int j0,
    const lu_int *begin,
    const lu_int *end,
    const lu_int *index,
    lu_int *jlist,
    lu_int *marked,
    lu_int *queue               /* size m workspace */
)
{
    lu_int j, k, pos, front, tail = 1, top = m, found = 0;

    queue[0] = j0;
    for (front = 0; front < tail && !found; front++)
    {
        j = queue[front];
        for (pos = begin[j]; pos < end[j]; pos++)
        {
            k = index[pos];
            if (k == j0)
            {
                found = 1;
                break;
            }
            if (marked[k] >= 0) /* not in queue yet */
            {
                marked[k] = FLIP(j); /* parent[k] = j  */
                queue[tail++] = k;   /* append to queue */
            }
        }
    }
    if (found)                  /* build path (j0,..,j) */
    {
        while (j != j0)
        {
            jlist[--top] = j;
            j = FLIP(marked[j]); /* go to parent */
            assert(j >= 0);
        }
        jlist[--top] = j0;
    }
    for (pos = 0; pos < tail; pos++)
        marked[queue[pos]] = 0; /* reset */

    return top;
}


/**
 * compress_packed()
 *
 * Compress matrix file to reuse memory gaps. Data line 0 <= i < m begins at
 * position begin[i] and ends before the first slot with index[slot] == GAP.
 * begin[m] points to the beginning of unused space at file end.
 *
 * An unused slot in the file must have index[slot] == GAP. All other slots
 * must have index[slot] > GAP. index[0] must be unused. On return
 * index[1..begin[m]-1] contains the data of nonempty lines. All empty lines
 * begin at slot 0. Each two subsequent lines are separated by one gap.
 *
 */
static lu_int compress_packed(const lu_int m, lu_int *begin, lu_int *index,
                              double *value)
{
    lu_int i, p, get, put, nz = 0;
    const lu_int end = begin[m];

    /* Mark the beginning of each nonempty line. */
    for (i = 0; i < m; i++)
    {
        p = begin[i];
        if (index[p] == GAP)
        {
            begin[i] = 0;
        }
        else
        {
            assert(index[p] > GAP);
            begin[i] = index[p]; /* temporarily store index here */
            index[p] = GAP-i-1;  /* mark beginning of line i */
        }
    }

    /* Compress nonempty lines. */
    assert(index[0] == GAP);
    i = -1;
    put = 1;
    for (get = 1; get < end; get++)
    {
        if (index[get] > GAP)   /* shift entry of line i */
        {
            assert(i >= 0);
            index[put] = index[get];
            value[put++] = value[get];
            nz++;
        }
        else if (index[get] < GAP) /* beginning of line i */
        {
            assert(i == -1);
            i = GAP - index[get] - 1;
            index[put] = begin[i]; /* store back */
            begin[i] = put;
            value[put++] = value[get];
            nz++;
        }
        else if (i >= 0)        /* line i ended at a gap */
        {
            i = -1;
            index[put++] = GAP;
        }
    }
    assert(i == -1);
    begin[m] = put;
    return nz;
}


/**
 * permute()
 *
 * Change row-column mappings for columns jlist[0..nswap]. When row i was
 * mapped to column jlist[n], then it will be mapped to column jlist[n+1].
 * When row i was mapped to column jlist[nswap], then it will be mapped to
 * column jlist[0].
 *
 * This requires to update pmap, qmap and the rowwise and columwise storage
 * of U. It also changes the pivot elements.
 *
 * Note: This is the most ugly part of the update code and looks horribly
 *       inefficient (in particular the list swaps). However, usually nswap
 *       is a small number (2, 3, 4, ..), so we don't need to give too much
 *       attention to it.
 *
 */
static void permute(struct lu *this, const lu_int *jlist, lu_int nswap)
{
    lu_int *pmap        = this->pmap;
    lu_int *qmap        = this->qmap;
    lu_int *Ubegin      = this->Ubegin;
    lu_int *Wbegin      = this->Wbegin;
    lu_int *Wend        = this->Wend;
    lu_int *Wflink      = this->Wflink;
    lu_int *Wblink      = this->Wblink;
    double *col_pivot   = this->col_pivot;
    double *row_pivot   = this->row_pivot;
    lu_int *Uindex      = this->Uindex;
    double *Uvalue      = this->Uvalue;
    lu_int *Windex      = this->Windex;
    double *Wvalue      = this->Wvalue;

    const lu_int j0 = jlist[0];
    const lu_int jn = jlist[nswap];
    const lu_int i0 = pmap[j0];
    const lu_int in = pmap[jn];

    lu_int begin, end, i, inext, j, jprev, n, where;
    double piv;

    assert(nswap >= 1);
    assert(qmap[i0] == j0);
    assert(qmap[in] == jn);
    assert(row_pivot[i0] == 0);
    assert(col_pivot[j0] == 0);

    /* --------------- */
    /* Update row file */
    /* --------------- */

    begin = Wbegin[jn];         /* keep for later */
    end = Wend[jn];
    piv = col_pivot[jn];

    for (n = nswap; n > 0; n--)
    {
        j = jlist[n];
        jprev = jlist[n-1];

        /* 
         * When row i was indexed by jprev in the row file before,
         * then it is indexed by j now.
         */
        Wbegin[j] = Wbegin[jprev];
        Wend[j] = Wend[jprev];
        lu_list_swap(Wflink, Wblink, j, jprev);

        /*
         * That row must have an entry in column j because (jprev,j) is an
         * edge in the augmenting path. This entry becomes a pivot element.
         * If jprev is not the first node in the path, then it has an entry
         * in the row (the old pivot) which becomes an off-diagonal entry now.
         */
        where = find(j, Windex, Wbegin[j], Wend[j]);
        assert(where < Wend[j]);
        if (n > 1)
        {
            assert(jprev != j0);
            Windex[where] = jprev;
            col_pivot[j] = Wvalue[where];
            assert(col_pivot[j]);
            Wvalue[where] = col_pivot[jprev];
        }
        else
        {
            assert(jprev == j0);
            col_pivot[j] = Wvalue[where];
            assert(col_pivot[j]);
            Wend[j]--;
            Windex[where] = Windex[Wend[j]];
            Wvalue[where] = Wvalue[Wend[j]];
        }
        this->min_pivot = fmin(this->min_pivot, fabs(col_pivot[j]));
        this->max_pivot = fmax(this->max_pivot, fabs(col_pivot[j]));
    }

    Wbegin[j0] = begin;
    Wend[j0] = end;
    where = find(j0, Windex, Wbegin[j0], Wend[j0]);
    assert(where < Wend[j0]);
    Windex[where] = jn;
    col_pivot[j0] = Wvalue[where];
    assert(col_pivot[j0]);
    Wvalue[where] = piv;
    this->min_pivot = fmin(this->min_pivot, fabs(col_pivot[j0]));
    this->max_pivot = fmax(this->max_pivot, fabs(col_pivot[j0]));

    /* ------------------ */
    /* Update column file */
    /* ------------------ */

    begin = Ubegin[i0];         /* keep for later */

    for (n = 0; n < nswap; n++)
    {
        i = pmap[jlist[n]];
        inext = pmap[jlist[n+1]];

        /*
         * When column j indexed by inext in the column file before,
         * then it is indexed by i now.
         */
        Ubegin[i] = Ubegin[inext];

        /*
         * That column must have an entry in row i because there is an
         * edge in the augmenting path. This entry becomes a pivot element.
         * There is also an entry in row inext (the old pivot), which now
         * becomes an off-diagonal entry.
         */
        where = find(i, Uindex, Ubegin[i], -1);
        assert(where >= 0);
        Uindex[where] = inext;
        row_pivot[i] = Uvalue[where];
        assert(row_pivot[i]);
        Uvalue[where] = row_pivot[inext];
    }

    Ubegin[in] = begin;
    where = find(in, Uindex, Ubegin[in], -1);
    assert(where >= 0);
    row_pivot[in] = Uvalue[where];
    assert(row_pivot[in]);
    for (end = where; Uindex[end] >= 0; end++) ;
    Uindex[where] = Uindex[end-1];
    Uvalue[where] = Uvalue[end-1];
    Uindex[end-1] = -1;

    /* -------------------------- */
    /* Update row-column mappings */
    /* -------------------------- */

    for (n = nswap; n > 0; n--)
    {
        j = jlist[n];
        i = pmap[jlist[n-1]];
        pmap[j] = i;
        qmap[i] = j;
    }
    pmap[j0] = in;
    qmap[in] = j0;

    #ifndef NDEBUG
    for (n = 0; n <= nswap; n++)
    {
        j = jlist[n];
        i = pmap[j];
        assert(row_pivot[i] == col_pivot[j]);
    }
    #endif
}


/**
 * check_consistency()
 *
 * Search for column file entry that is missing or different in row file,
 * and for row file entry that is missing or different in column file.
 * If such an entry is found, then its column and row are returned in *p_col
 * and *p_row. If no such entry exists, then *p_col and *p_row are negative.
 */
#ifdef DEBUG_EXTRA
static void check_consistency(struct lu *this, lu_int *p_col, lu_int *p_row)
{
    const lu_int m          = this->m;
    const lu_int *pmap      = this->pmap;
    const lu_int *qmap      = this->qmap;
    const lu_int *Ubegin    = this->Ubegin;
    const lu_int *Wbegin    = this->Wbegin;
    const lu_int *Wend      = this->Wend;
    const lu_int *Uindex    = this->Uindex;
    const double *Uvalue    = this->Uvalue;
    const lu_int *Windex    = this->Windex;
    const double *Wvalue    = this->Wvalue;
    lu_int i, ientry, j, jentry, pos, where, found;

    for (i = 0; i < m; i++)     /* process column file entries */
    {
        j = qmap[i];
        for (pos = Ubegin[i]; (ientry = Uindex[pos]) >= 0; pos++)
        {
            jentry = qmap[ientry];
            for (where = Wbegin[jentry];
                 where < Wend[jentry] && Windex[where] != j; where++)
                ;
            found = where < Wend[jentry] && Wvalue[where] == Uvalue[pos];
            if (!found)
            {
                *p_col = j;
                *p_row = ientry;
                return;
            }
        }
    }
    for (j = 0; j < m; j++)     /* process row file entries */
    {
        i = pmap[j];
        for (pos = Wbegin[j]; pos < Wend[j]; pos++)
        {
            jentry = Windex[pos];
            ientry = pmap[jentry];
            for (where = Ubegin[ientry];
                 Uindex[where] >= 0 && Uindex[where] != i; where++)
                ;
            found = Uindex[where] == i && Uvalue[where] == Wvalue[pos];
            if (!found)
            {
                *p_col = jentry;
                *p_row = i;
                return;
            }
        }
    }
    *p_col = -1;
    *p_row = -1;
}
#endif


/**
 * lu_update()
 *
 * Insert spike into U and restore triangularity. If the spiked matrix
 * is permuted triangular, then only permutations are updated. If the
 * spiked matrix is not permuted triangular, then the Forrest-Tomlin
 * update is used and the number of row eta matrices increases by 1.
 *
 * Return:
 *
 *  BASICLU_OK                      update successfully completed
 *  BASICLU_REALLOCATE              require more memory in W
 *  BASICLU_ERROR_singular_update   new pivot element is zero or < abstol
 *
 */
lu_int lu_update(struct lu *this, double xtbl)
{
    const lu_int m          = this->m;
    const lu_int nforrest   = this->nforrest;
    lu_int Unz              = this->Unz;
    const lu_int pad        = this->pad;
    const double stretch    = this->stretch;
    lu_int *pmap            = this->pmap;
    lu_int *qmap            = this->qmap;
    lu_int *pivotcol        = this->pivotcol;
    lu_int *pivotrow        = this->pivotrow;
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
    lu_int *marked          = this->marked;
    lu_int *iwork1          = this->iwork1;
    lu_int *iwork2          = iwork1 + m;
    double *work1           = this->work1;

    lu_int jpivot = this->btran_for_update;
    lu_int ipivot = pmap[jpivot];
    double oldpiv = col_pivot[jpivot];
    lu_int status = BASICLU_OK;

    lu_int i, j, jnext, n, nz, t, put, pos, end, where, room, grow, used,need,M;
    lu_int have_diag, intersect, istriangular, nz_roweta, nz_spike;
    lu_int nreach, *col_reach, *row_reach;
    double spike_diag, newpiv, piverr;
    double tic[2], elapsed;

    assert(nforrest < m);

    /*
     * Note: If the singularity test fails or memory is insufficient, then the
     *       update is aborted and the user may call this routine a second time.
     *       Changes made to data structures in the first call must not violate
     *       the logic in the second call.
     */

    /* ------- */
    /* Prepare */
    /* ------- */

    /* if present, move diagonal element to end of spike */
    spike_diag = 0.0;
    have_diag = 0;
    put = Ubegin[m];
    for (pos = put; (i = Uindex[pos]) >= 0; pos++)
    {
        if (i != ipivot)
        {
            Uindex[put] = i;
            Uvalue[put++] = Uvalue[pos];
        }
        else
        {
            spike_diag = Uvalue[pos];
            have_diag = 1;
        }
    }
    if (have_diag)
    {
        Uindex[put] = ipivot;
        Uvalue[put] = spike_diag;
    }
    nz_spike = put - Ubegin[m]; /* nz excluding diagonal */

    nz_roweta = Rbegin[nforrest+1] - Rbegin[nforrest];

    /* ------------- */
    /* Compute pivot */
    /* ------------- */

    /*
     * newpiv is the diagonal element in the spike column after the
     * Forrest-Tomlin update has been applied. It can be computed as
     *
     *    newpiv = spike_diag - dot(spike,row eta)                (1)
     * or
     *    newpiv = xtbl * oldpiv,                                 (2)
     *
     * where spike_diag is the diaognal element in the spike column
     * before the Forrest-Tomlin update and oldpiv was the pivot element
     * in column jpivot before inserting the spike. This routine uses
     * newpiv from (1) and reports the difference to (2) to the user
     * to monitor numerical stability.
     *
     * While computing (1), count intersection of patterns of spike and
     * row eta.
     */

    /* scatter row eta into work1 and mark positions */
    M = ++this->marker;
    for (pos = Rbegin[nforrest]; pos < Rbegin[nforrest+1]; pos++)
    {
        i = Lindex[pos];
        marked[i] = M;
        work1[i] = Lvalue[pos];
    }

    /* compute newpiv and count intersection */
    newpiv = spike_diag;
    intersect = 0;
    for (pos = Ubegin[m]; pos < Ubegin[m] + nz_spike; pos++)
    {
        i = Uindex[pos];
        assert(i != ipivot);
        if (marked[i] == M)
        {
            newpiv -= Uvalue[pos] * work1[i];
            intersect++;
        }
    }

    /* singularity test */
    if (newpiv == 0 || fabs(newpiv) < this->abstol)
    {
        status = BASICLU_ERROR_singular_update;
        return status;
    }

    /* stability measure */
    piverr = fabs(newpiv - xtbl*oldpiv);

    /* ------------ */
    /* Insert spike */
    /* ------------ */

    /* calculate bound on file growth */
    grow = 0;
    for (pos = Ubegin[m]; pos < Ubegin[m] + nz_spike; pos++)
    {
        i = Uindex[pos];
        assert(i != ipivot);
        j = qmap[i];
        jnext = Wflink[j];
        if (Wend[j] == Wbegin[jnext])
        {
            nz = Wend[j] - Wbegin[j];
            grow += nz+1;                 /* row including spike entry */
            grow += stretch*(nz+1) + pad; /* extra room */
        }
    }

    /* reallocate if necessary */
    room = Wend[m] - Wbegin[m];
    if (grow > room)
    {
        this->addmemW = grow-room;
        status = BASICLU_REALLOCATE;
        return status;
    }

    /* remove column jpivot from row file */
    nz = 0;
    for (pos = Ubegin[ipivot]; (i = Uindex[pos]) >= 0; pos++)
    {
        j = qmap[i];
        end = Wend[j]--;
        where = find(jpivot, Windex, Wbegin[j], end);
        assert(where < end);
        Windex[where] = Windex[end-1];
        Wvalue[where] = Wvalue[end-1];
        nz++;
    }
    Unz -= nz;

    /* erase column jpivot in column file */
    for (pos = Ubegin[ipivot]; Uindex[pos] >= 0; pos++)
    {
        Uindex[pos] = GAP;
    }

    /* set column pointers to spike, chop off diagonal */
    Ubegin[ipivot] = Ubegin[m];
    Ubegin[m] += nz_spike;
    Uindex[Ubegin[m]++] = GAP;

    /* insert spike into row file */
    for (pos = Ubegin[ipivot]; (i = Uindex[pos]) >= 0; pos++)
    {
        j = qmap[i];
        jnext = Wflink[j];
        if (Wend[j] == Wbegin[jnext])
        {
            nz = Wend[j] - Wbegin[j];
            room = 1 + stretch*(nz+1) + pad;
            lu_file_reappend(j, m, Wbegin, Wend, Wflink, Wblink, Windex,
                             Wvalue, room);
        }
        end = Wend[j]++;
        Windex[end] = jpivot;
        Wvalue[end] = Uvalue[pos];
    }
    Unz += nz_spike;

    /* insert diagonal */
    col_pivot[jpivot] = spike_diag;
    row_pivot[ipivot] = spike_diag;

    /* ------------------ */
    /* Test triangularity */
    /* ------------------ */

    if (have_diag)
    {
        /*
         * When the spike has a nonzero diagonal element, then the spiked matrix
         * is (symmetrically) permuted triangular if and only if reach(ipivot)
         * does not intersect with the spike pattern except for ipivot. Since
         * reach(ipivot) \ {ipivot} is the structural pattern of the row eta,
         * the matrix is permuted triangular iff the patterns of the row eta
         * and the spike do not intersect.
         *
         * To update the permutations below, we have to provide reach(ipivot)
         * and the associated column indices in topological order as arrays
         * row_reach[0..nreach-1] and col_reach[0..nreach-1]. Because the
         * pattern of the row eta was computed by a dfs, we obtain row_reach
         * simply by adding ipivot to the front. col_reach can then be obtained
         * through qmap.
         */
        istriangular = intersect == 0;
        if (istriangular)
        {
            this->min_pivot = fmin(this->min_pivot, fabs(newpiv));
            this->max_pivot = fmax(this->max_pivot, fabs(newpiv));

            /* build row_reach and col_reach in topological order */
            nreach = nz_roweta + 1;
            row_reach = iwork1;
            col_reach = iwork2;
            row_reach[0] = ipivot;
            col_reach[0] = jpivot;
            pos = Rbegin[nforrest];
            for (n = 1; n < nreach; n++)
            {
                i = Lindex[pos++];
                row_reach[n] = i;
                col_reach[n] = qmap[i];
            }
            this->nsymperm_total++;
        }
    }
    else
    {
        /*
         * The spike has a zero diagonal element, so the spiked matrix may only
         * be the *un*symmetric permutation of an upper triangular matrix.
         *
         * Part 1:
         *
         * Find an augmenting path in U[pmap,:] starting from jpivot.
         * An augmenting path is a sequence of column indices such that there
         * is an edge from each node to the next, and an edge from the final
         * node back to jpivot.
         *
         * bfs_path computes such a path in path[top..m-1].
         *
         * Because jpivot has no self-edge, the path must have at least two
         * nodes. The path must exist because otherwise the spiked matrix was
         * structurally singular and the singularity test above had failed.
         */
        lu_int *path = iwork1, top;
        lu_int *reach = iwork2, rtop;
        lu_int *pstack = (void *) work1;

        top = bfs_path(m, jpivot, Wbegin, Wend, Windex, path, marked, iwork2);
        assert(top < m-1);
        assert(path[top] == jpivot);

        /*
         * Part 2a:
         *
         * For each path index j (except the final one) mark the nodes in
         * reach(j), where the reach is computed in U[pmap,:] without the path
         * edges. If a path index is contained in the reach of an index that
         * comes before it in the path, then U is not permuted triangular.
         *
         * At the same time assemble the combined reach of all path nodes
         * (except the final one) in U[pmap_new,:], where pmap_new is the
         * column-row mapping after applying the permutation associated with
         * the augmenting path. We only have to replace each index where the
         * dfs starts by the next index in the path. The combined reach is
         * then assembled in topological order in
         *
         *    reach[rtop..m-1].
         */
        istriangular = 1;
        rtop = m;
        M = ++this->marker;
        for (t = top; t < m-1 && istriangular; t++)
        {
            j = path[t];
            jnext = path[t+1];
            where = find(jnext, Windex, Wbegin[j], Wend[j]);
            assert(where < Wend[j]);
            Windex[where] = j;  /* take out for a moment */
            rtop = lu_dfs(j, Wbegin, Wend, Windex, rtop, reach, pstack, marked,
                          M);
            assert(reach[rtop] == j);
            reach[rtop] = jnext;
            Windex[where] = jnext; /* restore */
            istriangular = marked[jnext] != M;
        }

        /*
         * Part 2b:
         *
         * If the matrix looks triangular so far, then also mark the reach of
         * the final path node, which is reach(jpivot) in U[pmap_new,:].
         * U is then permuted triangular iff the combined reach does not
         * intersect the spike pattern except in the final path index.
         */
        if (istriangular)
        {
            j = path[m-1];
            rtop = lu_dfs(j, Wbegin, Wend, Windex, rtop, reach, pstack, marked,
                          M);
            assert(reach[rtop] == j);
            reach[rtop] = jpivot;
            marked[j]--;        /* unmark for a moment */
            for (pos = Ubegin[ipivot]; (i = Uindex[pos]) >= 0; pos++)
            {
                if (marked[qmap[i]] == M) istriangular = 0;
            }
            marked[j]++;        /* restore */
        }

        /*
         * If U is permuted triangular, then permute to zero-free diagonal.
         * Set up row_reach[0..nreach-1] and col_reach[0..nreach-1] for
         * updating the permutations below. The column reach is the combined
         * reach of the path nodes. The row reach is is given through pmap.
         */
        if (istriangular)
        {
            lu_int nswap = m-top-1;
            permute(this, path+top, nswap);
            Unz--;
            assert(reach[rtop] == jpivot);
            col_reach = reach + rtop; /* stored in iwork2 */
            row_reach = iwork1 + rtop;
            nreach = m-rtop;
            for (n = 0; n < nreach; n++)
            {
                row_reach[n] = pmap[col_reach[n]];
            }
        }
    }

    /* --------------------- */
    /* Forrest-Tomlin update */
    /* --------------------- */

    if (!istriangular)
    {
        /* remove row ipivot from column file */
        for (pos = Wbegin[jpivot]; pos < Wend[jpivot]; pos++)
        {
            j = Windex[pos];
            assert(j != jpivot);
            where = -1;
            for (end = Ubegin[pmap[j]]; (i = Uindex[end]) >= 0; end++)
            {
                if (i == ipivot) where = end;
            }
            assert(where >= 0);
            Uindex[where] = Uindex[end-1];
            Uvalue[where] = Uvalue[end-1];
            Uindex[end-1] = -1;
            Unz--;
        }

        /* remove row ipivot from row file */
        Wend[jpivot] = Wbegin[jpivot];

        /* replace pivot */
        col_pivot[jpivot] = newpiv;
        row_pivot[ipivot] = newpiv;
        this->min_pivot = fmin(this->min_pivot, fabs(newpiv));
        this->max_pivot = fmax(this->max_pivot, fabs(newpiv));

        /* drop zeros from row eta; update max entry of row etas */
        nz = 0;
        put = Rbegin[nforrest];
        double max_eta = 0;
        for (pos = put; pos < Rbegin[nforrest+1]; pos++)
        {
            if (Lvalue[pos])
            {
                max_eta = fmax(max_eta, fabs(Lvalue[pos]));
                Lindex[put] = Lindex[pos];
                Lvalue[put++] = Lvalue[pos];
                nz++;
            }
        }
        Rbegin[nforrest+1] = put;
        this->Rnz += nz;
        this->max_eta = fmax(this->max_eta, max_eta);

        /* prepare permutation update */
        nreach = 1;
        row_reach = &ipivot;
        col_reach = &jpivot;
        this->nforrest++;
        this->nforrest_total++;
    }

    /* ------------------- */
    /* Update permutations */
    /* ------------------- */

    if (this->pivotlen + nreach > 2*m)
    {
        lu_garbage_perm(this);
    }

    /* append row indices row_reach[0..nreach-1] to end of pivot sequence */
    put = this->pivotlen;
    for (n = 0; n < nreach; n++) pivotrow[put++] = row_reach[n];

    /* append col indices col_reach[0..nreach-1] to end of pivot sequence */
    put = this->pivotlen;
    for (n = 0; n < nreach; n++) pivotcol[put++] = col_reach[n];

    this->pivotlen += nreach;

    /* -------- */
    /* Clean up */
    /* -------- */

    /* compress U if used memory is shrinked sufficiently */
    used = Ubegin[m];
    if (used-Unz-m > this->compress_thres * used)
    {
        nz = compress_packed(m, Ubegin, Uindex, Uvalue);
        assert(nz == Unz);
    }

    /* compress W if used memory is shrinked suficiently */
    used = Wbegin[m];
    need = Unz + stretch*Unz + m*pad;
    if ((used-need) > this->compress_thres * used)
    {
        nz = lu_file_compress(m, Wbegin, Wend, Wflink, Windex, Wvalue,
                              stretch, pad);
        assert(nz == Unz);
    }

    this->pivot_error = piverr / (1.0 + fabs(newpiv));
    this->Unz = Unz;
    this->btran_for_update = -1;
    this->ftran_for_update = -1;
    this->update_cost_numer += nz_roweta;
    this->nupdate++;
    this->nupdate_total++;

    #ifdef DEBUG_EXTRA
    {
        lu_int col, row;
        check_consistency(this, &col, &row);
        assert(col < 0 && row < 0);
    }
    #endif
    return status;
}
