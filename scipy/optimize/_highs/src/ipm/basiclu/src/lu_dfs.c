/*
 * lu_dfs.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * Depth first search in a graph
 *
 */

#include "lu_internal.h"

static lu_int dfs_end(
    lu_int i, const lu_int *begin, const lu_int *end, const lu_int *index,
    lu_int top, lu_int *xi, lu_int *pstack, lu_int *marked, const lu_int M);

static lu_int dfs(
    lu_int i, const lu_int *begin, const lu_int *index,
    lu_int top, lu_int *xi, lu_int *pstack, lu_int *marked, const lu_int M);


/*
 * lu_dfs() - compute reach(i) in a graph by depth first search
 *
 * @begin, @end, @index define the graph. When end is not NULL, then node j
 * has neighbours index[begin[j]..end[j]-1]. When end is NULL, then the
 * neighbour list is terminated by a negative index.
 *
 * On return xi[newtop..top-1] hold reach(i) in topological order; newtop is
 * the function return value. Nodes that were already marked are excluded from
 * the reach.
 *
 * @pstack is size m workspace (m the number of nodes in the graph); the
 * contents of pstack is undefined on entry/return.
 *
 * @marked is size m array. Node j is marked iff marked[j] == M.
 * On return nodes xi[newtop..top-1] are marked.
 *
 * If node i is marked on entry, the function does nothing.
 *
 */

lu_int lu_dfs(
    lu_int i, const lu_int *begin, const lu_int *end, const lu_int *index,
    lu_int top, lu_int *xi, lu_int *pstack, lu_int *marked, const lu_int M)
{
    if (marked[i] == M)
        return top;

    return end ?
        dfs_end(i, begin, end, index, top, xi, pstack, marked, M) :
        dfs(i, begin, index, top, xi, pstack, marked, M);
}


/* ==========================================================================
 * dfs_end
 *
 * adapted from T. Davis, CSPARSE
 * ========================================================================== */

static lu_int dfs_end(
    lu_int i, const lu_int *begin, const lu_int *end, const lu_int *index,
    lu_int top, lu_int *xi, lu_int *pstack, lu_int *marked, const lu_int M)
{
    lu_int inext, done, p, head = 0;
    assert(marked[i] != M);

    xi[0] = i;
    while (head >= 0)
    {
        i = xi[head];
        if (marked[i] != M)     /* node i has not been visited */
        {
            marked[i] = M;
            pstack[head] = begin[i];
        }
        done = 1;
        /* continue dfs at node i */
        for (p = pstack[head]; p < end[i]; p++)
        {
            inext = index[p];
            if (marked[inext] == M)
                continue;       /* skip visited node */
            pstack[head] = p+1;
            xi[++head] = inext; /* start dfs at node inext */
            done = 0;
            break;
        }
        if (done)               /* node i has no unvisited neighbours */
        {
            head--;
            xi[--top] = i;
        }
    }

    return top;
}


/* ==========================================================================
 * dfs
 *
 * adapted from T. Davis, CSPARSE
 * ========================================================================== */

static lu_int dfs(
    lu_int i, const lu_int *begin, const lu_int *index,
    lu_int top, lu_int *xi, lu_int *pstack, lu_int *marked, const lu_int M)
{
    lu_int inext, done, p, head = 0;
    assert(marked[i] != M);

    xi[0] = i;
    while (head >= 0)
    {
        i = xi[head];
        if (marked[i] != M)     /* node i has not been visited */
        {
            marked[i] = M;
            pstack[head] = begin[i];
        }
        done = 1;
        /* continue dfs at node i */
        for (p = pstack[head]; (inext = index[p]) >= 0; p++)
        {
            if (marked[inext] == M)
                continue;       /* skip visited node */
            pstack[head] = p+1;
            xi[++head] = inext; /* start dfs at node inext */
            done = 0;
            break;
        }
        if (done)               /* node i has no unvisited neighbours */
        {
            head--;
            xi[--top] = i;
        }
    }

    return top;
}
