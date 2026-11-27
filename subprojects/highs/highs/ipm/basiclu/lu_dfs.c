/*
 * lu_dfs.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 * Depth first search in a graph
 *
 */

#include "ipm/basiclu/lu_internal.h"

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
 * Method: iterative depth first search [1,2]. A node is appended to the
 * output when the search from it has finished. Writing the output from
 * xi[top-1] downwards then gives a topological order of the reached nodes [2].
 *
 * The simulated recursion uses two stacks indexed by depth d:
 *   pstack[d]  node being expanded at depth d,
 *   xi[d]      position in index[] of the next neighbour of pstack[d].
 * The second stack occupies xi[0..depth] while the output grows downward
 * from xi[top-1]. Every node on either stack is newly reached, and the
 * caller guarantees room for all newly reached nodes below top, so the two
 * never overlap.
 *
 * [1] R. Tarjan, "Depth-first search and linear graph algorithms",
 *     SIAM J. Comput. 1 (1972), pp. 146-160.
 * [2] J. Gilbert and T. Peierls, "Sparse partial pivoting in time
 *     proportional to arithmetic operations", SIAM J. Sci. Stat. Comput. 9
 *     (1988), pp. 862-874.
 */
lu_int lu_dfs(
    lu_int i, const lu_int *begin, const lu_int *end, const lu_int *index,
    lu_int top, lu_int *xi, lu_int *pstack, lu_int *marked, const lu_int M)
{
    lu_int depth, node, pos, child;

    if (marked[i] == M)
        return top;

    depth = 0;
    marked[i] = M;
    pstack[0] = i;
    xi[0] = begin[i];

    while (depth >= 0)
    {
        node = pstack[depth];
        pos = xi[depth];

        /* Skip neighbours that have been reached already. The neighbour list
         * ends at end[node], or at a negative entry when end is NULL. */
        if (end)
        {
            while (pos < end[node] && marked[index[pos]] == M)
                pos++;
            child = pos < end[node] ? index[pos] : -1;
        }
        else
        {
            while (index[pos] >= 0 && marked[index[pos]] == M)
                pos++;
            child = index[pos];
        }

        if (child >= 0)
        {
            /* Descend into child; resume node after it later. */
            xi[depth] = pos + 1;
            depth++;
            assert(depth < top);
            marked[child] = M;
            pstack[depth] = child;
            xi[depth] = begin[child];
        }
        else
        {
            /* All neighbours of node done. */
            depth--;
            xi[--top] = node;
        }
    }
    return top;
}
