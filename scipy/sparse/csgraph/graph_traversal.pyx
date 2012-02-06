"""
Routines for traversing graphs in compressed sparse format
"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) 2012

import numpy as np
cimport numpy as np

from scipy.sparse import csr_matrix, isspmatrix, isspmatrix_csr, isspmatrix_csc

cimport cython

from libc.stdlib cimport malloc, free

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t

# NULL_IDX is the index used in predecessor matrices to store a non-path
cdef ITYPE_t NULL_IDX = -9999


def cs_graph_breadth_first(csgraph, i_start,
                           directed=True, return_predecessors=True):
    """Return the breadth-first ordering starting with specified node.

    Parameters
    ----------
    csgraph: array-like or sparse matrix, shape=(N, N)
        compressed sparse graph.  Will be converted to csr format for
        the calculation.
    i_start: integer
        index of starting mode
    directed: bool (default=True)
        if True, then operate on a directed graph: only
        move from point i to point j along paths csgraph[i, j]
        if False, then find the shortest path on an undirected graph: the
        algorithm can progress from point i to j along csgraph[i, j] or
        csgraph[j, i]
    return_predecessors: bool (default=True)
        If True, return the predecesor array

    Returns
    -------
    node_array: np.ndarray, int, shape=(N_nodes,)
        breadth-first list of nodes, starting with specified node
    predecessors: np.ndarray, int, shape=(N_nodes,)
        returned only if return_predecessors == True
        list of predecessors of each node in a breadth-first tree    
    """
    global NULL_IDX

    # if csc matrix and the graph is nondirected, then we can convert to
    # csr using a transpose.
    if (not directed) and isspmatrix_csc(csgraph):
        csgraph = csgraph.T
    elif isspmatrix(csgraph):
        csgraph = csgraph.tocsr()
    else:
        csgraph = csr_matrix(csgraph)

    cdef int N = csgraph.shape[0]
    if csgraph.shape[1] != N:
        raise ValueError("csgraph must be a square matrix")

    cdef np.ndarray node_list = np.empty(N, dtype=ITYPE)
    cdef np.ndarray predecessors = np.empty(N, dtype=ITYPE)
    node_list.fill(NULL_IDX)
    predecessors.fill(NULL_IDX)

    if directed:
        length = _breadth_first(i_start,
                                csgraph.indices, csgraph.indptr,
                                node_list, predecessors)
    else:
        csgraph_T = csgraph.T.tocsr()
        length = _breadth_first_undirected(i_start,
                                           csgraph.indices, csgraph.indptr,
                                           csgraph_T.indices, csgraph_T.indptr,
                                           node_list, predecessors)

    if return_predecessors:
        return node_list[:length], predecessors
    else:
        return node_list[:length]
    

cdef unsigned int _breadth_first(
                           unsigned int head_node,
                           np.ndarray[ITYPE_t, ndim=1, mode='c'] indices,
                           np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr,
                           np.ndarray[ITYPE_t, ndim=1, mode='c'] node_list,
                           np.ndarray[ITYPE_t, ndim=1, mode='c'] predecessors):
    global NULL_IDX

    cdef unsigned int i, pnode, cnode
    cdef unsigned int i_nl, i_nl_end
    cdef unsigned int N = node_list.shape[0]

    node_list[0] = head_node
    i_nl = 0
    i_nl_end = 1

    while i_nl < i_nl_end:
        pnode = node_list[i_nl]

        for i from indptr[pnode] <= i < indptr[pnode + 1]:
            cnode = indices[i]
            if (cnode == head_node):
                continue
            elif (predecessors[cnode] == NULL_IDX):
                node_list[i_nl_end] = cnode
                predecessors[cnode] = pnode
                i_nl_end += 1

        i_nl += 1

    return i_nl
    

cdef unsigned int _breadth_first_undirected(
                           unsigned int head_node,
                           np.ndarray[ITYPE_t, ndim=1, mode='c'] indices1,
                           np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr1,
                           np.ndarray[ITYPE_t, ndim=1, mode='c'] indices2,
                           np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr2,
                           np.ndarray[ITYPE_t, ndim=1, mode='c'] node_list,
                           np.ndarray[ITYPE_t, ndim=1, mode='c'] predecessors):
    global NULL_IDX

    cdef unsigned int i, pnode, cnode
    cdef unsigned int i_nl, i_nl_end
    cdef unsigned int N = node_list.shape[0]

    node_list[0] = head_node
    i_nl = 0
    i_nl_end = 1

    while i_nl < i_nl_end:
        pnode = node_list[i_nl]

        for i from indptr1[pnode] <= i < indptr1[pnode + 1]:
            cnode = indices1[i]
            if (cnode == head_node):
                continue
            elif (predecessors[cnode] == NULL_IDX):
                node_list[i_nl_end] = cnode
                predecessors[cnode] = pnode
                i_nl_end += 1

        for i from indptr2[pnode] <= i < indptr2[pnode + 1]:
            cnode = indices2[i]
            if (cnode == head_node):
                continue
            elif (predecessors[cnode] == NULL_IDX):
                node_list[i_nl_end] = cnode
                predecessors[cnode] = pnode
                i_nl_end += 1

        i_nl += 1

    return i_nl
