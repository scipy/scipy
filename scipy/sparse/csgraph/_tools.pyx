"""
Tools and utilities for working with compressed sparse graphs
"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) 2012

import numpy as np
cimport numpy as np

from scipy.sparse import csr_matrix, isspmatrix, isspmatrix_csr, isspmatrix_csc
from _validation import validate_graph

cimport cython

from libc.stdlib cimport malloc, free

include 'parameters.pxi'


def reconstruct_path(csgraph, predecessors, directed=True):
    """Construct a tree from a graph and a predecessor list.

    Parameters
    ----------
    csgraph : array_like or sparse matrix
        The N x N matrix representing the directed or undirected graph
        from which the predecessors are drawn.
    predecessors : array_like, one dimension
        The length-N array of indices of predecessors for the tree.  The
        index of the parent of node i is given by predecessors[i].
    directed : bool, optional
        If True (default), then operate on a directed graph: only move from
        point i to point j along paths csgraph[i, j].
        If False, then operate on an undirected graph: the algorithm can
        progress from point i to j along csgraph[i, j] or csgraph[j, i].

    Returns
    -------
    cstree : csr matrix
        The N x N directed compressed-sparse representation of the tree drawn
        from csgraph which is encoded by the predecessor list.
    """
    csgraph = validate_graph(csgraph, directed, dense_output=False)

    N = csgraph.shape[0]

    nnull = (predecessors < 0).sum()

    indices = np.argsort(predecessors)[nnull:]
    pind = predecessors[indices]
    indptr = pind.searchsorted(np.arange(N + 1))

    if directed == True:
        data = csgraph[pind, indices]
    else:
        data1 = csgraph[pind, indices]
        data2 = csgraph[indices, pind]
        data1[data1 == 0] = np.inf
        data2[data2 == 0] = np.inf
        data = np.minimum(data1, data2)

    data = np.asarray(data).ravel()

    return csr_matrix((data, indices, indptr), shape=(N, N))


def construct_dist_matrix(graph,
                          predecessors,
                          directed=True):
    """Construct distance matrix from a predecessor matrix

    Parameters
    ----------
    graph: array_like or sparse
        The N x N matrix representation of a directed or undirected graph.
    predecessors: array_like
        The N x N matrix of predecessors of each node (see Notes below).
    directed: bool, optional
        If True (default), then operate on a directed graph: only move from
        point i to point j along paths csgraph[i, j].
        If False, then operate on an undirected graph: the algorithm can
        progress from point i to j along csgraph[i, j] or csgraph[j, i].

    Returns
    -------
    dist_matrix: ndarray
        The N x N matrix of distances between nodes along the path specified
        by the predecessor matrix.  If no path exists, the distance is zero.

    Notes
    -----
    The predecessor matrix is of the form returned by
    :func:`graph_shortest_path`.  Row i of the predecessor matrix contains
    information on the shortest paths from point i: each entry
    predecessors[i, j] gives the index of the previous node in the path from
    point i to point j.  If no path exists between point i and j, then
    predecessors[i, j] = -9999
    """
    graph = validate_graph(graph, directed, dtype=DTYPE,
                           csr_output=False, copy_if_dense=not directed)
    predecessors = np.asarray(predecessors)

    if predecessors.shape != graph.shape:
        raise ValueError("graph and predecessors "
                         "must have the same shape")

    dist_matrix = np.zeros(graph.shape, dtype=DTYPE)
    _construct_dist_matrix(graph, predecessors, dist_matrix, directed)
    
    return dist_matrix


cdef void _construct_dist_matrix(np.ndarray[DTYPE_t, ndim=2] graph,
                                 np.ndarray[ITYPE_t, ndim=2] pred,
                                 np.ndarray[DTYPE_t, ndim=2] dist,
                                 directed=True):
    # All matrices should be size N x N
    # note that graph will be modified if directed == False
    global NULL_IDX

    cdef int i, j, k1, k2, N
    N = graph.shape[0]

    #------------------------------------------
    # symmetrize matrix if necessary
    if not directed:
        graph[graph == 0] = np.inf
        for i from 0 <= i < N:
            for j from i + 1 <= j < N:
                if graph[j, i] <= graph[i, j]:
                    graph[i, j] = graph[j, i]
                else:
                    graph[j, i] = graph[i, j]
    #------------------------------------------

    for i from 0 <= i < N:
        for j from 0 <= j < N:
            k2 = j
            while k2 != i:
                k1 = pred[i, k2]
                if k1 == NULL_IDX:
                    break
                dist[i, j] += graph[k1, k2]
                k2 = k1
