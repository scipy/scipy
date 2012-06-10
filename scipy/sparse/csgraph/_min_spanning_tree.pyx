# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) 2011

import numpy as np
cimport numpy as np

from scipy.sparse import csr_matrix, isspmatrix_csc, isspmatrix
from scipy.sparse.csgraph._validation import validate_graph

include 'parameters.pxi'

def minimum_spanning_tree(csgraph, overwrite=False):
    r"""
    minimum_spanning_tree(csgraph, overwrite=False)

    Return a minimum spanning tree of an undirected graph

    A minimum spanning tree is a graph consisting of the subset of edges
    which together connect all connected nodes, while minimizing the total
    sum of weights on the edges.  This is computed using the Kruskal algorithm.

    Parameters
    ----------
    csgraph: array_like or sparse matrix, 2 dimensions
        The N x N matrix representing an undirected graph over N nodes
        (see notes below).
    overwrite: bool, optional
        if true, then parts of the input graph will be overwritten for
        efficiency.

    Returns
    -------
    span_tree: csr matrix
        The N x N compressed-sparse representation of the undirected minimum
        spanning tree over the input (see notes below).

    Notes
    -----
    This routine uses undirected graphs as input and output.  That is, if
    graph[i, j] and graph[j, i] are both zero, then nodes i and j do not
    have an edge connecting them.  If either is nonzero, then the two are
    connected by the minimum nonzero value of the two.

    Examples
    --------
    The following example shows the computation of a minimum spanning tree
    over a simple four-component graph::

         input graph             minimum spanning tree

             (0)                         (0)
            /   \                       /
           3     8                     3
          /       \                   /
        (3)---5---(1)               (3)---5---(1)
          \       /                           /
           6     2                           2
            \   /                           /
             (2)                         (2)

    It is easy to see from inspection that the minimum spanning tree involves
    removing the edges with weights 8 and 6.  In compressed sparse
    representation, the solution looks like this:

    >>> from scipy.sparse import csr_matrix
    >>> from scipy.sparse.csgraph import minimum_spanning_tree
    >>> X = csr_matrix([[0, 8, 0, 3],
    ...                 [0, 0, 2, 5],
    ...                 [0, 0, 0, 6],
    ...                 [0, 0, 0, 0]])
    >>> Tcsr = minimum_spanning_tree(X)
    >>> Tcsr.toarray().astype(int)
    array([[0, 0, 0, 3],
           [0, 0, 2, 5],
           [0, 0, 0, 0],
           [0, 0, 0, 0]])
    """
    global NULL_IDX

    csgraph = validate_graph(csgraph, True, DTYPE, dense_output=False,
                             copy_if_sparse=not overwrite)
    cdef int N = csgraph.shape[0]

    data = csgraph.data
    indices = csgraph.indices
    indptr = csgraph.indptr

    components = np.arange(N, dtype=ITYPE)
    predecessors = np.empty(N, dtype=ITYPE)
    predecessors.fill(NULL_IDX)

    i_sort = np.argsort(data).astype(ITYPE)
    row_indices = np.zeros(len(data), dtype=ITYPE)

    _min_spanning_tree(data, indices, indptr, i_sort,
                       row_indices, predecessors, components)

    sp_tree = csr_matrix((data, indices, indptr), (N, N))
    sp_tree.eliminate_zeros()

    return sp_tree


cdef _min_spanning_tree(np.ndarray[DTYPE_t, ndim=1, mode='c'] data,
                        np.ndarray[ITYPE_t, ndim=1, mode='c'] col_indices,
                        np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr,
                        np.ndarray[ITYPE_t, ndim=1, mode='c'] i_sort,
                        np.ndarray[ITYPE_t, ndim=1, mode='c'] row_indices,
                        np.ndarray[ITYPE_t, ndim=1, mode='c'] predecessors,
                        np.ndarray[ITYPE_t, ndim=1, mode='c'] components):
    # Work-horse routine for computing minimum spanning tree using
    #  Kruskal's algorithm.  By separating this code here, we get more
    #  efficient indexing.
    global NULL_IDX
    cdef unsigned int i, j, V1, V2
    cdef DTYPE_t E
    
    # Arrange `row_indices` to contain the row index of each value in `data`.
    # Note that the array `col_indices` already contains the column index.
    for i from 0 <= i < predecessors.shape[0]:
        for j from indptr[i] <= j < indptr[i + 1]:
            row_indices[j] = i
    
    # step through the edges from smallest to largest.
    #  V1 and V2 are the vertices, and E is the edge weight connecting them.
    for i from 0 <= i < i_sort.shape[0]:
        j = i_sort[i]
        V1 = row_indices[j]
        V2 = col_indices[j]
        E = data[j]

        # progress upward to the head node of each subtree
        while predecessors[V1] != NULL_IDX:
            V1 = predecessors[V1]
        while predecessors[V2] != NULL_IDX:
            V2 = predecessors[V2]

        # if the subtrees are different, then we connect them and keep the
        # edge.  Otherwise, we remove the edge: it duplicates one already
        # in the spanning tree.
        if components[V1] != components[V2]:
            predecessors[V2] = V1
        else:
            data[j] = 0
    
