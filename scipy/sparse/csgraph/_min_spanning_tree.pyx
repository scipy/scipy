# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) 2011

import numpy as np
cimport numpy as np
cimport cython

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph._validation import validate_graph
from scipy.sparse._sputils import is_pydata_spmatrix

np.import_array()

include 'parameters.pxi'

def minimum_spanning_tree(csgraph, overwrite=False):
    r"""
    minimum_spanning_tree(csgraph, overwrite=False)

    Return a minimum spanning tree of an undirected graph

    A minimum spanning tree is a graph consisting of the subset of edges
    which together connect all connected nodes, while minimizing the total
    sum of weights on the edges.  This is computed using the Kruskal algorithm.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : array_like or sparse matrix, 2 dimensions
        The N x N matrix representing an undirected graph over N nodes
        (see notes below).
    overwrite : bool, optional
        If true, then parts of the input graph will be overwritten for
        efficiency. Default is False.

    Returns
    -------
    span_tree : csr matrix
        The N x N compressed-sparse representation of the undirected minimum
        spanning tree over the input (see notes below).

    Notes
    -----
    This routine uses undirected graphs as input and output.  That is, if
    graph[i, j] and graph[j, i] are both zero, then nodes i and j do not
    have an edge connecting them.  If either is nonzero, then the two are
    connected by the minimum nonzero value of the two.

    This routine loses precision when users input a dense matrix.
    Small elements < 1E-8 of the dense matrix are rounded to zero.
    All users should input sparse matrices if possible to avoid it.

    If the graph is not connected, this routine returns the minimum spanning
    forest, i.e. the union of the minimum spanning trees on each connected
    component.

    If multiple valid solutions are possible, output may vary with SciPy and
    Python version.

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

    is_pydata_sparse = is_pydata_spmatrix(csgraph)
    if is_pydata_sparse:
        pydata_sparse_cls = csgraph.__class__
        pydata_sparse_fill_value = csgraph.fill_value
    csgraph = validate_graph(csgraph, True, DTYPE, dense_output=False,
                             copy_if_sparse=not overwrite)
    cdef int N = csgraph.shape[0]

    data = csgraph.data
    indices = csgraph.indices
    indptr = csgraph.indptr

    rank = np.zeros(N, dtype=ITYPE)
    predecessors = np.arange(N, dtype=ITYPE)

    # Stable sort is a necessary but not sufficient operation
    # to get to a canonical representation of solutions.
    i_sort = np.argsort(data, kind='stable').astype(ITYPE)
    row_indices = np.zeros(len(data), dtype=ITYPE)

    _min_spanning_tree(data, indices, indptr, i_sort,
                       row_indices, predecessors, rank)

    sp_tree = csr_matrix((data, indices, indptr), (N, N))
    sp_tree.eliminate_zeros()

    if is_pydata_sparse:
        # The `fill_value` keyword is new in PyData Sparse 0.15.4 (May 2024),
        # remove the `except` once the minimum supported version is >=0.15.4
        try:
            sp_tree = pydata_sparse_cls.from_scipy_sparse(
                sp_tree, fill_value=pydata_sparse_fill_value
            )
        except TypeError:
            sp_tree = pydata_sparse_cls.from_scipy_sparse(sp_tree)
    return sp_tree


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _min_spanning_tree(DTYPE_t[::1] data,
                             const ITYPE_t[::1] col_indices,
                             const ITYPE_t[::1] indptr,
                             const ITYPE_t[::1] i_sort,
                             ITYPE_t[::1] row_indices,
                             ITYPE_t[::1] predecessors,
                             ITYPE_t[::1] rank) noexcept nogil:
    # Work-horse routine for computing minimum spanning tree using
    #  Kruskal's algorithm.  By separating this code here, we get more
    #  efficient indexing.
    cdef unsigned int i, j, V1, V2, R1, R2, n_edges_in_mst, n_verts, n_data
    n_verts = predecessors.shape[0]
    n_data = i_sort.shape[0]

    # Arrange `row_indices` to contain the row index of each value in `data`.
    # Note that the array `col_indices` already contains the column index.
    for i in range(n_verts):
        for j in range(indptr[i], indptr[i + 1]):
            row_indices[j] = i

    # step through the edges from smallest to largest.
    #  V1 and V2 are connected vertices.
    n_edges_in_mst = 0
    i = 0
    while i < n_data and n_edges_in_mst < n_verts - 1:
        j = i_sort[i]
        V1 = row_indices[j]
        V2 = col_indices[j]

        # progress upward to the head node of each subtree
        R1 = V1
        while predecessors[R1] != R1:
            R1 = predecessors[R1]
        R2 = V2
        while predecessors[R2] != R2:
            R2 = predecessors[R2]

        # Compress both paths.
        while predecessors[V1] != R1:
            predecessors[V1] = R1
        while predecessors[V2] != R2:
            predecessors[V2] = R2

        # if the subtrees are different, then we connect them and keep the
        # edge.  Otherwise, we remove the edge: it duplicates one already
        # in the spanning tree.
        if R1 != R2:
            n_edges_in_mst += 1

            # Use approximate (because of path-compression) rank to try
            # to keep balanced trees.
            if rank[R1] > rank[R2]:
                predecessors[R2] = R1
            elif rank[R1] < rank[R2]:
                predecessors[R1] = R2
            else:
                predecessors[R2] = R1
                rank[R1] += 1
        else:
            data[j] = 0

        i += 1

    # We may have stopped early if we found a full-sized MST so zero out the rest
    while i < n_data:
        j = i_sort[i]
        data[j] = 0
        i += 1
