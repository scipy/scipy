import numpy as np

from scipy.sparse import csr_matrix, isspmatrix_csr

cimport cython
cimport numpy as np

include 'parameters.pxi'


class MaximumFlowResult:
    """Represents the result of a maximum flow calculation.

    Attributes
    ----------
    flow_value : int
        The value of the maximum flow.
    residual : csr_matrix
        The residual graph with respect to the maximum flow.
    """

    def __init__(self, flow_value, residual):
        self.flow_value = flow_value
        self.residual = residual

    def __repr__(self):
        return 'MaximumFlowResult with value of %d' % self.flow_value


def maximum_flow(csgraph, source, sink):
    r"""
    maximum_flow(csgraph, source, sink)

    Maximize the flow between two vertices in a graph.

    .. versionadded:: 1.4.0

    Parameters
    ----------
    csgraph : csr_matrix
        The square matrix representing a directed graph whose (i, j)'th entry
        is an integer representing the capacity of the edge between
        vertices i and j.
    source : int
        The source vertex from which the flow flows.
    sink : int
        The sink vertex to which the flow flows.

    Returns
    -------
    res : MaximumFlowResult
        A maximum flow represented by a ``MaximumFlowResult``
        which includes the value of the flow in ``flow_value``,
        and the residual graph in ``residual``.

    Raises
    ------
    TypeError:
        if the input graph is not in CSR format.

    ValueError:
        if the capacity values are not integers, or the source or sink are out
        of bounds.

    Notes
    -----
    This solves the maximum flow problem on a given directed weighted graph:
    A flow associates to every edge a value, also called a flow, less than the
    capacity of the edge, so that for every vertex (apart from the source and
    the sink vertices), the total incoming flow is equal to the total outgoing
    flow. The value of a flow is the sum of the flow of all edges leaving the
    source vertex, and the maximum flow problem consists of finding a flow
    whose value is maximal.

    By the max-flow min-cut theorem, the maximal value of the flow is also the
    total weight of the edges in a minimum cut.

    To solve the problem, we use the Edmonds--Karp algorithm. [1]_ This
    particular implementation strives to exploit sparsity. Its time complexity
    is :math:`O(VE^2)` and its space complexity is :math:`O(E)`.

    The maximum flow problem is usually defined with real valued capacities,
    but we require that all capacities are integral to ensure convergence. When
    dealing with rational capacities, or capacities belonging to
    :math:`x\mathbb{Q}` for some fixed :math:`x \in \mathbb{R}`, it is possible
    to reduce the problem to the integral case by scaling all capacities
    accordingly.

    Solving a maximum-flow problem can be used for example for graph cuts
    optimization in computer vision [3]_.

    References
    ----------
    .. [1] Edmonds, J. and Karp, R. M.
           Theoretical improvements in algorithmic efficiency for network flow
           problems. 1972. Journal of the ACM. 19 (2): pp. 248-264
    .. [2] Cormen, T. H. and Leiserson, C. E. and Rivest, R. L. and Stein C.
           Introduction to Algorithms. Second Edition. 2001. MIT Press.
    .. [3] https://en.wikipedia.org/wiki/Graph_cuts_in_computer_vision

    Examples
    --------
    Perhaps the simplest flow problem is that of a graph of only two vertices
    with an edge from source (0) to sink (1)::

        (0) --5--> (1)

    Here, the maximum flow is simply the capacity of the edge:

    >>> from scipy.sparse import csr_matrix
    >>> from scipy.sparse.csgraph import maximum_flow
    >>> graph = csr_matrix([[0, 5], [0, 0]])
    >>> maximum_flow(graph, 0, 1).flow_value
    5

    If, on the other hand, there is a bottleneck between source and sink, that
    can reduce the maximum flow::

        (0) --5--> (1) --3--> (2)

    >>> graph = csr_matrix([[0, 5, 0], [0, 0, 3], [0, 0, 0]])
    >>> maximum_flow(graph, 0, 2).flow_value
    3

    A less trivial example is given in [2]_, Chapter 26.1:

    >>> graph = csr_matrix([[0, 16, 13,  0,  0,  0],
    ...                     [0, 10,  0, 12,  0,  0],
    ...                     [0,  4,  0,  0, 14,  0],
    ...                     [0,  0,  9,  0,  0, 20],
    ...                     [0,  0,  0,  7,  0,  4],
    ...                     [0,  0,  0,  0,  0,  0]])
    >>> maximum_flow(graph, 0, 5).flow_value
    23

    It is possible to reduce the problem of finding a maximum matching in a
    bipartite graph to a maximum flow problem: Let :math:`G = ((U, V), E)` be a
    bipartite graph. Then, add to the graph a source vertex with edges to every
    vertex in :math:`U` and a sink vertex with edges from every vertex in
    :math:`V`. Finally, give every edge in the resulting graph a capacity of 1.
    Then, a maximum flow in the new graph gives a maximum matching in the
    original graph consisting of the edges in :math:`E` whose flow is positive.

    Assume that the edges are represented by a
    :math:`\lvert U \rvert \times \lvert V \rvert` matrix in CSR format whose
    :math:`(i, j)`'th entry is 1 if there is an edge from :math:`i \in U` to
    :math:`j \in V` and 0 otherwise; that is, the input is of the form required
    by :func:`maximum_bipartite_matching`. Then the CSR representation of the
    graph constructed above contains this matrix as a block. Here's an example:

    >>> graph = csr_matrix([[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 1, 0]])
    >>> print(graph.toarray())
    [[0 1 0 1]
     [1 0 1 0]
     [0 1 1 0]]
    >>> i, j = graph.shape
    >>> n = graph.nnz
    >>> indptr = np.concatenate([[0],
    ...                          graph.indptr + i,
    ...                          np.arange(n + i + 1, n + i + j + 1),
    ...                          [n + i + j]])
    >>> indices = np.concatenate([np.arange(1, i + 1),
    ...                           graph.indices + i + 1,
    ...                           np.repeat(i + j + 1, j)])
    >>> data = np.ones(n + i + j, dtype=int)
    >>>
    >>> graph_flow = csr_matrix((data, indices, indptr))
    >>> print(graph_flow.toarray())
    [[0 1 1 1 0 0 0 0 0]
     [0 0 0 0 0 1 0 1 0]
     [0 0 0 0 1 0 1 0 0]
     [0 0 0 0 0 1 1 0 0]
     [0 0 0 0 0 0 0 0 1]
     [0 0 0 0 0 0 0 0 1]
     [0 0 0 0 0 0 0 0 1]
     [0 0 0 0 0 0 0 0 1]
     [0 0 0 0 0 0 0 0 0]]

    At this point, we can find the maximum flow between the added sink and the
    added source and the desired matching can be obtained by restricting the
    residual graph to the block corresponding to the original graph:

    >>> flow = maximum_flow(graph_flow, 0, i+j+1)
    >>> matching = flow.residual[1:i+1, i+1:i+j+1]
    >>> print(matching.toarray())
    [[0 1 0 0]
     [1 0 0 0]
     [0 0 1 0]]

    This tells us that the first, second, and third vertex in :math:`U` are
    matched with the second, first, and third vertex in :math:`V` respectively.

    While this solves the maximum bipartite matching problem in general, note
    that algorithms specialized to that problem, such as
    :func:`maximum_bipartite_matching`, will generally perform better.

    This approach can also be used to solve various common generalizations of
    the maximum bipartite matching problem. If, for instance, some vertices can
    be matched with more than one other vertex, this may be handled by
    modifying the capacities of the new graph appropriately.

    """
    if not isspmatrix_csr(csgraph):
        raise TypeError("graph must be in CSR format")
    if not issubclass(csgraph.dtype.type, np.integer):
        raise ValueError("graph capacities must be integers")
    elif csgraph.dtype != ITYPE:
        csgraph = csgraph.astype(ITYPE)
    if source == sink:
        raise ValueError("source and sink vertices must differ")
    if csgraph.shape[0] != csgraph.shape[1]:
        raise ValueError("graph must be specified as a square matrix.")
    if source < 0 or source >= csgraph.shape[0]:
        raise ValueError('source value ({}) must be between '.format(source) +
                         '0 and {}'.format(csgraph.shape[0] - 1))
    if sink < 0 or sink >= csgraph.shape[0]:
        raise ValueError('sink value ({}) must be between '.format(sink) +
                         '0 and {}'.format(csgraph.shape[0] - 1))

    # Sorted indices are needed by both the _add_reverse_edges() and
    # the _make_edge_pointers() function.
    if not csgraph.has_sorted_indices:
      csgraph = csgraph.sorted_indices()

    # Our implementation of Edmonds--Karp assumes that edges always exist
    # in both directions, so we start by adding the reversed edges whenever
    # they are missing.
    m = _add_reverse_edges(csgraph)

    # Make edge pointers.
    rev_edge_ptr, tails = _make_edge_pointers(m)

    residual = _edmonds_karp(m.indptr, tails, m.indices,
                             m.data, rev_edge_ptr, source, sink)

    residual_array = np.asarray(residual)
    residual_matrix = csr_matrix((residual_array, m.indices, m.indptr),
                                 shape=m.shape)
    source_flow = residual_array[m.indptr[source]:m.indptr[source + 1]]
    return MaximumFlowResult(source_flow.sum(), residual_matrix)


@cython.boundscheck(False)
@cython.wraparound(False)
def _add_reverse_edges(a):
    """Add reversed edges to all edges in a graph.

    This adds to a given directed weighted graph all edges in the reverse
    direction and give them weight 0, unless they already exist.

    Parameters
    ----------
    a : csr_matrix
        The square matrix in CSR format representing a directed graph

    Returns
    -------
    res : csr_matrix
        A new matrix in CSR format in which the missing edges are represented
        by explicit zeros.

    """

    # Reference arrays of the input matrix.
    cdef int n = a.shape[0]
    cdef int[:] a_data_view = a.data
    cdef int[:] a_indices_view = a.indices
    cdef int[:] a_indptr_view = a.indptr

    # Create the transpose with the intent of using the resulting index
    # arrays for the addition of reverse edges with zero capacity.
    b = csr_matrix(a.transpose())
    cdef int[:] b_indices_view = b.indices
    cdef int[:] b_indptr_view = b.indptr

    # Create arrays for the result matrix with added reverse edges.
    c_data = np.zeros(2 * a.nnz, np.intc)
    cdef int[:] c_data_view = c_data
    c_indices = np.zeros(2 * a.nnz, np.intc)
    cdef int[:] c_indices_view = c_indices
    c_indptr = np.zeros(n + 1, np.intc)
    cdef int[:] c_indptr_view = c_indptr

    cdef int i = 0
    cdef int c_ptr = 0
    cdef int a_ptr, a_end, b_ptr, b_end
    while i != n:
        a_ptr = a_indptr_view[i]
        a_end = a_indptr_view[i + 1]
        b_ptr = b_indptr_view[i]
        b_end = b_indptr_view[i + 1]
        while a_ptr != a_end or b_ptr != b_end:
            if a_ptr != a_end and b_ptr != b_end:
                if a_indices_view[a_ptr] < b_indices_view[b_ptr]:
                    c_data_view[c_ptr] = a_data_view[a_ptr]
                    c_indices_view[c_ptr] = a_indices_view[a_ptr]
                    a_ptr += 1
                elif a_indices_view[a_ptr] > b_indices_view[b_ptr]:
                    c_data_view[c_ptr] = 0
                    c_indices_view[c_ptr] = b_indices_view[b_ptr]
                    b_ptr += 1
                else:
                    c_data_view[c_ptr] = a_data_view[a_ptr]
                    c_indices_view[c_ptr] = a_indices_view[a_ptr]
                    a_ptr += 1
                    b_ptr += 1
            elif a_ptr != a_end:
                c_data_view[c_ptr] = a_data_view[a_ptr]
                c_indices_view[c_ptr] = a_indices_view[a_ptr]
                a_ptr += 1
            elif b_ptr != b_end:
                c_indices_view[c_ptr] = b_indices_view[b_ptr]
                b_ptr += 1
            c_ptr += 1
        i += 1
        c_indptr_view[i] = c_ptr

    return csr_matrix((c_data, c_indices, c_indptr), shape=(n, n))


@cython.boundscheck(False)
@cython.wraparound(False)
def _make_edge_pointers(a):
    """Create for each edge pointers to its reverse and its tail."""
    cdef int n = a.shape[0]
    b_data = np.arange(a.data.shape[0], dtype=ITYPE)
    b_indices = a.indices.copy()
    b_indptr = a.indptr.copy()
    b = csr_matrix(
        (b_data, b_indices, b_indptr), shape=(n, n), dtype=ITYPE)
    b = csr_matrix(b.transpose())
    cdef int[:] b_indices_view = b_indices
    cdef int[:] b_indptr_view = b_indptr

    # Overwrite and reuse b_indices with the set of row indices for
    # each matrix entry, i.e., with the set of tail vertices of each edge.
    cdef int i, j
    for i in range(n):
        for j in range(b_indptr_view[i], b_indptr_view[i + 1]):
            b_indices_view[j] = i

    return b.data, b_indices


@cython.boundscheck(False)
@cython.wraparound(False)
cdef ITYPE_t[:] _edmonds_karp(
        ITYPE_t[:] edge_ptr,
        ITYPE_t[:] tails,
        ITYPE_t[:] heads,
        ITYPE_t[:] capacities,
        ITYPE_t[:] rev_edge_ptr,
        ITYPE_t source,
        ITYPE_t sink):
    """Solves the maximum flow problem using the Edmonds--Karp algorithm.

    This assumes that for every edge in the graph, the edge in the opposite
    direction is also in the graph (possibly with capacity 0).

    Parameters
    ----------
    edge_ptr : memoryview of length :math:`|V| + 1`
        For a given vertex v, the edges whose tail is ``v`` are those between
        ``edge_ptr[v]`` and ``edge_ptr[v + 1] - 1``.
    tails : memoryview of length :math:`|E|`
        For a given edge ``e``, ``tails[e]`` is the tail vertex of ``e``.
    heads : memoryview of length :math:`|E|`
        For a given edge ``e``, ``tails[e]`` is the head vertex of ``e``.
    capacities : memoryview of length :math:`|E|`
        For a given edge ``e``, ``capacities[e]`` is the capacity of ``e``.
    rev_edge_ptr : memoryview of length :math:`|E|`
        For a given edge ``e``, ``rev_edge_ptr[e]`` is the edge obtained by
        reversing ``e``. In particular, ``rev_edge_ptr[rev_edge_ptr[e]] == e``.
    source : int
        The source vertex.
    sink : int
        The sink vertex.

    Returns
    -------
    flow : memoryview of length :math:`|E|`
        The residual graph with respect to a maximum flow.

    """
    cdef ITYPE_t n_verts = edge_ptr.shape[0] - 1
    cdef ITYPE_t n_edges = capacities.shape[0]
    cdef ITYPE_t ITYPE_MAX = np.iinfo(ITYPE).max

    # Our result array will keep track of the flow along each edge
    cdef ITYPE_t[:] flow = np.zeros(n_edges, dtype=ITYPE)

    # Create a circular queue for breadth-first search. Elements are
    # popped dequeued at index start and queued at index end.
    cdef ITYPE_t[:] q = np.empty(n_verts, dtype=ITYPE)
    cdef ITYPE_t start, end

    # Create an array indexing predecessor edges
    cdef ITYPE_t[:] pred_edge = np.empty(n_verts, dtype=ITYPE)

    cdef bint path_found
    cdef ITYPE_t cur, df, t, e, edge, k

    # While augmenting paths from source to sink exist
    while True:
        for k in range(n_verts):
            pred_edge[k] = -1
        # Reset queue to consist only of source
        q[0] = source
        start = 0
        end = 1
        # While we have not found a path, and queue is not empty
        path_found = False
        while start != end and not path_found:
            # Pop queue
            cur = q[start]
            start += 1
            # Loop over all edges from the current vertex
            for e in range(edge_ptr[cur], edge_ptr[cur + 1]):
                t = heads[e]
                if pred_edge[t] == -1 and t != source and\
                        capacities[e] > flow[e]:
                    pred_edge[t] = e
                    if t == sink:
                        path_found = True
                        break
                    # Push to queue
                    q[end] = t
                    end += 1
        # Did we find an augmenting path?
        if path_found:
            df = ITYPE_MAX
            # Follow the path back from sink to source to find how
            # much flow can be pushed along the path.
            t = sink
            while t != source:
                e = pred_edge[t]
                df = min(df, capacities[e] - flow[e])
                t = tails[e]
            # Repeat the process, going from sink to source, but this
            # time push the flow that we found above.
            t = sink
            while t != source:
                e = pred_edge[t]
                flow[e] += df
                flow[rev_edge_ptr[e]] -= df
                t = tails[e]
        else:
            # If no augmenting path could be found, we're done.
            break
    return flow
