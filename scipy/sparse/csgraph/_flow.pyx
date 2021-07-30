# cython: wraparound=False, boundscheck=False

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

class MinCostFlowResult:

    def __init__(self, flow_value, residual):
        self.flow_value = flow_value
        self.residual = residual

    def __repr__(self):
        return 'MinCostFlowResult with value of %d' % self.flow_value


def maximum_flow(csgraph, source, sink, *, method='dinic'):
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
    method: {'edmonds_karp', 'dinic'}, optional
        The method/algorithm to be used for computing the maximum flow.
        Following methods are supported,

            * 'edmonds_karp': Edmonds Karp algorithm in [1]_.
            * 'dinic': Dinic's algorithm in [4]_.

        Default is 'dinic'.

        .. versionadded:: 1.8.0

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

    To solve the problem, we provide Edmonds--Karp [1]_ and Dinic's algorithm
    [4]_. The implementation of both algorithms strive to exploit sparsity.
    The time complexity of the former :math:`O(|V|\,|E|^2)` and its space
    complexity is :math:`O(|E|)`. The latter achieves its performance by
    building level graphs and finding blocking flows in them. Its time
    complexity is :math:`O(|V|^2\,|E|)` and its space complexity is
    :math:`O(|E|)`.

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
    .. [4] Dinic, Efim A.
           Algorithm for solution of a problem of maximum flow in networks with
           power estimation. In Soviet Math. Doklady, vol. 11, pp. 1277-1280.
           1970.

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
    >>> maximum_flow(graph, 0, 1, method='edmonds_karp').flow_value
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

    >>> flow = maximum_flow(graph_flow, 0, i+j+1, method='dinic')
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

    # Our maximum flow solvers assume that edges always exist
    # in both directions, so we start by adding the reversed edges whenever
    # they are missing.
    m = _add_reverse_edges(csgraph)
    rev_edge_ptr = _make_edge_pointers(m)
    if method == 'edmonds_karp':
        tails = _make_tails(m)
        residual = _edmonds_karp(m.indptr, tails, m.indices,
                                 m.data, rev_edge_ptr, source, sink)
    elif method == 'dinic':
        residual = _dinic(m.indptr, m.indices, m.data, rev_edge_ptr,
                          source, sink)
    else:
        raise ValueError('{} method is not supported yet.'.format(method))
    residual_array = np.asarray(residual)
    residual_matrix = csr_matrix((residual_array, m.indices, m.indptr),
                                 shape=m.shape)
    source_flow = residual_array[m.indptr[source]:m.indptr[source + 1]]
    return MaximumFlowResult(source_flow.sum(), residual_matrix)


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
    cdef ITYPE_t n = a.shape[0]
    cdef ITYPE_t[:] a_data_view = a.data
    cdef ITYPE_t[:] a_indices_view = a.indices
    cdef ITYPE_t[:] a_indptr_view = a.indptr

    # Create the transpose with the intent of using the resulting index
    # arrays for the addition of reverse edges with zero capacity. In
    # particular, we do not actually use the values in the transpose;
    # only the fact that the indices exist.
    at = csr_matrix(a.transpose())
    cdef ITYPE_t[:] at_indices_view = at.indices
    cdef ITYPE_t[:] at_indptr_view = at.indptr

    # Create arrays for the result matrix with added reverse edges. We
    # allocate twice the number of non-zeros in `a` for the data, which
    # will always be enough. It might be too many entries in case `a` has
    # some reverse edges already; in that case, over-allocating is not
    # a problem since csr_matrix implicitly truncates elements of data
    # and indices that go beyond the indices given by indptr.
    res_data = np.zeros(2 * a.nnz, ITYPE)
    cdef ITYPE_t[:] res_data_view = res_data
    res_indices = np.zeros(2 * a.nnz, ITYPE)
    cdef ITYPE_t[:] res_indices_view = res_indices
    res_indptr = np.zeros(n + 1, ITYPE)
    cdef ITYPE_t[:] res_indptr_view = res_indptr

    cdef ITYPE_t i = 0
    cdef ITYPE_t res_ptr = 0
    cdef ITYPE_t a_ptr, a_end, at_ptr, at_end
    cdef bint move_a, move_at
    # Loop over all rows
    while i != n:
        # For each row, to ensure that the resulting matrix has
        # sorted indices, we loop over the i'th rows in a and a.T
        # simultaneously, bumping the pointer in one matrix only
        # if that wouldn't break the sorting.
        a_ptr, a_end = a_indptr_view[i], a_indptr_view[i + 1]
        at_ptr, at_end = at_indptr_view[i], at_indptr_view[i + 1]
        while a_ptr != a_end or at_ptr != at_end:
            move_a = a_ptr != a_end \
                and (at_ptr == at_end
                     or a_indices_view[a_ptr] <= at_indices_view[at_ptr])
            move_at = at_ptr != at_end \
                and (a_ptr == a_end
                     or at_indices_view[at_ptr] <= a_indices_view[a_ptr])
            if move_a:
                # Note that it's possible that we move both pointers at once.
                # In that case, we explicitly want the value from the original
                # matrix.
                res_indices_view[res_ptr] = a_indices_view[a_ptr]
                res_data_view[res_ptr] = a_data_view[a_ptr]
                a_ptr += 1
            if move_at:
                res_indices_view[res_ptr] = at_indices_view[at_ptr]
                at_ptr += 1
            res_ptr += 1
        i += 1
        res_indptr_view[i] = res_ptr
    return csr_matrix((res_data, res_indices, res_indptr), shape=(n, n))


def _make_edge_pointers(a):
    """Create for each edge pointers to its reverse."""
    cdef int n = a.shape[0]
    b_data = np.arange(a.data.shape[0], dtype=ITYPE)
    b = csr_matrix(
        (b_data, a.indices, a.indptr), shape=(n, n), dtype=ITYPE)
    b = csr_matrix(b.transpose())
    return b.data


def _make_tails(a):
    """Create for each edge pointers to its tail."""
    cdef int n = a.shape[0]
    cdef ITYPE_t[:] tails = np.empty(a.data.shape[0], dtype=ITYPE)
    cdef ITYPE_t[:] a_indptr_view = a.indptr
    cdef ITYPE_t i, j
    for i in range(n):
        for j in range(a_indptr_view[i], a_indptr_view[i + 1]):
            tails[j] = i
    return tails


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
        For a given edge ``e``, ``heads[e]`` is the head vertex of ``e``.
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

cdef bint _build_level_graph(
        const ITYPE_t[:] edge_ptr,  # IN
        const ITYPE_t source,  # IN
        const ITYPE_t sink,  # IN
        const ITYPE_t[:] capacities,  # IN
        const ITYPE_t[:] heads,  # IN
        ITYPE_t[:] levels,  # IN/OUT
        ITYPE_t[:] q,  # IN/OUT
        ) nogil:
    """Builds layered graph from input graph using breadth first search.

    Parameters
    ----------
    edge_ptr : memoryview of length :math:`|V| + 1`
        For a given vertex ``v``, the edges whose tail is ``v`` are
        those between ``edge_ptr[v]`` and ``edge_ptr[v + 1] - 1``.
    source : int
        The source vertex.
    sink : int
        The sink vertex.
    capacities : memoryview of length :math:`|E|`
        For a given edge ``e``, ``capacities[e]`` is the capacity of ``e``.
    heads : memoryview of length :math:`|E|`
        For a given edge ``e``, ``heads[e]`` is the head vertex of ``e``.
    levels: memoryview of length :math:`|E|`
        For a given vertex ``v``, ``levels[v]`` is the level of ``v`` in
        the layered graph of input graph.
    q : memoryview of length :math:`|E|`
        Queue to be used in breadth first search. Passed to avoid repeated
        queue creation inside this function.

    Returns
    -------
    bool:
        ``True`` if the layered graph creation was successful,
        otherwise ``False``.

    """
    cdef ITYPE_t n_verts = edge_ptr.shape[0] - 1

    cdef ITYPE_t cur, start, end, dst_vertex, e

    q[0] = source
    start = 0
    end = 1
    levels[source] = 0

    while start != end:
        cur = q[start]
        start += 1
        if cur == sink:
            return 1
        for e in range(edge_ptr[cur], edge_ptr[cur + 1]):
            dst_vertex = heads[e]
            if capacities[e] > 0 and levels[dst_vertex] == -1:
                levels[dst_vertex] = levels[cur] + 1
                q[end] = dst_vertex
                end += 1
    return 0

cdef bint _augment_paths(
        const ITYPE_t[:] edge_ptr,  # IN
        const ITYPE_t source,  # IN
        const ITYPE_t sink,  # IN
        const ITYPE_t[:] levels,  # IN
        const ITYPE_t[:] heads,  # IN
        const ITYPE_t[:] rev_edge_ptr,  # IN
        ITYPE_t[:] capacities,  # IN/OUT
        ITYPE_t[:] progress,  # IN
        ITYPE_t[:] flows,  # OUT
        ITYPE_t[:, :] stack
        ) nogil:
    """Finds augmenting paths in layered graph using depth first search.

    Parameters
    ----------
    edge_ptr : memoryview of length :math:`|V| + 1`
        For a given vertex ``v``, the edges whose tail is ``v`` are
        those between ``edge_ptr[v]`` and ``edge_ptr[v + 1] - 1``.
    source : int
        The source vertex.
    sink : int
        The sink vertex.
    levels: memoryview of length :math:`|E|`
        For a given vertex ``v``, ``levels[v]`` is the level of ``v`` in
        the layered graph of input graph.
    heads : memoryview of length :math:`|E|`
        For a given edge ``e``, ``heads[e]`` is the head vertex of ``e``.
    rev_edge_ptr : memoryview of length :math:`|E|`
        For a given edge ``e``, ``rev_edge_ptr[e]`` is the edge obtained by
        reversing ``e``. In particular, ``rev_edge_ptr[rev_edge_ptr[e]] == e``.
    capacities : memoryview of length :math:`|E|`
        For a given edge ``e``, ``capacities[e]`` is the capacity of ``e``.
    progress: memoryview of length :math:`|E|`
        For a given vertex ``v``, ``progress[v]`` is the index of the next
        edge to be visited from ``v``.
    flows : memoryview of length :math:`|E|`
        The residual graph with respect to a maximum flow.
    stack : memoryview of length (:math:`|E|`, 2)
        Stack used during depth-first search.

    Returns
    -------
    bool
        True if and only if an augmenting path was found.

    """
    cdef ITYPE_t top, current, e, dst_vertex, current_flow, flow
    top = 0
    stack[top][0] = source
    stack[top][1] = 2147483647  # Max int

    while True:
        current = stack[top][0]
        flow = stack[top][1]
        e = progress[current]
        dst_vertex = heads[e]
        if (capacities[e] > 0 and
                levels[dst_vertex] == levels[current] + 1):
            current_flow = min(flow, capacities[e])
            if dst_vertex == sink:
                while top > -1:
                    e = progress[stack[top][0]]
                    capacities[e] -= current_flow
                    capacities[rev_edge_ptr[e]] += current_flow
                    flows[e] += current_flow
                    flows[rev_edge_ptr[e]] -= current_flow
                    top -= 1
                return True
            top += 1
            stack[top][0] = dst_vertex
            stack[top][1] = current_flow
        else:
            while progress[current] == edge_ptr[current + 1] - 1:
                top -= 1
                if top < 0: return False  # Did we pop the source?
                current = stack[top][0]
            progress[current] += 1

cdef ITYPE_t[:] _dinic(
        ITYPE_t[:] edge_ptr,
        ITYPE_t[:] heads,
        ITYPE_t[:] capacities,
        ITYPE_t[:] rev_edge_ptr,
        ITYPE_t source,
        ITYPE_t sink):
    """Solves the maximum flow problem using the Dinic's algorithm.

    This assumes that for every edge in the graph, the edge in the opposite
    direction is also in the graph (possibly with capacity 0).

    Parameters
    ----------
    edge_ptr : memoryview of length :math:`|V| + 1`
        For a given vertex ``v``, the edges whose tail is ``v`` are
        those between ``edge_ptr[v]`` and ``edge_ptr[v + 1] - 1``.
    heads : memoryview of length :math:`|E|`
        For a given edge ``e``, ``heads[e]`` is the head vertex of ``e``.
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
    flows : memoryview of length :math:`|E|`
        The residual graph with respect to a maximum flow.

    """
    cdef ITYPE_t n_verts = edge_ptr.shape[0] - 1
    cdef ITYPE_t n_edges = capacities.shape[0]

    cdef ITYPE_t[:] levels = np.empty(n_verts, dtype=ITYPE)
    cdef ITYPE_t[:] progress = np.empty(n_verts, dtype=ITYPE)
    cdef ITYPE_t[:] q = np.empty(n_verts, dtype=ITYPE)
    cdef ITYPE_t[:, :] stack = np.empty((n_verts, 2), dtype=ITYPE)
    cdef ITYPE_t[:] flows = np.zeros(n_edges, dtype=ITYPE)
    cdef ITYPE_t flow
    while True:
        for i in range(n_verts):
            levels[i] = -1
        if not _build_level_graph(edge_ptr, source, sink,
                                  capacities, heads, levels, q):
            break
        for i in range(n_verts):
            progress[i] = edge_ptr[i]
        while _augment_paths(edge_ptr, source, sink,
                             levels, heads, rev_edge_ptr,
                             capacities, progress, flows, stack):
            pass
    return flows

def minimum_cost_flow(csgraph, demand, cost):

    demand = np.asarray(demand)
    cost = np.asarray(cost)
    if not isspmatrix_csr(csgraph):
        raise TypeError("graph must be in CSR format")
    if not issubclass(csgraph.dtype.type, np.integer):
        raise ValueError("graph capacities must be integers")
    elif csgraph.dtype != ITYPE:
        csgraph = csgraph.astype(ITYPE)
    if not issubclass(demand.dtype.type, np.integer):
        raise ValueError("vertex demands must be integers")
    elif demand.dtype != ITYPE:
        demand = demand.astype(ITYPE)
    if not issubclass(cost.dtype.type, np.integer):
        raise ValueError("edge costs must be integers")
    elif cost.dtype != ITYPE:
        cost = cost.astype(ITYPE)

    _network_simplex_checks(csgraph.indptr, demand, cost, csgraph.indptr.shape[0] - 1)
    tails = _make_tails(csgraph)
    flow = _network_simplex(csgraph.indices, tails, csgraph.data,
                            demand, cost, csgraph.indptr.shape[0] - 1)
    flow_array = np.asarray(flow[csgraph.data.shape[0]:])
    flow_matrix = csr_matrix((flow_array, csgraph.indices, csgraph.indptr),
                                 shape=csgraph.shape)
    return MaximumFlowResult(flow_array.sum(), flow_matrix)

def _network_simplex_checks(
        ITYPE_t[:] capacities,
        ITYPE_t[:] demand,
        ITYPE_t[:] cost,
        ITYPE_t n_verts):
    cdef ITYPE_t n_edges = capacities.shape[0]
    cdef ITYPE_t demand_sum

    if n_verts == 0:
        raise ValueError("graph has no vertices")

    # print(n_verts, demand)
    demand_sum = 0
    for i in range(n_verts):
        if demand[i] == 1<<31 - 1 or demand[i] == -1<<31 :
            raise ValueError("vertex %d has infinite demand"%(i))
        demand_sum += demand[i]

    # TODO: Should we add support for selfloop edges?
    for i in range(n_edges):
        if cost[i] == 1<<31 - 1 or cost[i] == -1<<31 :
            raise ValueError("edge %d has infinite cost"%(i))
        if capacities[i] < 0:
            raise ValueError("edge %d has negative capacity"%(i))

    if demand_sum != 0:
        raise ValueError("sum of demands is not zero")

cdef void _initialize_spanning_tree(
    n_verts,
    n_edges,
    faux_inf,
    ITYPE_t[:] demand,
    ITYPE_t[:] edge_flow,
    ITYPE_t[:] vertex_potentials,
    ITYPE_t[:] parent,
    ITYPE_t[:] parent_edge,
    ITYPE_t[:] subtree_size,
    ITYPE_t[:] next_vertex_dft,
    ITYPE_t[:] prev_vertex_dft,
    ITYPE_t[:] last_descendent_dft):

    # edge_flow initialization
    for e in range(n_edges):
        edge_flow[e] = 0
    for v in range(n_verts):
        edge_flow[v + n_edges] = abs(demand[v])

    # vertex_potentials initialization
    for v in range(1, n_verts + 1):
        if demand[v - 1] <= 0:
            vertex_potentials[v] = faux_inf
        else:
            vertex_potentials[v] = -faux_inf

    # parent initialization
    for v in range(1, n_verts):
        parent[v] = 0
    parent[0] = -2

    # parent_edge initialization
    for v in range(1, n_verts + 1):
        parent_edge[v] = n_edges + v - 1

    # subtree_size initialization
    for v in range(1, n_verts + 1):
        subtree_size[v] = 1
    subtree_size[0] = n_verts + 1

    # next_vertex_dft initialization
    for v in range(1, n_verts + 1):
        next_vertex_dft[v] = v
    next_vertex_dft[n_verts] = 0
    next_vertex_dft[0] = 1

    # prev_vertex_dft initialization
    for i in range(1, n_verts + 1):
        prev_vertex_dft[i] = i - 1
    prev_vertex_dft[0] = n_verts

    # last_descendent_dft initialization
    for v in range(1, n_verts + 1):
        last_descendent_dft[v] = v
    last_descendent_dft[0] = n_verts

cdef ITYPE_t _reduced_cost(
    ITYPE_t e,
    ITYPE_t[:] edge_weights,
    ITYPE_t[:] vertex_potentials,
    ITYPE_t[:] edge_sources,
    ITYPE_t[:] edge_targets,
    ITYPE_t[:] edge_flow):
    # print("Reduced cost args: ", e, edge_weights[e], edge_sources[e], vertex_potentials[edge_sources[e]], edge_targets[e], vertex_potentials[edge_targets[e]])
    c = (edge_weights[e]
         - vertex_potentials[edge_sources[e]]
         + vertex_potentials[edge_targets[e]])
    return c if edge_flow[e] == 0 else -c

cdef bint _find_entering_edges(
    ITYPE_t n_edges,
    ITYPE_t[:] edge_weights,
    ITYPE_t[:] vertex_potentials,
    ITYPE_t[:] edge_flow,
    ITYPE_t[:] edge_sources,
    ITYPE_t[:] edge_targets,
    ITYPE_t[:] local_vars,
    bint prev_ret_value,
    ITYPE_t[:] result):
    cdef ITYPE_t l, i = -2, min_r_cost, p, q, c

    if not prev_ret_value:
        if n_edges == 0:
            return False

        local_vars[0] = np.ceil(np.sqrt(n_edges))  # B
        local_vars[1] = (n_edges + local_vars[0] - 1) // local_vars[0]  # M
        local_vars[2] = 0  # m
        # entering edges
        local_vars[3] = 0  # f
    # print(local_vars[0], local_vars[1], local_vars[2], local_vars[3])

    while local_vars[2] < local_vars[1]:
        # Determine the next block of edges.
        l = local_vars[3] + local_vars[0]
        min_r_cost = 1<<31 - 1
        # print("Initial min_r_cost: ", min_r_cost)
        if l <= n_edges:
            for e in range(local_vars[3], l):
                r_cost = _reduced_cost(e, edge_weights,
                                       vertex_potentials,
                                       edge_sources,
                                       edge_targets, edge_flow)
                # print("r_cost: ", r_cost, min_r_cost)
                if r_cost < min_r_cost:
                    i = e
                    min_r_cost = r_cost
        else:
            l -= n_edges
            for e in range(local_vars[3], n_edges):
                r_cost = _reduced_cost(e, edge_weights,
                                       vertex_potentials,
                                       edge_sources,
                                       edge_targets, edge_flow)
                if r_cost < min_r_cost:
                    i = e
                    min_r_cost = r_cost
            for e in range(l):
                r_cost = _reduced_cost(e, edge_weights,
                                       vertex_potentials,
                                       edge_sources,
                                       edge_targets, edge_flow)
                if r_cost < min_r_cost:
                    i = e
                    min_r_cost = r_cost
        local_vars[3] = l
        print("min_r_cost: ", min_r_cost)
        if min_r_cost >= 0:
            # No entering edge found in the current block.
            local_vars[2] += 1
        else:
            # print("i")
            # Entering edge found.
            if edge_flow[i] == 0:
                p = edge_sources[i]
                q = edge_targets[i]
            else:
                p = edge_targets[i]
                q = edge_sources[i]
            result[0] = i
            result[1] = p
            result[2] = q
            local_vars[2] = 0
            return True

    return False

cdef ITYPE_t _find_apex(
    ITYPE_t p,
    ITYPE_t q,
    ITYPE_t[:] subtree_size,
    ITYPE_t[:] parent,
    ITYPE_t n_Ver):
    cdef ITYPE_t size_p, size_q
    print("p, q: ", p, q)
    size_p = subtree_size[p]
    size_q = subtree_size[q]
    print("size_p, size_q: ", size_p, size_q)
    while True:
        while size_p < size_q:
            p = parent[p]
            size_p = subtree_size[p]
        while size_p > size_q:
            q = parent[q]
            size_q = subtree_size[q]
        if size_p == size_q:
            if p != q:
                p = parent[p]
                size_p = subtree_size[p]
                q = parent[q]
                size_q = subtree_size[q]
            else:
                return p

cdef ITYPE_t _trace_path(
    ITYPE_t p,
    ITYPE_t w,
    ITYPE_t[:] parent, 
    ITYPE_t[:] parent_edge,
    ITYPE_t[:] Wn,
    ITYPE_t[:] We,
    ITYPE_t n_verts):
    Wn[0] = p
    idx = 0
    while p != w:
        We[idx] = parent_edge[p]
        p = parent[p]
        Wn[idx + 1] = p
        idx += 1
    return idx

cdef void _find_cycle(
    ITYPE_t i,
    ITYPE_t p,
    ITYPE_t q,
    ITYPE_t[:] subtree_size,
    ITYPE_t[:] parent,
    ITYPE_t[:] parent_edge,
    ITYPE_t n_edges,
    ITYPE_t n_verts,
    ITYPE_t[:] Wn,
    ITYPE_t[:] We,
    ITYPE_t[:] Wne_len):
    cdef ITYPE_t num_edges, num_edges_R
    cdef ITYPE_t offset, Wn_len, We_len
    cdef ITYPE_t[:] WnR, WeR
    w = _find_apex(p, q, subtree_size, parent, n_verts + 1)
    print("w: ", w)
    print("After _find_apex")
    num_edges = _trace_path(p, w, parent, 
                            parent_edge,
                            Wn, We, n_verts)
    print("Wn: ", end="")
    _print_array(Wn)
    print("We: ", end="")
    _print_array(We)
    print("Wne_len: ", num_edges)
    print("After _trace_path")
    Wn_len = num_edges + 1
    We_len = num_edges
    print("Wn_len: ", Wn_len)
    print("We_len: ", We_len)

    WnR = np.empty(n_verts, dtype=ITYPE)
    WeR = np.empty(n_edges, dtype=ITYPE)
    num_edges_R = _trace_path(q, w, parent, 
                              parent_edge,
                              WnR, WeR,
                              n_verts)
    for idx in range(We_len//2):
        tmp = We[idx]
        We[idx] = We[We_len - idx - 1]
        We[We_len - idx - 1] = tmp
    for idx in range(Wn_len//2):
        tmp = Wn[idx]
        Wn[idx] = Wn[Wn_len - idx - 1]
        Wn[Wn_len - idx - 1] = tmp
    if ((We_len == 1 and We[0] != i) or
        We_len != 1):
        We[We_len] = i
        We_len += 1
    for idx in range(num_edges_R):
        Wn[idx + Wn_len] = WnR[idx]
    Wn_len += num_edges_R
    for idx in range(num_edges_R):
        We[idx + We_len] = WeR[idx]
    We_len += num_edges_R
    Wne_len[0] = Wn_len
    Wne_len[1] = We_len

cdef ITYPE_t _residual_capacity(
    ITYPE_t i,
    ITYPE_t p,
    ITYPE_t[:] edge_sources,
    ITYPE_t[:] edge_capacities,
    ITYPE_t[:] edge_flow):
    if edge_sources[i] == p:
        return edge_capacities[i] - edge_flow[i]
    else:
        return edge_flow[i]

cdef void _find_leaving_edge(
    ITYPE_t[:] Wn,
    ITYPE_t[:] We,
    ITYPE_t[:] Wne_len,
    ITYPE_t[:] edge_sources,
    ITYPE_t[:] edge_targets,
    ITYPE_t[:] edge_capacities,
    ITYPE_t[:] edge_flow,
    ITYPE_t[:] ret_values):
    cdef ITYPE_t min_res_cap, res_cap
    cdef ITYPE_t i = -2, p, j = -2, s = -2, t, idx_e, idx_n
    min_res_cap = 1<<31 - 1
    idx_n = Wne_len[0] - 1
    idx_e = Wne_len[1] - 1
    while idx_n >= 0 and idx_e >= 0:
        i = We[idx_e]
        p = Wn[idx_n]
        res_cap = _residual_capacity(i, p,
                                     edge_sources,
                                     edge_capacities,
                                     edge_flow)
        if res_cap < min_res_cap:
            min_res_cap = res_cap
            j = i
            s = p
        idx_n -= 1
        idx_e -= 1
    if edge_sources[j] == s:
        t = edge_targets[j]
    else:
        t = edge_sources[j]
    ret_values[0] = j
    ret_values[1] = s
    ret_values[2] = t

cdef ITYPE_t _index(
    ITYPE_t[:] W,
    ITYPE_t W_len,
    ITYPE_t i):
    cdef ITYPE_t idx
    for idx in range(W_len):
        if W[idx] == i:
            return idx
    return -1

cdef void _augment_flow(
    ITYPE_t[:] Wn,
    ITYPE_t[:] We,
    ITYPE_t[:] Wne_len,
    ITYPE_t f,
    ITYPE_t[:] edge_sources,
    ITYPE_t[:] edge_flow):
    cdef ITYPE_t i, p, idx_e, idx_n
    idx_e = 0
    idx_n = 0
    while (idx_n < Wne_len[0] and
           idx_e < Wne_len[1]):
        i = We[idx_n]
        p = Wn[idx_e]
        if edge_sources[i] == p:
            edge_flow[i] += f
        else:
            edge_flow[i] -= f
        idx_e += 1
        idx_n += 1

cdef void _remove_edge(
    ITYPE_t s,
    ITYPE_t t,
    ITYPE_t[:] subtree_size,
    ITYPE_t[:] prev_vertex_dft,
    ITYPE_t[:] last_descendent_dft,
    ITYPE_t[:] next_vertex_dft,
    ITYPE_t[:] parent,
    ITYPE_t[:] parent_edge,
    ITYPE_t n_verts):
    cdef ITYPE_t size_t, prev_t, last_t
    cdef ITYPE_t next_last_t
    size_t = subtree_size[t]
    prev_t = prev_vertex_dft[t]
    last_t = last_descendent_dft[t]
    next_last_t = next_vertex_dft[last_t]
    # Remove (s, t).
    parent[t] = -2
    parent_edge[t] = -2
    # Remove the subtree rooted at t from the depth-first thread.
    next_vertex_dft[prev_t] = next_last_t
    prev_vertex_dft[next_last_t] = prev_t
    next_vertex_dft[last_t] = t
    prev_vertex_dft[t] = last_t
    # Update the subtree sizes and last descendants of the (old) acenstors
    # of t.
    while s != -2:
        subtree_size[s] -= size_t
        if last_descendent_dft[s] == last_t:
            last_descendent_dft[s] = prev_t
        s = parent[s]
        # print("s: ", s)

cdef void _make_root(
    ITYPE_t n_verts,
    ITYPE_t q,
    ITYPE_t[:] subtree_size,
    ITYPE_t[:] last_descendent_dft,
    ITYPE_t[:] prev_vertex_dft,
    ITYPE_t[:] next_vertex_dft,
    ITYPE_t[:] parent,
    ITYPE_t[:] parent_edge):
    cdef ITYPE_t idx, path_len, p
    cdef ITYPE_t size_p, last_p, last_q, prev_q
    cdef ITYPE_t next_last_q
    cdef ITYPE_t[:] ancestors = np.empty(n_verts, dtype=ITYPE)
    idx = 0
    while q != -2:
        ancestors[idx] = q
        q = parent[q]
        idx += 1
    path_len = idx
    idx = path_len - 1
    while idx > 0:
        p = ancestors[idx]
        q = ancestors[idx - 1]
        size_p = subtree_size[p]
        last_p = last_descendent_dft[p]
        prev_q = prev_vertex_dft[q]
        last_q = last_descendent_dft[q]
        next_last_q = next_vertex_dft[last_q]
        # Make p a child of q.
        parent[p] = q
        parent[q] = -2
        parent_edge[p] = parent_edge[q]
        parent_edge[q] = -2
        subtree_size[p] = size_p - subtree_size[q]
        subtree_size[q] = size_p
        # Remove the subtree rooted at q from the depth-first thread.
        next_vertex_dft[prev_q] = next_last_q
        prev_vertex_dft[next_last_q] = prev_q
        next_vertex_dft[last_q] = q
        prev_vertex_dft[q] = last_q
        if last_p == last_q:
            last_descendent_dft[p] = prev_q
            last_p = prev_q
        # Add the remaining parts of the subtree rooted at p as a subtree
        # of q in the depth-first thread.
        prev_vertex_dft[p] = last_q
        next_vertex_dft[last_q] = p
        next_vertex_dft[last_p] = q
        prev_vertex_dft[q] = last_p
        last_descendent_dft[q] = last_p
        idx -= 1

cdef void _add_edge(
    ITYPE_t i,
    ITYPE_t p,
    ITYPE_t q,
    ITYPE_t[:] last_descendent_dft,
    ITYPE_t[:] next_vertex_dft,
    ITYPE_t[:] prev_vertex_dft,
    ITYPE_t[:] subtree_size,
    ITYPE_t[:] parent,
    ITYPE_t[:] parent_edge,
    ITYPE_t n_verts):
    cdef ITYPE_t last_p, next_last_p
    cdef ITYPE_t size_q, last_q

    last_p = last_descendent_dft[p]
    next_last_p = next_vertex_dft[last_p]
    size_q = subtree_size[q]
    last_q = last_descendent_dft[q]
    # Make q a child of p.
    parent[q] = p
    parent_edge[q] = i
    # Insert the subtree rooted at q into the depth-first thread.
    next_vertex_dft[last_p] = q
    prev_vertex_dft[q] = last_p
    prev_vertex_dft[next_last_p] = last_q
    next_vertex_dft[last_q] = next_last_p
    # Update the subtree sizes and last descendants of the (new) ancestors
    # of q.
    while p != -2:
        subtree_size[p] += size_q
        if last_descendent_dft[p] == last_p:
            last_descendent_dft[p] = last_q
        p = parent[p]

cdef bint _trace_subtree(
    ITYPE_t p,
    ITYPE_t[:] last_descendent_dft,
    ITYPE_t[:] next_vertex_dft,
    ITYPE_t[:] subtree):
    cdef ITYPE_t idx, l

    subtree[0] = p
    idx = 1
    l = last_descendent_dft[p]
    while p != l:
        p = next_vertex_dft[p]
        subtree[idx] = p
        idx += 1
    return idx

cdef void _print_array(ITYPE_t[:] arr):
    cdef ITYPE_t i
    for i in range(arr.shape[0]):
        print(arr[i], end=", ")
    print()

cdef void _update_potentials(
    ITYPE_t i,
    ITYPE_t p,
    ITYPE_t q,
    ITYPE_t n_verts,
    ITYPE_t[:] edge_targets,
    ITYPE_t[:] edge_weights,
    ITYPE_t[:] last_descendent_dft,
    ITYPE_t[:] next_vertex_dft,
    ITYPE_t[:] vertex_potentials):
    cdef ITYPE_t idx, d, path_len
    cdef ITYPE_t[:] subtree = np.empty(n_verts, dtype=ITYPE)
    if q == edge_targets[i]:
        d = vertex_potentials[p] - edge_weights[i] - vertex_potentials[q]
    else:
        d = vertex_potentials[p] + edge_weights[i] - vertex_potentials[q]
    path_len = _trace_subtree(p, last_descendent_dft,
                              next_vertex_dft, subtree)
    for idx in range(path_len):
        q = subtree[idx]
        vertex_potentials[q] += d

cdef ITYPE_t[:] _network_simplex(
        ITYPE_t[:] heads,
        ITYPE_t[:] tails,
        ITYPE_t[:] capacities,
        ITYPE_t[:] demand,
        ITYPE_t[:] cost,
        ITYPE_t n_verts):
    cdef ITYPE_t n_edges = capacities.shape[0]
    cdef ITYPE_t idx, faux_inf, capacities_sum, cost_sum
    cdef ITYPE_t j, s, t, i, tmp, p, q
    cdef ITYPE_t[:] edge_flow, vertex_potentials, Wne_len
    cdef ITYPE_t[:] parent, parent_edge, subtree_size
    cdef ITYPE_t[:] next_vertex_dft, prev_vertex_dft
    cdef ITYPE_t[:] last_descendent_dft
    cdef ITYPE_t[:] local_vars, result, Wn, We
    cdef bint prev_ret_value = False

    cdef ITYPE_t[:] edge_sources = np.empty(n_edges + n_verts + 1, dtype=ITYPE)
    cdef ITYPE_t[:] edge_targets = np.empty(n_edges + n_verts + 1, dtype=ITYPE)
    cdef ITYPE_t[:] edge_capacities = np.empty(n_edges + n_verts + 1, dtype=ITYPE)
    cdef ITYPE_t[:] edge_weights = np.empty(n_edges + n_verts + 1, dtype=ITYPE)
    idx = 0
    for e in range(n_edges):
        if capacities[e] != 0:
            edge_sources[idx] = tails[e] + 1
            edge_targets[idx] = heads[e] + 1
            idx += 1
    # print("edge_sources: ", end=" ")
    _print_array(edge_sources)
    # print("edge_targets: ", end=" ")
    _print_array(edge_targets)
    for v in range(n_verts):
        if demand[v] > 0:
            edge_sources[idx] = 0
            edge_targets[idx] = v + 1
        else:
            edge_sources[idx] = v + 1
            edge_targets[idx] = 0
        idx += 1

    capacities_sum, cost_sum = 0, 0
    for e in range(n_edges):
        if capacities[e] < 1<<31 - 1:
            capacities_sum += capacities[e]
        cost_sum += abs(cost[e])
    faux_inf = max(cost_sum, capacities_sum)
    for v in range(n_verts):
        faux_inf = max(faux_inf, demand[v])
    faux_inf = 3*faux_inf or 1
    # print("faux_inf: ", faux_inf)
    for e in range(n_edges):
        edge_capacities[e] = capacities[e]
        edge_weights[e] = cost[e]
    for v in range(n_verts):
        edge_capacities[v + n_edges] = faux_inf
        edge_weights[v + n_edges] = faux_inf

    edge_flow = np.zeros(n_edges + n_verts, dtype=ITYPE)
    vertex_potentials = np.empty(n_verts + 1, dtype=ITYPE)
    parent = np.empty(n_verts + 1, dtype=ITYPE)
    parent_edge = np.empty(n_edges + n_verts, dtype=ITYPE)
    subtree_size = np.empty(n_verts + 1, dtype=ITYPE)
    next_vertex_dft = np.empty(n_verts + 2, dtype=ITYPE)
    prev_vertex_dft = np.empty(n_verts + 1, dtype=ITYPE)
    last_descendent_dft = np.empty(n_verts + 1, dtype=ITYPE)
    _initialize_spanning_tree(n_verts, n_edges, faux_inf,
                              demand, edge_flow, vertex_potentials,
                              parent, parent_edge, subtree_size,
                              next_vertex_dft, prev_vertex_dft,
                              last_descendent_dft)

    local_vars = np.empty(4, dtype=ITYPE)
    result = np.empty(3, dtype=ITYPE)
    prev_ret_value = _find_entering_edges(n_edges, edge_weights,
                                          vertex_potentials, edge_flow,
                                          edge_sources, edge_targets,
                                          local_vars, prev_ret_value,
                                          result)
    itr = 0
    Wn = np.empty(n_edges + 1, dtype=ITYPE)
    We = np.empty(n_edges + 1, dtype=ITYPE)
    while prev_ret_value and itr < 100:
        Wne_len = np.empty(2, dtype=ITYPE)
        i, p, q = result[0], result[1], result[2]
        print("i, p, q: ", i, p, q)
        _find_cycle(i, p, q,
                    subtree_size,
                    parent, parent_edge,
                    n_edges, n_verts,
                    Wn, We, Wne_len)
        print("Wn: ", end="")
        _print_array(Wn)
        print("We: ", end="")
        _print_array(We)
        print("Wne_len: ", Wne_len[0], Wne_len[1])
        print("After _find_cycle")
        # return edge_flow
        _find_leaving_edge(Wn, We, Wne_len,
                           edge_sources, edge_targets,
                           edge_capacities, edge_flow,
                           result)
        print("After _find_leaving_edge")
        j = result[0]
        s = result[1]
        t = result[2]
        print("j, s, t: ", j, s, t)
        _augment_flow(Wn, We, Wne_len,
                      _residual_capacity(j, s, edge_sources,
                                         edge_capacities,
                                         edge_flow),
                      edge_sources, edge_flow)
        print("After _augment_flow")
        print("i, j: ", i, j)
        # Do nothing more if the entering edge is the same as the leaving edge.
        if i != j:
            if parent[t] != s:
                # Ensure that s is the parent of t.
                tmp = t
                t = s
                s = tmp
            if _index(We, Wne_len[1], i) > _index(We, Wne_len[1], j):
                # Ensure that q is in the subtree rooted at t.
                tmp = q
                q = p
                p = tmp
            print("_index_i, _index_j: ", _index(We, Wne_len[1], i),
                                          _index(We, Wne_len[1], j))
            print("parent[-1]: ", parent[n_verts])
            _remove_edge(s, t, subtree_size,
                         prev_vertex_dft,
                         last_descendent_dft,
                         next_vertex_dft,
                         parent, parent_edge,
                         n_verts)
            print("After _remove_edge")
            # return edge_flow
            _make_root(n_verts, q,
                       subtree_size,
                       last_descendent_dft,
                       prev_vertex_dft,
                       next_vertex_dft,
                       parent, parent_edge)
            print("After _make_root")
            _add_edge(i, p, q,
                      last_descendent_dft,
                      next_vertex_dft,
                      prev_vertex_dft,
                      subtree_size,
                      parent, parent_edge,
                      n_verts)
            print("After _add_edge")
            _update_potentials(i, p, q, n_verts,
                               edge_targets, edge_weights,
                               last_descendent_dft,
                               next_vertex_dft,
                               vertex_potentials)
            print("After _update_potentials")
        prev_ret_value = _find_entering_edges(n_edges, edge_weights,
                                              vertex_potentials, edge_flow,
                                              edge_sources, edge_targets,
                                              local_vars, prev_ret_value,
                                              result)
        print("After _find_entering_edges")
        itr += 1

    return edge_flow