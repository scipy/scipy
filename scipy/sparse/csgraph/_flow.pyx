# cython: wraparound=False, boundscheck=False

import numpy as np

from scipy.sparse import csr_matrix, issparse
from scipy.sparse._sputils import convert_pydata_sparse_to_scipy, is_pydata_spmatrix

cimport numpy as np

include 'parameters.pxi'

np.import_array()


class MaximumFlowResult:
    """Represents the result of a maximum flow calculation.

    Attributes
    ----------
    flow_value : int
        The value of the maximum flow.
    flow : csr_matrix
        The maximum flow.
    """

    def __init__(self, flow_value, flow):
        self.flow_value = flow_value
        self.flow = flow

    def __repr__(self):
        return 'MaximumFlowResult with value of %d' % self.flow_value


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
        and the flow graph in ``flow``.

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

    >>> import numpy as np
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
    ...                     [0,  0, 10, 12,  0,  0],
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
    flow function to the block corresponding to the original graph:

    >>> result = maximum_flow(graph_flow, 0, i+j+1, method='dinic')
    >>> matching = result.flow[1:i+1, i+1:i+j+1]
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
    is_pydata_sparse = is_pydata_spmatrix(csgraph)
    if is_pydata_sparse:
        pydata_sparse_cls = csgraph.__class__
    csgraph = convert_pydata_sparse_to_scipy(csgraph, target_format="csr")
    if not (issparse(csgraph) and csgraph.format == "csr"):
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
        flow = _edmonds_karp(m.indptr, tails, m.indices,
                             m.data, rev_edge_ptr, source, sink)
    elif method == 'dinic':
        flow = _dinic(m.indptr, m.indices, m.data, rev_edge_ptr,
                      source, sink)
    else:
        raise ValueError('{} method is not supported yet.'.format(method))
    flow_array = np.asarray(flow)
    flow_matrix = csr_matrix((flow_array, m.indices, m.indptr),
                             shape=m.shape)
    if is_pydata_sparse:
        flow_matrix = pydata_sparse_cls.from_scipy_sparse(flow_matrix)
    source_flow = flow_array[m.indptr[source]:m.indptr[source + 1]]
    return MaximumFlowResult(source_flow.sum(), flow_matrix)


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
        ITYPE_t sink) noexcept:
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
        The flow graph with respect to a maximum flow.

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
    cdef ITYPE_t cur, df, t, e, k

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
        ) noexcept nogil:
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
        ) noexcept nogil:
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
        The flow graph with respect to a maximum flow.
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
        ITYPE_t sink) noexcept:
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
        The flow graph with respect to a maximum flow.

    """
    cdef ITYPE_t n_verts = edge_ptr.shape[0] - 1
    cdef ITYPE_t n_edges = capacities.shape[0]

    cdef ITYPE_t[:] levels = np.empty(n_verts, dtype=ITYPE)
    cdef ITYPE_t[:] progress = np.empty(n_verts, dtype=ITYPE)
    cdef ITYPE_t[:] q = np.empty(n_verts, dtype=ITYPE)
    cdef ITYPE_t[:, :] stack = np.empty((n_verts, 2), dtype=ITYPE)
    cdef ITYPE_t[:] flows = np.zeros(n_edges, dtype=ITYPE)
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
