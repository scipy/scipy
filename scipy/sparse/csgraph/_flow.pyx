# cython: wraparound=False, boundscheck=False

import numpy as np

from scipy.sparse import csr_matrix, issparse
from scipy.sparse._sputils import convert_pydata_sparse_to_scipy, is_pydata_spmatrix

cimport numpy as np
from libc.math cimport ceil, sqrt
from libc.stdlib cimport abs, malloc, free
from libc.limits cimport INT_MAX, INT_MIN

include 'parameters.pxi'

np.import_array()

DEF NO_PARENT_PLACEHOLDER = -2


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

    @property
    def residual(self):
        warnings.warn(
            "The attribute `residual` has been renamed to `flow`"
            " and will be removed in SciPy 1.11.",
            DeprecationWarning, stacklevel=2
        )
        return self.flow

ctypedef struct edge_result:
    # A structure to store edges found out
    # by _find_entering_edges function.
    ITYPE_t i
    ITYPE_t p
    ITYPE_t q

ctypedef struct return_struct:
    # A structure to store the return value
    # and current state of _find_entering_edges
    # function.
    ITYPE_t B, M, m, f
    edge_result i_p_q
    bint prev_ret_value

class MinimumCostFlowResult:

    def __init__(self, flow_, flow_cost):
        self.flow = flow_
        self.flow_cost = flow_cost

    def __repr__(self):
        return 'MinimumCostFlowResult with a cost of %d' % self.flow_cost


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
        const ITYPE_t[:] edge_ptr,
        const ITYPE_t[:] tails,
        const ITYPE_t[:] heads,
        const ITYPE_t[:] capacities,
        const ITYPE_t[:] rev_edge_ptr,
        const ITYPE_t source,
        const ITYPE_t sink) noexcept:
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
        const ITYPE_t[:] edge_ptr,
        const ITYPE_t[:] heads,
        ITYPE_t[:] capacities,
        const ITYPE_t[:] rev_edge_ptr,
        const ITYPE_t source,
        const ITYPE_t sink) noexcept:
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


def minimum_cost_flow(csgraph, demand, cost):
    r"""
    minimum_cost_flow(csgraph, demand, cost)

    Finds a flow which minimizes the total cost while satisfying all
    vertices' demands in the given directed graph.

    .. versionadded:: 1.9.0

    Parameters
    ----------
    csgraph : csr_matrix
        The square matrix representing a directed graph whose (i, j)'th entry
        is an integer representing the capacity of the edge between
        vertices i and j.
    demand : array_like
        One dimensional array whose v'th entry is an integer representing the demand
        of vertex v.
    cost : array_like
        One dimensional array whose e'th entry is an integer representing the cost
        of per unit flow passing through edge e.

    Returns
    -------
    res : MinimumCostFlowResult
        A minimum cost flow represented by a ``MinimumCostFlowResult``
        which includes the cost of the flow in ``flow_cost`` and
        the flow graph in ``flow``.

    Raises
    ------
    TypeError:
        if the input graph is not in CSR format.

    ValueError:
        If any of the following conditions is met:
          - the capacity values, vertex demands or edge costs are not integers
          - an edge cost or vertex demand are infinite in absolute value
          - an edge has negative capacity
          - the sum of demands is not zero
          - there doesn't exist a flow which satisfies all vertex demands.

    Notes
    -----
    This solves the minimum cost flow problem [1]_ on a given directed weighted graph:
    a flow associates to every edge a value, also called a flow, less than the
    capacity of the edge, and for provided cost of the edge so that for every
    vertex (apart from the source and the sink vertices), the difference between
    the total incoming flow and the total outgoing flow is equal to the demand
    of that vertex. The value of a flow is the sum of the flow on all edges and
    the minimum cost flow problem consists of finding a flow whose cost (which is
    defined as the scalar product of the edges' flows and costs) is minimal.

    To solve the problem, we use Network Simplex [2]_.
    The time complexity of the former :math:`O(|V|^2\,|E|\,log(|V|\,C))`.

    The minimum cost flow problem is usually defined with real valued capacities,
    edge costs and vertex demands but we require that all capacities are integral.
    When dealing with rational capacities, or capacities belonging to
    :math:`x\mathbb{Q}` for some fixed :math:`x \in \mathbb{R}`, it is possible
    to reduce the problem to the integral case by scaling all capacities
    accordingly.

    Solving a minimum cost flow problem has many practical uses such as
    in assignment problem, covers and matching in bipartite graphs [3]_.

    We have ported the Python implementation of NetworkX [4]_ to Cython.

    References
    ----------

    .. [1] https://en.wikipedia.org/wiki/Minimum-cost_flow_problem#Definition
    .. [2] Orlin, J.B. A polynomial time primal network simplex algorithm
           for minimum cost flows. Mathematical Programming 78,
           109-129 (1997).
    .. [3] https://en.wikipedia.org/wiki/Network_simplex_algorithm#Applications
    .. [4] A Hagberg, D Schult, P Swart, Exploring Network Structure, Dynamics,
           and Function using NetworkX in Proceedings of the 7th Python in
           Science conference (SciPy 2008), G Varoquaux, T Vaught, J Millman
           (Eds.), pp. 11-15

    Examples
    --------
    >>> from scipy.sparse import csr_matrix
    >>> from scipy.sparse.csgraph import minimum_cost_flow
    >>> A = csr_matrix([[0, 4, 10, 0],
    ...                 [0, 0, 0, 9],
    ...                 [0, 0, 0, 5],
    ...                 [0, 0, 0, 0]])
    >>> demand = [-5, 0, 0, 5]
    >>> cost = [3, 6, 1, 2]
    >>> result = minimum_cost_flow(A, demand, cost)
    >>> result.flow_cost
    24
    >>> result.flow.toarray()
    array([[0, 4, 1, 0],
        [-4, 0, 0, 4],
        [-1, 0, 0, 1],
        [0, -4, -1, 0]], dtype=int32)
    """
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

    n_verts = csgraph.indptr.shape[0] - 1
    n_edges = cost.shape[0]

    _network_simplex_checks(csgraph.data, demand, cost)
    cdef ITYPE_t[:] tails = _make_tails(csgraph)
    cdef ITYPE_t[:] row = np.empty((2 * n_edges,), dtype=ITYPE)
    cdef ITYPE_t[:] col = np.empty((2 * n_edges,), dtype=ITYPE)
    cdef ITYPE_t[:] flow_data = np.empty((2 * n_edges,), dtype=ITYPE)
    ns_result = _network_simplex(csgraph.indices, tails, csgraph.data,
                                 demand, cost, n_verts,
                                 row, col, flow_data)
    if not ns_result.is_correct:
        raise ValueError("no flow satisfies all node demands")
    flow_matrix = csr_matrix((flow_data[0:ns_result.size],
                             (row[0:ns_result.size], col[0:ns_result.size])),
                             shape=(n_verts, n_verts))
    return MinimumCostFlowResult(flow_matrix, ns_result.flow_cost)


def _network_simplex_checks(
    capacities,           # IN
    demand,               # IN
    cost,                 # IN
):
    """
    Performs checks on input parameter to make sure that they
    comply with the values of the network simplex algorithms.
    """
    if np.any(demand == INT_MAX) or np.any(demand == INT_MIN):
        raise ValueError("one of the vertices has infinite demand")

    if np.any(cost == INT_MAX) or np.any(cost == INT_MIN):
        raise ValueError("one of the edges has infinite cost")

    if np.any(capacities < 0):
        raise ValueError("one of the edges has negative capacity")

    demand_sum = demand.sum()

    if demand_sum != 0:
        raise ValueError("sum of demands is not zero")


cdef void _initialize_spanning_tree(
    ITYPE_t n_verts,                 # IN
    ITYPE_t n_non_zero_edges,        # IN
    ITYPE_t faux_inf,                # IN
    ITYPE_t[:] demand,               # IN
    ITYPE_t[:] edge_flow,            # IN/OUT
    ITYPE_t[:] vertex_potentials,    # IN/OUT
    ITYPE_t[:] parent,               # IN/OUT
    ITYPE_t[:] parent_edge,          # IN/OUT
    ITYPE_t[:] subtree_size,         # IN/OUT
    ITYPE_t[:] next_vertex_dft,      # IN/OUT
    ITYPE_t[:] prev_vertex_dft,      # IN/OUT
    ITYPE_t[:] last_descendent_dft   # IN/OUT
) nogil:
    """Initialises a spanning tree which can be a feasible solution to the
    input minimum cost flow problem.

    The _network_simplex function then performs the algorithm to iteratively
    update this tree until an optimal solution is achieved.
    """
    cdef ITYPE_t e, v

    # edge_flow initialization
    for e in range(n_non_zero_edges):
        edge_flow[e] = 0
    for v in range(n_verts):
        edge_flow[v + n_non_zero_edges] = abs(demand[v])

    # vertex_potentials initialization
    # 0-th index is not needed in the algorithm
    # so it is left un-initialized
    for v in range(1, n_verts + 1):
        if demand[v - 1] <= 0:
            vertex_potentials[v] = faux_inf
        else:
            vertex_potentials[v] = -faux_inf

    # parent initialization
    parent[0] = NO_PARENT_PLACEHOLDER
    subtree_size[0] = n_verts + 1
    prev_vertex_dft[0] = n_verts
    last_descendent_dft[0] = n_verts
    for v in range(1, n_verts + 1):
        parent[v] = 0
        parent_edge[v] = n_non_zero_edges + v - 1
        subtree_size[v] = 1
        prev_vertex_dft[v] = v - 1
        last_descendent_dft[v] = v

    # next_vertex_dft initialization
    for v in range(n_verts):
        next_vertex_dft[v] = v + 1
    next_vertex_dft[n_verts] = 0


cdef inline ITYPE_t _reduced_cost(
    ITYPE_t e,                       # IN
    ITYPE_t[:] edge_weights,         # IN
    ITYPE_t[:] vertex_potentials,    # IN
    ITYPE_t[:] edge_sources,         # IN
    ITYPE_t[:] edge_targets,         # IN
    ITYPE_t[:] edge_flow             # IN
) nogil:
    """Returns the reduced cost of an edge computed using
    potentials of the vertices of the given edge and its weight.

    The reduced cost is needed to determine whether adding the
    edge to current spanning tree improves circulation or not.
    """
    cdef ITYPE_t c = (edge_weights[e]
         - vertex_potentials[edge_sources[e]]
         + vertex_potentials[edge_targets[e]]
    )
    return c if edge_flow[e] == 0 else -c

cdef class network_simplex_result:
    """
    Cython class to store the network
    simplex result.
    """
    cdef readonly ITYPE_t flow_cost
    cdef readonly ITYPE_t size
    cdef readonly bint is_correct

cdef return_struct _find_entering_edges(
    ITYPE_t n_edges,                     # IN
    ITYPE_t[:] edge_weights,             # IN
    ITYPE_t[:] vertex_potentials,        # IN
    ITYPE_t[:] edge_flow,                # IN
    ITYPE_t[:] edge_sources,             # IN
    ITYPE_t[:] edge_targets,             # IN
    return_struct r_struct,              # IN
) nogil:
    """Pivoting strategy of the network simplex algorithm.

    We have used the NetworkX implementation of this function.
    However, we don't use yield and instead manually maintain
    the state of this function using return_struct defined above.

    Following is the explanation from ``networkx.algorithms.flow.networksimplex``:

        "Entering edges are found by combining Dantzig's rule and Bland's
        rule. The edges are cyclically grouped into blocks of size B. Within
        each block, Dantzig's rule is applied to find an entering edge. The
        blocks to search is determined following Bland's rule."
    """
    cdef:
        ITYPE_t l, i = NO_PARENT_PLACEHOLDER, min_r_cost, p, q, c
        ITYPE_t e, r_cost, B

    if not r_struct.prev_ret_value:
        if n_edges == 0:
            r_struct.prev_ret_value = False
            return r_struct

        r_struct.B = < ITYPE_t > ceil(sqrt(n_edges))  # B
        B = r_struct.B
        r_struct.M = (n_edges + B - 1) // B  # M
        r_struct.m = 0  # m
        # entering edges
        r_struct.f = 0  # f

    while r_struct.m < r_struct.M:
        # Determine the next block of edges.
        l = r_struct.f + r_struct.B
        min_r_cost = INT_MAX
        if l <= n_edges:
            for e in range(r_struct.f, l):
                r_cost = _reduced_cost(e, edge_weights,
                                       vertex_potentials,
                                       edge_sources,
                                       edge_targets, edge_flow)
                if r_cost < min_r_cost:
                    i = e
                    min_r_cost = r_cost
        else:
            l -= n_edges
            for e in range(r_struct.f, n_edges):
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
        r_struct.f = l
        if min_r_cost >= 0:
            # No entering edge found in the current block.
            r_struct.m += 1
        else:
            # Entering edge found.
            if edge_flow[i] == 0:
                p = edge_sources[i]
                q = edge_targets[i]
            else:
                p = edge_targets[i]
                q = edge_sources[i]
            r_struct.i_p_q.i = i
            r_struct.i_p_q.p = p
            r_struct.i_p_q.q = q
            r_struct.m = 0
            r_struct.prev_ret_value = True
            return r_struct

    r_struct.prev_ret_value = False
    return r_struct

cdef ITYPE_t _find_apex(
    ITYPE_t p,                           # IN
    ITYPE_t q,                           # IN
    ITYPE_t[:] subtree_size,             # IN
    ITYPE_t[:] parent                    # IN
) nogil:
    """Find the lowest common ancestor of nodes p and q in the spanning tree."""
    cdef:
        ITYPE_t size_p = subtree_size[p]
        ITYPE_t size_q = subtree_size[q]
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
    ITYPE_t p,                           # IN
    ITYPE_t w,                           # IN
    ITYPE_t[:] parent,                   # IN
    ITYPE_t[:] parent_edge,              # IN
    ITYPE_t[:] Wn,                       # IN/OUT
    ITYPE_t[:] We                        # IN/OUT
) nogil:
    """Return the nodes and edges on the path from node p to its ancestor w."""
    cdef:
        ITYPE_t idx = 0
    Wn[0] = p
    while p != w:
        We[idx] = parent_edge[p]
        p = parent[p]
        Wn[idx + 1] = p
        idx += 1
    return idx

cdef void _find_cycle(
    ITYPE_t i,                     # IN
    ITYPE_t p,                     # IN
    ITYPE_t q,                     # IN
    ITYPE_t[:] subtree_size,       # IN
    ITYPE_t[:] parent,             # IN
    ITYPE_t[:] parent_edge,        # IN
    ITYPE_t n_edges,               # IN
    ITYPE_t n_verts,               # IN
    ITYPE_t[:] WnR,                # IN
    ITYPE_t[:] WeR,                # IN
    ITYPE_t[:] Wn,                 # IN/OUT
    ITYPE_t[:] We,                 # IN/OUT
    ITYPE_t[:] Wne_len             # IN/OUT
) nogil:
    """Return the nodes and edges on the cycle containing edge i == (p, q) when the
    latter is added to the spanning tree.

    The cycle is oriented in the direction from p to q.
    """
    cdef:
        ITYPE_t w = _find_apex(p, q, subtree_size, parent)
        ITYPE_t num_edges = _trace_path(p, w, parent, parent_edge, Wn, We)
        ITYPE_t Wn_len = num_edges + 1
        ITYPE_t We_len = num_edges
        ITYPE_t num_edges_R = _trace_path(q, w, parent, parent_edge, WnR, WeR)
        ITYPE_t idx, offset, tmp

    for idx in range(We_len//2):
        We[idx], We[We_len - idx - 1] = We[We_len - idx - 1], We[idx]
    for idx in range(Wn_len//2):
        Wn[idx], Wn[Wn_len - idx - 1] = Wn[Wn_len - idx - 1], Wn[idx]
    if (We_len == 1 and We[0] != i) or We_len != 1:
        We[We_len] = i
        We_len += 1
    for idx in range(num_edges_R):
        Wn[idx + Wn_len] = WnR[idx]
        We[idx + We_len] = WeR[idx]
    Wn_len += num_edges_R
    We_len += num_edges_R
    Wne_len[0] = Wn_len
    Wne_len[1] = We_len

cdef inline ITYPE_t _residual_capacity(
    ITYPE_t i,                     # IN
    ITYPE_t p,                     # IN
    ITYPE_t[:] edge_sources,       # IN
    ITYPE_t[:] edge_capacities,    # IN
    ITYPE_t[:] edge_flow           # IN
) nogil:
    """Return the residual capacity of an edge i in the direction away from
    its endpoint p."""
    if edge_sources[i] == p:
        return edge_capacities[i] - edge_flow[i]
    else:
        return edge_flow[i]

cdef edge_result _find_leaving_edge(
    ITYPE_t[:] Wn,                 # IN
    ITYPE_t[:] We,                 # IN
    ITYPE_t[:] Wne_len,            # IN
    ITYPE_t[:] edge_sources,       # IN
    ITYPE_t[:] edge_targets,       # IN
    ITYPE_t[:] edge_capacities,    # IN
    ITYPE_t[:] edge_flow,          # IN
    edge_result ret_values         # IN/OUT
) nogil:
    """Return the leaving edge in a cycle represented by Wn and We."""
    cdef:
        ITYPE_t min_res_cap = INT_MAX
        ITYPE_t res_cap
        ITYPE_t i = NO_PARENT_PLACEHOLDER
        ITYPE_t j = NO_PARENT_PLACEHOLDER
        ITYPE_t s = NO_PARENT_PLACEHOLDER
        ITYPE_t p, t
        ITYPE_t idx_e =  Wne_len[1] - 1
        ITYPE_t idx_n = Wne_len[0] - 1

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
    ret_values.i = j
    ret_values.p = s
    ret_values.q = t
    return ret_values

cdef inline ITYPE_t _index(
    ITYPE_t[:] W,      # IN
    ITYPE_t W_len,     # IN
    ITYPE_t i          # IN
) nogil:
    """Returns the index of value in the given array."""
    cdef ITYPE_t idx
    for idx in range(W_len):
        if W[idx] == i:
            return idx
    return -1

cdef inline void _augment_flow(
    ITYPE_t[:] Wn,                # IN
    ITYPE_t[:] We,                # IN
    ITYPE_t[:] Wne_len,           # IN
    ITYPE_t f,                    # IN
    ITYPE_t[:] edge_sources,      # IN
    ITYPE_t[:] edge_flow          # IN/OUT
) nogil:
    """Augment f units of flow along a cycle represented by Wn and We."""
    cdef:
        ITYPE_t i, p
        ITYPE_t idx_e = 0
        ITYPE_t idx_n = 0

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
    ITYPE_t s,                        # IN/OUT
    ITYPE_t t,                        # IN/OUT
    ITYPE_t[:] subtree_size,          # IN/OUT
    ITYPE_t[:] prev_vertex_dft,       # IN/OUT
    ITYPE_t[:] last_descendent_dft,   # IN/OUT
    ITYPE_t[:] next_vertex_dft,       # IN/OUT
    ITYPE_t[:] parent,                # IN/OUT
    ITYPE_t[:] parent_edge,           # IN/OUT
    ITYPE_t n_verts                   # IN/OUT
) nogil:
    """Remove an edge (s, t) where parent[t] == s from the spanning tree."""
    cdef:
        ITYPE_t size_t = subtree_size[t]
        ITYPE_t prev_t = prev_vertex_dft[t]
        ITYPE_t last_t = last_descendent_dft[t]
        ITYPE_t next_last_t = next_vertex_dft[last_t]

    # Remove (s, t).
    parent[t] = NO_PARENT_PLACEHOLDER
    parent_edge[t] = NO_PARENT_PLACEHOLDER
    # Remove the subtree rooted at t from the depth-first thread.
    next_vertex_dft[prev_t] = next_last_t
    prev_vertex_dft[next_last_t] = prev_t
    next_vertex_dft[last_t] = t
    prev_vertex_dft[t] = last_t
    # Update the subtree sizes and last descendants of the (old) ancestors
    # of t.
    while s != NO_PARENT_PLACEHOLDER:
        subtree_size[s] -= size_t
        if last_descendent_dft[s] == last_t:
            last_descendent_dft[s] = prev_t
        s = parent[s]

cdef void _make_root(
    ITYPE_t n_verts,                   # IN
    ITYPE_t q,                         # IN
    ITYPE_t[:] subtree_size,           # IN/OUT
    ITYPE_t[:] last_descendent_dft,    # IN/OUT
    ITYPE_t[:] prev_vertex_dft,        # IN/OUT
    ITYPE_t[:] next_vertex_dft,        # IN/OUT
    ITYPE_t[:] parent,                 # IN/OUT
    ITYPE_t[:] parent_edge,            # IN/OUT
    ITYPE_t[:] ancestors,              # IN/OUT
) nogil:
    """Make a node q the root of its containing subtree."""
    cdef:
        ITYPE_t idx = 0
        ITYPE_t path_len, p
        ITYPE_t size_p, last_p, last_q, prev_q
        ITYPE_t next_last_q

    while q != NO_PARENT_PLACEHOLDER:
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
        parent[q] = NO_PARENT_PLACEHOLDER
        parent_edge[p] = parent_edge[q]
        parent_edge[q] = NO_PARENT_PLACEHOLDER
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
    ITYPE_t i,                         # IN
    ITYPE_t p,                         # IN
    ITYPE_t q,                         # IN
    ITYPE_t[:] last_descendent_dft,    # IN/OUT
    ITYPE_t[:] next_vertex_dft,        # IN/OUT
    ITYPE_t[:] prev_vertex_dft,        # IN/OUT
    ITYPE_t[:] subtree_size,           # IN/OUT
    ITYPE_t[:] parent,                 # IN/OUT
    ITYPE_t[:] parent_edge,            # IN/OUT
) nogil:
    """Add an edge (p, q) to the spanning tree where q is the root of a subtree."""
    cdef:
        ITYPE_t last_p = last_descendent_dft[p]
        ITYPE_t next_last_p = next_vertex_dft[last_p]
        ITYPE_t size_q = subtree_size[q]
        ITYPE_t last_q = last_descendent_dft[q]

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
    while p != NO_PARENT_PLACEHOLDER:
        subtree_size[p] += size_q
        if last_descendent_dft[p] == last_p:
            last_descendent_dft[p] = last_q
        p = parent[p]

cdef inline bint _trace_subtree(
    ITYPE_t p,                         # IN
    ITYPE_t[:] last_descendent_dft,    # IN
    ITYPE_t[:] next_vertex_dft,        # IN
    ITYPE_t[:] subtree                 # IN/OUT
) nogil:
    """Return the nodes in the subtree rooted at a node p."""
    cdef:
        ITYPE_t idx = 1, l
    subtree[0] = p

    l = last_descendent_dft[p]
    while p != l:
        p = next_vertex_dft[p]
        subtree[idx] = p
        idx += 1
    return idx

cdef void _update_potentials(
    ITYPE_t i,                         # IN
    ITYPE_t p,                         # IN
    ITYPE_t q,                         # IN
    ITYPE_t[:] edge_targets,           # IN
    ITYPE_t[:] edge_weights,           # IN
    ITYPE_t[:] last_descendent_dft,    # IN
    ITYPE_t[:] next_vertex_dft,        # IN
    ITYPE_t[:] vertex_potentials,      # IN/OUT
    ITYPE_t[:] subtree,                # IN/OUT
) nogil:
    """
    Update the potentials of the nodes in the subtree rooted at a node
    q connected to its parent p by an edge i.
    """
    cdef ITYPE_t idx, d, path_len

    if q == edge_targets[i]:
        d = vertex_potentials[p] - edge_weights[i] - vertex_potentials[q]
    else:
        d = vertex_potentials[p] + edge_weights[i] - vertex_potentials[q]
    path_len = _trace_subtree(q, last_descendent_dft,
                              next_vertex_dft, subtree)
    for idx in range(path_len):
        q = subtree[idx]
        vertex_potentials[q] += d

cdef inline void _add_entry(
    ITYPE_t edge_source,               # IN
    ITYPE_t edge_target,               # IN
    ITYPE_t flow,                      # IN
    ITYPE_t idx,                       # IN
    ITYPE_t n_edges,
    ITYPE_t[:] row,                    # IN/OUT
    ITYPE_t[:] col,                    # IN/OUT
    ITYPE_t[:] flow_data               # IN/OUT
) nogil:
    """Add the result values as in entry in row, col and flow_data arrays."""
    row[idx] = edge_source - 1
    row[idx + 1] = edge_target - 1
    col[idx] = edge_target - 1
    col[idx + 1] = edge_source - 1
    flow_data[idx] = flow
    flow_data[idx + 1] = -flow

cdef network_simplex_result _network_simplex(
    ITYPE_t[:] heads,                  # IN
    ITYPE_t[:] tails,                  # IN
    ITYPE_t[:] capacities,             # IN
    ITYPE_t[:] demand,                 # IN
    ITYPE_t[:] cost,                   # IN
    ITYPE_t n_verts,                   # IN
    ITYPE_t[:] row,                    # IN/OUT
    ITYPE_t[:] col,                    # IN/OUT
    ITYPE_t[:] flow_data               # IN/OUT
):
    cdef:
        ITYPE_t n_edges = capacities.shape[0]
        ITYPE_t idx, faux_inf, capacities_sum, cost_sum
        ITYPE_t j, s, t, i, tmp, p, q, n_non_zero_edges, v
        ITYPE_t flow_cost, flow_value
        ITYPE_t[:] edge_flow = np.zeros(n_edges + n_verts, dtype=ITYPE)
        ITYPE_t[:] vertex_potentials = np.empty(n_verts + 1, dtype=ITYPE)
        ITYPE_t[:] Wne_len = np.empty(2, dtype=ITYPE)
        ITYPE_t[:] parent = np.empty(n_verts + 1, dtype=ITYPE)
        ITYPE_t[:] parent_edge = np.empty(n_edges + n_verts, dtype=ITYPE)
        ITYPE_t[:] subtree_size = np.empty(n_verts + 1, dtype=ITYPE)
        ITYPE_t[:] next_vertex_dft = np.empty(n_verts + 1, dtype=ITYPE)
        ITYPE_t[:] prev_vertex_dft = np.empty(n_verts + 1, dtype=ITYPE)
        ITYPE_t[:] last_descendent_dft = np.empty(n_verts + 1, dtype=ITYPE)
        ITYPE_t[:] Wn = np.empty(n_edges + 1, dtype=ITYPE)
        ITYPE_t[:] We = np.empty(n_edges + 1, dtype=ITYPE)
        ITYPE_t[:] WnR = np.empty(n_verts, dtype=ITYPE)
        ITYPE_t[:] WeR = np.empty(n_edges, dtype=ITYPE)
        ITYPE_t[:] ancestors = np.empty(n_verts + 1, dtype=ITYPE)
        ITYPE_t[:] subtree = np.empty(2*n_verts, dtype=ITYPE)
        ITYPE_t[:] edge_sources = np.empty(n_edges + n_verts + 1, dtype=ITYPE)
        ITYPE_t[:] edge_targets = np.empty(n_edges + n_verts + 1, dtype=ITYPE)
        ITYPE_t[:] edge_capacities = np.empty(n_edges + n_verts + 1, dtype=ITYPE)
        ITYPE_t[:] edge_weights = np.empty(n_edges + n_verts + 1, dtype=ITYPE)

        return_struct r_struct
        edge_result i_p_q
        network_simplex_result ns_result = network_simplex_result()

    r_struct.prev_ret_value = False

    with nogil:
        idx = 0
        n_non_zero_edges = 0
        for e in range(n_edges):
            if capacities[e] != 0:
                n_non_zero_edges += 1
                edge_sources[idx] = tails[e] + 1
                edge_targets[idx] = heads[e] + 1
                idx += 1
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
            if capacities[e] < INT_MAX:
                capacities_sum += capacities[e]
            cost_sum += abs(cost[e])
        faux_inf = max(cost_sum, capacities_sum)
        for v in range(n_verts):
            faux_inf = max(faux_inf, demand[v])
        faux_inf = 3*faux_inf or 1
        idx = 0
        for e in range(n_edges):
            if capacities[e] != 0:
                edge_capacities[idx] = capacities[e]
                edge_weights[idx] = cost[e]
                idx += 1
        for v in range(n_verts):
            edge_capacities[idx] = faux_inf
            edge_weights[idx] = faux_inf
            idx += 1

        _initialize_spanning_tree(n_verts, n_non_zero_edges, faux_inf,
                                  demand, edge_flow, vertex_potentials,
                                  parent, parent_edge, subtree_size,
                                  next_vertex_dft, prev_vertex_dft,
                                  last_descendent_dft)

        r_struct = _find_entering_edges(n_non_zero_edges, edge_weights,
                                        vertex_potentials, edge_flow,
                                        edge_sources, edge_targets,
                                        r_struct)

        while r_struct.prev_ret_value:
            i, p, q = r_struct.i_p_q.i, r_struct.i_p_q.p, r_struct.i_p_q.q
            _find_cycle(i, p, q,
                        subtree_size,
                        parent, parent_edge,
                        n_edges, n_verts,
                        WnR, WeR, Wn, We,
                        Wne_len)
            i_p_q = _find_leaving_edge(Wn, We, Wne_len,
                                       edge_sources, edge_targets,
                                       edge_capacities, edge_flow,
                                       i_p_q)
            j = i_p_q.i
            s = i_p_q.p
            t = i_p_q.q
            _augment_flow(Wn, We, Wne_len,
                          _residual_capacity(j, s, edge_sources,
                                             edge_capacities,
                                             edge_flow),
                          edge_sources, edge_flow)
            # Do nothing more if the entering edge is the
            # same as the leaving edge.
            if i != j:
                if parent[t] != s:
                    # Ensure that s is the parent of t.
                    t, s = s, t
                if _index(We, Wne_len[1], i) > _index(We, Wne_len[1], j):
                    # Ensure that q is in the subtree rooted at t.
                    q, p = p, q
                _remove_edge(s, t, subtree_size,
                             prev_vertex_dft,
                             last_descendent_dft,
                             next_vertex_dft,
                             parent, parent_edge,
                             n_verts)
                _make_root(n_verts, q,
                           subtree_size,
                           last_descendent_dft,
                           prev_vertex_dft,
                           next_vertex_dft,
                           parent, parent_edge,
                           ancestors)
                _add_edge(i, p, q,
                          last_descendent_dft,
                          next_vertex_dft,
                          prev_vertex_dft,
                          subtree_size,
                          parent, parent_edge)
                _update_potentials(i, p, q,
                                   edge_targets, edge_weights,
                                   last_descendent_dft,
                                   next_vertex_dft,
                                   vertex_potentials, subtree)
            r_struct = _find_entering_edges(n_non_zero_edges,
                                            edge_weights,
                                            vertex_potentials, edge_flow,
                                            edge_sources, edge_targets,
                                            r_struct)

        ns_result.is_correct = True
        for v in range(n_verts):
            if edge_flow[v + n_non_zero_edges] != 0:
                ns_result.is_correct = False
                break

        if ns_result.is_correct:
            flow_cost = 0
            for i in range(n_edges):
                flow_cost += edge_weights[i]*edge_flow[i]
            idx = 0
            for i in range(n_edges):
                if edge_flow[i] != 0 and edge_capacities[i] != 0:
                    _add_entry(edge_sources[i], edge_targets[i],
                               edge_flow[i], idx, 0,
                               row, col, flow_data)
                    idx += 2
            ns_result.flow_cost = flow_cost
            ns_result.size = idx

    # end: with nogil
    return ns_result
