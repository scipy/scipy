import warnings

cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport INFINITY


from scipy.sparse import issparse

np.import_array()

include "parameters.pxi"


def maximum_bipartite_matching(graph, perm_type='row'):
    r"""
    maximum_bipartite_matching(graph, perm_type='row')

    Returns a matching of a bipartite graph whose cardinality is as least that
    of any given matching of the graph.

    Parameters
    ----------
    graph : sparse matrix
        Input sparse in CSR format whose rows represent one partition of the
        graph and whose columns represent the other partition. An edge between
        two vertices is indicated by the corresponding entry in the matrix
        existing in its sparse representation.
    perm_type : str, {'row', 'column'}
        Which partition to return the matching in terms of: If ``'row'``, the
        function produces an array whose length is the number of columns in the
        input, and whose :math:`j`'th element is the row matched to the
        :math:`j`'th column. Conversely, if ``perm_type`` is ``'column'``, this
        returns the columns matched to each row.

    Returns
    -------
    perm : ndarray
        A matching of the vertices in one of the two partitions. Unmatched
        vertices are represented by a ``-1`` in the result.

    Notes
    -----
    This function implements the Hopcroft--Karp algorithm [1]_. Its time
    complexity is :math:`O(\lvert E \rvert \sqrt{\lvert V \rvert})`, and its
    space complexity is linear in the number of rows. In practice, this
    asymmetry between rows and columns means that it can be more efficient to
    transpose the input if it contains more columns than rows.

    By Konig's theorem, the cardinality of the matching is also the number of
    vertices appearing in a minimum vertex cover of the graph.

    Note that if the sparse representation contains explicit zeros, these are
    still counted as edges.

    The implementation was changed in SciPy 1.4.0 to allow matching of general
    bipartite graphs, where previous versions would assume that a perfect
    matching existed. As such, code written against 1.4.0 will not necessarily
    work on older versions.

    References
    ----------
    .. [1] John E. Hopcroft and Richard M. Karp. "An n^{5 / 2} Algorithm for
           Maximum Matchings in Bipartite Graphs" In: SIAM Journal of Computing
           2.4 (1973), pp. 225--231. :doi:`10.1137/0202019`

    Examples
    --------
    >>> from scipy.sparse import csr_matrix
    >>> from scipy.sparse.csgraph import maximum_bipartite_matching

    As a simple example, consider a bipartite graph in which the partitions
    contain 2 and 3 elements respectively. Suppose that one partition contains
    vertices labelled 0 and 1, and that the other partition contains vertices
    labelled A, B, and C. Suppose that there are edges connecting 0 and C,
    1 and A, and 1 and B. This graph would then be represented by the following
    sparse matrix:

    >>> graph = csr_matrix([[0, 0, 1], [1, 1, 0]])

    Here, the 1s could be anything, as long as they end up being stored as
    elements in the sparse matrix. We can now calculate maximum matchings as
    follows:

    >>> print(maximum_bipartite_matching(graph, perm_type='column'))
    [2 0]
    >>> print(maximum_bipartite_matching(graph, perm_type='row'))
    [ 1 -1  0]

    The first output tells us that 1 and 2 are matched with C and A
    respectively, and the second output tells us that A, B, and C are matched
    with 1, nothing, and 0 respectively.

    Note that explicit zeros are still converted to edges. This means that a
    different way to represent the above graph is by using the CSR structure
    directly as follows:

    >>> data = [0, 0, 0]
    >>> indices = [2, 0, 1]
    >>> indptr = [0, 1, 3]
    >>> graph = csr_matrix((data, indices, indptr))
    >>> print(maximum_bipartite_matching(graph, perm_type='column'))
    [2 0]
    >>> print(maximum_bipartite_matching(graph, perm_type='row'))
    [ 1 -1  0]

    When one or both of the partitions are empty, the matching is empty as
    well:

    >>> graph = csr_matrix((2, 0))
    >>> print(maximum_bipartite_matching(graph, perm_type='column'))
    [-1 -1]
    >>> print(maximum_bipartite_matching(graph, perm_type='row'))
    []

    When the input matrix is square, and the graph is known to admit a perfect
    matching, i.e. a matching with the property that every vertex in the graph
    belongs to some edge in the matching, then one can view the output as the
    permutation of rows (or columns) turning the input matrix into one with the
    property that all diagonal elements are non-empty:

    >>> a = [[0, 1, 2, 0], [1, 0, 0, 1], [2, 0, 0, 3], [0, 1, 3, 0]]
    >>> graph = csr_matrix(a)
    >>> perm = maximum_bipartite_matching(graph, perm_type='row')
    >>> print(graph[perm].toarray())
    [[1 0 0 1]
     [0 1 2 0]
     [0 1 3 0]
     [2 0 0 3]]

    """
    if not issparse(graph):
        raise TypeError("graph must be sparse")
    if graph.format not in ("csr", "csc", "coo"):
        raise TypeError("graph must be in CSC, CSR, or COO format.")
    graph = graph.tocsr()
    i, j = graph.shape
    x, y = _hopcroft_karp(graph.indices, graph.indptr, i, j)
    return np.asarray(x if perm_type == 'column' else y)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef tuple _hopcroft_karp(ITYPE_t[:] indices, ITYPE_t[:] indptr,
                          ITYPE_t i, ITYPE_t j):
    cdef ITYPE_t INF = np.iinfo(ITYPE).max
    # x will end up containing the matchings of rows to columns, while
    # y will contain the matchings of columns to rows.
    cdef ITYPE_t[:] x = np.empty(i, dtype=ITYPE)
    cdef ITYPE_t[:] y = np.empty(j, dtype=ITYPE)

    # During the BFS step, dist will keep track of the level of the search. We
    # only keep track of this for the vertices in the left partition. Note that
    # the array is one longer than the cardinality of the left partition:
    # throughout the algorithm, we make use of the standard trick of adding an
    # auxiliary vertex whose index will be i, and whose semantics are that
    # every unmatched column will be matched with this vertex.
    cdef ITYPE_t[:] dist = np.empty(i + 1, dtype=ITYPE)

    cdef ITYPE_t k, v, w, up, u, yu, u_old

    # At the end of the day, unmatched vertices will have a value of -1. As
    # mentioned above, unmatched vertices in the right partition will be
    # matched to the auxiliary vertex i, and we will set their value to -1
    # only when we have found a maximum matching.
    for k in range(i):
        x[k] = -1

    for k in range(j):
        y[k] = i

    # Set up three structures for use in our searches: q will represent a queue
    # of vertices for use during the BFS step. As above, we will only keep
    # track of the vertices in the left partition, as well as the auxiliary
    # vertex; as no vertex is added more than once, the capacity of the queue
    # will be i + 1. As in a circular buffer, we use two pointers, head and
    # tail to keep track of the ends of the queue: Elements are dequeued from
    # head and queued at tail.
    cdef ITYPE_t[:] q = np.empty(i + 1, dtype=ITYPE)
    cdef ITYPE_t head, tail

    # Similarly, we use a stack for our depth-first search. As above, we only
    # represent vertices in the left partition, and since no augmenting path
    # will visit more than i of thse before encountering an unmatched vertex
    # (as represented by i), the stack capacity can be limited to i + 1.
    # Elements will be pushed to stack_head and popped from stack_head - 1.
    cdef ITYPE_t[:] stack = np.empty(i + 1, dtype=ITYPE)
    cdef ITYPE_t stack_head

    # Finally, during our depth-first search, we keep track of the path along
    # which we move. This will simplify the updates to the matching that occur
    # when an augmenting path is found.
    cdef ITYPE_t[:] parents = np.empty(i, dtype=ITYPE)

    # The breadth-first search part of the algorithm. This will terminate when
    # we are unable to find a path to an unassigned vertex, which boils down to
    # not being able to find a path to the auxiliary vertex i.
    while True:
        # Empty the queue by resetting the two pointers.
        head = 0
        tail = 0
        for v in range(i):
            if x[v] < 0:
                dist[v] = 0
                # Enqueue v. Note that in an ordinary circular buffer, we would
                # avoid overflows by wrapping around indices, but since we will
                # never enqueue more than i + 1 different elements at a single
                # iteration, we can avoid doing so.
                q[tail] = v
                tail += 1
            else:
                dist[v] = INF
        dist[i] = INF
        # Check if the queue is empty. Note than in an ordinary circular buffer
        # some care needs to be taken with this check when the buffer is full,
        # which we can avoid as above.
        while head < tail:
            v = q[head]
            head += 1
            if dist[v] < dist[i]:
                for up in range(indptr[v], indptr[v+1]):
                    u = indices[up]
                    if dist[y[u]] == INF:
                        dist[y[u]] = dist[v] + 1
                        q[tail] = y[u]
                        tail += 1
        # Vertices not encountered during the BFS will have a dist of INF. In
        # particular, dist[i] will be INF exactly when we did not encounter any
        # unmatched vertices.
        if dist[i] == INF:
            break

        # Depth-first search starting from every unmatched vertex.
        for w in range(i):
            if x[w] < 0:
                done = False
                # Initialize stack to contain only w and reset path.
                stack[0] = w
                stack_head = 1
                while stack_head != 0:
                    # Pop v from stack.
                    stack_head -= 1
                    v = stack[stack_head]
                    could_augment = False
                    for up in range(indptr[v], indptr[v + 1]):
                        u = indices[up]
                        yu = y[u]
                        if dist[yu] == dist[v] + 1:
                            could_augment = True
                            # If yu is unmatched, we have found an augmenting
                            # path. We update the matching and move on to the
                            # next unmatched vertex.
                            if yu == i:
                                done = True
                                # Unwind and follow the path back to the root.
                                while True:
                                    # Mark v as visited to ensure that it
                                    # features in only one augmenting path in
                                    # this sequence of DFS runs.
                                    dist[v] = INF
                                    u_old = x[v]
                                    y[u] = v
                                    x[v] = u
                                    u = u_old
                                    if v == w:
                                        break
                                    v = parents[v]
                                break
                            stack[stack_head] = yu
                            stack_head += 1
                            parents[yu] = v
                    if done:
                        break

    for k in range(j):
        if y[k] == i:
            y[k] = -1

    return x, y


def min_weight_full_bipartite_matching(biadjacency_matrix, maximize=False):
    r"""
    min_weight_full_bipartite_matching(biadjacency_matrix, maximize=False)

    Returns the minimum weight full matching of a bipartite graph.

    .. versionadded:: 1.6.0

    Parameters
    ----------
    biadjacency_matrix : sparse matrix
        Biadjacency matrix of the bipartite graph: A sparse matrix in CSR, CSC,
        or COO format whose rows represent one partition of the graph and whose
        columns represent the other partition. An edge between two vertices is
        indicated by the corresponding entry in the matrix, and the weight of
        the edge is given by the value of that entry. This should not be
        confused with the full adjacency matrix of the graph, as we only need
        the submatrix defining the bipartite structure.

    maximize : bool (default: False)
        Calculates a maximum weight matching if true.

    Returns
    -------
    row_ind, col_ind : array
        An array of row indices and one of corresponding column indices giving
        the optimal matching. The total weight of the matching can be computed
        as ``graph[row_ind, col_ind].sum()``. The row indices will be
        sorted; in the case of a square matrix they will be equal to
        ``numpy.arange(graph.shape[0])``.

    Notes
    -----

    Let :math:`G = ((U, V), E)` be a weighted bipartite graph with non-zero
    weights :math:`w : E \to \mathbb{R} \setminus \{0\}`. This function then
    produces a matching :math:`M \subseteq E` with cardinality

    .. math::
       \lvert M \rvert = \min(\lvert U \rvert, \lvert V \rvert),

    which minimizes the sum of the weights of the edges included in the
    matching, :math:`\sum_{e \in M} w(e)`, or raises an error if no such
    matching exists.

    When :math:`\lvert U \rvert = \lvert V \rvert`, this is commonly
    referred to as a perfect matching; here, since we allow
    :math:`\lvert U \rvert` and :math:`\lvert V \rvert` to differ, we
    follow Karp [1]_ and refer to the matching as *full*.

    This function implements the LAPJVsp algorithm [2]_, short for "Linear
    assignment problem, Jonker--Volgenant, sparse".

    The problem it solves is equivalent to the rectangular linear assignment
    problem. [3]_ As such, this function can be used to solve the same problems
    as :func:`scipy.optimize.linear_sum_assignment`. That function may perform
    better when the input is dense, or for certain particular types of inputs,
    such as those for which the :math:`(i, j)`'th entry is the distance between
    two points in Euclidean space.

    If no full matching exists, this function raises a ``ValueError``. For
    determining the size of the largest matching in the graph, see
    :func:`maximum_bipartite_matching`.

    We require that weights are non-zero only to avoid issues with the handling
    of explicit zeros when converting between different sparse representations.
    Zero weights can be handled by adding a constant to all weights, so that
    the resulting matrix contains no zeros.

    References
    ----------
    .. [1] Richard Manning Karp:
       An algorithm to Solve the m x n Assignment Problem in Expected Time
       O(mn log n).
       Networks, 10(2):143-152, 1980.
    .. [2] Roy Jonker and Anton Volgenant:
       A Shortest Augmenting Path Algorithm for Dense and Sparse Linear
       Assignment Problems.
       Computing 38:325-340, 1987.
    .. [3] https://en.wikipedia.org/wiki/Assignment_problem

    Examples
    --------
    >>> from scipy.sparse import csr_matrix
    >>> from scipy.sparse.csgraph import min_weight_full_bipartite_matching

    Let us first consider an example in which all weights are equal:

    >>> biadjacency_matrix = csr_matrix([[1, 1, 1], [1, 0, 0], [0, 1, 0]])

    Here, all we get is a perfect matching of the graph:

    >>> print(min_weight_full_bipartite_matching(biadjacency_matrix)[1])
    [2 0 1]

    That is, the first, second, and third rows are matched with the third,
    first, and second column respectively. Note that in this example, the 0
    in the input matrix does *not* correspond to an edge with weight 0, but
    rather a pair of vertices not paired by an edge.

    Note also that in this case, the output matches the result of applying
    :func:`maximum_bipartite_matching`:

    >>> from scipy.sparse.csgraph import maximum_bipartite_matching
    >>> biadjacency = csr_matrix([[1, 1, 1], [1, 0, 0], [0, 1, 0]])
    >>> print(maximum_bipartite_matching(biadjacency, perm_type='column'))
    [2 0 1]

    When multiple edges are available, the ones with lowest weights are
    preferred:

    >>> biadjacency = csr_matrix([[3, 3, 6], [4, 3, 5], [10, 1, 8]])
    >>> row_ind, col_ind = min_weight_full_bipartite_matching(biadjacency)
    >>> print(col_ind)
    [0 2 1]

    The total weight in this case is :math:`3 + 5 + 1 = 9`:

    >>> print(biadjacency[row_ind, col_ind].sum())
    9

    When the matrix is not square, i.e. when the two partitions have different
    cardinalities, the matching is as large as the smaller of the two
    partitions:

    >>> biadjacency = csr_matrix([[0, 1, 1], [0, 2, 3]])
    >>> row_ind, col_ind = min_weight_full_bipartite_matching(biadjacency)
    >>> print(row_ind, col_ind)
    [0 1] [2 1]
    >>> biadjacency = csr_matrix([[0, 1], [3, 1], [1, 4]])
    >>> row_ind, col_ind = min_weight_full_bipartite_matching(biadjacency)
    >>> print(row_ind, col_ind)
    [0 2] [1 0]

    When one or both of the partitions are empty, the matching is empty as
    well:

    >>> biadjacency = csr_matrix((2, 0))
    >>> row_ind, col_ind = min_weight_full_bipartite_matching(biadjacency)
    >>> print(row_ind, col_ind)
    [] []

    In general, we will always reach the same sum of weights as if we had used
    :func:`scipy.optimize.linear_sum_assignment` but note that for that one,
    missing edges are represented by a matrix entry of ``float('inf')``. Let us
    generate a random sparse matrix with integer entries between 1 and 10:

    >>> import numpy as np
    >>> from scipy.sparse import random
    >>> from scipy.optimize import linear_sum_assignment
    >>> sparse = random(10, 10, random_state=42, density=.5, format='coo') * 10
    >>> sparse.data = np.ceil(sparse.data)
    >>> dense = sparse.toarray()
    >>> dense = np.full(sparse.shape, np.inf)
    >>> dense[sparse.row, sparse.col] = sparse.data
    >>> sparse = sparse.tocsr()
    >>> row_ind, col_ind = linear_sum_assignment(dense)
    >>> print(dense[row_ind, col_ind].sum())
    28.0
    >>> row_ind, col_ind = min_weight_full_bipartite_matching(sparse)
    >>> print(sparse[row_ind, col_ind].sum())
    28.0

    """
    if not issparse(biadjacency_matrix):
        raise TypeError("graph must be sparse")
    if biadjacency_matrix.format not in ("csr", "csc", "coo"):
        raise TypeError("graph must be in CSC, CSR, or COO format.")

    if not (np.issubdtype(biadjacency_matrix.dtype, np.number) or
            biadjacency_matrix.dtype == np.dtype(np.bool_)):
        raise ValueError("expected a matrix containing numerical entries, " +
                         "got %s" % (biadjacency_matrix.dtype,))

    biadjacency_matrix = biadjacency_matrix.astype(np.double)

    if maximize:
        biadjacency_matrix = -biadjacency_matrix

    # Change all infinities to zeros, then remove those zeros, but warn the
    # user if any zeros were present in the first place.
    if not np.all(biadjacency_matrix.data):
        warnings.warn('explicit zero weights are removed before matching')

    biadjacency_matrix.data[np.isposinf(biadjacency_matrix.data)] = 0
    biadjacency_matrix.eliminate_zeros()

    i, j = biadjacency_matrix.shape

    a = np.arange(np.min(biadjacency_matrix.shape))

    # The algorithm expects more columns than rows in the graph, so
    # we use the transpose if that is not already the case. We also
    # ensure that we have a full matching. In principle, it should be
    # possible to avoid this check for a performance improvement, by
    # checking for infeasibility during the execution of _lapjvsp below
    # instead, but some cases are not yet handled there.
    if j < i:
        biadjacency_matrix_t = biadjacency_matrix.T
        if biadjacency_matrix_t.format != "csr":
            biadjacency_matrix_t = biadjacency_matrix_t.tocsr()
        matching, _ = _hopcroft_karp(biadjacency_matrix_t.indices,
                                     biadjacency_matrix_t.indptr,
                                     j, i)
        matching = np.asarray(matching)
        if np.sum(matching != -1) != min(i, j):
            raise ValueError('no full matching exists')
        b = np.asarray(_lapjvsp(biadjacency_matrix_t.indptr,
                                biadjacency_matrix_t.indices,
                                biadjacency_matrix_t.data,
                                j, i))
        indices = np.argsort(b)
        return (b[indices], a[indices])
    else:
        if biadjacency_matrix.format != "csr":
            biadjacency_matrix = biadjacency_matrix.tocsr()
        matching, _ = _hopcroft_karp(biadjacency_matrix.indices,
                                     biadjacency_matrix.indptr,
                                     i, j)
        matching = np.asarray(matching)
        if np.sum(matching != -1) != min(i, j):
            raise ValueError('no full matching exists')
        b = np.asarray(_lapjvsp(biadjacency_matrix.indptr,
                                biadjacency_matrix.indices,
                                biadjacency_matrix.data,
                                i, j))
        return (a, b)


# We will use uint8 to represent booleans to simplify arrays of booleans below.
BTYPE = np.uint8
ctypedef np.uint8_t BTYPE_t


@cython.boundscheck(False)
@cython.wraparound(False)
cdef ITYPE_t[:] _lapjvsp(ITYPE_t[:] first,
                         ITYPE_t[:] kk,
                         DTYPE_t[:] cc,
                         ITYPE_t nr,
                         ITYPE_t nc) noexcept:
    """Solves the minimum weight bipartite matching problem using LAPJVsp.

    The implementation at hand is a straightforward port of the original Pascal
    implementation by Anton Volgenant. The code is currently available as
    LAPJVS.P on [1]_ and is licensed under the BSD-3 license, copyright
    A. Volgenant/Amsterdam School of Economics, University of Amsterdam.

    The original code is uncommented, and is easiest to follow by referring to
    [2]_ rather than the code itself. As such, in our Cython port, we strive to
    stay as close to the original implementation as possible, and prefer code
    that is easy to recognize from its Pascal version over optimized code. When
    we do deviate from the original implementation in non-obvious ways, we
    explicitly spell out how and why. In general, we have updated all indices
    (since Pascal is 1-indexed), introduced a few more variables to highlight
    the scope of some variables, converted the goto statements to functions,
    and use double infinity instead of large integers.

    Parameters
    ----------
    first : memoryview of length :math:`|U| + 1`
        Corresponds to the ``indptr`` attribute of the graph in CSR format.
    kk : memoryview
        Corresponds to the ``indices`` attribute of the graph in CSR format.
    cc : memoryview
        Corresponds to the ``data`` attribute of the graph in CSR format.
    nr : int
        The number of rows of the graph in CSR format.
    nc : int
        The number of columns of the graph in CSR format. Must be at least as
        large as ``nr``.

    Returns
    -------
    row_ind : memoryview of length :math:`|U|`
        For each row, the column matched with that row.

    References
    ----------
    .. [1] http://www.assignmentproblems.com/LAPJV.htm
    .. [2] Roy Jonker and Anton Volgenant:
       A Shortest Augmenting Path Algorithm for Dense and Sparse Linear
       Assignment Problems.
       Computing 38:325-340, 1987.

    """
    cdef ITYPE_t l0, jp, t, i, lp, j1, tp, j0p, j1p, l0p, h, i0, td1
    cdef DTYPE_t min_diff, v0, vj, dj
    cdef DTYPE_t[:] v = np.zeros(nc, dtype=DTYPE)
    cdef ITYPE_t[:] x = np.empty(nr, dtype=ITYPE)
    for i in range(nr):
        x[i] = -1
    cdef ITYPE_t[:] y = np.empty(nc, dtype=ITYPE)
    for j in range(nc):
        y[j] = -1
    cdef DTYPE_t[:] u = np.zeros(nr, dtype=DTYPE)
    cdef DTYPE_t[:] d = np.zeros(nc, dtype=DTYPE)
    cdef BTYPE_t[:] ok = np.zeros(nc, dtype=BTYPE)
    cdef BTYPE_t[:] xinv = np.zeros(nr, dtype=BTYPE)
    cdef ITYPE_t[:] free = np.empty(nr, dtype=ITYPE)
    for i in range(nr):
        free[i] = -1
    cdef ITYPE_t[:] todo = np.empty(nc, dtype=ITYPE)
    for j in range(nc):
        todo[j] = -1
    cdef ITYPE_t[:] lab = np.zeros(nc, dtype=ITYPE)

    # We skip the initialization entirely in the non-square case and instead
    # fill all of `free` explicitly.
    if nr == nc:
        # From here, proceed from line 55 in the original Pascal code, LAPJVS.P
        for jp in range(nc):
            v[jp] = INFINITY
        for i in range(nr):
            for t in range(first[i], first[i + 1]):
                jp = kk[t]
                if cc[t] < v[jp]:
                    v[jp] = cc[t]
                    y[jp] = i
        for jp in range(nc - 1, -1, -1):
            i = y[jp]
            # If no row has been matched with column jp at this point, that
            # can only mean that the column has no incident rows at all.
            if i == -1:
                raise ValueError('no full matching exists')
            if x[i] == -1:
                x[i] = jp
            else:
                y[jp] = -1
                # Here, the original Pascal code simply inverts the sign of
                # x[i]; as that doesn't play too well with zero-indexing, we
                # explicitly keep track of uniqueness instead.
                xinv[i] = 1
        lp = 0
        for i in range(nr):
            if xinv[i] == 1:
                continue
            if x[i] != -1:
                min_diff = INFINITY
                j1 = x[i]
                for t in range(first[i], first[i + 1]):
                    jp = kk[t]
                    if jp != j1:
                        if cc[t] - v[jp] < min_diff:
                            min_diff = cc[t] - v[jp]
                u[i] = min_diff
                tp = first[i]
                while kk[tp] != j1:
                    tp += 1
                v[j1] = cc[tp] - min_diff
            else:
                # The following two lines are swapped in the Pascal code: This
                # works because of the implicit index correction, where we
                # initialize lp to 0 as opposed to the -1 that would otherwise
                # correspond to the 0 initialization in Pascal (as we recall
                # that all indices are shifted by one).
                free[lp] = i
                lp += 1
        for _ in range(2):
            h = 0
            l0p = lp
            lp = 0
            while h < l0p:
                i = free[h]
                h += 1
                # Note: In the original Pascal code, the indices of the lowest
                # and second-lowest reduced costs are never reset. This can
                # cause issues for infeasible problems; see
                # https://stackoverflow.com/q/62875232/5085211
                j0p = -1
                j1p = -1
                v0 = INFINITY
                vj = INFINITY
                for t in range(first[i], first[i + 1]):
                    jp = kk[t]
                    dj = cc[t] - v[jp]
                    if dj < vj:
                        if dj >= v0:
                            vj = dj
                            j1p = jp
                        else:
                            vj = v0
                            v0 = dj
                            j1p = j0p
                            j0p = jp
                # If the index of the column with the largest reduced cost has
                # not been set, no assignment is possible for this row.
                if j0p < 0:
                    raise ValueError('no full matching exists')
                i0 = y[j0p]
                u[i] = vj
                if v0 < vj:
                    v[j0p] += v0 - vj
                elif i0 != -1:
                    j0p = j1p
                    i0 = y[j0p]
                x[i] = j0p
                y[j0p] = i
                if i0 != -1:
                    if v0 < vj:
                        h -= 1
                        free[h] = i0
                    else:
                        free[lp] = i0
                        lp += 1
        l0 = lp
    else:
        l0 = nr
        for i in range(nr):
            free[i] = i

    # In the original Pascal code, the solution for each l is inlined,
    # but to avoid goto statements, and convoluted logic to break out
    # of nested loops, we extract the body of the loop instead. The
    # function thus corresponds to lines 109--154 in the Pascal code,
    # with lines 155 and 156 seperated into two separate functions below.
    td1 = -1
    for l in range(l0):
        td1 = _lapjvsp_single_l(l, nc, d, ok, free, first, kk, cc, v, lab,
                                todo, y, x, td1)
    return x


@cython.boundscheck(False)
@cython.wraparound(False)
cdef ITYPE_t _lapjvsp_single_l(ITYPE_t l, ITYPE_t nc, DTYPE_t[:] d,
                               BTYPE_t[:] ok, ITYPE_t[:] free,
                               ITYPE_t[:] first, ITYPE_t[:] kk,
                               DTYPE_t[:] cc, DTYPE_t[:] v, ITYPE_t[:] lab,
                               ITYPE_t[:] todo, ITYPE_t[:] y, ITYPE_t[:] x,
                               ITYPE_t td1) except -2:
    cdef ITYPE_t jp, i0, j, t, td2, hp, last, j0, i, tp
    cdef DTYPE_t min_diff, dj, h, vj

    for jp in range(nc):
        d[jp] = INFINITY
        ok[jp] = 0
    min_diff = INFINITY
    i0 = free[l]

    for t in range(first[i0], first[i0 + 1]):
        j = kk[t]
        dj = cc[t] - v[j]
        d[j] = dj
        lab[j] = i0
        if dj <= min_diff:
            if dj < min_diff:
                td1 = -1
                min_diff = dj
            td1 += 1
            todo[td1] = j

    for hp in range(td1 + 1):
        j = todo[hp]
        if y[j] == -1:
            _lapjvsp_update_assignments(lab, y, x, j, i0)
            return td1
        ok[j] = 1
    td2 = nc - 1
    last = nc

    while True:
        # If td1 is negative at this point, that corresponds to no assignments
        # having been made in the previous run, so no full matching exists
        # and we error out.
        if td1 < 0:
            raise ValueError('no full matching exists')
        j0 = todo[td1]
        td1 -= 1
        i = y[j0]
        todo[td2] = j0
        td2 -= 1
        tp = first[i]
        while kk[tp] != j0:
            tp += 1
        h = cc[tp] - v[j0] - min_diff

        for t in range(first[i], first[i + 1]):
            j = kk[t]
            if ok[j] == 0:
                vj = cc[t] - v[j] - h
                if vj < d[j]:
                    d[j] = vj
                    lab[j] = i
                    if vj == min_diff:
                        if y[j] == -1:
                            _lapjvsp_update_dual(nc, d, v, todo,
                                                 last, min_diff)
                            _lapjvsp_update_assignments(lab, y, x, j, i0)
                            return td1
                        td1 += 1
                        todo[td1] = j
                        ok[j] = 1

        if td1 == -1:
            # The original Pascal code uses finite numbers instead of INFINITY
            # so we need to adjust slightly here.
            min_diff = INFINITY
            last = td2 + 1

            for jp in range(nc):
                if d[jp] != INFINITY and d[jp] <= min_diff and ok[jp] == 0:
                    if d[jp] < min_diff:
                        td1 = -1
                        min_diff = d[jp]
                    td1 += 1
                    todo[td1] = jp

            for hp in range(td1 + 1):
                j = todo[hp]
                if y[j] == -1:
                    _lapjvsp_update_dual(nc, d, v, todo, last, min_diff)
                    _lapjvsp_update_assignments(lab, y, x, j, i0)
                    return td1
                ok[j] = 1


@cython.boundscheck(False)
@cython.wraparound(False)
cdef _lapjvsp_update_dual(ITYPE_t nc, DTYPE_t[:] d, DTYPE_t[:] v,
                          ITYPE_t[:] todo, ITYPE_t last, DTYPE_t min_diff):
    cdef ITYPE_t j0
    for k in range(last, nc):
        j0 = todo[k]
        v[j0] += d[j0] - min_diff


@cython.boundscheck(False)
@cython.wraparound(False)
cdef _lapjvsp_update_assignments(ITYPE_t[:] lab, ITYPE_t[:] y, ITYPE_t[:] x,
                                 ITYPE_t j, ITYPE_t i0):
    cdef ITYPE_t i, k
    while True:
        i = lab[j]
        y[j] = i
        k = j
        j = x[i]
        x[i] = k
        if i == i0:
            return
