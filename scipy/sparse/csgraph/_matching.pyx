cimport cython
import numpy as np
cimport numpy as np


from scipy.sparse import (csr_matrix,
                          isspmatrix_coo, isspmatrix_csc, isspmatrix_csr)

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
           2.4 (1973), pp. 225--231. <https://dx.doi.org/10.1137/0202019>.

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
    if isspmatrix_csc(graph) or isspmatrix_coo(graph):
        graph = graph.tocsr()
    elif not isspmatrix_csr(graph):
        raise TypeError("graph must be in CSC, CSR, or COO format.")
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
