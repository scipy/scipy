"""
Routines for performing shortest-path graph searches

The main interface is in the function :func:`shortest_path`.  This
calls cython routines that compute the shortest path using
the Floyd-Warshall algorithm, Dijkstra's algorithm with priority queue,
the Bellman-Ford algorithm, or Johnson's Algorithm.

Yen's k-Shortest Path Algorithm is available for
finding the k-shortest paths between two nodes in a graph.
"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) 2011
import warnings

import numpy as np
cimport numpy as np

from scipy.sparse import csr_array, issparse
from scipy.sparse.csgraph._validation import validate_graph
from scipy.sparse._sputils import (convert_pydata_sparse_to_scipy,
                                   safely_cast_index_arrays)

cimport cython

from libc.math cimport INFINITY

from libcpp.algorithm cimport lower_bound
from libcpp.queue cimport priority_queue
from libcpp.pair cimport pair
from libcpp.vector cimport vector

np.import_array()

include 'parameters.pxi'

# EPS is the precision of DTYPE (float64, from parameters.pxi)
DEF DTYPE_EPS = 1E-15


class NegativeCycleError(Exception):
    def __init__(self, message=''):
        Exception.__init__(self, message)


def shortest_path(csgraph, method='auto',
                  directed=True,
                  return_predecessors=False,
                  unweighted=False,
                  overwrite=False,
                  indices=None):
    """
    shortest_path(csgraph, method='auto', directed=True, return_predecessors=False,
                  unweighted=False, overwrite=False, indices=None)

    Perform a shortest-path graph search on a positive directed or
    undirected graph.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : array_like, or sparse array or matrix, 2 dimensions
        The N x N array of distances representing the input graph.
    method : string ['auto'|'FW'|'D'], optional
        Algorithm to use for shortest paths.  Options are:

           'auto' -- (default) select the best among 'FW', 'D', 'BF', or 'J'
                     based on the input data.

           'FW'   -- Floyd-Warshall algorithm.
                     Computational cost is approximately ``O[N^3]``.
                     The input csgraph will be converted to a dense representation.

           'D'    -- Dijkstra's algorithm with priority queue.
                     Computational cost is approximately ``O[I * (E + N) * log(N)]``,
                     where ``E`` is the number of edges in the graph,
                     and ``I = len(indices)`` if ``indices`` is passed. Otherwise,
                     ``I = N``.
                     The input csgraph will be converted to a csr representation.

           'BF'   -- Bellman-Ford algorithm.
                     This algorithm can be used when weights are negative.
                     If a negative cycle is encountered, an error will be raised.
                     Computational cost is approximately ``O[N(N^2 k)]``, where
                     ``k`` is the average number of connected edges per node.
                     The input csgraph will be converted to a csr representation.

           'J'    -- Johnson's algorithm.
                     Like the Bellman-Ford algorithm, Johnson's algorithm is
                     designed for use when the weights are negative. It combines
                     the Bellman-Ford algorithm with Dijkstra's algorithm for
                     faster computation.

    directed : bool, optional
        If True (default), then find the shortest path on a directed graph:
        only move from point i to point j along paths csgraph[i, j].
        If False, then find the shortest path on an undirected graph: the
        algorithm can progress from point i to j along csgraph[i, j] or
        csgraph[j, i]
    return_predecessors : bool, optional
        If True, return the size (N, N) predecessor matrix.
    unweighted : bool, optional
        If True, then find unweighted distances.  That is, rather than finding
        the path between each point such that the sum of weights is minimized,
        find the path such that the number of edges is minimized.
    overwrite : bool, optional
        If True, overwrite csgraph with the result.  This applies only if
        method == 'FW' and csgraph is a dense, c-ordered array with
        dtype=float64.
    indices : array_like or int, optional
        If specified, only compute the paths from the points at the given
        indices. Incompatible with method == 'FW'.

    Returns
    -------
    dist_matrix : ndarray
        The N x N matrix of distances between graph nodes. dist_matrix[i,j]
        gives the shortest distance from point i to point j along the graph.
    predecessors : ndarray, shape (n_indices, n_nodes,)
        Returned only if return_predecessors == True.
        If `indices` is None then ``n_indices = n_nodes`` and the shape of
        the matrix becomes ``(n_nodes, n_nodes)``.
        The matrix of predecessors, which can be used to reconstruct
        the shortest paths.  Row i of the predecessor matrix contains
        information on the shortest paths from point i: each entry
        predecessors[i, j] gives the index of the previous node in the
        path from point i to point j.  If no path exists between point
        i and j, then predecessors[i, j] = -9999

    Raises
    ------
    NegativeCycleError:
        if there are negative cycles in the graph

    See Also
    --------
    :ref:`word-ladders-example` : An illustratation of the ``shortest_path`` API with a meaninful example.
                                  It also reconstructs the shortest path by using predecessors matrix returned
                                  by this function.

    Notes
    -----
    As currently implemented, Dijkstra's algorithm and Johnson's algorithm
    do not work for graphs with direction-dependent distances when
    directed == False.  i.e., if csgraph[i,j] and csgraph[j,i] are non-equal
    edges, method='D' may yield an incorrect result.

    If multiple valid solutions are possible, output may vary with SciPy and
    Python version.

    Examples
    --------
    >>> from scipy.sparse import csr_array
    >>> from scipy.sparse.csgraph import shortest_path

    >>> graph = [
    ... [0, 0, 7, 0],
    ... [0, 0, 8, 5],
    ... [7, 8, 0, 0],
    ... [0, 5, 0, 0]
    ... ]
    >>> graph = csr_array(graph)
    >>> print(graph)
    <Compressed Sparse Row sparse array of dtype 'int64'
    	with 6 stored elements and shape (4, 4)>
    	Coords	Values
    	(0, 2)	7
    	(1, 2)	8
    	(1, 3)	5
    	(2, 0)	7
    	(2, 1)	8
    	(3, 1)	5

    >>> sources = [0, 2]
    >>> dist_matrix, predecessors = shortest_path(csgraph=graph, directed=False, indices=sources, return_predecessors=True)
    >>> dist_matrix
    array([[ 0., 15.,  7., 20.],
           [ 7.,  8.,  0., 13.]])
    >>> predecessors
    array([[-9999,     2,     0,     1],
           [    2,     2, -9999,     1]], dtype=int32)

    Reconstructing shortest paths from sources to all the nodes of the graph.

    >>> shortest_paths = {}
    >>> for idx in range(len(sources)):
    ...     for node in range(4):
    ...         curr_node = node # start from the destination node
    ...         path = []
    ...         while curr_node != -9999: # no previous node available, exit the loop
    ...             path = [curr_node] + path # prefix the previous node obtained from the last iteration
    ...             curr_node = int(predecessors[idx][curr_node]) # set current node to previous node
    ...         shortest_paths[(sources[idx], node)] = path
    ...

    Computing the length of the shortest path from node 0 to node 3
    of the graph. It can be observed that computed length and the
    ``dist_matrix`` value are exactly same.

    >>> shortest_paths[(0, 3)]
    [0, 2, 1, 3]
    >>> path03 = shortest_paths[(0, 3)]
    >>> sum([graph[path03[0], path03[1]], graph[path03[1], path03[2]], graph[path03[2], path03[3]]])
    np.int64(20)
    >>> dist_matrix[0][3]
    np.float64(20.0)

    Another example of computing shortest path length from node 2 to node 3.
    Here, ``dist_matrix[1][3]`` is used to get the length of the path returned by
    ``shortest_path``. This is because node 2 is the second source, so the
    lengths of the path from it to other nodes in the graph will be at index 1
    in ``dist_matrix``.

    >>> shortest_paths[(2, 3)]
    [2, 1, 3]
    >>> path23 = shortest_paths[(2, 3)]
    >>> sum([graph[path23[0], path23[1]], graph[path23[1], path23[2]]])
    np.int64(13)
    >>> dist_matrix[1][3]
    np.float64(13.0)

    """
    csgraph = convert_pydata_sparse_to_scipy(csgraph, accept_fv=[0, np.inf, np.nan])

    # validate here to catch errors early but don't store the result;
    # we'll validate again later
    validate_graph(csgraph, directed, DTYPE,
                   copy_if_dense=(not overwrite),
                   copy_if_sparse=(not overwrite))

    cdef bint is_sparse
    cdef ssize_t N      # XXX cdef ssize_t Nk fails in Python 3 (?)

    if method == 'auto':
        # guess fastest method based on number of nodes and edges
        N = csgraph.shape[0]
        csgraph = convert_pydata_sparse_to_scipy(csgraph)
        is_sparse = issparse(csgraph)
        if is_sparse:
            Nk = csgraph.nnz
            if csgraph.format in ('csr', 'csc', 'coo'):
                edges = csgraph.data
            else:
                edges = csgraph.tocoo().data
        elif np.ma.isMaskedArray(csgraph):
            Nk = csgraph.count()
            edges = csgraph.compressed()
        else:
            edges = csgraph[np.isfinite(csgraph)]
            edges = edges[edges != 0]
            Nk = edges.size

        if indices is not None or Nk < N * N / 4:
            if np.any(edges < 0):
                method = 'J'
            else:
                method = 'D'
        else:
            method = 'FW'

    if method == 'FW':
        if indices is not None:
            raise ValueError("Cannot specify indices with method == 'FW'.")
        return floyd_warshall(csgraph, directed,
                              return_predecessors=return_predecessors,
                              unweighted=unweighted,
                              overwrite=overwrite)

    elif method == 'D':
        return dijkstra(csgraph, directed,
                        return_predecessors=return_predecessors,
                        unweighted=unweighted, indices=indices)

    elif method == 'BF':
        return bellman_ford(csgraph, directed,
                            return_predecessors=return_predecessors,
                            unweighted=unweighted, indices=indices)

    elif method == 'J':
        return johnson(csgraph, directed,
                       return_predecessors=return_predecessors,
                       unweighted=unweighted, indices=indices)

    else:
        raise ValueError("unrecognized method '%s'" % method)


def floyd_warshall(csgraph, directed=True,
                   return_predecessors=False,
                   unweighted=False,
                   overwrite=False):
    """
    floyd_warshall(csgraph, directed=True, return_predecessors=False,
                   unweighted=False, overwrite=False)

    Compute the shortest path lengths using the Floyd-Warshall algorithm

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : array_like, or sparse array or matrix, 2 dimensions
        The N x N array of distances representing the input graph.
    directed : bool, optional
        If True (default), then find the shortest path on a directed graph:
        only move from point i to point j along paths csgraph[i, j].
        If False, then find the shortest path on an undirected graph: the
        algorithm can progress from point i to j along csgraph[i, j] or
        csgraph[j, i]
    return_predecessors : bool, optional
        If True, return the size (N, N) predecessor matrix.
    unweighted : bool, optional
        If True, then find unweighted distances.  That is, rather than finding
        the path between each point such that the sum of weights is minimized,
        find the path such that the number of edges is minimized.
    overwrite : bool, optional
        If True, overwrite csgraph with the result.  This applies only if
        csgraph is a dense, c-ordered array with dtype=float64.

    Returns
    -------
    dist_matrix : ndarray
        The N x N matrix of distances between graph nodes. dist_matrix[i,j]
        gives the shortest distance from point i to point j along the graph.

    predecessors : ndarray
        Returned only if return_predecessors == True.
        The N x N matrix of predecessors, which can be used to reconstruct
        the shortest paths.  Row i of the predecessor matrix contains
        information on the shortest paths from point i: each entry
        predecessors[i, j] gives the index of the previous node in the
        path from point i to point j.  If no path exists between point
        i and j, then predecessors[i, j] = -9999

    Raises
    ------
    NegativeCycleError:
        if there are negative cycles in the graph

    Notes
    -----
    If multiple valid solutions are possible, output may vary with SciPy and
    Python version.

    Examples
    --------
    >>> from scipy.sparse import csr_array
    >>> from scipy.sparse.csgraph import floyd_warshall

    >>> graph = [
    ... [0, 1, 2, 0],
    ... [0, 0, 0, 1],
    ... [2, 0, 0, 3],
    ... [0, 0, 0, 0]
    ... ]
    >>> graph = csr_array(graph)
    >>> print(graph)
    <Compressed Sparse Row sparse array of dtype 'int64'
    	with 5 stored elements and shape (4, 4)>
    	Coords	Values
    	(0, 1)	1
    	(0, 2)	2
    	(1, 3)	1
    	(2, 0)	2
    	(2, 3)	3

    >>> dist_matrix, predecessors = floyd_warshall(csgraph=graph, directed=False, return_predecessors=True)
    >>> dist_matrix
    array([[0., 1., 2., 2.],
           [1., 0., 3., 1.],
           [2., 3., 0., 3.],
           [2., 1., 3., 0.]])
    >>> predecessors
    array([[-9999,     0,     0,     1],
           [    1, -9999,     0,     1],
           [    2,     0, -9999,     2],
           [    1,     3,     3, -9999]], dtype=int32)

    """
    dist_matrix = validate_graph(csgraph, directed, DTYPE,
                                 csr_output=False,
                                 copy_if_dense=not overwrite)
    if not issparse(csgraph):
        # for dense array input, zero entries represent non-edge
        dist_matrix[dist_matrix == 0] = INFINITY

    if unweighted:
        dist_matrix[~np.isinf(dist_matrix)] = 1

    if return_predecessors:
        predecessor_matrix = np.empty(dist_matrix.shape,
                                      dtype=ITYPE, order='C')
    else:
        predecessor_matrix = np.empty((0, 0), dtype=ITYPE)

    _floyd_warshall(dist_matrix,
                    predecessor_matrix,
                    int(directed))

    if np.any(dist_matrix.diagonal() < 0):
        raise NegativeCycleError("Negative cycle in nodes %s"
                                 % np.where(dist_matrix.diagonal() < 0)[0])

    if return_predecessors:
        return dist_matrix, predecessor_matrix
    else:
        return dist_matrix


@cython.boundscheck(False)
cdef void _floyd_warshall(
               np.ndarray[DTYPE_t, ndim=2, mode='c'] dist_matrix,
               np.ndarray[ITYPE_t, ndim=2, mode='c'] predecessor_matrix,
               int directed=0) noexcept:
    # dist_matrix : in/out
    #    on input, the graph
    #    on output, the matrix of shortest paths
    # dist_matrix should be a [N,N] matrix, such that dist_matrix[i, j]
    # is the distance from point i to point j.  Zero-distances imply that
    # the points are not connected.
    cdef unsigned int N = dist_matrix.shape[0]
    assert dist_matrix.shape[1] == N

    cdef unsigned int i, j, k

    cdef DTYPE_t d_ijk

    # ----------------------------------------------------------------------
    #  Initialize distance matrix
    #   - set diagonal to zero
    #   - symmetrize matrix if non-directed graph is desired
    dist_matrix.flat[::N + 1] = 0
    if not directed:
        for i in range(N):
            for j in range(i + 1, N):
                if dist_matrix[j, i] <= dist_matrix[i, j]:
                    dist_matrix[i, j] = dist_matrix[j, i]
                else:
                    dist_matrix[j, i] = dist_matrix[i, j]

    #----------------------------------------------------------------------
    #  Initialize predecessor matrix
    #   - check matrix size
    #   - initialize diagonal and all non-edges to NULL
    #   - initialize all edges to the row index
    cdef int store_predecessors = False

    if predecessor_matrix.size > 0:
        store_predecessors = True
        assert predecessor_matrix.shape[0] == N
        assert predecessor_matrix.shape[1] == N
        predecessor_matrix.fill(NULL_IDX)
        i_edge = np.where(~np.isinf(dist_matrix))
        predecessor_matrix[i_edge] = i_edge[0]
        predecessor_matrix.flat[::N + 1] = NULL_IDX

    # Now perform the Floyd-Warshall algorithm.
    # In each loop, this finds the shortest path from point i
    #  to point j using intermediate nodes 0 ... k
    if store_predecessors:
        for k in range(N):
            for i in range(N):
                if dist_matrix[i, k] == INFINITY:
                    continue
                for j in range(N):
                    d_ijk = dist_matrix[i, k] + dist_matrix[k, j]
                    if d_ijk < dist_matrix[i, j]:
                        dist_matrix[i, j] = d_ijk
                        predecessor_matrix[i, j] = predecessor_matrix[k, j]
    else:
        for k in range(N):
            for i in range(N):
                if dist_matrix[i, k] == INFINITY:
                    continue
                for j in range(N):
                    d_ijk = dist_matrix[i, k] + dist_matrix[k, j]
                    if d_ijk < dist_matrix[i, j]:
                        dist_matrix[i, j] = d_ijk


def dijkstra(csgraph, directed=True, indices=None,
             return_predecessors=False,
             unweighted=False, limit=np.inf,
             bint min_only=False):
    """
    dijkstra(csgraph, directed=True, indices=None, return_predecessors=False,
             unweighted=False, limit=np.inf, min_only=False)

    Dijkstra algorithm using priority queue

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : array_like, or sparse array or matrix, 2 dimensions
        The N x N array of non-negative distances representing the input graph.
    directed : bool, optional
        If True (default), then find the shortest path on a directed graph:
        only move from point i to point j along paths csgraph[i, j] and from
        point j to i along paths csgraph[j, i].
        If False, then find the shortest path on an undirected graph: the
        algorithm can progress from point i to j or j to i along either
        csgraph[i, j] or csgraph[j, i].

        .. warning:: Refer the notes below while using with ``directed=False``.
    indices : array_like or int, optional
        if specified, only compute the paths from the points at the given
        indices.
    return_predecessors : bool, optional
        If True, return the size (N, N) predecessor matrix.
    unweighted : bool, optional
        If True, then find unweighted distances.  That is, rather than finding
        the path between each point such that the sum of weights is minimized,
        find the path such that the number of edges is minimized.
    limit : float, optional
        The maximum distance to calculate, must be >= 0. Using a smaller limit
        will decrease computation time by aborting calculations between pairs
        that are separated by a distance > limit. For such pairs, the distance
        will be equal to np.inf (i.e., not connected).

        .. versionadded:: 0.14.0
    min_only : bool, optional
        If False (default), for every node in the graph, find the shortest path
        from every node in indices.
        If True, for every node in the graph, find the shortest path from any
        of the nodes in indices (which can be substantially faster).

        .. versionadded:: 1.3.0

    Returns
    -------
    dist_matrix : ndarray, shape ([n_indices, ]n_nodes,)
        The matrix of distances between graph nodes. If min_only=False,
        dist_matrix has shape (n_indices, n_nodes) and dist_matrix[i, j]
        gives the shortest distance from point i to point j along the graph.
        If min_only=True, dist_matrix has shape (n_nodes,) and contains for
        a given node the shortest path to that node from any of the nodes
        in indices.
    predecessors : ndarray, shape ([n_indices, ]n_nodes,)
        If ``min_only=False``, this has shape ``(n_indices, n_nodes)``,
        otherwise it has shape ``(n_nodes,)``.
        If `indices` is None and ``min_only=False`` then ``n_indices = n_nodes``
        and the shape of the matrix becomes ``(n_nodes, n_nodes)``.
        Returned only if return_predecessors == True.
        The matrix of predecessors, which can be used to reconstruct
        the shortest paths.  Row i of the predecessor matrix contains
        information on the shortest paths from point i: each entry
        predecessors[i, j] gives the index of the previous node in the
        path from point i to point j.  If no path exists between point
        i and j, then predecessors[i, j] = -9999

    sources : ndarray, shape (n_nodes,)
        Returned only if min_only=True and return_predecessors=True.
        Contains the index of the source which had the shortest path
        to each target.  If no path exists within the limit,
        this will contain -9999.  The value at the indices passed
        will be equal to that index (i.e. the fastest way to reach
        node i, is to start on node i).

    Notes
    -----
    As currently implemented, Dijkstra's algorithm does not work for
    graphs with direction-dependent distances when directed == False.
    i.e., if csgraph[i,j] and csgraph[j,i] are not equal and
    both are nonzero, setting directed=False will not yield the correct
    result.

    Also, this routine does not work for graphs with negative
    distances.  Negative distances can lead to infinite cycles that must
    be handled by specialized algorithms such as Bellman-Ford's algorithm
    or Johnson's algorithm.

    If multiple valid solutions are possible, output may vary with SciPy and
    Python version.

    Examples
    --------
    >>> from scipy.sparse import csr_array
    >>> from scipy.sparse.csgraph import dijkstra

    >>> graph = [
    ... [0, 1, 2, 0],
    ... [0, 0, 0, 1],
    ... [0, 0, 0, 3],
    ... [0, 0, 0, 0]
    ... ]
    >>> graph = csr_array(graph)
    >>> print(graph)
    <Compressed Sparse Row sparse array of dtype 'int64'
    	with 4 stored elements and shape (4, 4)>
    	Coords	Values
    	(0, 1)	1
    	(0, 2)	2
    	(1, 3)	1
    	(2, 3)	3

    >>> dist_matrix, predecessors = dijkstra(csgraph=graph, directed=False, indices=0, return_predecessors=True)
    >>> dist_matrix
    array([0., 1., 2., 2.])
    >>> predecessors
    array([-9999,     0,     0,     1], dtype=int32)

    """
    #------------------------------
    # validate csgraph and convert to csr
    csgraph = validate_graph(csgraph, directed, DTYPE, dense_output=False)

    if np.any(csgraph.data < 0):
        warnings.warn("Graph has negative weights: dijkstra will give "
                      "inaccurate results if the graph contains negative "
                      "cycles. Consider johnson or bellman_ford.")

    N = csgraph.shape[0]

    #------------------------------
    # initialize/validate indices
    if indices is None:
        indices = np.arange(N, dtype=ITYPE)
        if min_only:
            return_shape = (N,)
        else:
            return_shape = indices.shape + (N,)
    else:
        indices = np.array(indices, order='C', dtype=ITYPE, copy=True)
        if min_only:
            return_shape = (N,)
        else:
            return_shape = indices.shape + (N,)
        indices = np.atleast_1d(indices).reshape(-1)
        indices[indices < 0] += N
        if np.any(indices < 0) or np.any(indices >= N):
            raise ValueError("indices out of range 0...N")

    cdef DTYPE_t limitf = limit
    if limitf < 0:
        raise ValueError('limit must be >= 0')

    #------------------------------
    # initialize dist_matrix for output
    if min_only:
        dist_matrix = np.full(N, np.inf, dtype=DTYPE)
        dist_matrix[indices] = 0
    else:
        dist_matrix = np.full((len(indices), N), np.inf, dtype=DTYPE)
        dist_matrix[np.arange(len(indices)), indices] = 0

    #------------------------------
    # initialize predecessors for output
    if return_predecessors:
        if min_only:
            predecessor_matrix = np.empty((N), dtype=ITYPE)
            predecessor_matrix.fill(NULL_IDX)
            source_matrix = np.empty((N), dtype=ITYPE)
            source_matrix.fill(NULL_IDX)
        else:
            predecessor_matrix = np.empty((len(indices), N), dtype=ITYPE)
            predecessor_matrix.fill(NULL_IDX)
            source_matrix = np.empty((len(indices), 0), dtype=ITYPE) # unused
    else:
        if min_only:
            predecessor_matrix = np.empty(0, dtype=ITYPE)
            source_matrix = np.empty(0, dtype=ITYPE) # unused
        else:
            predecessor_matrix = np.empty((len(indices), 0), dtype=ITYPE)
            source_matrix = np.empty((len(indices), 0), dtype=ITYPE) # unused

    if unweighted:
        csr_data = np.ones(csgraph.data.shape)
    else:
        csr_data = csgraph.data
    csr_indices, csr_indptr = safely_cast_index_arrays(csgraph, ITYPE, msg="csgraph")

    if directed:
        # for null transposed CSR
        dummy_double_array = np.empty(0, dtype=DTYPE)
        dummy_int_array = np.empty(0, dtype=ITYPE)
        if min_only:
            _dijkstra(indices,
                      csr_data, csr_indices, csr_indptr,
                      dummy_double_array, dummy_int_array, dummy_int_array,
                      dist_matrix, predecessor_matrix, source_matrix,
                      limitf)
        else:
            _dijkstra_multi_separate(
                      indices,
                      csr_data, csr_indices, csr_indptr,
                      dummy_double_array, dummy_int_array, dummy_int_array,
                      dist_matrix, predecessor_matrix, source_matrix,
                      limitf)

    else:
        csrT = csgraph.T.tocsr()
        csrT_indices, csrT_indptr = safely_cast_index_arrays(csrT, ITYPE, msg="csgraph")
        if unweighted:
            csrT_data = csr_data
        else:
            csrT_data = csrT.data
        if min_only:
            _dijkstra(indices,
                    csr_data, csr_indices, csr_indptr,
                    csrT_data, csrT_indices, csrT_indptr,
                    dist_matrix, predecessor_matrix, source_matrix,
                    limitf)
        else:
            _dijkstra_multi_separate(
                                 indices,
                                 csr_data, csr_indices, csr_indptr,
                                 csrT_data, csrT_indices, csrT_indptr,
                                 dist_matrix, predecessor_matrix, source_matrix,
                                 limitf)

    if return_predecessors:
        if min_only:
            return (dist_matrix.reshape(return_shape),
                    predecessor_matrix.reshape(return_shape),
                    source_matrix.reshape(return_shape))
        else:
            return (dist_matrix.reshape(return_shape),
                    predecessor_matrix.reshape(return_shape))
    else:
        return dist_matrix.reshape(return_shape)


ctypedef unsigned int uint_t
ctypedef pair[DTYPE_t, uint_t] dist_index_pair_t
ctypedef priority_queue[dist_index_pair_t] dijkstra_queue_t

@cython.boundscheck(False)
cdef void _dijkstra_scan_heap(dijkstra_queue_t &heap,
                         dist_index_pair_t v,
                         const double[:] csr_weights,
                         const int[:] csr_indices,
                         const int[:] csr_indptr,
                         double[:] dist_matrix,
                         int[:] pred,
                         int return_pred,
                         int[:] sources,
                         int return_source,
                         DTYPE_t limit) noexcept nogil:
    cdef:
        ITYPE_t j
        unsigned int j_current
        DTYPE_t next_val

    # v is a dist_index_pair_t poped from the queue
    # v.first: the distance of the vertex
    # v.second: index of the vertex
    for j in range(csr_indptr[v.second], csr_indptr[v.second + 1]):
        j_current = csr_indices[j]
        next_val = v.first + csr_weights[j]
        if next_val <= limit:
            if dist_matrix[j_current] > next_val:
                dist_matrix[j_current] = next_val
                # The same vertex may be pushed multiple times to the queue, but
                # anything with suboptimal distance is ignored when poped
                heap.push(dist_index_pair_t(-next_val, j_current))
                if return_pred:
                    pred[j_current] = v.second
                if return_source:
                    sources[j_current] = sources[v.second]


@cython.boundscheck(False)
cdef int _dijkstra(
            const int[:] source_indices,
            const double[:] csr_weights,
            const int[:] csr_indices,
            const int[:] csr_indptr,
            const double[:] csrT_weights,
            const int[:] csrT_indices,
            const int[:] csrT_indptr,
            double[:] dist_matrix,
            int[:] pred,
            int[:] sources,
            DTYPE_t limit) except -1 nogil:
    cdef:
        unsigned int Nind = source_indices.shape[0]
        unsigned int N = dist_matrix.shape[0]
        unsigned int i, j_source
        bint return_pred = (pred.shape[0] > 0)
        bint return_sources = (sources.shape[0] > 0)
        bint directed = (csrT_weights.shape[0] == 0)

        # pairs of {-distance, vertex index} will be pushed
        # to treat it as a min-heap instead of max-heap
        dijkstra_queue_t heap = dijkstra_queue_t()
        dist_index_pair_t v

    if return_pred and pred.shape[0] != N:
        raise RuntimeError(
            f"Invalid predecessors array shape {pred.shape}. Expected {(N,)}."
        )
    if return_sources and sources.shape[0] != N:
        raise RuntimeError(
            f"Invalid sources array shape {sources.shape}. Expected {(N,)}."
        )

    for i in range(Nind):
        j_source = source_indices[i]
        dist_matrix[j_source] = 0
        heap.push(dist_index_pair_t(-dist_matrix[j_source], j_source))
        if return_sources:
            sources[j_source] = j_source


    while heap.size():
        v = heap.top()
        heap.pop()
        v.first = -v.first
        # Do not process v if its distance has been updated
        # after v was pushed to the queue, in which case
        # _dijkstra_scan_heap should have already been called with
        # the vertex v.second
        # This assures _dijkstra_scan_heap is only called once per vertex
        # and the total complexity is O(Mlog(M)) per source
        if dist_matrix[v.second] < v.first :
            continue

        _dijkstra_scan_heap(heap, v, csr_weights, csr_indices, csr_indptr,
                            dist_matrix, pred, return_pred,
                            sources, return_sources, limit)
        if not directed:
            _dijkstra_scan_heap(heap, v,
                                csrT_weights, csrT_indices, csrT_indptr,
                                dist_matrix, pred, return_pred,
                                sources, return_sources, limit)
    return 0

@cython.boundscheck(False)
cdef int _dijkstra_multi_separate(
            const int[:] source_indices,
            const double[:] csr_weights,
            const int[:] csr_indices,
            const int[:] csr_indptr,
            const double[:] csrT_weights,
            const int[:] csrT_indices,
            const int[:] csrT_indptr,
            double[:, :] dist_matrix,
            int[:, :] pred,
            int[:, :] sources,
            DTYPE_t limit) except -1:
    cdef:
        unsigned int Nind = source_indices.shape[0]
        unsigned int i
        int source_list[1]

    if dist_matrix.shape[0] != Nind:
        raise RuntimeError(
            f"Not enough rows in distances matrix. Got {dist_matrix.shape[0]}, expected {Nind}."
        )
    if pred.shape[0] != Nind:
        raise RuntimeError(
            f"Not enough rows in predecessors matrix. Got {pred.shape[0]}, expected {Nind}."
        )
    if sources.shape[0] != Nind:
        raise RuntimeError(
            f"Not enough rows in sources matrix. Got {sources.shape[0]}, expected {Nind}."
        )

    for i in range(Nind):
        source_list[0] = source_indices[i]
        _dijkstra(source_list,
                  csr_weights, csr_indices, csr_indptr,
                  csrT_weights, csrT_indices, csrT_indptr,
                  dist_matrix[i], pred[i], sources[i], limit)


def bellman_ford(csgraph, directed=True, indices=None,
                 return_predecessors=False,
                 unweighted=False):
    """
    bellman_ford(csgraph, directed=True, indices=None, return_predecessors=False,
                 unweighted=False)

    Compute the shortest path lengths using the Bellman-Ford algorithm.

    The Bellman-Ford algorithm can robustly deal with graphs with negative
    weights.  If a negative cycle is detected, an error is raised.  For
    graphs without negative edge weights, Dijkstra's algorithm may be faster.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : array_like, or sparse array or matrix, 2 dimensions
        The N x N array of distances representing the input graph.
    directed : bool, optional
        If True (default), then find the shortest path on a directed graph:
        only move from point i to point j along paths csgraph[i, j].
        If False, then find the shortest path on an undirected graph: the
        algorithm can progress from point i to j along csgraph[i, j] or
        csgraph[j, i]
    indices : array_like or int, optional
        if specified, only compute the paths from the points at the given
        indices.
    return_predecessors : bool, optional
        If True, return the size (N, N) predecessor matrix.
    unweighted : bool, optional
        If True, then find unweighted distances.  That is, rather than finding
        the path between each point such that the sum of weights is minimized,
        find the path such that the number of edges is minimized.

    Returns
    -------
    dist_matrix : ndarray
        The N x N matrix of distances between graph nodes. dist_matrix[i,j]
        gives the shortest distance from point i to point j along the graph.

    predecessors : ndarray, shape (n_indices, n_nodes,)
        Returned only if ``return_predecessors=True``.
        If `indices` is None then ``n_indices = n_nodes`` and the shape of
        the matrix becomes ``(n_nodes, n_nodes)``.
        The matrix of predecessors, which can be used to reconstruct
        the shortest paths.  Row i of the predecessor matrix contains
        information on the shortest paths from point i: each entry
        predecessors[i, j] gives the index of the previous node in the
        path from point i to point j.  If no path exists between point
        i and j, then predecessors[i, j] = -9999

    Raises
    ------
    NegativeCycleError:
        if there are negative cycles in the graph

    Notes
    -----
    This routine is specially designed for graphs with negative edge weights.
    If all edge weights are positive, then Dijkstra's algorithm is a better
    choice.

    If multiple valid solutions are possible, output may vary with SciPy and
    Python version.

    Examples
    --------
    >>> from scipy.sparse import csr_array
    >>> from scipy.sparse.csgraph import bellman_ford

    >>> graph = [
    ... [0, 1 ,2, 0],
    ... [0, 0, 0, 1],
    ... [2, 0, 0, 3],
    ... [0, 0, 0, 0]
    ... ]
    >>> graph = csr_array(graph)
    >>> print(graph)
    <Compressed Sparse Row sparse array of dtype 'int64'
    	with 5 stored elements and shape (4, 4)>
    	Coords	Values
    	(0, 1)	1
    	(0, 2)	2
    	(1, 3)	1
    	(2, 0)	2
    	(2, 3)	3

    >>> dist_matrix, predecessors = bellman_ford(csgraph=graph, directed=False, indices=0, return_predecessors=True)
    >>> dist_matrix
    array([0., 1., 2., 2.])
    >>> predecessors
    array([-9999,     0,     0,     1], dtype=int32)

    """
    # ------------------------------
    # validate csgraph and convert to csr matrix
    csgraph = validate_graph(csgraph, directed, DTYPE,
                             dense_output=False)
    N = csgraph.shape[0]

    # ------------------------------
    # initialize/validate indices
    if indices is None:
        indices = np.arange(N, dtype=ITYPE)
    else:
        indices = np.array(indices, order='C', dtype=ITYPE)
        indices[indices < 0] += N
        if np.any(indices < 0) or np.any(indices >= N):
            raise ValueError("indices out of range 0...N")
    return_shape = indices.shape + (N,)
    indices = np.atleast_1d(indices).reshape(-1)

    # ------------------------------
    # initialize dist_matrix for output
    dist_matrix = np.empty((len(indices), N), dtype=DTYPE)
    dist_matrix.fill(np.inf)
    dist_matrix[np.arange(len(indices)), indices] = 0

    # ------------------------------
    # initialize predecessors for output
    if return_predecessors:
        predecessor_matrix = np.empty((len(indices), N), dtype=ITYPE)
        predecessor_matrix.fill(NULL_IDX)
    else:
        predecessor_matrix = np.empty((0, N), dtype=ITYPE)

    if unweighted:
        csr_data = np.ones(csgraph.data.shape)
    else:
        csr_data = csgraph.data
    csr_indices, csr_indptr = safely_cast_index_arrays(csgraph, ITYPE, msg="csgraph")

    if directed:
        ret = _bellman_ford_directed(indices,
                                     csr_data, csr_indices, csr_indptr,
                                     dist_matrix, predecessor_matrix)
    else:
        ret = _bellman_ford_undirected(indices,
                                       csr_data, csr_indices, csr_indptr,
                                       dist_matrix, predecessor_matrix)

    if ret >= 0:
        raise NegativeCycleError("Negative cycle detected on node %i" % ret)

    if return_predecessors:
        return (dist_matrix.reshape(return_shape),
                predecessor_matrix.reshape(return_shape))
    else:
        return dist_matrix.reshape(return_shape)


cdef int _bellman_ford_directed(
            const int[:] source_indices,
            const double[:] csr_weights,
            const int[:] csr_indices,
            const int[:] csr_indptr,
            double[:, :] dist_matrix,
            int[:, :] pred) noexcept:
    cdef:
        unsigned int Nind = dist_matrix.shape[0]
        unsigned int N = dist_matrix.shape[1]
        unsigned int i, j, j_source, count
        ITYPE_t k
        DTYPE_t d1, d2, w12
        int return_pred = (pred.size > 0)

    for i in range(Nind):
        j_source = source_indices[i]

        # relax all edges N-1 times
        for count in range(N - 1):
            for j in range(N):
                d1 = dist_matrix[i, j]
                for k in range(csr_indptr[j], csr_indptr[j + 1]):
                    w12 = csr_weights[k]
                    d2 = dist_matrix[i, csr_indices[k]]
                    if d1 + w12 < d2:
                        dist_matrix[i, csr_indices[k]] = d1 + w12
                        if return_pred:
                            pred[i, csr_indices[k]] = j

        # check for negative-weight cycles
        for j in range(N):
            d1 = dist_matrix[i, j]
            for k in range(csr_indptr[j], csr_indptr[j + 1]):
                w12 = csr_weights[k]
                d2 = dist_matrix[i, csr_indices[k]]
                if d1 + w12 + DTYPE_EPS < d2:
                    return j_source

    return -1


cdef int _bellman_ford_undirected(
            const int[:] source_indices,
            const double[:] csr_weights,
            const int[:] csr_indices,
            const int[:] csr_indptr,
            double[:, :] dist_matrix,
            int[:, :] pred) noexcept:
    cdef:
        unsigned int Nind = dist_matrix.shape[0]
        unsigned int N = dist_matrix.shape[1]
        unsigned int i, j, j_source, ind_k, count
        ITYPE_t k
        DTYPE_t d1, d2, w12
        int return_pred = (pred.size > 0)

    for i in range(Nind):
        j_source = source_indices[i]

        # relax all edges N-1 times
        for count in range(N - 1):
            for j in range(N):
                d1 = dist_matrix[i, j]
                for k in range(csr_indptr[j], csr_indptr[j + 1]):
                    w12 = csr_weights[k]
                    ind_k = csr_indices[k]
                    d2 = dist_matrix[i, ind_k]
                    if d1 + w12 < d2:
                        dist_matrix[i, ind_k] = d2 = d1 + w12
                        if return_pred:
                            pred[i, ind_k] = j
                    if d2 + w12 < d1:
                        dist_matrix[i, j] = d1 = d2 + w12
                        if return_pred:
                            pred[i, j] = ind_k

        # check for negative-weight cycles
        for j in range(N):
            d1 = dist_matrix[i, j]
            for k in range(csr_indptr[j], csr_indptr[j + 1]):
                w12 = csr_weights[k]
                d2 = dist_matrix[i, csr_indices[k]]
                if abs(d2 - d1) > w12 + DTYPE_EPS:
                    return j_source

    return -1


def johnson(csgraph, directed=True, indices=None,
            return_predecessors=False,
            unweighted=False):
    """
    johnson(csgraph, directed=True, indices=None, return_predecessors=False,
            unweighted=False)

    Compute the shortest path lengths using Johnson's algorithm.

    Johnson's algorithm combines the Bellman-Ford algorithm and Dijkstra's
    algorithm to quickly find shortest paths in a way that is robust to
    the presence of negative cycles.  If a negative cycle is detected,
    an error is raised.  For graphs without negative edge weights,
    dijkstra may be faster.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : array_like, or sparse array or matrix, 2 dimensions
        The N x N array of distances representing the input graph.
    directed : bool, optional
        If True (default), then find the shortest path on a directed graph:
        only move from point i to point j along paths csgraph[i, j].
        If False, then find the shortest path on an undirected graph: the
        algorithm can progress from point i to j along csgraph[i, j] or
        csgraph[j, i]
    indices : array_like or int, optional
        if specified, only compute the paths from the points at the given
        indices.
    return_predecessors : bool, optional
        If True, return the size (N, N) predecessor matrix.
    unweighted : bool, optional
        If True, then find unweighted distances.  That is, rather than finding
        the path between each point such that the sum of weights is minimized,
        find the path such that the number of edges is minimized.

    Returns
    -------
    dist_matrix : ndarray
        The N x N matrix of distances between graph nodes. dist_matrix[i,j]
        gives the shortest distance from point i to point j along the graph.

    predecessors : ndarray, shape (n_indices, n_nodes,)
        Returned only if return_predecessors == True.
        If `indices` is None then ``n_indices = n_nodes`` and the shape of
        the matrix becomes ``(n_nodes, n_nodes)``.
        The matrix of predecessors, which can be used to reconstruct
        the shortest paths.  Row i of the predecessor matrix contains
        information on the shortest paths from point i: each entry
        predecessors[i, j] gives the index of the previous node in the
        path from point i to point j.  If no path exists between point
        i and j, then predecessors[i, j] = -9999

    Raises
    ------
    NegativeCycleError:
        if there are negative cycles in the graph

    Notes
    -----
    This routine is specially designed for graphs with negative edge weights.
    If all edge weights are positive, then Dijkstra's algorithm is a better
    choice.

    If multiple valid solutions are possible, output may vary with SciPy and
    Python version.

    Examples
    --------
    >>> from scipy.sparse import csr_array
    >>> from scipy.sparse.csgraph import johnson

    >>> graph = [
    ... [0, 1, 2, 0],
    ... [0, 0, 0, 1],
    ... [2, 0, 0, 3],
    ... [0, 0, 0, 0]
    ... ]
    >>> graph = csr_array(graph)
    >>> print(graph)
    <Compressed Sparse Row sparse array of dtype 'int64'
    	with 5 stored elements and shape (4, 4)>
    	Coords	Values
    	(0, 1)	1
    	(0, 2)	2
    	(1, 3)	1
    	(2, 0)	2
    	(2, 3)	3

    >>> dist_matrix, predecessors = johnson(csgraph=graph, directed=False, indices=0, return_predecessors=True)
    >>> dist_matrix
    array([0., 1., 2., 2.])
    >>> predecessors
    array([-9999,     0,     0,     1], dtype=int32)

    """
    # ------------------------------
    # if unweighted, there are no negative weights: we just use dijkstra
    if unweighted:
        return dijkstra(csgraph, directed, indices,
                        return_predecessors, unweighted)

    # ------------------------------
    # validate csgraph and convert to csr matrix
    csgraph = validate_graph(csgraph, directed, DTYPE,
                             dense_output=False)
    N = csgraph.shape[0]

    # ------------------------------
    # initialize/validate indices
    if indices is None:
        indices = np.arange(N, dtype=ITYPE)
        return_shape = indices.shape + (N,)
    else:
        indices = np.array(indices, order='C', dtype=ITYPE)
        return_shape = indices.shape + (N,)
        indices = np.atleast_1d(indices).reshape(-1)
        indices[indices < 0] += N
        if np.any(indices < 0) or np.any(indices >= N):
            raise ValueError("indices out of range 0...N")

    #------------------------------
    # initialize dist_matrix for output
    dist_matrix = np.empty((len(indices), N), dtype=DTYPE)
    dist_matrix.fill(np.inf)
    dist_matrix[np.arange(len(indices)), indices] = 0

    #------------------------------
    # initialize predecessors for output
    if return_predecessors:
        predecessor_matrix = np.empty((len(indices), N), dtype=ITYPE)
        predecessor_matrix.fill(NULL_IDX)
    else:
        predecessor_matrix = np.empty((len(indices), 0), dtype=ITYPE)

    #------------------------------
    # initialize distance array
    dist_array = np.zeros(N, dtype=DTYPE)

    csr_data = csgraph.data.copy()
    csr_indices, csr_indptr = safely_cast_index_arrays(csgraph, ITYPE, msg="csgraph")

    #------------------------------
    # here we first add a single node to the graph, connected by a
    # directed edge of weight zero to each node, and perform bellman-ford
    if directed:
        ret = _johnson_directed(csr_data, csr_indices, csr_indptr, dist_array)
    else:
        ret = _johnson_undirected(csr_data, csr_indices, csr_indptr, dist_array)

    if ret >= 0:
        raise NegativeCycleError("Negative cycle detected on node %i" % ret)

    #------------------------------
    # add the bellman-ford weights to the data
    _johnson_add_weights(csr_data, csr_indices, csr_indptr, dist_array)

    dummy_source_matrix = np.empty((len(indices), 0), dtype=ITYPE)
    if directed:
        # for null transposed CSR
        dummy_double_array = np.empty(0, dtype=DTYPE)
        dummy_int_array = np.empty(0, dtype=ITYPE)
        _dijkstra_multi_separate(
            indices,
            csr_data, csr_indices, csr_indptr,
            dummy_double_array, dummy_int_array, dummy_int_array,
            dist_matrix, predecessor_matrix, dummy_source_matrix, np.inf)
    else:
        csgraphT = csr_array((csr_data, csr_indices, csr_indptr),
                               csgraph.shape).T.tocsr()
        _johnson_add_weights(csgraphT.data, csgraphT.indices,
                             csgraphT.indptr, dist_array)
        _dijkstra_multi_separate(
            indices,
            csr_data,csr_indices, csr_indptr,
            csgraphT.data, csgraphT.indices, csgraphT.indptr,
            dist_matrix, predecessor_matrix, dummy_source_matrix,
            np.inf)

    # ------------------------------
    # correct the distance matrix for the bellman-ford weights
    dist_matrix += dist_array
    dist_matrix -= dist_array[:, None][indices]

    if return_predecessors:
        return (dist_matrix.reshape(return_shape),
                predecessor_matrix.reshape(return_shape))
    else:
        return dist_matrix.reshape(return_shape)


cdef void _johnson_add_weights(
            double[:] csr_weights,
            const int[:] csr_indices,
            const int[:] csr_indptr,
            const double[:] dist_array) noexcept:
    # let w(u, v) = w(u, v) + h(u) - h(v)
    cdef:
        unsigned int j, N = dist_array.shape[0]
        ITYPE_t k

    for j in range(N):
        for k in range(csr_indptr[j], csr_indptr[j + 1]):
            csr_weights[k] += dist_array[j]
            csr_weights[k] -= dist_array[csr_indices[k]]


cdef int _johnson_directed(
            const double[:] csr_weights,
            const int[:] csr_indices,
            const int[:] csr_indptr,
            double[:] dist_array) noexcept:
    # Note: The contents of dist_array must be initialized to zero on entry
    cdef:
        unsigned int N = dist_array.shape[0]
        unsigned int j, count
        ITYPE_t k
        DTYPE_t d1, d2, w12

    # relax all edges (N+1) - 1 times
    for count in range(N):
        for j in range(N):
            d1 = dist_array[j]
            for k in range(csr_indptr[j], csr_indptr[j + 1]):
                w12 = csr_weights[k]
                d2 = dist_array[csr_indices[k]]
                if d1 + w12 < d2:
                    dist_array[csr_indices[k]] = d1 + w12

    # check for negative-weight cycles
    for j in range(N):
        d1 = dist_array[j]
        for k in range(csr_indptr[j], csr_indptr[j + 1]):
            w12 = csr_weights[k]
            d2 = dist_array[csr_indices[k]]
            if d1 + w12 + DTYPE_EPS < d2:
                return j

    return -1


cdef int _johnson_undirected(
            const double[:] csr_weights,
            const int[:] csr_indices,
            const int[:] csr_indptr,
            double[:] dist_array) noexcept:
    # Note: The contents of dist_array must be initialized to zero on entry
    cdef:
        unsigned int N = dist_array.shape[0]
        unsigned int j, ind_k, count
        ITYPE_t k
        DTYPE_t d1, d2, w12

    # relax all edges (N+1) - 1 times
    for count in range(N):
        for j in range(N):
            d1 = dist_array[j]
            for k in range(csr_indptr[j], csr_indptr[j + 1]):
                w12 = csr_weights[k]
                ind_k = csr_indices[k]
                d2 = dist_array[ind_k]
                if d1 + w12 < d2:
                    dist_array[ind_k] = d1 + w12
                if d2 + w12 < d1:
                    dist_array[j] = d1 = d2 + w12

    # check for negative-weight cycles
    for j in range(N):
        d1 = dist_array[j]
        for k in range(csr_indptr[j], csr_indptr[j + 1]):
            w12 = csr_weights[k]
            d2 = dist_array[csr_indices[k]]
            if abs(d2 - d1) > w12 + DTYPE_EPS:
                return j

    return -1


# Author: Tomer Sery  -- <tomersery28@gmail.com>
# License: BSD 3-clause ("New BSD License"), (C) 2024

def yen(
    csgraph,
    source,
    sink,
    K,
    *,
    directed=True,
    return_predecessors=False,
    unweighted=False,
):
    """
    yen(csgraph, source, sink, K, *, directed=True, return_predecessors=False,
        unweighted=False)

    Yen's K-Shortest Paths algorithm on a directed or undirected graph.

    .. versionadded:: 1.14.0

    Parameters
    ----------
    csgraph : array_like, or sparse array or matrix, 2 dimensions
        The N x N array of distances representing the input graph.
    source : int
        The index of the starting node for the paths.
    sink : int
        The index of the ending node for the paths.
    K : int
        The number of shortest paths to find.
    directed : bool, optional
        If ``True`` (default), then find the shortest path on a directed graph:
        only move from point ``i`` to point ``j`` along paths ``csgraph[i, j]``.
        If False, then find the shortest path on an undirected graph: the
        algorithm can progress from point i to j along ``csgraph[i, j]`` or
        ``csgraph[j, i]``.
    return_predecessors : bool, optional
        If ``True``, return the size ``(M, N)`` predecessor matrix. Default: ``False``.
    unweighted : bool, optional
        If ``True``, then find unweighted distances. That is, rather than finding
        the path between each point such that the sum of weights is minimized,
        find the path such that the number of edges is minimized. Default: ``False``.

    Returns
    -------
    dist_array : ndarray
        Array of size ``M`` of shortest distances between the source and sink nodes.
        ``dist_array[i]`` gives the i-th shortest distance from the source to the sink
        along the graph. ``M`` is the number of shortest paths found, which is less than or
        equal to `K`.
    predecessors : ndarray
        Returned only if ``return_predecessors == True``.
        The M x N matrix of predecessors, which can be used to reconstruct
        the shortest paths.
        ``M`` is the number of shortest paths found, which is less than or equal to `K`.
        Row ``i`` of the predecessor matrix contains
        information on the ``i``-th shortest path from the source to the sink: each
        entry ``predecessors[i, j]`` gives the index of the previous node in the
        path from the source to node ``j``.  If the path does not pass via node ``j``,
        then ``predecessors[i, j] = -9999``.

    Raises
    ------
    NegativeCycleError:
        If there are negative cycles in the graph

    Notes
    -----
    Yen's algorithm is a graph search algorithm that finds single-source `K`-shortest
    loopless paths for a graph with nonnegative edge cost. The algorithm was published
    by Jin Y. Yen in 1971 and employs any shortest path algorithm to find the best path,
    then proceeds to find ``K - 1`` deviations of the best path.

    The algorithm is based on Dijsktra's algorithm for finding each shortest path.
    In case there are negative edges in the graph, Johnson's algorithm is applied.

    If multiple valid solutions are possible, output may vary with SciPy and
    Python version.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Yen%27s_algorithm
    .. [2] https://www.ams.org/journals/qam/1970-27-04/S0033-569X-1970-0253822-7/

    Examples
    --------
    >>> from scipy.sparse import csr_array
    >>> from scipy.sparse.csgraph import yen

    >>> graph = [
    ... [0, 1, 2, 0],
    ... [0, 0, 0, 1],
    ... [2, 0, 0, 3],
    ... [0, 0, 0, 0]
    ... ]
    >>> graph = csr_array(graph)
    >>> print(graph)
    <Compressed Sparse Row sparse array of dtype 'int64'
    	with 5 stored elements and shape (4, 4)>
    	Coords	Values
    	(0, 1)	1
    	(0, 2)	2
    	(1, 3)	1
    	(2, 0)	2
    	(2, 3)	3

    >>> dist_array, predecessors = yen(csgraph=graph, source=0, sink=3, K=2,
    ...                                directed=False, return_predecessors=True)
    >>> dist_array
    array([2., 5.])
    >>> predecessors
    array([[-9999,     0, -9999,     1],
        [-9999, -9999,     0,     2]], dtype=int32)

    """

    csgraph = validate_graph(csgraph, directed, DTYPE, dense_output=False)

    cdef int N = csgraph.shape[0]
    cdef int has_negative_weights = False
    dist_array = np.full(K, INFINITY, dtype=DTYPE)

    predecessor_matrix = np.full((K, N), NULL_IDX, dtype=ITYPE)

    if unweighted:
        csr_data = np.ones(csgraph.data.shape)
    else:
        csr_data = csgraph.data.copy()
        if np.any(csr_data < 0):
            # Use Johnson's algorithm to handle negative weights
            has_negative_weights = True
            johnson_dist_array = np.zeros(N, dtype=DTYPE)
            if directed:
                ret = _johnson_directed(csr_data, csgraph.indices,
                                        csgraph.indptr, johnson_dist_array)
            else:
                ret = _johnson_undirected(csr_data, csgraph.indices,
                                          csgraph.indptr, johnson_dist_array)
            if ret >= 0:
                raise NegativeCycleError("Negative cycle detected on node %i" % ret)
    if has_negative_weights:
        _johnson_add_weights(csr_data, csgraph.indices, csgraph.indptr,
                             johnson_dist_array)

    if directed:
        csgraphT = csgraph
        csrT_data = np.empty(0, dtype=DTYPE)
    else:
        csgraphT = csgraph.T.tocsr()
        if unweighted:
            csrT_data = csr_data
        else:
            if has_negative_weights:
                _johnson_add_weights(csgraphT.data, csgraphT.indices,
                                     csgraphT.indptr, johnson_dist_array)
            csrT_data = csgraphT.data

    _yen(
        source, sink,
        csr_data, csgraph.indices, csgraph.indptr,
        csrT_data, csgraphT.indices, csgraphT.indptr,
        dist_array, predecessor_matrix,
    )
    if has_negative_weights:
        dist_array += johnson_dist_array[sink] - johnson_dist_array[source]


    num_paths_found = sum(dist_array < INFINITY)
    return_shape = (num_paths_found, N)
    if return_predecessors:
        return (dist_array[:num_paths_found].reshape((num_paths_found,)),
                predecessor_matrix[:num_paths_found].reshape(return_shape))
    return dist_array[:num_paths_found].reshape((num_paths_found,))


ctypedef vector[int] yen_path_t

cdef struct YenDistancePathStruct:
    double distance
    yen_path_t path
    int spur_node


cdef inline bint _yen_compare_distance(YenDistancePathStruct a, YenDistancePathStruct b):
    return a.distance < b.distance


cdef class _YenCandidatePaths:
    cdef:
        vector[YenDistancePathStruct] _distances_and_paths
        size_t _required_paths

    def __cinit__(self, K: int):
        self._distances_and_paths = vector[YenDistancePathStruct]()
        self._required_paths = K

    @cython.boundscheck(False)
    cdef void insert_path(
        self,
        const double distance,
        const int[:] source_to_spur_path,
        const int[:] spur_to_sink_path,
        const int spur_node,
        const int sink,
    ):
        cdef:
            yen_path_t path_to_insert
            int idx = sink
            YenDistancePathStruct new_element
            vector[YenDistancePathStruct].iterator it

        if self._distances_and_paths.size() >= self._required_paths and distance >= self.max_distance():
            # The new path is longer than the longest path in the vector - return
            return

        # Store the path in reverse order, from sink to source
        # path_to_insert[0] = sink -> vec[1] -> ... -> path_to_insert[N-1] = source
        while idx != spur_node:
            path_to_insert.push_back(idx)
            idx = spur_to_sink_path[idx]
        while idx != NULL_IDX:
            path_to_insert.push_back(idx)
            idx = source_to_spur_path[idx]

        new_element.distance = distance
        new_element.path = path_to_insert
        new_element.spur_node = spur_node
        it = lower_bound(self._distances_and_paths.begin(), self._distances_and_paths.end(), new_element,
                         _yen_compare_distance)
        self._distances_and_paths.insert(it, new_element)

        # Reduce the number of paths to amount required
        while self._distances_and_paths.size() > self._required_paths:
            self._distances_and_paths.pop_back()

    cdef double min_distance(self):
        if self.empty():
            return INFINITY
        return self._distances_and_paths[0].distance

    cdef double max_distance(self):
        if self.empty():
            return -INFINITY
        return self._distances_and_paths.back().distance

    @cython.boundscheck(False)
    cdef int pop_path_to_memory_view(
        self,
        int[:] target,
    ):
        cdef:
            yen_path_t shortest_path
            size_t idx
            int spur_node

        if self.empty():
            raise RuntimeError("No paths to pop")

        shortest_path = self._distances_and_paths[0].path
        spur_node = self._distances_and_paths[0].spur_node
        self._distances_and_paths.erase(self._distances_and_paths.begin())

        # Restore the path in the correct order
        for idx in range(shortest_path.size() - 1):
            target[shortest_path[idx]] = shortest_path[idx + 1]

        self._required_paths -= 1

        return spur_node

    cdef bint empty(self):
        return self._distances_and_paths.empty()


@cython.boundscheck(False)
cdef void _yen(
    const int source,
    const int sink,
    const double[:] original_weights, const int[:] csr_indices, const int[:] csr_indptr,
    const double[:] originalT_weights, const int[:] csrT_indices, const int[:] csrT_indptr,
    double[:] shortest_distances,
    int[:, :] shortest_paths_predecessors,
):
    cdef:
        int K = shortest_paths_predecessors.shape[0] # Number of paths to find
        int N = shortest_paths_predecessors.shape[1] # Number of nodes in graph
        bint directed = originalT_weights.size == 0

        # Dijkstra's operands and results arrays
        int[:] indice_node_arr = np.array([source], dtype=ITYPE)
        int[:] predecessor_matrix = np.full((N), NULL_IDX, dtype=ITYPE)
        double[:] dist_matrix = np.full((N), np.inf, dtype=DTYPE)
        int[:] dummy_source_matrix = np.empty((0), dtype=ITYPE) # unused
        _YenCandidatePaths candidate_paths = _YenCandidatePaths(K)

    dist_matrix[source] = 0

    # ---------------------------------------------------
    # Compute and store the shortest path
    _dijkstra(
        indice_node_arr,
        original_weights, csr_indices, csr_indptr,
        originalT_weights, csrT_indices, csrT_indptr,
        dist_matrix, predecessor_matrix, dummy_source_matrix,
        INFINITY,
    )

    shortest_distances[0] = dist_matrix[sink]
    if shortest_distances[0] == INFINITY:
        # No paths between source and sink
        return

    cdef:
        # Store the original graph weights for restoring the graph
        double[:] csr_weights = original_weights.copy()
        double[:] csrT_weights

        int k, i, spur_node, node, short_path_idx, tmp_i
        int spur_node_k_minus_1
        double root_path_distance, total_distance, tmp_d

    # Avoid copying a size 0 memory view
    if directed:
        csrT_weights = np.empty(0, dtype=DTYPE)
    else:
        csrT_weights = originalT_weights.copy()

    # Copy shortest path to shortest_paths_predecessors
    node = sink
    while node != NULL_IDX:
        shortest_paths_predecessors[0, node] = predecessor_matrix[node]
        node = predecessor_matrix[node]


    # ---------------------------------------------------
    # Compute and store the K-1 shortest paths
    spur_node_k_minus_1 = source
    for k in range(1, K):
        # Set spur node as sink
        spur_node = sink
        # Set the original path distance as the previous shortest distance
        root_path_distance = shortest_distances[k-1]

        # ---------------------------------------------------
        # For each spur_node in the previous k-shortest path
        # Search for a new short path from it to the sink
        while spur_node != spur_node_k_minus_1:
            # Decrease the root path distance by the distance of it's final edge and
            # set the source of the final edge as the new spur node
            tmp_i = shortest_paths_predecessors[k-1][spur_node] # previous node
            tmp_d = INFINITY # last edge distance
            for i in range(csr_indptr[tmp_i], csr_indptr[tmp_i + 1]):
                if csr_indices[i] == spur_node:
                    tmp_d = csr_weights[i]
                    break
            if not directed:
                for i in range(csrT_indptr[tmp_i], csrT_indptr[tmp_i + 1]):
                    if csrT_indices[i] == spur_node and csrT_weights[i] < tmp_d:
                        tmp_d = csrT_weights[i]
                        break
            if tmp_d == INFINITY:
                raise ValueError(f"No edge between nodes {tmp_i} and {spur_node}")
            root_path_distance -= tmp_d
            spur_node = tmp_i

            # ---------------------------------------------------
            # Avoid following found shortest paths
            for short_path_idx in range(k):
                # For each shortest path
                # Remove the edge {spur_node -> next node} in shortest path
                # If the original path coincides with the current shortest path up to
                # spur node.
                node = spur_node
                while (
                    shortest_paths_predecessors[short_path_idx][node]
                    == shortest_paths_predecessors[k-1][node]
                ):
                    if node == source:
                        # Remove edge spur_node -> next node
                        for i in range(
                            csr_indptr[spur_node], csr_indptr[spur_node + 1]
                        ):
                            if (
                                spur_node
                                == shortest_paths_predecessors[short_path_idx][
                                    csr_indices[i]
                                ]
                            ):
                                csr_weights[i] = INFINITY
                        if not directed:
                            for i in range(
                                csrT_indptr[spur_node], csrT_indptr[spur_node + 1]
                            ):
                                if (
                                    spur_node
                                    == shortest_paths_predecessors[short_path_idx][
                                        csrT_indices[i]
                                    ]
                                ):
                                    csrT_weights[i] = INFINITY

                        break
                    node = shortest_paths_predecessors[short_path_idx][node]

            # ---------------------------------------------------
            # Avoid loops in paths by removing all nodes of the root path from the graph
            # except for the spur node.
            # A node is removed from the graph by setting all its out-edges to infinity
            node = shortest_paths_predecessors[k-1][spur_node]
            while node != NULL_IDX:
                csr_weights[csr_indptr[node]: csr_indptr[node + 1]] = INFINITY
                if not directed:
                    csrT_weights[csrT_indptr[node]: csrT_indptr[node + 1]] = INFINITY
                node = shortest_paths_predecessors[k-1][node]

            # ---------------------------------------------------
            # Search for the shortest path from spur_node to sink

            # Reset the distance and predecessor matrix
            predecessor_matrix[:] = NULL_IDX
            dist_matrix[:] = INFINITY
            dist_matrix[source] = 0
            # Search only for paths starting for spur_node
            indice_node_arr[0] = spur_node
            _dijkstra(
                indice_node_arr,
                csr_weights, csr_indices, csr_indptr,
                csrT_weights, csrT_indices, csrT_indptr,
                dist_matrix, predecessor_matrix, dummy_source_matrix,
                INFINITY,
            )

            # Compute the total distance of the found path
            total_distance = dist_matrix[sink] + root_path_distance

            # ---------------------------------------------------
            # Add the found path to arrays of candidates
            if total_distance != INFINITY:
                candidate_paths.insert_path(
                    total_distance,
                    shortest_paths_predecessors[k-1],
                    predecessor_matrix,
                    spur_node,
                    sink,
                )

           # ---------------------------------------------------
            # Restore graph weights
            node = spur_node
            while node != NULL_IDX:
                csr_weights[csr_indptr[node]: csr_indptr[node + 1]] = (
                    original_weights[csr_indptr[node]: csr_indptr[node + 1]]
                )
                if not directed:
                    csrT_weights[csrT_indptr[node]: csrT_indptr[node + 1]] = (
                        originalT_weights[csrT_indptr[node]: csrT_indptr[node + 1]]
                    )
                node = shortest_paths_predecessors[k-1, node]


        # ---------------------------------------------------
        # Find shortest path in candidates and add to result arrays
        total_distance = candidate_paths.min_distance()
        if total_distance == INFINITY:
            # There are no more paths
            break
        else:
            shortest_distances[k] = total_distance
            spur_node_k_minus_1 = candidate_paths.pop_path_to_memory_view(shortest_paths_predecessors[k])
