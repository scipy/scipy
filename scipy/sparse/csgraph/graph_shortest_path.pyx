"""
Routines for performing shortest-path graph searches

The main interface is in the function `cs_graph_shortest_path`.  This
calls cython routines that compute the shortest path using either
the Floyd-Warshall algorithm, or Dykstra's algorithm with Fibonacci Heaps.
"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) 2011

import numpy as np
cimport numpy as np

from scipy.sparse import csr_matrix, isspmatrix, isspmatrix_csr, isspmatrix_csc

cimport cython

from libc.stdlib cimport malloc, free

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t
    

def cs_graph_shortest_path(csgraph, directed=True,
                           method='auto', overwrite=True):
    """
    Perform a shortest-path graph search on a positive directed or
    undirected graph.

    Parameters
    ----------
    csgraph : array, matrix, or sparse matrix, shape = (N,N)
        Array of non-negative distances.
        If vertex i is connected to vertex j, then csgraph[i,j] gives
        the distance between the vertices.
        If vertex i is not connected to vertex j, then csgraph[i,j] = 0
    directed : boolean
        if True, then find the shortest path on a directed graph: only
        progress from a point to its neighbors, not the other way around.
        if False, then find the shortest path on an undirected graph: the
        algorithm can progress from a point to its neighbors and vice versa.
    method : string ['auto'|'FW'|'D']
        method to use.  Options are
        'auto' : attempt to choose the best method for the current problem
        'FW' : Floyd-Warshall algorithm.  O[N^3]
        'D' : Dijkstra's algorithm with Fibonacci heaps.  O[(k+log(N))N^2]
    overwrite : bool, default=True
        Overwrite csgraph with the result.  This applies only if
        csgraph is a dense, c-ordered array with dtype=float64.
        Otherwise, a copy will be made.

    Returns
    -------
    dist_matrix : np.ndarray, float, shape = [N,N]
        dist_matrix[i,j] gives the shortest distance from point i to point j
        along the graph.

    Notes
    -----
    As currently implemented, Dijkstra's algorithm does not work for
    graphs with direction-dependent distances when directed == False.
    i.e., if csgraph[i,j] and csgraph[j,i] are not equal and
    both are nonzero, method='D' will not necessarily yield the correct
    result.

    Also, these routines have not been tested for graphs with negative
    distances.  Negative distances can lead to infinite cycles that must
    be handled by specialized algorithms.
    """
    if not isspmatrix(csgraph):
        csgraph = csr_matrix(csgraph)

    N = csgraph.shape[0]
    Nk = len(csgraph.data)

    if method == 'auto':
        if Nk < N * N / 4:
            method = 'D'
        else:
            method = 'FW'

    if method == 'FW':
        dist_matrix = floyd_warshall(csgraph, directed, overwrite)
    elif method == 'D':
        dist_matrix = dijkstra(csgraph, directed)
    else:
        raise ValueError("unrecognized method '%s'" % method)

    return dist_matrix


def floyd_warshall(csgraph, directed=False, overwrite=True):
    """
    Compute the shortest path lengths using the Floyd-Warshall algorithm

    Parameters
    ----------
    csgraph : array, matrix, or sparse matrix, shape=(N, N)
        Array of positive distances.
        If vertex i is connected to vertex j, then csgraph[i,j] gives
        the distance between the vertices.
        If vertex i is not connected to vertex j, then csgraph[i,j] = 0
    directed : bool, default = False
        if True, then find the shortest path on a directed graph: only
        progress from a point to its neighbors, not the other way around.
        if False, then find the shortest path on an undirected graph: the
        algorithm can progress from a point to its neighbors and vice versa.
    overwrite : bool, default=True
        Overwrite csgraph with the result.  This applies only if
        csgraph is a dense, c-ordered array with dtype=float64.
        Otherwise, a copy will be made.

    Returns
    -------
    dist_matrix : ndarray, shape=[N, N]
        the matrix of shortest paths between points.
        If no path exists, the path length is zero

    Notes
    -----
    Thes routine has not been tested for graphs with negative
    distances.  Negative distances can lead to unanticipated results.
    """
    # graph needs to be a dense, C-ordered copy of csgraph with the
    # correct dtype.  Addisionally, if overwrite is False, we need to
    # assure that a copy is made.
    if isspmatrix(csgraph):
        dist_matrix = np.asarray(csgraph.toarray(), dtype=DTYPE, order='C')
    else:
        if overwrite:
            dist_matrix = np.asarray(csgraph, dtype=DTYPE, order='C')
        else:
            dist_matrix = np.array(csgraph, dtype=DTYPE, order='C')

    if dist_matrix.ndim != 2 or dist_matrix.shape[0] != dist_matrix.shape[1]:
        raise ValueError("csgraph must have shape (N, N)")

    return _floyd_warshall(dist_matrix, int(directed))


@cython.boundscheck(False)
cdef np.ndarray _floyd_warshall(
                   np.ndarray[DTYPE_t, ndim=2, mode='c'] dist_matrix,
                   int directed=0):
    # dist_matrix should be a [N,N] matrix, such that dist_matrix[i, j]
    # is the distance from point i to point j.  Zero-distances imply that
    # the points are not connected.
    cdef int N = dist_matrix.shape[0]
    assert dist_matrix.shape[1] == N

    cdef unsigned int i, j, k, m

    cdef DTYPE_t infinity = np.inf
    cdef DTYPE_t sum_ijk

    #initialize all distances to infinity
    dist_matrix[np.where(dist_matrix == 0)] = infinity

    # ensure dist_matrix[i,i] is zero
    dist_matrix.flat[::N + 1] = 0

    # for a non-directed dist_matrix, we need to symmetrize the distances
    if not directed:
        for i from 0 <= i < N:
            for j from i + 1 <= j < N:
                if dist_matrix[j, i] <= dist_matrix[i, j]:
                    dist_matrix[i, j] = dist_matrix[j, i]
                else:
                    dist_matrix[j, i] = dist_matrix[i, j]

    #now perform the Floyd-Warshall algorithm
    for k from 0 <= k < N:
        for i from 0 <= i < N:
            if dist_matrix[i, k] == infinity:
                continue
            for j from 0 <= j < N:
                sum_ijk = dist_matrix[i, k] + dist_matrix[k, j]
                if sum_ijk < dist_matrix[i, j]:
                    dist_matrix[i, j] = sum_ijk

    dist_matrix[np.where(np.isinf(dist_matrix))] = 0

    return dist_matrix


def dijkstra(csgraph, directed=False, indices=None):
    """
    Dijkstra algorithm using Fibonacci Heaps

    Parameters
    ----------
    csgraph : array, matrix, or sparse matrix, shape=(N, N)
        Array of positive distances: this will be internally converted to
        a csr_matrix with dtype=np.float64.
        If vertex i is connected to vertex j, then csgraph[i,j] gives
        the distance between the vertices.
        If vertex i is not connected to vertex j, then csgraph[i,j] = 0
    directed : bool, default = False
        if True, then find the shortest path on a directed graph: only
        progress from a point to its neighbors, not the other way around.
        if False, then find the shortest path on an undirected graph: the
        algorithm can progress from a point to its neighbors and vice versa.
        If directed == False, then csgraph must be of a certain form:
        see the notes below.
    indices : 1d array or None
        if specified, only compute the paths for the points at the given
        indices.

    Returns
    -------
    dist_matrix : array, shape = (Nind, N)
        the matrix of shortest paths between points.
        If no path exists, the path length is zero.
        If indices == None, then Nind = N.
        If indices is specified, then Nind = len(indices)

    Notes
    -----
    As currently implemented, Dijkstra's algorithm does not work for
    graphs with direction-dependent distances when directed == False.
    i.e., if csgraph[i,j] and csgraph[j,i] are not equal and
    both are nonzero, setting directed=False will not yield the correct
    result.

    Also, this routine does not work for graphs with negative
    distances.  Negative distances can lead to infinite cycles that must
    be handled by specialized algorithms.
    """
    # We need to convert csgraph to a csr_matrix, and create
    # dense c-ordered float64 matrix to contain the result.
    if isspmatrix_csr(csgraph):
        pass
    elif (not directed) and isspmatrix_csc(csgraph):
        # if directed, then it's safe to assume that csgraph
        # and csgraph.T can be treated equivalently.  We
        # can get it into csr_matrix format very efficiently.
        csgraph = csgraph.T
    else:
        csgraph = csr_matrix(csgraph, dtype=DTYPE)

    if np.any(csgraph.data < 0):
        raise ValueError("Negative distances are not supported")

    if csgraph.ndim != 2 or csgraph.shape[0] != csgraph.shape[1]:
        raise ValueError("csgraph must have shape (N, N)")

    N = csgraph.shape[0]

    if indices is None:
        indices = np.arange(N, dtype=ITYPE)
    else:
        indices = np.asarray(indices, order='C', dtype=ITYPE)
    
    return_shape = indices.shape + (N,)
    indices = np.atleast_1d(indices).reshape(-1)
    dist_matrix = np.zeros((len(indices), N), dtype=DTYPE)

    dist_matrix = _dijkstra(csgraph, dist_matrix, indices, directed)

    return dist_matrix.reshape(return_shape)

@cython.boundscheck(False)
cdef np.ndarray _dijkstra(csgraph,
                          np.ndarray[DTYPE_t, ndim=2, mode='c'] dist_matrix,
                          np.ndarray[ITYPE_t, ndim=1, mode='c'] compute_ind,
                          int directed=0):
    # csgraph is a square csr_matrix, or object with attributes `data`,
    # `indices`, and `indptr` which store a matrix in csr format.
    # `graph` is an uninitialized array which will store the output.
    # `compute_ind` gives the indices of paths to compute
    # `graph` is assumed to be shape [Nind, N], where csgraph is shape
    #   (N, N) and compute_ind is length Nind.as csgraph.  If this is
    #   not the case, then a memory error/segfault could result.
    # if directed is false, we convert the csr matrix to a csc matrix
    #   in order to find bi-directional distances.
    cdef unsigned int Nind = dist_matrix.shape[0]
    cdef unsigned int N = dist_matrix.shape[1]
    cdef unsigned int i

    cdef FibonacciHeap heap

    cdef FibonacciNode* nodes = <FibonacciNode*> malloc(N *
                                                        sizeof(FibonacciNode))

    cdef np.ndarray distances, neighbors, indptr
    cdef np.ndarray distances2, neighbors2, indptr2

    if not isspmatrix_csr(csgraph):
        csgraph = csr_matrix(csgraph)

    distances = np.asarray(csgraph.data, dtype=DTYPE, order='C')
    neighbors = np.asarray(csgraph.indices, dtype=ITYPE, order='C')
    indptr = np.asarray(csgraph.indptr, dtype=ITYPE, order='C')

    for i from 0 <= i < N:
        initialize_node(&nodes[i], i)

    heap.min_node = NULL

    if directed:
        for i from 0 <= i < Nind:
            _dijkstra_directed_one_row(compute_ind[i],
                                       neighbors, distances, indptr,
                                       dist_matrix, &heap, nodes)
    else:
        #use the csr -> csc sparse matrix conversion to quickly get
        # both directions of neigbors
        csgraph_T = csgraph.T.tocsr()

        distances2 = np.asarray(csgraph_T.data,
                                dtype=DTYPE, order='C')
        neighbors2 = np.asarray(csgraph_T.indices,
                                dtype=ITYPE, order='C')
        indptr2 = np.asarray(csgraph_T.indptr,
                             dtype=ITYPE, order='C')

        for i from 0 <= i < Nind:
            _dijkstra_one_row(compute_ind[i],
                              neighbors, distances, indptr,
                              neighbors2, distances2, indptr2,
                              dist_matrix, &heap, nodes)

    free(nodes)

    return dist_matrix


######################################################################
# FibonacciNode structure
#  This structure and the operations on it are the nodes of the
#  Fibonacci heap.
#
cdef enum FibonacciState:
    SCANNED
    NOT_IN_HEAP
    IN_HEAP

cdef struct FibonacciNode:
    unsigned int index
    unsigned int rank
    FibonacciState state
    DTYPE_t val
    FibonacciNode* parent
    FibonacciNode* left_sibling
    FibonacciNode* right_sibling
    FibonacciNode* children


cdef void initialize_node(FibonacciNode* node,
                          unsigned int index,
                          DTYPE_t val=0):
    # Assumptions: - node is a valid pointer
    #              - node is not currently part of a heap
    node.index = index
    node.val = val
    node.rank = 0
    node.state = NOT_IN_HEAP

    node.parent = NULL
    node.left_sibling = NULL
    node.right_sibling = NULL
    node.children = NULL


cdef FibonacciNode* rightmost_sibling(FibonacciNode* node):
    # Assumptions: - node is a valid pointer
    cdef FibonacciNode* temp = node
    while(temp.right_sibling):
        temp = temp.right_sibling
    return temp


cdef FibonacciNode* leftmost_sibling(FibonacciNode* node):
    # Assumptions: - node is a valid pointer
    cdef FibonacciNode* temp = node
    while(temp.left_sibling):
        temp = temp.left_sibling
    return temp


cdef void add_child(FibonacciNode* node, FibonacciNode* new_child):
    # Assumptions: - node is a valid pointer
    #              - new_child is a valid pointer
    #              - new_child is not the sibling or child of another node
    new_child.parent = node

    if node.children:
        add_sibling(node.children, new_child)
    else:
        node.children = new_child
        new_child.right_sibling = NULL
        new_child.left_sibling = NULL
        node.rank = 1


cdef void add_sibling(FibonacciNode* node, FibonacciNode* new_sibling):
    # Assumptions: - node is a valid pointer
    #              - new_sibling is a valid pointer
    #              - new_sibling is not the child or sibling of another node
    cdef FibonacciNode* temp = rightmost_sibling(node)
    temp.right_sibling = new_sibling
    new_sibling.left_sibling = temp
    new_sibling.right_sibling = NULL
    new_sibling.parent = node.parent
    if new_sibling.parent:
        new_sibling.parent.rank += 1


cdef void remove(FibonacciNode* node):
    # Assumptions: - node is a valid pointer
    if node.parent:
        node.parent.rank -= 1
        if node.left_sibling:
            node.parent.children = node.left_sibling
        elif node.right_sibling:
            node.parent.children = node.right_sibling
        else:
            node.parent.children = NULL

    if node.left_sibling:
        node.left_sibling.right_sibling = node.right_sibling
    if node.right_sibling:
        node.right_sibling.left_sibling = node.left_sibling

    node.left_sibling = NULL
    node.right_sibling = NULL
    node.parent = NULL


######################################################################
# FibonacciHeap structure
#  This structure and operations on it use the FibonacciNode
#  routines to implement a Fibonacci heap

ctypedef FibonacciNode* pFibonacciNode


cdef struct FibonacciHeap:
    FibonacciNode* min_node
    pFibonacciNode[100] roots_by_rank  # maximum number of nodes is ~2^100.


cdef void insert_node(FibonacciHeap* heap,
                      FibonacciNode* node):
    # Assumptions: - heap is a valid pointer
    #              - node is a valid pointer
    #              - node is not the child or sibling of another node
    if heap.min_node:
        add_sibling(heap.min_node, node)
        if node.val < heap.min_node.val:
            heap.min_node = node
    else:
        heap.min_node = node


cdef void decrease_val(FibonacciHeap* heap,
                       FibonacciNode* node,
                       DTYPE_t newval):
    # Assumptions: - heap is a valid pointer
    #              - newval <= node.val
    #              - node is a valid pointer
    #              - node is not the child or sibling of another node
    node.val = newval
    if node.parent and (node.parent.val >= newval):
        remove(node)
        insert_node(heap, node)


cdef void link(FibonacciHeap* heap, FibonacciNode* node):
    # Assumptions: - heap is a valid pointer
    #              - node is a valid pointer
    #              - node is already within heap

    cdef FibonacciNode *linknode, *parent, *child

    if heap.roots_by_rank[node.rank] == NULL:
        heap.roots_by_rank[node.rank] = node
    else:
        linknode = heap.roots_by_rank[node.rank]
        heap.roots_by_rank[node.rank] = NULL

        if node.val < linknode.val or node == heap.min_node:
            remove(linknode)
            add_child(node, linknode)
            link(heap, node)
        else:
            remove(node)
            add_child(linknode, node)
            link(heap, linknode)


cdef FibonacciNode* remove_min(FibonacciHeap* heap):
    # Assumptions: - heap is a valid pointer
    #              - heap.min_node is a valid pointer
    cdef FibonacciNode *temp, *temp_right, *out
    cdef unsigned int i

    # make all min_node children into root nodes
    if heap.min_node.children:
        temp = leftmost_sibling(heap.min_node.children)
        temp_right = NULL

        while temp:
            temp_right = temp.right_sibling
            remove(temp)
            add_sibling(heap.min_node, temp)
            temp = temp_right

        heap.min_node.children = NULL

    # choose a root node other than min_node
    temp = leftmost_sibling(heap.min_node)
    if temp == heap.min_node:
        if heap.min_node.right_sibling:
            temp = heap.min_node.right_sibling
        else:
            out = heap.min_node
            heap.min_node = NULL
            return out

    # remove min_node, and point heap to the new min
    out = heap.min_node
    remove(heap.min_node)
    heap.min_node = temp

    # re-link the heap
    for i from 0 <= i < 100:
        heap.roots_by_rank[i] = NULL

    while temp:
        if temp.val < heap.min_node.val:
            heap.min_node = temp
        temp_right = temp.right_sibling
        link(heap, temp)
        temp = temp_right

    return out


######################################################################
# Debugging: Functions for printing the Fibonacci heap
#
#cdef void print_node(FibonacciNode* node, int level=0):
#    print '%s(%i,%i) %i' % (level*'   ', node.index, node.val, node.rank)
#    if node.children:
#        print_node(leftmost_sibling(node.children), level+1)
#    if node.right_sibling:
#        print_node(node.right_sibling, level)
#
#
#cdef void print_heap(FibonacciHeap* heap):
#    print "---------------------------------"
#    if heap.min_node:
#        print_node(leftmost_sibling(heap.min_node))
#    else:
#        print "[empty heap]"


@cython.boundscheck(False)
cdef void _dijkstra_directed_one_row(
                          unsigned int i_node,
                          np.ndarray[ITYPE_t, ndim=1, mode='c'] neighbors,
                          np.ndarray[DTYPE_t, ndim=1, mode='c'] distances,
                          np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr,
                          np.ndarray[DTYPE_t, ndim=2, mode='c'] graph,
                          FibonacciHeap* heap,
                          FibonacciNode* nodes):
    # Calculate distances from a single point to all targets using a
    # directed graph.
    #
    # Parameters
    # ----------
    # i_node : index of source point
    # neighbors : array, shape = [N,]
    #     indices of neighbors for each point
    # distances : array, shape = [N,]
    #     lengths of edges to each neighbor
    # indptr : array, shape = (N+1,)
    #     the neighbors of point i are given by
    #     neighbors[indptr[i]:indptr[i+1]]
    # graph : array, shape = (Nind,N)
    #     on return, graph[i_node] contains the path lengths from
    #     i_node to each target
    # heap: the Fibonacci heap object to use
    # nodes : the array of nodes to use
    cdef unsigned int N = graph.shape[1]
    cdef unsigned int i
    cdef FibonacciNode *v, *current_neighbor
    cdef DTYPE_t dist

    # initialize nodes
    for i from 0 <= i < N:
        initialize_node(&nodes[i], i)

    heap.min_node = NULL
    insert_node(heap, &nodes[i_node])

    while heap.min_node:
        v = remove_min(heap)
        v.state = SCANNED

        for i from indptr[v.index] <= i < indptr[v.index + 1]:
            current_neighbor = &nodes[neighbors[i]]
            if current_neighbor.state != SCANNED:
                dist = distances[i]
                if current_neighbor.state == NOT_IN_HEAP:
                    current_neighbor.state = IN_HEAP
                    current_neighbor.val = v.val + dist
                    insert_node(heap, current_neighbor)
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)

        #v has now been scanned: add the distance to the results
        graph[i_node, v.index] = v.val


@cython.boundscheck(False)
cdef void _dijkstra_one_row(unsigned int i_node,
                            np.ndarray[ITYPE_t, ndim=1, mode='c'] neighbors1,
                            np.ndarray[DTYPE_t, ndim=1, mode='c'] distances1,
                            np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr1,
                            np.ndarray[ITYPE_t, ndim=1, mode='c'] neighbors2,
                            np.ndarray[DTYPE_t, ndim=1, mode='c'] distances2,
                            np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr2,
                            np.ndarray[DTYPE_t, ndim=2, mode='c'] graph,
                            FibonacciHeap* heap,
                            FibonacciNode* nodes):
    # Calculate distances from a single point to all targets using an
    # undirected graph.
    #
    # Parameters
    # ----------
    # i_node : index of source point
    # neighbors[1,2] : array, shape = [N,]
    #     indices of neighbors for each point
    # distances[1,2] : array, shape = [N,]
    #     lengths of edges to each neighbor
    # indptr[1,2] : array, shape = (N+1,)
    #     the neighbors of point i are given by
    #     neighbors1[indptr1[i]:indptr1[i+1]] and
    #     neighbors2[indptr2[i]:indptr2[i+1]]
    # graph : array, shape = (Nind, N)
    #     on return, graph[i_node] contains the path lengths from
    #     i_node to each target
    # heap: the Fibonacci heap object to use
    # nodes : the array of nodes to use
    cdef unsigned int N = graph.shape[1]
    cdef unsigned int i
    cdef FibonacciNode *v, *current_neighbor
    cdef DTYPE_t dist

    # re-initialize nodes
    # children, parent, left_sibling, right_sibling should already be NULL
    # rank should already be 0, index will already be set
    # we just need to re-set state and val
    for i from 0 <= i < N:
        nodes[i].state = NOT_IN_HEAP
        nodes[i].val = 0

    insert_node(heap, &nodes[i_node])

    while heap.min_node:
        v = remove_min(heap)
        v.state = SCANNED

        for i from indptr1[v.index] <= i < indptr1[v.index + 1]:
            current_neighbor = &nodes[neighbors1[i]]
            if current_neighbor.state != SCANNED:
                dist = distances1[i]
                if current_neighbor.state == NOT_IN_HEAP:
                    current_neighbor.state = IN_HEAP
                    current_neighbor.val = v.val + dist
                    insert_node(heap, current_neighbor)
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)

        for i from indptr2[v.index] <= i < indptr2[v.index + 1]:
            current_neighbor = &nodes[neighbors2[i]]
            if current_neighbor.state != SCANNED:
                dist = distances2[i]
                if current_neighbor.state == NOT_IN_HEAP:
                    current_neighbor.state = IN_HEAP
                    current_neighbor.val = v.val + dist
                    insert_node(heap, current_neighbor)
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)

        #v has now been scanned: add the distance to the results
        graph[i_node, v.index] = v.val
