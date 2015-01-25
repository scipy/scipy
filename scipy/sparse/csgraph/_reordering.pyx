# Author: Paul Nation  -- <pnation@korea.ac.kr>
# Original Source: QuTiP: Quantum Toolbox in Python (qutip.org)
# License: New BSD, (C) 2014

import numpy as np
cimport numpy as np

from scipy.sparse import isspmatrix_csr, isspmatrix_csc, isspmatrix_coo

include 'parameters.pxi'

def reverse_cuthill_mckee(graph, symmetric_mode=False):
    """
    Returns the permutation array that orders a sparse CSR or CSC matrix
    in Reverse-Cuthill McKee ordering.  
    
    It is assumed by default, ``symmetric_mode=False``, that the input matrix 
    is not symmetric and works on the matrix ``A+A.T``. If you are 
    guaranteed that the matrix is symmetric in structure (values of matrix 
    elements do not matter) then set ``symmetric_mode=True``.
    
    Parameters
    ----------
    graph : sparse matrix
        Input sparse in CSC or CSR sparse matrix format.
    symmetric_mode : bool, optional
        Is input matrix guaranteed to be symmetric.

    Returns
    -------
    perm : ndarray
        Array of permuted row and column indices.
 
    References
    ----------
    E. Cuthill and J. McKee, "Reducing the Bandwidth of Sparse Symmetric Matrices",
    ACM '69 Proceedings of the 1969 24th national conference, (1969).
    
    """
    if not (isspmatrix_csc(graph) or isspmatrix_csr(graph)):
        raise TypeError('Input must be in CSC or CSR sparse matrix format.')
    nrows = graph.shape[0]
    if not symmetric_mode:
        graph = graph+graph.transpose()
    return _reverse_cuthill_mckee(graph.indices, graph.indptr, nrows)


def maximum_bipartite_matching(graph, perm_type='row'):
    """
    Returns an array of row or column permutations that makes
    the diagonal of a nonsingular square CSC sparse matrix zero free.  
    
    Such a permutation is always possible provided that the matrix 
    is nonsingular. This function looks at the structure of the matrix 
    only. The input matrix will be converted to CSC matrix format if
    necessary.

    Parameters
    ----------
    graph : sparse matrix
        Input sparse in CSC format
    perm_type : str, {'row', 'column'}
        Type of permutation to generate.

    Returns
    -------
    perm : ndarray
        Array of row or column permutations.

    Notes
    -----
    This function relies on a maximum cardinality bipartite matching 
    algorithm based on a breadth-first search (BFS) of the underlying 
    graph.

    References
    ----------
    I. S. Duff, K. Kaya, and B. Ucar, "Design, Implementation, and 
    Analysis of Maximum Transversal Algorithms", ACM Trans. Math. Softw.
    38, no. 2, (2011).

    """
    cdef np.npy_intp nrows = graph.shape[0]
    if nrows != graph.shape[1]:
        raise ValueError('Maximum bipartite matching requires a square matrix.')
    if isspmatrix_csr(graph) or isspmatrix_coo(graph):
        graph = graph.tocsc()
    elif not isspmatrix_csc(graph):
        raise TypeError("graph must be in CSC, CSR, or COO format.")
    if perm_type == 'column':
        graph = graph.transpose().tocsc()
    perm = _maximum_bipartite_matching(graph.indices, graph.indptr, nrows)
    if np.any(perm==-1):
        raise Exception('Possibly singular input matrix.')
    return perm


cdef _node_degrees(
        np.ndarray[int32_or_int64, ndim=1, mode="c"] ind,
        np.ndarray[int32_or_int64, ndim=1, mode="c"] ptr,
        np.npy_intp num_rows):
    """
    Find the degree of each node (matrix row) in a graph represented
    by a sparse CSR or CSC matrix.
    """
    cdef np.npy_intp ii, jj
    cdef np.ndarray[int32_or_int64] degree = np.zeros(num_rows, dtype=ind.dtype)
    
    for ii in range(num_rows):
        degree[ii] = ptr[ii + 1] - ptr[ii]
        for jj in range(ptr[ii], ptr[ii + 1]):
            if ind[jj] == ii:
                # add one if the diagonal is in row ii
                degree[ii] += 1
                break
    return degree
    

def _reverse_cuthill_mckee(np.ndarray[int32_or_int64, ndim=1, mode="c"] ind,
        np.ndarray[int32_or_int64, ndim=1, mode="c"] ptr,
        np.npy_intp num_rows):
    """
    Reverse Cuthill-McKee ordering of a sparse symmetric CSR or CSC matrix.  
    We follow the original Cuthill-McKee paper and always start the routine
    at a node of lowest degree for each connected component.
    """
    cdef np.npy_intp N = 0, N_old, level_start, level_end, temp
    cdef np.npy_intp zz, ii, jj, kk, ll, level_len
    cdef np.ndarray[int32_or_int64] order = np.zeros(num_rows, dtype=ind.dtype)
    cdef np.ndarray[int32_or_int64] degree = _node_degrees(ind, ptr, num_rows)
    cdef np.ndarray[np.npy_intp] inds = np.argsort(degree)
    cdef np.ndarray[np.npy_intp] rev_inds = np.argsort(inds)
    cdef np.ndarray[ITYPE_t] temp_degrees = np.zeros(np.max(degree), dtype=ITYPE)
    cdef int32_or_int64 i, j, seed, temp2
    
    # loop over zz takes into account possible disconnected graph.
    for zz in range(num_rows):
        if inds[zz] != -1:   # Do BFS with seed=inds[zz]
            seed = inds[zz]
            order[N] = seed
            N += 1
            inds[rev_inds[seed]] = -1
            level_start = N - 1
            level_end = N

            while level_start < level_end:
                for ii in range(level_start, level_end):
                    i = order[ii]
                    N_old = N

                    # add unvisited neighbors
                    for jj in range(ptr[i], ptr[i + 1]):
                        # j is node number connected to i
                        j = ind[jj]
                        if inds[rev_inds[j]] != -1:
                            inds[rev_inds[j]] = -1
                            order[N] = j
                            N += 1

                    # Add values to temp_degrees array for insertion sort
                    level_len = 0
                    for kk in range(N_old, N):
                        temp_degrees[level_len] = degree[order[kk]]
                        level_len += 1
                
                    # Do insertion sort for nodes from lowest to highest degree
                    for kk in range(1,level_len):
                        temp = temp_degrees[kk]
                        temp2 = order[N_old+kk]
                        ll = kk
                        while (ll > 0) and (temp < temp_degrees[ll-1]):
                            temp_degrees[ll] = temp_degrees[ll-1]
                            order[N_old+ll] = order[N_old+ll-1]
                            ll -= 1
                        temp_degrees[ll] = temp
                        order[N_old+ll] = temp2
                
                # set next level start and end ranges
                level_start = level_end
                level_end = N

        if N == num_rows:
            break

    # return reversed order for RCM ordering
    return order[::-1]


def _maximum_bipartite_matching(
        np.ndarray[int32_or_int64, ndim=1, mode="c"] inds,
        np.ndarray[int32_or_int64, ndim=1, mode="c"] ptrs,
        np.npy_intp n):
    """
    Maximum bipartite matching of a graph in CSC format.
    """
    cdef np.ndarray[int32_or_int64] visited = np.zeros(n, dtype=inds.dtype)
    cdef np.ndarray[ITYPE_t] queue = np.zeros(n, dtype=ITYPE)
    cdef np.ndarray[ITYPE_t] previous = np.zeros(n, dtype=ITYPE)
    cdef np.ndarray[int32_or_int64] match = np.empty(n, dtype=inds.dtype)
    cdef np.ndarray[ITYPE_t] row_match = np.empty(n, dtype=ITYPE)
    cdef np.npy_intp queue_ptr, queue_col, ptr, i, j, queue_size
    cdef np.npy_intp col, next_num = 1
    cdef int32_or_int64 row, temp, eptr

    for i in range(n):
        match[i] = -1
        row_match[i] = -1

    for i in range(n):
        if match[i] == -1 and (ptrs[i] != ptrs[i + 1]):
            queue[0] = i
            queue_ptr = 0
            queue_size = 1
            while (queue_size > queue_ptr):
                queue_col = queue[queue_ptr]
                queue_ptr += 1
                eptr = ptrs[queue_col + 1]
                for ptr in range(ptrs[queue_col], eptr):
                    row = inds[ptr]
                    temp = visited[row]
                    if (temp != next_num and temp != -1):
                        previous[row] = queue_col
                        visited[row] = next_num
                        col = row_match[row]
                        if (col == -1):
                            while (row != -1):
                                col = previous[row]
                                temp = match[col]
                                match[col] = row
                                row_match[row] = col
                                row = temp
                            next_num += 1
                            queue_size = 0
                            break
                        else:
                            queue[queue_size] = col
                            queue_size += 1

            if match[i] == -1:
                for j in range(1, queue_size):
                    visited[match[queue[j]]] = -1

    return match
