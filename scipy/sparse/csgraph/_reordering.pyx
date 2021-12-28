# Author: Paul Nation  -- <nonhermitian@gmail.com>
# Original Source: QuTiP: Quantum Toolbox in Python (qutip.org)
# License: New BSD, (C) 2014

import numpy as np
cimport numpy as np
from warnings import warn
from scipy.sparse import (csc_matrix, csr_matrix, isspmatrix, isspmatrix_coo,
                          isspmatrix_csc, isspmatrix_csr,
                          SparseEfficiencyWarning)
from . import maximum_bipartite_matching

np.import_array()

include 'parameters.pxi'

def reverse_cuthill_mckee(graph, symmetric_mode=False):
    """
    reverse_cuthill_mckee(graph, symmetric_mode=False)
    
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
 
    Notes
    -----
    .. versionadded:: 0.15.0

    References
    ----------
    E. Cuthill and J. McKee, "Reducing the Bandwidth of Sparse Symmetric Matrices",
    ACM '69 Proceedings of the 1969 24th national conference, (1969).

    Examples
    --------
    >>> from scipy.sparse import csr_matrix
    >>> from scipy.sparse.csgraph import reverse_cuthill_mckee

    >>> graph = [
    ... [0, 1, 2, 0],
    ... [0, 0, 0, 1],
    ... [2, 0, 0, 3],
    ... [0, 0, 0, 0]
    ... ]
    >>> graph = csr_matrix(graph)
    >>> print(graph)
      (0, 1)	1
      (0, 2)	2
      (1, 3)	1
      (2, 0)	2
      (2, 3)	3

    >>> reverse_cuthill_mckee(graph)
    array([3, 2, 1, 0], dtype=int32)
    
    """
    if not (isspmatrix_csc(graph) or isspmatrix_csr(graph)):
        raise TypeError('Input must be in CSC or CSR sparse matrix format.')
    nrows = graph.shape[0]
    if not symmetric_mode:
        graph = graph+graph.transpose()
    return _reverse_cuthill_mckee(graph.indices, graph.indptr, nrows)


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


def structural_rank(graph):
    """
    structural_rank(graph)
    
    Compute the structural rank of a graph (matrix) with a given 
    sparsity pattern.

    The structural rank of a matrix is the number of entries in the maximum 
    transversal of the corresponding bipartite graph, and is an upper bound 
    on the numerical rank of the matrix. A graph has full structural rank 
    if it is possible to permute the elements to make the diagonal zero-free.

    .. versionadded:: 0.19.0

    Parameters
    ----------
    graph : sparse matrix
        Input sparse matrix.

    Returns
    -------
    rank : int
        The structural rank of the sparse graph.
    
    References
    ----------
    .. [1] I. S. Duff, "Computing the Structural Index", SIAM J. Alg. Disc. 
            Meth., Vol. 7, 594 (1986).
    
    .. [2] http://www.cise.ufl.edu/research/sparse/matrices/legend.html

    Examples
    --------
    >>> from scipy.sparse import csr_matrix
    >>> from scipy.sparse.csgraph import structural_rank

    >>> graph = [
    ... [0, 1, 2, 0],
    ... [1, 0, 0, 1],
    ... [2, 0, 0, 3],
    ... [0, 1, 3, 0]
    ... ]
    >>> graph = csr_matrix(graph)
    >>> print(graph)
      (0, 1)	1
      (0, 2)	2
      (1, 0)	1
      (1, 3)	1
      (2, 0)	2
      (2, 3)	3
      (3, 1)	1
      (3, 2)	3

    >>> structural_rank(graph)
    4

    """
    if not isspmatrix:
        raise TypeError('Input must be a sparse matrix')
    if not isspmatrix_csr(graph):
        if not (isspmatrix_csc(graph) or isspmatrix_coo(graph)):
            warn('Input matrix should be in CSC, CSR, or COO matrix format',
                    SparseEfficiencyWarning)
        graph = csr_matrix(graph)
    # If A is a tall matrix, then transpose.
    if graph.shape[0] > graph.shape[1]:
        graph = graph.T.tocsr()
    rank = np.sum(maximum_bipartite_matching(graph) >= 0)
    return rank
