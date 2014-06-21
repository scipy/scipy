# Author: Paul Nation  -- <pnation@korea.ac.kr>
# Original Source: QuTiP: Quantum Toolbox in Python (qutip.org)
# License: New BSD, (C) 2014

import numpy as np
cimport numpy as np

from scipy.sparse import isspmatrix_csr, isspmatrix_csc

include 'parameters.pxi'

def symrcm(graph, symmetric_mode=False):
    """
    Returns the permutation array that orders a sparse CSR or CSC matrix
    in Reverse-Cuthill McKee ordering.  Since the input matrix must be 
    symmetric, this routine works on the matrix A+Trans(A) if the 
    symmetric_mode flag is set to False (default).

    It is assumed by default that the input matrix is not symmetric.  
    This is because it is faster to do A+Trans(A) than it is to check for 
    symmetry for a generic matrix.  If you are guaranteed that the matrix 
    is symmetric in structure (values of matrix element do not matter) then 
    set symmetric_mode=True.
    
    Parameters
    ----------
    graph : sparse matrix
        Input sparse CSC or CSR sparse matrix format.
    symmetric_mode : bool {False, True}
        Is input matrix guaranteed to be symmetric.

    Returns
    -------
    perm : array_like
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
    

def _node_degrees(
        np.ndarray[ITYPE_t, ndim=1, mode="c"] ind,
        np.ndarray[ITYPE_t, ndim=1, mode="c"] ptr,
        int num_rows):
    """
    Find the degree of each node (matrix row) in a graph represented
    by a sparse CSR or CSC matrix.
    """
    cdef unsigned int ii, jj
    cdef np.ndarray[ITYPE_t] degree = np.zeros(num_rows, dtype=ITYPE)
    
    for ii in range(num_rows):
        degree[ii] = ptr[ii + 1] - ptr[ii]
        for jj in range(ptr[ii], ptr[ii + 1]):
            if ind[jj] == ii:
                # add one if the diagonal is in row ii
                degree[ii] += 1
                break
    return degree
    

def _reverse_cuthill_mckee(np.ndarray[ITYPE_t, ndim=1, mode="c"] ind,
        np.ndarray[ITYPE_t, ndim=1, mode="c"] ptr,
        int num_rows):
    """
    Reverse Cuthill-McKee ordering of a sparse symmetric CSR or CSC matrix.  
    We follow the original Cuthill-McKee paper and always start the routine
    at a node of lowest degree for each connected component.
    """
    cdef unsigned int N = 0, N_old, seed, level_start, level_end, temp, temp2
    cdef unsigned int zz, i, j, ii, jj, kk, ll
    cdef np.ndarray[ITYPE_t] order = np.zeros(num_rows, dtype=ITYPE)
    cdef np.ndarray[ITYPE_t] degree = _node_degrees(ind, ptr, num_rows).astype(ITYPE)
    cdef np.ndarray[ITYPE_t] inds = np.argsort(degree).astype(ITYPE)
    cdef np.ndarray[ITYPE_t] rev_inds = np.argsort(inds).astype(ITYPE)
    cdef np.ndarray[ITYPE_t] temp_degrees = np.zeros(num_rows, dtype=ITYPE)

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
                    temp = 0
                    for kk in range(N_old, N):
                        temp_degrees[temp] = degree[order[kk]]
                        temp += 1

                    # Do insertion sort for nodes from lowest to highest degree
                    for kk in range(N_old, N - 1):
                        temp = temp_degrees[kk - N_old]
                        temp2 = order[kk]
                        ll = kk - 1
                        while (ll >= N_old) and (temp_degrees[ll] > temp):
                            temp_degrees[ll + 1 - N_old] = temp_degrees[ll - N_old]
                            order[ll + 1] = order[ll]
                            ll -= 1
                        temp_degrees[ll + 1 - N_old] = temp
                        order[ll + 1] = temp2

                # set next level start and end ranges
                level_start = level_end
                level_end = N

        if N == num_rows:
            break

    # return reveresed order for RCM ordering
    return order[::-1]  