"""
Laplacian of a compressed-sparse graph
"""

# Authors: Aric Hagberg <hagberg@lanl.gov>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Jake Vanderplas <vanderplas@astro.washington.edu>
# License: BSD

import numpy as np
from scipy.sparse import isspmatrix

###############################################################################
# Graph laplacian
def _graph_laplacian_sparse(graph, normed=False, return_diag=False):
    n_nodes = graph.shape[0]
    if not graph.format == 'coo':
        lap = (-graph).tocoo()
    else:
        lap = -graph.copy()
    diag_mask = (lap.row == lap.col)
    if not diag_mask.sum() == n_nodes:
        # The sparsity pattern of the matrix has holes on the diagonal,
        # we need to fix that
        diag_idx = lap.row[diag_mask]

        lap = lap.tolil()

        diagonal_holes = list(set(range(n_nodes)).difference(
                                diag_idx))
        lap[diagonal_holes, diagonal_holes] = 1
        lap = lap.tocoo()
        diag_mask = (lap.row == lap.col)
    lap.data[diag_mask] = 0
    w = -np.asarray(lap.sum(axis=1)).squeeze()
    if normed:
        w = np.sqrt(w)
        w_zeros = w == 0
        w[w_zeros] = 1
        lap.data /= w[lap.row]
        lap.data /= w[lap.col]
        lap.data[diag_mask] = (1 - w_zeros).astype(lap.data.dtype)
    else:
        lap.data[diag_mask] = w[lap.row[diag_mask]]
    if return_diag:
        return lap, w
    return lap


def _graph_laplacian_dense(graph, normed=False, return_diag=False):
    n_nodes = graph.shape[0]
    lap = -graph.copy()
    lap.flat[::n_nodes + 1] = 0
    w = -lap.sum(axis=0)
    if normed:
        w = np.sqrt(w)
        w_zeros = w == 0
        w[w_zeros] = 1
        lap /= w
        lap /= w[:, np.newaxis]
        lap.flat[::n_nodes + 1] = 1 - w_zeros
    else:
        lap.flat[::n_nodes + 1] = w
    if return_diag:
        return lap, w
    return lap


def laplacian(graph, normed=False, return_diag=False):
    """ Return the Laplacian matrix of the given graph.

    Parameters
    ----------
    graph: array-like or sparse matrix, shape=(N, N)
        directed compressed-sparse graph
    normed: boolean (optional)
        if True, then compute normalized Laplacian
    return_diag: boolean (optional)
        if True, then return diagonal as well as laplacian

    Returns
    -------
    lap: ndarray, shape=(N, N)
        the laplacian matrix of graph
    
    diag: ndarray, size=N [if return_diag == True]
        the diagonal of the laplacian matrix

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csgraph
    >>> G = np.arange(5) * np.arange(5)[:, np.newaxis]
    >>> G
    array([[ 0,  0,  0,  0,  0],
           [ 0,  1,  2,  3,  4],
           [ 0,  2,  4,  6,  8],
           [ 0,  3,  6,  9, 12],
           [ 0,  4,  8, 12, 16]])
    >>> csgraph.laplacian(G, normed=False)
    array([[  0,   0,   0,   0,   0],
           [  0,   9,  -2,  -3,  -4],
           [  0,  -2,  16,  -6,  -8],
           [  0,  -3,  -6,  21, -12],
           [  0,  -4,  -8, -12,  24]])

    Notes
    -----
    The Laplacian matrix of a graph is sometimes referred to as the
    "Kirchoff matrix" or the "admittance matrix", and is useful in many
    parts of spectral graph theory.  In particular, the eigen-decomposition
    of the laplacian matrix can give insight into many properties of the
    graph.
    """
    if normed and (np.issubdtype(graph.dtype, np.int)
                    or np.issubdtype(graph.dtype, np.uint)):
        graph = graph.astype(np.float)
    if isspmatrix(graph):
        return _graph_laplacian_sparse(graph, normed=normed,
                                       return_diag=return_diag)
    else:
        # We have a numpy array
        return _graph_laplacian_dense(graph, normed=normed,
                                       return_diag=return_diag)
