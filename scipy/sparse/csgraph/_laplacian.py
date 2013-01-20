"""
Laplacian of a compressed-sparse graph
"""

# Authors: Aric Hagberg <hagberg@lanl.gov>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Jake Vanderplas <vanderplas@astro.washington.edu>
# License: BSD

from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.sparse import isspmatrix, coo_matrix


###############################################################################
# Graph laplacian
def laplacian(csgraph, normed=False, return_diag=False):
    """ Return the Laplacian matrix of a directed graph.

    For non-symmetric graphs the out-degree is used in the computation.

    Parameters
    ----------
    csgraph : array_like or sparse matrix, 2 dimensions
        compressed-sparse graph, with shape (N, N).
    normed : bool, optional
        If True, then compute normalized Laplacian.
    return_diag : bool, optional
        If True, then return diagonal as well as laplacian.

    Returns
    -------
    lap : ndarray
        The N x N laplacian matrix of graph.
    diag : ndarray
        The length-N diagonal of the laplacian matrix.
        diag is returned only if return_diag is True.

    Notes
    -----
    The Laplacian matrix of a graph is sometimes referred to as the
    "Kirchoff matrix" or the "admittance matrix", and is useful in many
    parts of spectral graph theory.  In particular, the eigen-decomposition
    of the laplacian matrix can give insight into many properties of the graph.

    For non-symmetric directed graphs, the laplacian is computed using the
    out-degree of each node.

    Examples
    --------
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
    """
    if csgraph.ndim != 2 or csgraph.shape[0] != csgraph.shape[1]:
        raise ValueError('csgraph must be a square matrix or array')

    if normed and (np.issubdtype(csgraph.dtype, np.int)
                   or np.issubdtype(csgraph.dtype, np.uint)):
        csgraph = csgraph.astype(np.float)

    if isspmatrix(csgraph):
        return _laplacian_sparse(csgraph, normed=normed,
                                 return_diag=return_diag)
    else:
        return _laplacian_dense(csgraph, normed=normed,
                                return_diag=return_diag)


def _laplacian_sparse(graph, normed=False, return_diag=False):
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
        diagonal_holes = list(set(range(n_nodes)).difference(
                                diag_idx))
        new_data = np.concatenate([lap.data, np.ones(len(diagonal_holes))])
        new_row = np.concatenate([lap.row, diagonal_holes])
        new_col = np.concatenate([lap.col, diagonal_holes])
        lap = coo_matrix((new_data, (new_row, new_col)), shape=lap.shape)
        diag_mask = (lap.row == lap.col)

    lap.data[diag_mask] = 0
    w = -np.asarray(lap.sum(axis=1)).squeeze()
    if normed:
        w = np.sqrt(w)
        w_zeros = (w == 0)
        w[w_zeros] = 1
        lap.data /= w[lap.row]
        lap.data /= w[lap.col]
        lap.data[diag_mask] = (1 - w_zeros[lap.row[diag_mask]]).astype(lap.data.dtype)
    else:
        lap.data[diag_mask] = w[lap.row[diag_mask]]

    if return_diag:
        return lap, w
    return lap


def _laplacian_dense(graph, normed=False, return_diag=False):
    n_nodes = graph.shape[0]
    lap = -np.asarray(graph)  # minus sign leads to a copy

    # set diagonal to zero
    lap.flat[::n_nodes + 1] = 0
    w = -lap.sum(axis=0)
    if normed:
        w = np.sqrt(w)
        w_zeros = (w == 0)
        w[w_zeros] = 1
        lap /= w
        lap /= w[:, np.newaxis]
        lap.flat[::n_nodes + 1] = 1 - w_zeros
    else:
        lap.flat[::n_nodes + 1] = w

    if return_diag:
        return lap, w
    return lap
