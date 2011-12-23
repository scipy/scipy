"""
Graph utilities and algorithms

Graphs are represented with their adjacency matrices, preferably using
sparse matrices.
"""

# Authors: Aric Hagberg <hagberg@lanl.gov>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD

import numpy as np
from scipy import sparse

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


def graph_laplacian(graph, normed=False, return_diag=False):
    """ Return the Laplacian of the given graph.
    """
    if normed and (np.issubdtype(graph.dtype, np.int)
                    or np.issubdtype(graph.dtype, np.uint)):
        graph = graph.astype(np.float)
    if sparse.isspmatrix(graph):
        return _graph_laplacian_sparse(graph, normed=normed,
                                       return_diag=return_diag)
    else:
        # We have a numpy array
        return _graph_laplacian_dense(graph, normed=normed,
                                       return_diag=return_diag)
