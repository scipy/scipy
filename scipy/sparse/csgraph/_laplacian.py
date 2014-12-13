"""
Laplacian of a compressed-sparse graph
"""

# Authors: Aric Hagberg <hagberg@lanl.gov>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Jake Vanderplas <vanderplas@astro.washington.edu>
# License: BSD

from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.sparse import isspmatrix


###############################################################################
# Graph laplacian
def laplacian(csgraph, normed=False, return_diag=False, use_out_degree=False):
    """
    Return the Laplacian matrix of a directed graph.

    Parameters
    ----------
    csgraph : array_like or sparse matrix, 2 dimensions
        compressed-sparse graph, with shape (N, N).
    normed : bool, optional
        If True, then compute normalized Laplacian.
    return_diag : bool, optional
        If True, then also return an array related to vertex degrees.
    use_out_degree : bool, optional
        If True, then use out-degree instead of in-degree.
        This distinction matters only if the graph is asymmetric.
        Default: False.

    Returns
    -------
    lap : ndarray
        The N x N laplacian matrix of graph.
    diag : ndarray, optional
        The length-N diagonal of the Laplacian matrix.
        For the normalized Laplacian, this is the array of square roots
        of vertex degrees or 1 if the degree is zero.

    Notes
    -----
    The Laplacian matrix of a graph is sometimes referred to as the
    "Kirchoff matrix" or the "admittance matrix", and is useful in many
    parts of spectral graph theory.  In particular, the eigen-decomposition
    of the laplacian matrix can give insight into many properties of the graph.

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

    create_lap = _laplacian_sparse if isspmatrix(csgraph) else _laplacian_dense
    degree_axis = 1 if use_out_degree else 0
    lap, d = create_lap(csgraph, normed=normed, axis=degree_axis)
    if return_diag:
        return lap, d
    return lap


def _setdiag_dense(A, d):
    A.flat[::len(d)+1] = d


def _laplacian_sparse(graph, normed=False, axis=0):
    if graph.format == 'coo':
        m = graph.copy()
    else:
        m = graph.tocoo()
    w = m.sum(axis=axis).getA1() - m.diagonal()
    if normed:
        w = np.sqrt(w)
        isolated_node_mask = (w == 0)
        w[isolated_node_mask] = 1
        m.data /= w[m.row]
        m.data /= w[m.col]
        m.data *= -1
        m.setdiag(1 - isolated_node_mask)
    else:
        m.data *= -1
        m.setdiag(w)
    return m, w


def _laplacian_dense(graph, normed=False, axis=0):
    m = np.array(graph)
    np.fill_diagonal(m, 0)
    w = m.sum(axis=axis)
    if normed:
        w = np.sqrt(w)
        isolated_node_mask = (w == 0)
        w[isolated_node_mask] = 1
        m /= w
        m /= w[:, np.newaxis]
        m *= -1
        _setdiag_dense(m, 1 - isolated_node_mask)
    else:
        m *= -1
        _setdiag_dense(m, w)
    return m, w
