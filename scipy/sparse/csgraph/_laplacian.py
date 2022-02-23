"""
Laplacian of a compressed-sparse graph
"""

# Authors: Aric Hagberg <hagberg@lanl.gov>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Jake Vanderplas <vanderplas@astro.washington.edu>
#          Andrew Knyazev <Andrew.Knyazev@ucdenver.edu>
# License: BSD

import numpy as np
from scipy.sparse import isspmatrix


###############################################################################
# Graph laplacian
def laplacian(csgraph, normed=False, return_diag=False, use_out_degree=False,
              *, copy=True):
    """
    Return the Laplacian matrix of a directed graph.

    Parameters
    ----------
    csgraph : array_like or sparse matrix, 2 dimensions
        compressed-sparse graph, with shape (N, N).
    normed : bool, optional
        If True, then compute symmetric normalized Laplacian.
        Default: False.
    return_diag : bool, optional
        If True, then also return an array related to vertex degrees.
        Default: False.
    use_out_degree : bool, optional
        If True, then use out-degree instead of in-degree.
        This distinction matters only if the graph is asymmetric.
        Default: False.
    copy: bool, optional
        If False, then change `csgraph` in place if possible,
        avoiding doubling the memory use.
        Default: True, for backward compatibility.

    Returns
    -------
    lap : ndarray or sparse matrix
        The N x N laplacian matrix of csgraph. It will be a NumPy array (dense)
        if the input was dense, or a sparse matrix otherwise.
    diag : ndarray, optional
        The length-N diagonal of the Laplacian matrix.
        For the normalized Laplacian, this is the array of square roots
        of vertex degrees or 1 if the degree is zero.

    Notes
    -----
    The Laplacian matrix of a graph is sometimes referred to as the
    "Kirchoff matrix" or just the "Laplacian", and is useful in many
    parts of spectral graph theory. In particular, the eigen-decomposition
    of the Laplacian can give insight into many properties of the graph, e.g.,
    is commonly used for spectal data enmedding and clustering.

    The constructed Laplacian doubles the memory use if ``copy=True``,
    which is the default. Choosing ``copy=False`` has no effect unless
    the matrix is sparse in the ``coo`` format, or dense array, except
    for the integer input with ``normed=True`` that forces the float output.

    Sparse input is reformatted into ``coo``.

    If the input adjacency matrix is not symmetic, the Laplacian is also
    non-symmetric and may need to be symmetrized; e.g., ``lap += lap.T``,
    before the eigen-decomposition.

    Diagonal entries of the input adjacency matrix are ignored and
    replaced with zeros for the purpose of normalization where ``normed=True``.
    The normalization uses the inverse square roots of row-sums of the input
    adjacency matrix, and thus may fail if the row-sums contain
    zeros, negative, or complex with a non-zero imaginary part values.

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

    if normed and (np.issubdtype(csgraph.dtype, np.signedinteger)
                   or np.issubdtype(csgraph.dtype, np.uint)):
        csgraph = csgraph.astype(float)

    create_lap = _laplacian_sparse if isspmatrix(csgraph) else _laplacian_dense
    degree_axis = 1 if use_out_degree else 0
    lap, d = create_lap(csgraph, normed=normed, axis=degree_axis,
                        copy=copy)
    if return_diag:
        return lap, d
    return lap


def _setdiag_dense(A, d):
    A.flat[::len(d)+1] = d


def _laplacian_sparse(graph, normed=False, axis=0,
                      copy=True):
    needs_copy = False
    if graph.format in ('lil', 'dok'):
        m = graph.tocoo()
    else:
        m = graph
        if copy:
            needs_copy = True
    w = m.sum(axis=axis).getA1() - m.diagonal()
    if normed:
        m = m.tocoo(copy=needs_copy)
        isolated_node_mask = (w == 0)
        w = np.where(isolated_node_mask, 1, np.sqrt(w))
        m.data /= w[m.row]
        m.data /= w[m.col]
        m.data *= -1
        m.setdiag(1 - isolated_node_mask)
    else:
        if m.format == 'dia':
            m = m.copy()
        else:
            m = m.tocoo(copy=needs_copy)
        m.data *= -1
        m.setdiag(w)
    return m, w


def _laplacian_dense(graph, normed=False, axis=0,
                     copy=True):
    if copy:
        m = np.array(graph)
    else:
        m = np.asarray(graph)

    np.fill_diagonal(m, 0)
    w = m.sum(axis=axis)
    if normed:
        isolated_node_mask = (w == 0)
        w = np.where(isolated_node_mask, 1, np.sqrt(w))
        m /= w
        m /= w[:, np.newaxis]
        m *= -1
        _setdiag_dense(m, 1 - isolated_node_mask)
    else:
        m *= -1
        _setdiag_dense(m, w)
    return m, w
