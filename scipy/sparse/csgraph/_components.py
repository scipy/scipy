from __future__ import division, print_function, absolute_import

import numpy as np

from scipy.sparse.sparsetools import cs_graph_components as _cs_graph_components

from scipy.sparse.csr import csr_matrix
from scipy.sparse.base import isspmatrix

_msg0 = 'x must be a symmetric square matrix!'
_msg1 = _msg0 + '(has shape %s)'


def cs_graph_components(x):
    """
    Determine connected components of a graph stored as a compressed
    sparse row or column matrix.

    For speed reasons, the symmetry of the matrix x is not checked. A
    nonzero at index `(i, j)` means that node `i` is connected to node
    `j` by an edge. The number of rows/columns of the matrix thus
    corresponds to the number of nodes in the graph.

    Parameters
    -----------
    x : array_like or sparse matrix, 2 dimensions
        The adjacency matrix of the graph. Only the upper triangular part
        is used.

    Returns
    --------
    n_comp : int
        The number of connected components.
    label : ndarray (ints, 1 dimension):
        The label array of each connected component (-2 is used to
        indicate empty rows in the matrix: 0 everywhere, including
        diagonal). This array has the length of the number of nodes,
        i.e. one label for each node of the graph. Nodes having the same
        label belong to the same connected component.

    Notes
    ------
    The matrix is assumed to be symmetric and the upper triangular part
    of the matrix is used. The matrix is converted to a CSR matrix unless
    it is already a CSR.

    Examples
    --------
    >>> from scipy.sparse.csgraph import connected_components
    >>> D = np.eye(4)
    >>> D[0,1] = D[1,0] = 1
    >>> cs_graph_components(D)
    (3, array([0, 0, 1, 2]))
    >>> from scipy.sparse import dok_matrix
    >>> cs_graph_components(dok_matrix(D))
    (3, array([0, 0, 1, 2]))

    """
    try:
        shape = x.shape
    except AttributeError:
        raise ValueError(_msg0)

    if not ((len(x.shape) == 2) and (x.shape[0] == x.shape[1])):
        raise ValueError(_msg1 % x.shape)

    if isspmatrix(x):
        x = x.tocsr()
    else:
        x = csr_matrix(x)

    label = np.empty((shape[0],), dtype=x.indptr.dtype)

    n_comp = _cs_graph_components(shape[0], x.indptr, x.indices, label)

    return n_comp, label
