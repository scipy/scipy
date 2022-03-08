"""
Laplacian of a compressed-sparse graph
"""

import numpy as np
from scipy.sparse import isspmatrix
from scipy.sparse.linalg import LinearOperator


###############################################################################
# Graph laplacian
def laplacian(
    csgraph,
    normed=False,
    return_diag=False,
    use_out_degree=False,
    *,
    copy=True,
    form="array",
    dtype=None,
    symmetrized=False
):
    """
    Return the Laplacian of a directed graph.

    Parameters
    ----------
    csgraph : array_like or sparse matrix, 2 dimensions
        compressed-sparse graph, with shape (N, N).
    normed : bool, optional
        If True, then compute symmetrically normalized Laplacian.
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
    form: `array`, or `function`, or `lo`
        Determines the format of the output Laplacian:
        `array` is a numpy array;
        `function` is a pointer to evaluating the Laplacian-vector
        or Laplacian-matrix product;
        `lo` results in the format of the `LinearOperator`.
        Choosing `function` or `lo` always avoids doubling
        the memory use, ignoring `copy` value.
        Default: `array`, for backward compatibility.
    dtype: None or one of numeric numpy dtypes, optional
        The dtype of the output. If `dtype=None`, the dtype of the
        output matches the dtype of the input csgraph, except for
        the case `normed=True` and integer-like csgraph, where
        the output dtype is 'float' allowing accurate normalization,
        but dramatically increasing the memory use.
        Default: None, for backward compatibility.
    symmetrized: bool, optional
        If True, then the output Laplacian is symmetric/Hermitian.
        The symmetrization is done by `csgraph + csgraph.T.conj`
        without dividing by 2 to preserve integer dtypes if possible
        prior to the construction of the Laplacian.
        The symmetrization will increase the memory footprint of
        sparse matrices unless the sparsity pattern is symmetric or
        `form` is `function` or `lo`.
        Default: False, for backward compatibility.

    Returns
    -------
    lap : ndarray, or sparse matrix, or `LinearOperator`
        The N x N Laplacian of csgraph. It will be a NumPy array (dense)
        if the input was dense, or a sparse matrix otherwise, or
        the format of a function or `LinearOperator` if
        `form` equals `function` or `lo`, correspondingly.
    diag : ndarray, optional
        The length-N main diagonal of the Laplacian matrix.
        For the normalized Laplacian, this is the array of square roots
        of vertex degrees or 1 if the degree is zero.

    Notes
    -----
    The Laplacian matrix of a graph is sometimes referred to as the
    "Kirchhoff matrix" or just the "Laplacian", and is useful in many
    parts of spectral graph theory.
    In particular, the eigen-decomposition of the Laplacian can give
    insight into many properties of the graph, e.g.,
    is commonly used for spectral data embedding and clustering.

    The constructed Laplacian doubles the memory use if ``copy=True`` and
    `form=`array`` which is the default.
    Choosing ``copy=False`` has no effect unless `form=`array``
    or the matrix is sparse in the ``coo`` format, or dense array, except
    for the integer input with ``normed=True`` that forces the float output.

    Sparse input is reformatted into ``coo`` if `form=`array``,
    which is the default.

    If the input adjacency matrix is not symmetic, the Laplacian is also
    non-symmetric and may need to be symmetrized; e.g., ``lap += lap.T``,
    before the eigen-decomposition.

    Diagonal entries of the input adjacency matrix are ignored and
    replaced with zeros for the purpose of normalization where ``normed=True``.
    The normalization uses the inverse square roots of row-sums of the input
    adjacency matrix, and thus may fail if the row-sums contain
    zeros, negative, or complex with a non-zero imaginary part values.

    The normalization is symmetric, making the normalized Laplacian also
    symmetric if the input csgraph was symmetric.

    References
    ----------
    .. [1] Laplacian matrix. https://en.wikipedia.org/wiki/Laplacian_matrix

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
    >>> csgraph.laplacian(G)
    array([[  0,   0,   0,   0,   0],
           [  0,   9,  -2,  -3,  -4],
           [  0,  -2,  16,  -6,  -8],
           [  0,  -3,  -6,  21, -12],
           [  0,  -4,  -8, -12,  24]])
    >>> G = np.arange(9).reshape(3, 3)
    >>> G
    array([[0, 1, 2],
           [3, 4, 5],
           [6, 7, 8]])
    >>> L_in_degree = csgraph.laplacian(G)
    >>> L_in_degree
    array([[ 9, -1, -2],
           [-3,  8, -5],
           [-6, -7,  7]])
    >>> L_out_degree = csgraph.laplacian(G, use_out_degree=True)
    >>> L_out_degree
    array([[ 3, -1, -2],
           [-3,  8, -5],
           [-6, -7, 13]])
    >>> L_in_degree + L_out_degree.T
    array([[ 12,  -4,  -8],
            [ -4,  16, -12],
            [ -8, -12,  20]])
    >>> csgraph.laplacian(G, symmetrized=True)
    array([[ 12,  -4,  -8],
           [ -4,  16, -12],
           [ -8, -12,  20]])
    >>> csgraph.laplacian(G + G.T)
    array([[ 12,  -4,  -8],
           [ -4,  16, -12],
           [ -8, -12,  20]])
    """
    if csgraph.ndim != 2 or csgraph.shape[0] != csgraph.shape[1]:
        raise ValueError("csgraph must be a square matrix or array")

    if normed and (
        np.issubdtype(csgraph.dtype, np.signedinteger)
        or np.issubdtype(csgraph.dtype, np.uint)
    ):
        csgraph = csgraph.astype(np.float64)

    if form == "array":
        create_lap = (
            _laplacian_sparse if isspmatrix(csgraph) else _laplacian_dense
        )
    else:
        create_lap = (
            _laplacian_sparse_flo
            if isspmatrix(csgraph)
            else _laplacian_dense_flo
        )

    degree_axis = 1 if use_out_degree else 0

    lap, d = create_lap(
        csgraph,
        normed=normed,
        axis=degree_axis,
        copy=copy,
        form=form,
        dtype=dtype,
        symmetrized=symmetrized,
    )
    if return_diag:
        return lap, d
    return lap


def _setdiag_dense(m, d):
    step = len(d) + 1
    m.flat[::step] = d


def _laplace(m, d):
    return lambda v: v * d[np.newaxis, :] - m @ v


def _laplace_normed(m, d, nd):
    md_sym = _laplace(m, d)
    return lambda v: nd * md_sym(v) * nd[:, np.newaxis]


def _laplace_sym(m, d):
    return (
        lambda v: v * d[np.newaxis, :]
        - m @ v
        - np.transpose(np.conjugate(np.transpose(np.conjugate(v)) @ m))
    )


def _laplace_normed_sym(m, d, nd):
    md_sym = _laplace_sym(m, d)
    return lambda v: nd * md_sym(v) * nd[:, np.newaxis]


def _linearoperator(mv, shape, dtype):
    return LinearOperator(matvec=mv, matmat=mv, shape=shape, dtype=dtype)


def _laplacian_sparse_flo(graph, normed, axis, copy, form, dtype, symmetrized):
    # The keyword argument ``copy`` is unused and has no effect here.

    if dtype is None:
        dtype = graph.dtype

    graph_sum = graph.sum(axis=axis).getA1()
    graph_diagonal = graph.diagonal()
    diag = graph_sum - graph_diagonal
    if symmetrized:
        graph_sum += graph.sum(axis=1 - axis).getA1()
        diag = graph_sum - graph_diagonal - graph_diagonal

    if normed:
        isolated_node_mask = diag == 0
        w = np.where(isolated_node_mask, 1, np.sqrt(diag))
        if symmetrized:
            md = _laplace_normed_sym(graph, graph_sum, 1.0 / w)
        else:
            md = _laplace_normed(graph, graph_sum, 1.0 / w)
        if form == "function":
            return md, w.astype(dtype, copy=False)
        elif form == "lo":
            m = _linearoperator(md, shape=graph.shape, dtype=dtype)
            return m, w.astype(dtype, copy=False)
        else:
            raise ValueError(f"Invalid form: {repr(form)}")
    else:
        if symmetrized:
            md = _laplace_sym(graph, graph_sum)
        else:
            md = _laplace(graph, graph_sum)
        if form == "function":
            return md, diag.astype(dtype, copy=False)
        elif form == "lo":
            m = _linearoperator(md, shape=graph.shape, dtype=dtype)
            return m, diag.astype(dtype, copy=False)
        else:
            raise ValueError(f"Invalid form: {repr(form)}")


def _laplacian_sparse(graph, normed, axis, copy, form, dtype, symmetrized):

    if form != "array":
        raise ValueError('{form} must be "array"'.format(form=repr(form)))

    if dtype is None:
        dtype = graph.dtype

    needs_copy = False
    if graph.format in ("lil", "dok"):
        m = graph.tocoo()
    else:
        m = graph
        if copy:
            needs_copy = True

    if symmetrized:
        m += m.T.conj()

    w = m.sum(axis=axis).getA1() - m.diagonal()
    if normed:
        m = m.tocoo(copy=needs_copy)
        isolated_node_mask = w == 0
        w = np.where(isolated_node_mask, 1, np.sqrt(w))
        m.data /= w[m.row]
        m.data /= w[m.col]
        m.data *= -1
        m.setdiag(1 - isolated_node_mask)
    else:
        if m.format == "dia":
            m = m.copy()
        else:
            m = m.tocoo(copy=needs_copy)
        m.data *= -1
        m.setdiag(w)

    return m.astype(dtype, copy=False), w.astype(dtype)


def _laplacian_dense_flo(graph, normed, axis, copy, form, dtype, symmetrized):

    if dtype is None:
        dtype = graph.dtype

    if copy:
        m = np.array(graph)
    else:
        m = np.asarray(graph)

    if dtype is None:
        dtype = m.dtype

    graph_sum = m.sum(axis=axis)
    graph_diagonal = m.diagonal()
    diag = graph_sum - graph_diagonal
    if symmetrized:
        graph_sum += m.sum(axis=1 - axis)
        diag = graph_sum - graph_diagonal - graph_diagonal

    if normed:
        isolated_node_mask = diag == 0
        w = np.where(isolated_node_mask, 1, np.sqrt(diag))
        if symmetrized:
            md = _laplace_normed_sym(m, graph_sum, 1.0 / w)
        else:
            md = _laplace_normed(m, graph_sum, 1.0 / w)
        if form == "function":
            return md, w.astype(dtype, copy=False)
        elif form == "lo":
            m = _linearoperator(md, shape=graph.shape, dtype=dtype)
            return m, w.astype(dtype, copy=False)
        else:
            raise ValueError(f"Invalid form: {repr(form)}")
    else:
        if symmetrized:
            md = _laplace_sym(m, graph_sum)
        else:
            md = _laplace(m, graph_sum)
        if form == "function":
            return md, diag.astype(dtype, copy=False)
        elif form == "lo":
            m = _linearoperator(md, shape=graph.shape, dtype=dtype)
            return m, diag.astype(dtype, copy=False)
        else:
            raise ValueError(f"Invalid form: {repr(form)}")


def _laplacian_dense(graph, normed, axis, copy, form, dtype, symmetrized):

    if form != "array":
        raise ValueError('{form} must be "array"'.format(form=repr(form)))

    if dtype is None:
        dtype = graph.dtype

    if copy:
        m = np.array(graph)
    else:
        m = np.asarray(graph)

    if dtype is None:
        dtype = m.dtype

    if symmetrized:
        m += m.T.conj()
    np.fill_diagonal(m, 0)
    w = m.sum(axis=axis)
    if normed:
        isolated_node_mask = w == 0
        w = np.where(isolated_node_mask, 1, np.sqrt(w))
        m /= w
        m /= w[:, np.newaxis]
        m *= -1
        _setdiag_dense(m, 1 - isolated_node_mask)
    else:
        m *= -1
        _setdiag_dense(m, w)

    return m.astype(dtype, copy=False), w.astype(dtype, copy=False)
