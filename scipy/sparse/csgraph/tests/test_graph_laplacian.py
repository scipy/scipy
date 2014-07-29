# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
#         Jake Vanderplas <vanderplas@astro.washington.edu>
# License: BSD
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (run_module_suite, assert_allclose,
        assert_array_almost_equal)
from scipy import sparse

from scipy.sparse import csgraph


def _explicit_laplacian(x, normed=False):
    if sparse.issparse(x):
        x = x.todense()
    x = np.asarray(x)
    y = -1.0 * x
    for j in range(y.shape[0]):
        y[j,j] = x[j,j+1:].sum() + x[j,:j].sum()
    if normed:
        d = np.diag(y).copy()
        d[d == 0] = 1.0
        y /= d[:,None]**.5
        y /= d[None,:]**.5
    return y


def _check_graph_laplacian(mat, normed):
    if not hasattr(mat, 'shape'):
        mat = eval(mat, dict(np=np, sparse=sparse))

    if sparse.issparse(mat):
        sp_mat = mat
        mat = sp_mat.todense()
    else:
        sp_mat = sparse.csr_matrix(mat)

    laplacian = csgraph.laplacian(mat, normed=normed)
    n_nodes = mat.shape[0]
    if not normed:
        assert_array_almost_equal(laplacian.sum(axis=0), np.zeros(n_nodes))
    assert_array_almost_equal(laplacian.T, laplacian)
    assert_array_almost_equal(laplacian,
        csgraph.laplacian(sp_mat, normed=normed).todense())

    assert_array_almost_equal(laplacian,
        _explicit_laplacian(mat, normed=normed))


def test_graph_laplacian():
    mats = ('np.arange(10) * np.arange(10)[:, np.newaxis]',
            'np.ones((7, 7))',
            'np.eye(19)',
            'sparse.diags([1, 1], [-1, 1], shape=(4,4))',
            'sparse.diags([1, 1], [-1, 1], shape=(4,4)).todense()',
            'np.asarray(sparse.diags([1, 1], [-1, 1], shape=(4,4)).todense())',
            'np.vander(np.arange(4)) + np.vander(np.arange(4)).T',
            )

    for mat_str in mats:
        for normed in (True, False):
            yield _check_graph_laplacian, mat_str, normed


def _assert_allclose_sparse(a, b, **kwargs):
    # helper function that can deal with sparse matrices
    if sparse.issparse(a):
        a = a.toarray()
    if sparse.issparse(b):
        b = a.toarray()
    assert_allclose(a, b, **kwargs)


def test_asymmetric_laplacian():
    # adjacency matrix, laplacian, normalized laplacian
    A = [[0, 1, 0],
         [4, 2, 0],
         [0, 0, 0]]
    L = [[1, -1, 0],
         [-4, 4, 0],
         [0, 0, 0]]
    Ld = [1, 4, 0]
    M = [[1, -0.5, 0],
         [-2, 1, 0],
         [0, 0, 0]]
    Md = [1, 2, 1]
    for arr_type in np.array, sparse.csr_matrix, sparse.coo_matrix:
        for t in int, float, complex:
            adj = arr_type(A, dtype=t)
            # check laplacian
            m = csgraph.laplacian(adj, normed=False, return_diag=False)
            _assert_allclose_sparse(m, L, atol=1e-12)
            m, d = csgraph.laplacian(adj, normed=False, return_diag=True)
            _assert_allclose_sparse(m, L, atol=1e-12)
            _assert_allclose_sparse(d, Ld, atol=1e-12)
            # check normalized laplacian
            m = csgraph.laplacian(adj, normed=True, return_diag=False)
            _assert_allclose_sparse(m, M, atol=1e-12)
            m, d = csgraph.laplacian(adj, normed=True, return_diag=True)
            _assert_allclose_sparse(m, M, atol=1e-12)
            _assert_allclose_sparse(d, Md, atol=1e-12)


if __name__ == '__main__':
    run_module_suite()
