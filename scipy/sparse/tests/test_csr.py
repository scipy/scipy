from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_array_almost_equal, run_module_suite, assert_
from scipy.sparse import csr_matrix


def _check_csr_rowslice(i, sl, X, Xcsr):
    np_slice = X[i, sl]
    csr_slice = Xcsr[i, sl]
    assert_array_almost_equal(np_slice, csr_slice.toarray()[0])
    assert_(type(csr_slice) is csr_matrix)


def test_csr_rowslice():
    N = 10
    np.random.seed(0)
    X = np.random.random((N, N))
    X[X > 0.7] = 0
    Xcsr = csr_matrix(X)

    slices = [slice(None, None, None),
              slice(None, None, -1),
              slice(1, -2, 2),
              slice(-2, 1, -2)]

    for i in range(N):
        for sl in slices:
            yield _check_csr_rowslice, i, sl, X, Xcsr

def test_csr_getrowsequence():
    N = 3
    np.random.seed(0)
    X = np.random.random((N, N))

    # test over all 2 ** (N ** 2) sparsity patterns on NxN matrix
    for i in xrange(2 ** (N ** 2)):
        bits = [int(c) for c in bin(i)[2:]]
        if len(bits) < (N ** 2):
            bits = ([0] * (N ** 2 - len(bits))) + bits

        Z = X * np.array(bits, dtype=np.int).reshape((N, N))
        Zcsr = csr_matrix(Z)

        assert_array_almost_equal(Z[[0], :], Zcsr[[0], :].toarray())
        assert_array_almost_equal(Z[[2], :], Zcsr[[2], :].toarray())
        assert_array_almost_equal(Z[[0, 2], :], Zcsr[[0, 2], :].toarray())
        assert_array_almost_equal(Z[[1, 2], :], Zcsr[[1, 2], :].toarray())
        assert_array_almost_equal(Z[[2, 0], :], Zcsr[[2, 1], :].toarray())

def test_csr_getrow():
    N = 10
    np.random.seed(0)
    X = np.random.random((N, N))
    X[X > 0.7] = 0
    Xcsr = csr_matrix(X)

    for i in range(N):
        arr_row = X[i:i + 1, :]
        csr_row = Xcsr.getrow(i)

        assert_array_almost_equal(arr_row, csr_row.toarray())
        assert_(type(csr_row) is csr_matrix)


def test_csr_getcol():
    N = 10
    np.random.seed(0)
    X = np.random.random((N, N))
    X[X > 0.7] = 0
    Xcsr = csr_matrix(X)

    for i in range(N):
        arr_col = X[:, i:i + 1]
        csr_col = Xcsr.getcol(i)

        assert_array_almost_equal(arr_col, csr_col.toarray())
        assert_(type(csr_col) is csr_matrix)


if __name__ == "__main__":
    run_module_suite()
