import numpy as np
from numpy.testing import assert_array_almost_equal, run_module_suite
from scipy.sparse import csr_matrix


def _check_csr_rowslice(i, sl, X, Xcsr):
    np_slice = X[i, sl]
    csr_slice = Xcsr[i, sl]
    assert_array_almost_equal(np_slice, csr_slice.toarray()[0])
    assert type(csr_slice) is csr_matrix


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


def _check_csr_getrow(i, X, Xcsr):
    arr_row = X[i:i + 1, :]
    csr_row = Xcsr.getrow(i)

    assert_array_almost_equal(arr_row, csr_row.toarray())
    assert type(csr_row) is csr_matrix


def _check_csr_getcol(i, X, Xcsr):
    arr_col = X[:, i:i + 1]
    csr_col = Xcsr.getcol(i)

    assert_array_almost_equal(arr_col, csr_col.toarray())
    assert type(csr_col) is csr_matrix


def test_csr_getrow_getcol():
    N = 10
    np.random.seed(0)
    X = np.random.random((N, N))
    X[X > 0.7] = 0
    Xcsr = csr_matrix(X)

    for i in range(N):
        yield _check_csr_getrow, i, X, Xcsr
        yield _check_csr_getcol, i, X, Xcsr


if __name__ == "__main__":
    run_module_suite()
