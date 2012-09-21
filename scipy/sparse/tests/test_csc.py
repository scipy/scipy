import numpy as np
from numpy.testing import assert_array_almost_equal, run_module_suite
from scipy.sparse import csr_matrix, csc_matrix


def _check_csc_getrow(i, X, Xcsc):
    arr_row = X[i:i + 1, :]
    csc_row = Xcsc.getrow(i)

    assert_array_almost_equal(arr_row, csc_row.toarray())
    assert type(csc_row) is csr_matrix


def _check_csc_getcol(i, X, Xcsc):
    arr_col = X[:, i:i + 1]
    csc_col = Xcsc.getcol(i)

    assert_array_almost_equal(arr_col, csc_col.toarray())
    assert type(csc_col) is csc_matrix


def test_csc_getrow_getcol():
    N = 10
    np.random.seed(0)
    X = np.random.random((N, N))
    X[X > 0.7] = 0
    Xcsc = csc_matrix(X)

    for i in range(N):
        yield _check_csc_getrow, i, X, Xcsc
        yield _check_csc_getcol, i, X, Xcsc


if __name__ == "__main__":
    run_module_suite()
