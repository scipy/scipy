import numpy as np
from numpy.testing import assert_array_almost_equal, run_module_suite
from scipy.sparse import csr_matrix


def check_csr_rowslice(X, Xcsr, i, sl):
    np_slice = X[i, sl]
    csr_slice = Xcsr[i, sl]
    assert_array_almost_equal(np_slice, csr_slice.toarray()[0])
    

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
            yield (check_csr_rowslice, X, Xcsr, i, sl)


def test_csr_getrow_getcol():
    N = 10
    np.random.seed(0)
    X = np.random.random((N, N))
    X[X > 0.7] = 0
    Xcsr = csr_matrix(X)
    
    for i in range(N):
        # check getrow()
        yield (assert_array_almost_equal,
               X[i:i + 1, :], Xcsr.getrow(i).toarray())
        # check getcol()
        yield (assert_array_almost_equal,
               X[:, i:i + 1], Xcsr.getcol(i).toarray())


if __name__ == "__main__":
    run_module_suite()
