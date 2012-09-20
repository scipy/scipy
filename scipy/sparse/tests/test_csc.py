import numpy as np
from numpy.testing import assert_array_almost_equal, run_module_suite
from scipy.sparse import csc_matrix


def test_csc_getrow_getcol():
    N = 10
    np.random.seed(0)
    X = np.random.random((N, N))
    X[X > 0.7] = 0
    Xcsr = csc_matrix(X)
    
    for i in range(N):
        # check getrow()
        yield (assert_array_almost_equal,
               X[i:i + 1, :], Xcsr.getrow(i).toarray())
        # check getcol()
        yield (assert_array_almost_equal,
               X[:, i:i + 1], Xcsr.getcol(i).toarray())


if __name__ == "__main__":
    run_module_suite()
