import numpy as np
from numpy.testing import assert_array_almost_equal
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import csgraph_from_dense, csgraph_to_dense


def test_csgraph_from_dense():
    G = np.random.random((10, 10))
    some_nulls = (G < 0.4)
    all_nulls = (G < 0.8)

    # test with single null value
    G[all_nulls] = 0
    G_csr = csgraph_from_dense(G, null_values=0)
    assert_array_almost_equal(G, G_csr.toarray())

    # test with multiple null values
    G[some_nulls] = 10
    G_csr = csgraph_from_dense(G, null_values=[0, 10])
    G[all_nulls] = 0
    assert_array_almost_equal(G, G_csr.toarray())

    # test with infinities
    G[all_nulls] = 10
    G[some_nulls] = np.inf
    G_csr = csgraph_from_dense(G, null_values=[0, 10])
    G[all_nulls] = 0
    assert_array_almost_equal(G, G_csr.toarray())


def test_csgraph_to_dense():
    G = np.random.random((10, 10))
    nulls = (G < 0.8)
    G[nulls] = np.inf
    
    G_csr = csgraph_from_dense(G)

    for null_value in [0, 10, -np.inf, np.inf]:
        G[nulls] = null_value
        assert_array_almost_equal(G, csgraph_to_dense(G_csr, null_value))


def test_multiple_edges():
    # create a random sqare matrix with an even number of elements
    X = np.random.random((10, 10))
    Xcsr = csr_matrix(X)

    # now double-up every other column
    Xcsr.indices[::2] = Xcsr.indices[1::2]

    # normal sparse toarray() will sum the duplicated edges
    Xdense = Xcsr.toarray()
    assert_array_almost_equal(Xdense[:, 1::2],
                              X[:, ::2] + X[:, 1::2])

    # csgraph_to_dense chooses the minimum of each duplicated edge
    Xdense = csgraph_to_dense(Xcsr)
    assert_array_almost_equal(Xdense[:, 1::2],
                              np.minimum(X[:, ::2], X[:, 1::2]))
