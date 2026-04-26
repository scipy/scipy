import numpy as np
from numpy.testing import assert_allclose

from scipy.sparse.csgraph import euclidean_minimum_spanning_tree

np.random.seed(0)


def test_euclidean_min_spanning_tree():
    n_samples = 100
    n_dim = 2

    # uniform distribution
    X = np.random.rand(n_samples, n_dim)
    mst_pairwise = euclidean_minimum_spanning_tree(X, method="pairwise")
    mst_delaunay = euclidean_minimum_spanning_tree(X, method="delaunay")
    assert_allclose(mst_pairwise.toarray(), mst_delaunay.toarray())

    # normal distribution
    X = np.random.randn(n_samples, n_dim)
    mst_pairwise = euclidean_minimum_spanning_tree(X, method="pairwise")
    mst_delaunay = euclidean_minimum_spanning_tree(X, method="delaunay")
    assert_allclose(mst_pairwise.toarray(), mst_delaunay.toarray())
