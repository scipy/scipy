"""Test the minimum spanning tree function"""
import numpy as np
from numpy.testing import assert_
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree

def construct_sparse_graph(N):
    np.random.seed(1234)
    csgraph = np.random.random((N, N))
    ind = np.random.randint(N, size=(2, N * N / 2))
    csgraph[ind[0], ind[1]] = 0
    return csr_matrix(csgraph)

def test_minimum_spanning_tree():
    for N in (5, 10, 15, 20):
        csgraph = construct_sparse_graph(N)
        mintree = minimum_spanning_tree(csgraph)
        assert_(mintree.nnz < N)
