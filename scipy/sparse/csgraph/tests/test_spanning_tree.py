"""Test the minimum spanning tree function"""
import numpy as np
from numpy.testing import assert_
import numpy.testing as npt
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree

def test_minimum_spanning_tree():

    np.random.seed(1234)
    
    for N in (5, 10, 15, 20):
        
        # Create a graph as above where every edge is >= 3 except at
        # most N that are 2.
        graph = 3 + np.random.random((N, N))
        ind = np.random.randint(N, size=(2, N * N / 2))
        graph[ind[0], ind[1]] = 2
        
        csgraph = csr_matrix(graph)

        # We thus expect at most N+1 in the spanning tree.        
        mintree = minimum_spanning_tree(csgraph)
        assert_(mintree.nnz < N)
        
        # Now set the sub diagonal to 1, to create a known spanning tree.
        idx = np.arange(N-1)
        graph[idx,idx+1] = 1
        csgraph = csr_matrix(graph)
        mintree = minimum_spanning_tree(csgraph)
        
        # We expect to see this pattern in the spanning tree and otherwise
        # have this zero.
        expected = np.zeros((N,N))
        expected[idx,idx+1] = 1
        
        npt.assert_array_equal(mintree.todense(), expected, 'Incorrect spanning tree found.')