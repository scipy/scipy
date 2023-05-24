"""Test the minimum spanning tree function"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree


def test_minimum_spanning_tree():
    # Create a graph with two connected components.
    graph = np.array([[0,1,0,0,0],
                      [1,0,0,0,0],
                      [0,0,0,8,5],
                      [0,0,8,0,1],
                      [0,0,5,1,0]])

    # Create the expected spanning tree(s). With two connected components,
    # we have four possible results (as each tree can be transposed).
    expected_results = (
        np.array([[0,1,0,0,0],
                  [0,0,0,0,0],
                  [0,0,0,0,5],
                  [0,0,0,0,1],
                  [0,0,0,0,0]], dtype=float),
        np.array([[0,0,0,0,0],
                  [1,0,0,0,0],
                  [0,0,0,0,5],
                  [0,0,0,0,1],
                  [0,0,0,0,0]], dtype=float),
        np.array([[0,1,0,0,0],
                  [0,0,0,0,0],
                  [0,0,0,0,0],
                  [0,0,0,0,0],
                  [0,0,5,1,0]], dtype=float),
        np.array([[0,0,0,0,0],
                  [1,0,0,0,0],
                  [0,0,0,0,0],
                  [0,0,0,0,0],
                  [0,0,5,1,0]], dtype=float),
    )

    # Ensure minimum spanning tree code gives one of the expected outputs.
    csgraph = csr_matrix(graph)
    mintree = minimum_spanning_tree(csgraph).toarray()
    assert any(np.array_equal(mintree, tree) for tree in expected_results)

    # Ensure that the original graph was not modified.
    npt.assert_array_equal(csgraph.toarray(), graph,
        'Original graph was modified.')

    # Now let the algorithm modify the csgraph in place.
    mintree2 = minimum_spanning_tree(csgraph, overwrite=True).toarray()
    assert any(np.array_equal(mintree2, tree) for tree in expected_results)


@pytest.mark.parametrize('N', [5, 10, 15, 20])
def test_minimum_spanning_tree_random_graphs(N: int):
    np.random.seed(1234)

    # Create a random graph.
    graph = 3 + np.random.random((N, N))
    csgraph = csr_matrix(graph)

    # The spanning tree has at most N - 1 edges.
    mintree = minimum_spanning_tree(csgraph)
    assert mintree.nnz < N

    # Set the sub diagonal to 1 to create a known spanning tree.
    idx = np.arange(N - 1)
    graph[idx, idx+1] = 1
    csgraph = csr_matrix(graph)
    mintree = minimum_spanning_tree(csgraph).toarray()

    # We expect to see this pattern in the spanning tree and otherwise
    # have this zero.
    expected = np.zeros((N, N))
    expected[idx, idx+1] = 1

    npt.assert_array_equal(mintree, expected, 'Incorrect spanning tree found.')
