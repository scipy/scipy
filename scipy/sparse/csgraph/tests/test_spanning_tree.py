"""Test the minimum spanning tree function"""
import numpy as np
import pytest
from numpy.testing import assert_
import numpy.testing as npt
from scipy.sparse import csr_array, coo_array
from scipy.sparse.csgraph import minimum_spanning_tree


def test_minimum_spanning_tree():

    # Create a graph with two connected components.
    graph = [[0,1,0,0,0],
             [1,0,0,0,0],
             [0,0,0,8,5],
             [0,0,8,0,1],
             [0,0,5,1,0]]
    graph = np.asarray(graph)

    # Create the expected spanning tree.
    expected = [[0,1,0,0,0],
                [0,0,0,0,0],
                [0,0,0,0,5],
                [0,0,0,0,1],
                [0,0,0,0,0]]
    expected = np.asarray(expected)

    # Ensure minimum spanning tree code gives this expected output.
    csgraph = csr_array(graph)
    mintree = minimum_spanning_tree(csgraph)
    mintree_array = mintree.toarray()
    npt.assert_array_equal(mintree_array, expected,
                           'Incorrect spanning tree found.')

    # Ensure that the original graph was not modified.
    npt.assert_array_equal(csgraph.toarray(), graph,
        'Original graph was modified.')

    # Now let the algorithm modify the csgraph in place.
    mintree = minimum_spanning_tree(csgraph, overwrite=True)
    npt.assert_array_equal(mintree.toarray(), expected,
        'Graph was not properly modified to contain MST.')

    np.random.seed(1234)
    for N in (5, 10, 15, 20):

        # Create a random graph.
        graph = 3 + np.random.random((N, N))
        csgraph = csr_array(graph)

        # The spanning tree has at most N - 1 edges.
        mintree = minimum_spanning_tree(csgraph)
        assert_(mintree.nnz < N)

        # Set the sub diagonal to 1 to create a known spanning tree.
        idx = np.arange(N-1)
        graph[idx,idx+1] = 1
        csgraph = csr_array(graph)
        mintree = minimum_spanning_tree(csgraph)

        # We expect to see this pattern in the spanning tree and otherwise
        # have this zero.
        expected = np.zeros((N, N))
        expected[idx, idx+1] = 1

        npt.assert_array_equal(mintree.toarray(), expected,
            'Incorrect spanning tree found.')

@pytest.mark.parametrize("dtype", [np.int32, np.int64])
def test_mst_with_various_index_dtypes(dtype):
    #CSR array
    # Row indices
    indptr = np.array([0, 2, 4, 5], dtype=dtype)
    indices = np.array([1, 2, 0, 2, 1], dtype=dtype)
    data = np.array([2, 0, 2, 3, 3], dtype=float)

    graph = csr_array((data, indices, indptr), shape=(3, 3))

    # Check whether the dtype of indices is as expected
    assert graph.indices.dtype == dtype
    assert graph.indptr.dtype == dtype

    # Compute MST
    mst = minimum_spanning_tree(graph)

    # Check whether the dtype of indices is as expected
    assert mst.indices.dtype == dtype
    assert mst.indptr.dtype == dtype

    # Expected MST has 2 edges: (0->1, weight 2)
    expected = np.array([
        [0, 2, 0],
        [0, 0, 0],
        [0, 0, 0]
    ], dtype=float)

    npt.assert_array_almost_equal(mst.toarray(), expected)

    # COO array
    data = np.array([1.0, 4.0, 3.0, 2.0, 5.0], dtype=float)
    row = np.array([0, 1, 2, 3, 4], dtype=dtype)
    col = np.array([1, 2, 3, 4, 5], dtype=dtype)

    # Build COO graph with (9, 9) shape
    adj = coo_array((data, (row, col)), shape=(9, 9))

    # Check whether the dtype of indices is as expected
    assert adj.row.dtype == dtype
    assert adj.col.dtype == dtype

    # Compute MST
    mst = minimum_spanning_tree(adj)

    # Check whether the dtype of indices is as expected
    assert mst.indices.dtype == dtype
    assert mst.indptr.dtype == dtype

    # Expected MST will include all 5 edges since graph is disconnected elsewhere
    expected = np.zeros((9, 9), dtype=float)
    expected[0, 1] = 1.0
    expected[1, 2] = 4.0
    expected[2, 3] = 3.0
    expected[3, 4] = 2.0
    expected[4, 5] = 5.0

    npt.assert_array_almost_equal(mst.toarray(), expected)
