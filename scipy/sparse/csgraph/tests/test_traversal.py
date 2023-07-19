import numpy as np
from numpy.testing import assert_array_almost_equal
from scipy.sparse import csr_array
from scipy.sparse.csgraph import (breadth_first_tree, depth_first_tree,
    csgraph_to_dense, csgraph_from_dense)


def test_graph_breadth_first():
    csgraph = np.array([[0, 1, 2, 0, 0],
                        [1, 0, 0, 0, 3],
                        [2, 0, 0, 7, 0],
                        [0, 0, 7, 0, 1],
                        [0, 3, 0, 1, 0]])
    csgraph = csgraph_from_dense(csgraph, null_value=0)

    bfirst = np.array([[0, 1, 2, 0, 0],
                       [0, 0, 0, 0, 3],
                       [0, 0, 0, 7, 0],
                       [0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0]])

    for directed in [True, False]:
        bfirst_test = breadth_first_tree(csgraph, 0, directed)
        assert_array_almost_equal(csgraph_to_dense(bfirst_test),
                                  bfirst)


def test_graph_depth_first():
    csgraph = np.array([[0, 1, 2, 0, 0],
                        [1, 0, 0, 0, 3],
                        [2, 0, 0, 7, 0],
                        [0, 0, 7, 0, 1],
                        [0, 3, 0, 1, 0]])
    csgraph = csgraph_from_dense(csgraph, null_value=0)

    dfirst = np.array([[0, 1, 0, 0, 0],
                       [0, 0, 0, 0, 3],
                       [0, 0, 0, 0, 0],
                       [0, 0, 7, 0, 0],
                       [0, 0, 0, 1, 0]])

    for directed in [True, False]:
        dfirst_test = depth_first_tree(csgraph, 0, directed)
        assert_array_almost_equal(csgraph_to_dense(dfirst_test),
                                  dfirst)


def test_graph_breadth_first_trivial_graph():
    csgraph = np.array([[0]])
    csgraph = csgraph_from_dense(csgraph, null_value=0)

    bfirst = np.array([[0]])

    for directed in [True, False]:
        bfirst_test = breadth_first_tree(csgraph, 0, directed)
        assert_array_almost_equal(csgraph_to_dense(bfirst_test),
                                  bfirst)


def test_graph_depth_first_trivial_graph():
    csgraph = np.array([[0]])
    csgraph = csgraph_from_dense(csgraph, null_value=0)

    bfirst = np.array([[0]])

    for directed in [True, False]:
        bfirst_test = depth_first_tree(csgraph, 0, directed)
        assert_array_almost_equal(csgraph_to_dense(bfirst_test),
                                  bfirst)


def test_breadth_first_int64_indices():
    # See https://github.com/scipy/scipy/issues/18716
    g = csr_array((np.array([1]), np.array([[0], [1]])), shape=(2, 2))
    assert g.indices.dtype == np.int64
    for directed in [True, False]:
        bf = breadth_first_tree(g, 0, directed=directed)
        assert_array_almost_equal(csgraph_to_dense(bf), [[0, 1], [0, 0]])


def test_depth_first_int64_indices():
    # See https://github.com/scipy/scipy/issues/18716
    g = csr_array((np.array([1]), np.array([[0], [1]])), shape=(2, 2))
    assert g.indices.dtype == np.int64
    for directed in [True, False]:
        df = depth_first_tree(g, 0, directed=directed)
        assert_array_almost_equal(csgraph_to_dense(df), [[0, 1], [0, 0]])
