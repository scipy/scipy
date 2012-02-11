import numpy as np
from numpy.testing import assert_array_almost_equal
from scipy.sparse.csgraph import breadth_first_tree, depth_first_tree


def test_graph_breadth_first():
    csgraph = np.array([[ 0, 1, 2, 0, 0],
                        [ 1, 0, 0, 0, 3],
                        [ 2, 0, 0, 7, 0],
                        [ 0, 0, 7, 0, 1],
                        [ 0, 3, 0, 1, 0]])
    
    bfirst = np.array([[ 0, 1, 2, 0, 0],
                       [ 0, 0, 0, 0, 3],
                       [ 0, 0, 0, 7, 0],
                       [ 0, 0, 0, 0, 0],
                       [ 0, 0, 0, 0, 0]])

    for directed in [True, False]:
        bfirst_test = breadth_first_tree(csgraph, 0, directed)
        assert_array_almost_equal(bfirst_test.toarray(),
                                  bfirst)
        

def test_graph_depth_first():
    csgraph = np.array([[ 0, 1, 2, 0, 0],
                        [ 1, 0, 0, 0, 3],
                        [ 2, 0, 0, 7, 0],
                        [ 0, 0, 7, 0, 1],
                        [ 0, 3, 0, 1, 0]])

    dfirst = np.array([[ 0, 1, 0, 0, 0],
                       [ 0, 0, 0, 0, 3],
                       [ 0, 0, 0, 0, 0],
                       [ 0, 0, 7, 0, 0],
                       [ 0, 0, 0, 1, 0]])

    for directed in [True, False]:
        dfirst_test = depth_first_tree(csgraph, 0, directed)
        assert_array_almost_equal(dfirst_test.toarray(),
                                  dfirst)
