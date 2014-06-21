from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_equal
from scipy.sparse.csgraph import symrcm
from scipy.sparse import csr_matrix

def test_graph_symrcm():
    A = np.array([[1, 0, 0, 0, 1, 0, 0, 0],
                [0, 1, 1, 0, 0, 1, 0, 1],
                [0, 1, 1, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 1, 0],
                [1, 0, 1, 0, 1, 0, 0, 0],
                [0, 1, 0, 0, 0, 1, 0, 1],
                [0, 0, 0, 1, 0, 0, 1, 0],
                [0, 1, 0, 0, 0, 1, 0, 1]], dtype=np.int32)
    
    graph = csr_matrix(A)
    perm = symrcm(graph)
    correct_perm = np.array([6, 3, 7, 5, 1, 2, 4, 0])
    assert_equal((perm - correct_perm).all(), 0)