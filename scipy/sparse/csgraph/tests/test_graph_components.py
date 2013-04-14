from __future__ import division, print_function, absolute_import

import warnings

import numpy as np
from numpy.testing import assert_, assert_equal
from numpy.testing.utils import WarningManager
from scipy.sparse import csr_matrix, csgraph


def test_cs_graph_components():
    D = np.eye(4, dtype=np.bool)

    warn_ctx = WarningManager()
    warn_ctx.__enter__()
    try:
        warnings.filterwarnings("ignore",
                    message="`cs_graph_components` is deprecated")

        n_comp, flag = csgraph.cs_graph_components(csr_matrix(D))
        assert_(n_comp == 4)
        assert_equal(flag, [0, 1, 2, 3])

        D[0, 1] = D[1, 0] = 1

        n_comp, flag = csgraph.cs_graph_components(csr_matrix(D))
        assert_(n_comp == 3)
        assert_equal(flag, [0, 0, 1, 2])

        # A pathological case...
        D[2, 2] = 0
        n_comp, flag = csgraph.cs_graph_components(csr_matrix(D))
        assert_(n_comp == 2)
        assert_equal(flag, [0, 0, -2, 1])
    finally:
        warn_ctx.__exit__()
