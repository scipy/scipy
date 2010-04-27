import numpy as np
from numpy.testing import run_module_suite, assert_almost_equal

import scipy.sparse as sp
import scipy.sparse.linalg as spla

def test_gmres_basic():
    A = np.vander(np.arange(10) + 1)[:, ::-1]
    b = np.zeros(10)
    b[0] = 1
    x = np.linalg.solve(A, b)

    x_gm, err = spla.gmres(A, b, restart=5, maxiter=1)

    assert_almost_equal(x_gm[0], 0.359, decimal=2)

if __name__ == "__main__":
    run_module_suite()
