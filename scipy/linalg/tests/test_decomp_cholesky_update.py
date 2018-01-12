from __future__ import division, print_function, absolute_import
import numpy as np
from numpy.testing import assert_, assert_allclose, assert_equal
from scipy import linalg
from scipy.linalg._decomp_cholesky_update import *
from scipy.linalg.tests.test_decomp_update import assert_upper_tri

def gen_input(n=1000):
    A = np.random.rand(n, n)
    A = A + A.T + n * np.eye(n)
    R = linalg.cholesky(A)
    z = np.random.rand(n)
    # Uncomment the following to make it occasionally trip up
    #z *= 0.95 * np.sqrt(n)
    return R, z

class CholeskyUpdateTest(object):
    def test_update(self):
        # Test rank-1 update
        n = 1000
        R, z = gen_input(n)
        U = cholesky_update(R, z)
        tol = n * np.spacing(100.)
        assert_allclose(U.T.dot(U) - (R.T.dot(R) + np.outer(z, z)),
                        np.zeros([n, n]), atol=tol)

    def test_downdate(self):
        # Test rank-1 downdate
        n = 1000
        R, z = gen_input(n)
        D = cholesky_update(R, z, mod_type="-")
        tol = n * np.spacing(100.)
        assert_allclose(D.T.dot(D) - (R.T.dot(R) - np.outer(z, z)),
                        np.zeros([n, n]), atol=tol)
