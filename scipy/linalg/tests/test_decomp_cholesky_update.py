from __future__ import division, print_function, absolute_import
import numpy as np
from numpy.testing import assert_allclose, assert_raises, assert_array_equal
from scipy import linalg
from scipy.linalg._decomp_cholesky_update import *
from numpy.random import RandomState

def gen_input(n=1000, seed=72479):
    rand_state = RandomState(seed)
    A = rand_state.rand(n, n)
    A = A + A.T + n * np.eye(n)
    R = linalg.cholesky(A)
    z = rand_state.rand(n)
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
        D = cholesky_update(R, z, downdate=True)
        tol = n * np.spacing(100.)
        assert_allclose(D.T.dot(D) - (R.T.dot(R) - np.outer(z, z)),
                        np.zeros([n, n]), atol=tol)

    def test_list_input(self):
        # Test list as input
        n = 1000
        R, z = gen_input(n)
        D = cholesky_update(list(R), list(z))
        tol = n * np.spacing(100.)
        assert_allclose(D.T.dot(D) - (R.T.dot(R) - np.outer(z, z)),
                        np.zeros([n, n]), atol=tol)

    def test_finite(self):
        # Test with nan and inf in the input
        n = 10
        R, z = gen_input(n)

        # Test exception if nan present in array
        R[0, 0] = np.nan
        assert_raises(cholesky_update, R, z)

        # Test exception if inf present in array
        R[0, 0] = np.inf
        assert_raises(cholesky_update, R, z)

        # Test exception if nan present in update vector
        R, z = gen_input(n)

        z[0] = np.nan
        assert_raises(cholesky_update, R, z)

        # Test exception if inf present in update vector
        z[0] = np.inf
        assert_raises(cholesky_update, R, z)

    def test_dimensions(self):
        # Tests related to dimension of both array and vector
        n = 10
        R, z = gen_input(n)

        # Array: Too many dims
        assert_raises(cholesky_update, np.expand_dims(R, 0), z)

        # Array: Too few dims
        assert_raises(cholesky_update, z, z)

        # Vector: Too many dims
        assert_raises(cholesky_update, R, R)

        # Rectangular array
        assert_raises(cholesky_update, R[:, :-1], z)

        # Incompatible dims. between R and z
        assert_raises(cholesky_update, R[:-1, :-1], z)

    def test_param_overwrite(self):
        # Test param overwrite
        n = 100

        # Test with overwrite_R=False and overwrite_z=False
        R, z = gen_input(n)
        R_copy = R.copy()
        z_copy = z.copy()
        cholesky_update(R, z)
        assert_array_equal(R, R_copy)
        assert_array_equal(z, z_copy)

        # Test with overwrite_R=True and overwrite_z=True
        R, z = gen_input(n)
        R_copy = R.copy()
        z_copy = z.copy()
        cholesky_update(R, z, overwrite_R=True, overwrite_z=True)
        assert_raises(AssertionError, assert_array_equal, R, R_copy)
        assert_raises(AssertionError, assert_array_equal, z, z_copy)

    def test_param_lower(self):
        # Test the 'lower' parameter
        n = 100
        R, z = gen_input(n)
        L = R.T
        assert_allclose(cholesky_update(R, z),
                        cholesky_update(L, z, lower=True).T)

    def test_param_eps(self):
        # Test the 'eps' parameter
        n = 100
        R, z = gen_input(n)
        eps_val = n * np.spacing(1.)
        assert_allclose(cholesky_update(R, z),
                        cholesky_update(R, z, eps=eps_val))
