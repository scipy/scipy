""" Test functions in the scipy.linalg._matrix_norms module.

"""
from __future__ import division, print_function, absolute_import

import numpy as np
import numpy.linalg

from numpy.testing import (
        TestCase, assert_, assert_equal, assert_raises, assert_array_equal,
        assert_almost_equal, assert_allclose, run_module_suite, dec)

import scipy
import scipy.linalg
from scipy.linalg._matrix_norms import (
        elementwise_norm, frobenius_norm, nuclear_norm, spectral_norm,
        schatten_norm, induced_norm, ky_fan_norm)

old_assert_almost_equal = assert_almost_equal


def assert_almost_equal(a, b, **kw):
    if np.asarray(a).dtype.type in (np.single, np.csingle):
        decimal = 6
    else:
        decimal = 12
    old_assert_almost_equal(a, b, decimal=decimal, **kw)


class _TestMatrixNorms(object):
    # Test matrix norms without using the 'axis' or 'keepdims' keyword args.
    # These tests are copied from or inspired by numpy matrix norm tests.

    def test_matrix_2x2(self):
        for arraytype in np.matrix, np.array:
            A = arraytype([[1, 3], [5, 7]], dtype=self.dt)

            # Frobenius norm.
            desired = 84**0.5
            assert_almost_equal(frobenius_norm(A), desired)
            assert_almost_equal(elementwise_norm(A, 2), desired)
            assert_almost_equal(schatten_norm(A, 2), desired)
            assert_almost_equal(scipy.linalg.norm(A), desired)
            assert_almost_equal(scipy.linalg.norm(A, 'fro'), desired)
            assert_almost_equal(np.linalg.norm(A), desired)
            assert_almost_equal(np.linalg.norm(A, 'fro'), desired)

            # Spectral norm.
            desired = 9.1231056256176615
            assert_almost_equal(spectral_norm(A), desired)
            assert_almost_equal(induced_norm(A, 2), desired)
            assert_almost_equal(schatten_norm(A, np.inf), desired)
            assert_almost_equal(ky_fan_norm(A, 1), desired)
            assert_almost_equal(scipy.linalg.norm(A, 2), desired)
            assert_almost_equal(np.linalg.norm(A, 2), desired)

            # Nuclear norm.
            desired = 10.0
            assert_almost_equal(nuclear_norm(A), desired)
            assert_almost_equal(schatten_norm(A, 1), desired)
            assert_almost_equal(ky_fan_norm(A, min(A.shape)), desired)

            # Maximum absolute row sum norm.
            desired = 12.0
            assert_almost_equal(induced_norm(A, np.inf), desired)
            assert_almost_equal(scipy.linalg.norm(A, np.inf), desired)
            assert_almost_equal(np.linalg.norm(A, np.inf), desired)

            # Maximum absolute column sum norm.
            desired = 10.0
            assert_almost_equal(induced_norm(A, 1), desired)
            assert_almost_equal(scipy.linalg.norm(A, 1), desired)
            assert_almost_equal(np.linalg.norm(A, 1), desired)


    """
    def test_matrix_3x3(self):
        # This test has been added because the 2x2 example
        # happened to have equal nuclear norm and induced 1-norm.
        # The 1/10 scaling factor accommodates the absolute tolerance
        # used in assert_almost_equal.
        A = (1/10) * np.array([[1, 2, 3], [6, 0, 5], [3, 2, 1]], dtype=self.dt)
        assert_almost_equal(norm(A), (1/10) * 89**0.5)
        assert_almost_equal(norm(A, 'fro'), (1/10) * 89**0.5)
        assert_almost_equal(norm(A, 'nuc'), 1.3366836911774836)
        assert_almost_equal(norm(A, inf), 1.1)
        assert_almost_equal(norm(A, -inf), 0.6)
        assert_almost_equal(norm(A, 1), 1.0)
        assert_almost_equal(norm(A, -1), 0.4)
        assert_almost_equal(norm(A, 2), 0.88722940323461277)
        assert_almost_equal(norm(A, -2), 0.19456584790481812)
    """


class TestNormDouble(_TestMatrixNorms):
    dt = np.double
    dec = 12


class TestNormSingle(_TestMatrixNorms):
    dt = np.float32
    dec = 6


class TestNormInt64(_TestMatrixNorms):
    dt = np.int64
    dec = 12


if __name__ == "__main__":
    run_module_suite()
