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


class _TestMatrixNorms(object):
    # Test matrix norms without using the 'axis' or 'keepdims' keyword args.
    # These tests are copied from or inspired by numpy matrix norm tests.

    def _assert_allclose(self, a, b):
        # Customize the relative tolerance according to dtype.
        assert_allclose(a, b, rtol=self.rtol, atol=0)

    def _check_frobenius_norm(self, A, desired):
        self._assert_allclose(frobenius_norm(A), desired)
        self._assert_allclose(elementwise_norm(A, 2), desired)
        self._assert_allclose(schatten_norm(A, 2), desired)
        self._assert_allclose(scipy.linalg.norm(A), desired)
        self._assert_allclose(scipy.linalg.norm(A, 'fro'), desired)
        self._assert_allclose(np.linalg.norm(A), desired)
        self._assert_allclose(np.linalg.norm(A, 'fro'), desired)

    def _check_spectral_norm(self, A, desired):
        self._assert_allclose(spectral_norm(A), desired)
        self._assert_allclose(induced_norm(A, 2), desired)
        self._assert_allclose(schatten_norm(A, np.inf), desired)
        self._assert_allclose(ky_fan_norm(A, 1), desired)
        self._assert_allclose(scipy.linalg.norm(A, 2), desired)
        self._assert_allclose(np.linalg.norm(A, 2), desired)

    def _check_nuclear_norm(self, A, desired):
        self._assert_allclose(nuclear_norm(A), desired)
        self._assert_allclose(schatten_norm(A, 1), desired)
        self._assert_allclose(ky_fan_norm(A, min(A.shape)), desired)

    def _check_max_absolute_row_sum_norm(self, A, desired):
        self._assert_allclose(induced_norm(A, np.inf), desired)
        self._assert_allclose(scipy.linalg.norm(A, np.inf), desired)
        self._assert_allclose(np.linalg.norm(A, np.inf), desired)

    def _check_max_absolute_column_sum_norm(self, A, desired):
        self._assert_allclose(induced_norm(A, 1), desired)
        self._assert_allclose(scipy.linalg.norm(A, 1), desired)
        self._assert_allclose(np.linalg.norm(A, 1), desired)

    def _check_uninteresting_schatten_norms(self, A):
        for p in 1.2, 3:
            s = scipy.linalg.svdvals(A)
            desired = np.power(s, p).sum()**(1/p)
            self._assert_allclose(schatten_norm(A, p), desired)

    def _check_uninteresting_elementwise_norms(self, A):
        for p in 1.2, 3:
            desired = np.power(A, p).sum()**(1/p)
            v = np.asarray(A).ravel()
            self._assert_allclose(elementwise_norm(A, p), desired)
            self._assert_allclose(elementwise_norm(v, p), desired)
            self._assert_allclose(scipy.linalg.norm(v, p), desired)
            self._assert_allclose(np.linalg.norm(v, p), desired)

    def _check_uninteresting_ky_fan_norms(self, A):
        for k in range(2, min(A.shape)):
            s = scipy.linalg.svdvals(A)
            desired = s[:k].sum()
            self._assert_allclose(ky_fan_norm(A, k), desired)

    def _check_parameterized_norms(self, A):
        self._check_uninteresting_schatten_norms(A)
        self._check_uninteresting_elementwise_norms(A)
        self._check_uninteresting_ky_fan_norms(A)

    def test_2x2(self):
        A = self.arraytype([[1, 3], [5, 7]], dtype=self.dt)

        # Check special norms and parameterized norms.
        self._check_frobenius_norm(A, np.sqrt(84))
        self._check_spectral_norm(A, 9.1231056256176615)
        self._check_nuclear_norm(A, 10.0)
        self._check_max_absolute_row_sum_norm(A, 12.0)
        self._check_max_absolute_column_sum_norm(A, 10.0)
        self._check_parameterized_norms(A)

    def test_3x3(self):
        # This test has been added because the 2x2 example
        # happened to have equal nuclear norm and induced 1-norm.
        # Also, for matrices smaller than 3x3 every ky-fan norm
        # is either a spectral norm or a nuclear norm.
        A = np.array([[1, 2, 3], [6, 0, 5], [3, 2, 1]], dtype=self.dt)

        # Check special norms and parameterized norms.
        self._check_frobenius_norm(A, np.sqrt(89))
        self._check_spectral_norm(A, 8.8722940323461277)
        self._check_nuclear_norm(A, 13.366836911774836)
        self._check_max_absolute_row_sum_norm(A, 11.0)
        self._check_max_absolute_column_sum_norm(A, 10.0)
        self._check_parameterized_norms(A)

    def test_2x3(self):
        # Check norms of a rectangular matrix.
        A = np.array([[1, 2, 3], [6, 0, 5]], dtype=self.dt)

        # Check special norms and parameterized norms.
        self._check_frobenius_norm(A, np.sqrt(75))
        self._check_spectral_norm(A, 8.3075790106768661)
        self._check_nuclear_norm(A, 10.753827358936125)
        self._check_max_absolute_row_sum_norm(A, 11.0)
        self._check_max_absolute_column_sum_norm(A, 8.0)
        self._check_parameterized_norms(A)


class TestNormDoubleArray(_TestMatrixNorms):
    arraytype = np.array
    dt = np.double
    rtol = 1e-12


class TestNormSingleArray(_TestMatrixNorms):
    arraytype = np.array
    dt = np.float32
    rtol = 1e-6


class TestNormInt64Array(_TestMatrixNorms):
    arraytype = np.array
    dt = np.int64
    rtol = 1e-12


class TestNormDoubleMatrix(_TestMatrixNorms):
    arraytype = np.matrix
    dt = np.double
    rtol = 1e-12


class TestNormSingleMatrix(_TestMatrixNorms):
    arraytype = np.matrix
    dt = np.float32
    rtol = 1e-6


class TestNormInt64Matrix(_TestMatrixNorms):
    arraytype = np.matrix
    dt = np.int64
    rtol = 1e-12


if __name__ == "__main__":
    run_module_suite()
