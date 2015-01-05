""" Test functions in the scipy.linalg._matrix_norms module.

"""
from __future__ import division, print_function, absolute_import

import itertools

import numpy as np
import numpy.linalg
from numpy.testing import (assert_equal, assert_raises, assert_allclose,
        assert_array_less, run_module_suite, dec)

import scipy
import scipy.linalg
from scipy.lib._version import NumpyVersion
from scipy.linalg._matrix_norms import (
        elementwise_norm, frobenius_norm, nuclear_norm, spectral_norm,
        schatten_norm, induced_norm, ky_fan_norm)


def _specialized_norm(mynorm, param):
    # functools.partial does not quite work here because of argument order.
    def _inner(A, **kwargs):
        return mynorm(A, param, **kwargs)
    return _inner


# These matrix norms raise errors on arrays with ndim < 2.
_pure_matrix_norms = [
        frobenius_norm,
        nuclear_norm,
        spectral_norm,
        _specialized_norm(induced_norm, 1),
        _specialized_norm(induced_norm, 2),
        _specialized_norm(induced_norm, np.inf),
        _specialized_norm(schatten_norm, 1),
        _specialized_norm(schatten_norm, 1.2),
        _specialized_norm(schatten_norm, 2),
        _specialized_norm(schatten_norm, 3),
        _specialized_norm(schatten_norm, np.inf),
        _specialized_norm(ky_fan_norm, 1),
        _specialized_norm(ky_fan_norm, 2)]

# Elementwise norms can be matrix norms but also generalize to other shapes.
_matrix_norms = _pure_matrix_norms + [
        _specialized_norm(elementwise_norm, 1.2),
        _specialized_norm(elementwise_norm, 1),
        _specialized_norm(elementwise_norm, 2),
        _specialized_norm(elementwise_norm, 3),
        _specialized_norm(elementwise_norm, np.inf)]

# Classify norms according to unitary invariance for testing.
_unitarily_invariant_norms = [
        frobenius_norm,
        nuclear_norm,
        spectral_norm,
        _specialized_norm(induced_norm, 2),
        _specialized_norm(schatten_norm, 1),
        _specialized_norm(schatten_norm, 1.2),
        _specialized_norm(schatten_norm, 2),
        _specialized_norm(schatten_norm, 3),
        _specialized_norm(schatten_norm, np.inf),
        _specialized_norm(ky_fan_norm, 1),
        _specialized_norm(ky_fan_norm, 2),
        _specialized_norm(elementwise_norm, 2)]

_not_unitarily_invariant_norms = [
        _specialized_norm(induced_norm, 1),
        _specialized_norm(induced_norm, np.inf),
        _specialized_norm(elementwise_norm, 1.2),
        _specialized_norm(elementwise_norm, 1),
        _specialized_norm(elementwise_norm, 3),
        _specialized_norm(elementwise_norm, np.inf)]

# List some submultiplicative norms.
_ordinary_submultiplicative_norms = _unitarily_invariant_norms + [
        _specialized_norm(induced_norm, 1),
        _specialized_norm(induced_norm, np.inf)]

# These specializations cannot be computed for whatever reason.
# Some of them may not actually correspond to norms.
# Some of them may be just too difficult to compute.
_bad_matrix_norms = [
        _specialized_norm(induced_norm, -2),
        _specialized_norm(induced_norm, -1),
        _specialized_norm(induced_norm, 0),
        _specialized_norm(induced_norm, 0.5),
        _specialized_norm(induced_norm, 1.5),
        _specialized_norm(induced_norm, 3),
        _specialized_norm(ky_fan_norm, 0),
        _specialized_norm(ky_fan_norm, -1),
        _specialized_norm(ky_fan_norm, 0.5),
        _specialized_norm(ky_fan_norm, 1.5),
        _specialized_norm(ky_fan_norm, np.inf),
        _specialized_norm(schatten_norm, -1),
        _specialized_norm(schatten_norm, 0),
        _specialized_norm(schatten_norm, 0.5),
        _specialized_norm(elementwise_norm, -1),
        _specialized_norm(elementwise_norm, 0),
        _specialized_norm(elementwise_norm, 0.5)]


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

    @dec.knownfailureif(True, 'empty norms are not yet allowed')
    def test_empty(self):
        A = self.arraytype([[]])

        # Check special norms and parameterized norms.
        self._check_frobenius_norm(A, 0)
        self._check_spectral_norm(A, 0)
        self._check_nuclear_norm(A, 0)
        self._check_max_absolute_row_sum_norm(A, 0)
        self._check_max_absolute_column_sum_norm(A, 0)
        self._check_parameterized_norms(A)

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

    def test_3x2(self):
        # Check norms of a rectangular matrix.
        A = np.array([[1, 6], [2, 0], [3, 5]], dtype=self.dt)

        # Check special norms and parameterized norms.
        self._check_frobenius_norm(A, np.sqrt(75))
        self._check_spectral_norm(A, 8.3075790106768661)
        self._check_nuclear_norm(A, 10.753827358936125)
        self._check_max_absolute_row_sum_norm(A, 8.0)
        self._check_max_absolute_column_sum_norm(A, 11.0)
        self._check_parameterized_norms(A)

    @dec.knownfailureif(NumpyVersion(np.__version__) < '1.8.0.dev')
    def test_axis(self):
        B = np.arange(1, 25, dtype=self.dt).reshape(2, 3, 4)

        for norm in _matrix_norms:

            n = norm(B, axis=(1, 2))
            desired = [norm(B[k]) for k in range(B.shape[0])]
            self._assert_allclose(n, desired)

            n = norm(B, axis=(2, 1))
            desired = [norm(B[k].T) for k in range(B.shape[0])]
            self._assert_allclose(n, desired)

            n = norm(B, axis=(0, 2))
            desired = [norm(B[:, k,:]) for k in range(B.shape[1])]
            self._assert_allclose(n, desired)

            n = norm(B, axis=(0, 1))
            desired = [norm(B[:,:, k]) for k in range(B.shape[2])]
            self._assert_allclose(n, desired)

    @dec.knownfailureif(NumpyVersion(np.__version__) < '1.10.0.dev')
    def test_keepdims(self):
        A = np.arange(1,25, dtype=self.dt).reshape(2,3,4)

        allclose_err = 'axis = {0}'
        shape_err = 'Shape mismatch found {0}, expected {1}, axis={2}'
        for norm in _matrix_norms:
            for k in itertools.permutations(range(A.ndim), 2):
                expected = norm(A, axis=k)
                found = norm(A, axis=k, keepdims=True)
                assert_allclose(np.squeeze(found), expected,
                        err_msg=allclose_err.format(k))
                expected_shape = list(A.shape)
                expected_shape[k[0]] = 1
                expected_shape[k[1]] = 1
                expected_shape = tuple(expected_shape)
                assert_equal(found.shape, expected_shape,
                    shape_err.format(found.shape, expected_shape, k))

    @dec.knownfailureif(NumpyVersion(np.__version__) < '1.8.0.dev')
    def test_bad_args(self):
        # Check that bad arguments raise the appropriate exceptions.

        A = np.array([[1, 2, 3], [4, 5, 6]], dtype=self.dt)
        B = np.arange(1, 25, dtype=self.dt).reshape(2, 3, 4)

        # Using `axis=<integer>` or passing in a 1-D array implies vector
        # norms are being computed, so pure matrix norms raise a ValueError.
        for norm in _pure_matrix_norms:
            assert_raises(ValueError, norm, A, axis=0)
            assert_raises(ValueError, norm, [3, 4], axis=None)

        # Bad matrix norm, ok axis.
        for norm in _bad_matrix_norms:
            assert_raises(ValueError, norm, A, axis=None)
            assert_raises(ValueError, norm, A, axis=(0, 1))
            assert_raises(ValueError, norm, B, axis=(1, 2))

        # Invalid axis.
        for norm in _pure_matrix_norms:
            assert_raises(ValueError, norm, B, axis=3)
            assert_raises(ValueError, norm, B, axis=(2, 3))
            assert_raises(ValueError, norm, B, axis=(0, 1, 2))


def test_unitary_invariance():
    np.random.seed(1234)
    n = 3
    A = np.array([[1, 2, 3], [6, 0, 5], [3, 2, 1]], dtype=float)
    U = scipy.linalg.orth(np.random.randn(n, n) + 1j * np.random.randn(n, n))
    assert_allclose(U.dot(U.T.conj()), np.identity(n), atol=1e-12)
    B = U.dot(A).dot(U.T.conj())
    for norm in _unitarily_invariant_norms:
        assert_allclose(norm(A), norm(B))
    for norm in _not_unitarily_invariant_norms:
        assert_raises(AssertionError, assert_allclose, norm(A), norm(B))


def test_ordinary_submultiplicativity():
    # Most of the matrix norms are submultiplicative.
    np.random.seed(1234)
    n = 4
    for i in range(10):
        A = np.random.randn(n, n)
        B = np.random.randn(n, n)
        for norm in _ordinary_submultiplicative_norms:
            assert_array_less(norm(A.dot(B)), norm(A)*norm(B))


def test_hadamard_submultiplicativity():
    # Among unitarily invariant norms, ordinary and hadamard
    # submultiplicativity are equivalent.
    # Theorem 5.5.7
    # Topics in Matrix Analysis
    # Roger Horn
    np.random.seed(1234)
    n = 4
    for i in range(10):
        A = np.random.randn(n, n)
        B = np.random.randn(n, n)
        for norm in _unitarily_invariant_norms:
            assert_array_less(norm(A * B), norm(A)*norm(B))


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
