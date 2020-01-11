"""Test functions for linalg.matmul_toeplitz function
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.linalg import toeplitz, matmul_toeplitz

import pytest
from pytest import raises as assert_raises
from numpy.testing import assert_allclose
np.random.seed(42)

class ToeplitzTestCase:
    def test_real(self):
        cases = []

        n = 1
        c = np.random.normal(size=n)
        r = np.random.normal(size=n)
        b = np.random.normal(size=(n, 1))
        cases.append((b, c, r, False))

        n = 2
        c = np.random.normal(size=n)
        r = np.random.normal(size=n)
        b = np.random.normal(size=(n, 1))
        cases.append((b, c, r, False))

        n = 10
        c = np.random.normal(size=n)
        r = np.random.normal(size=n)
        b = np.random.normal(size=(n, 1))
        cases.append((b, c, r, True))

        n = 1000
        c = np.random.normal(size=n)
        r = np.random.normal(size=n)
        b = np.random.normal(size=(n, 1))
        cases.append((b, c, r, False))

        n = 100
        c = np.random.normal(size=n)
        r = np.random.normal(size=n)
        b = np.random.normal(size=(n, np.random.randint(1, 10)))
        cases.append((b, c, r, False))

        n = 100
        c = np.random.normal(size=(n, 1))
        r = np.random.normal(size=(n, 1))
        b = np.random.normal(size=(n, np.random.randint(1, 10)))
        cases.append((b, c, r, True))

        n = 100
        c = np.random.normal(size=(n, 1))
        r = None
        b = np.random.normal(size=(n, np.random.randint(1, 10)))
        cases.append((b, c, r, True))

        n = 100
        c = np.random.normal(size=(n, 1))
        r = None
        b = np.random.normal(size=n)
        cases.append((b, c, r, False))

        [self.do(*i) for i in cases]

    def test_complex(self):
        n = 127
        c = np.random.normal(size=(n, 1)) + np.random.normal(size=(n, 1))*1j
        r = np.random.normal(size=(n, 1)) + np.random.normal(size=(n, 1))*1j
        b = np.random.normal(size=(n, 3)) + np.random.normal(size=(n, 3))*1j
        self.do(b, c, r, False)

    def test_exceptions(self):

        n = 100
        c = np.random.normal(size=2*n)
        r = np.random.normal(size=n)
        b = np.random.normal(size=n)
        assert_raises(ValueError, matmul_toeplitz, (c, r), b, True)

        n = 100
        c = np.random.normal(size=n)
        r = np.random.normal(size=2*n)
        b = np.random.normal(size=n)
        assert_raises(ValueError, matmul_toeplitz, (c, r), b, True)

        n = 100
        c = np.random.normal(size=n)
        r = np.random.normal(size=n)
        b = np.random.normal(size=n-1)
        assert_raises(ValueError, matmul_toeplitz, (c, r), b, True)


class TestMul(ToeplitzTestCase):
    # For toeplitz matrices, matmul_toeplitz() should be equivalent to @.
    def do(self, b, c, r=None, check_finite=False):
        if r is None:
            actual = matmul_toeplitz(c, b, check_finite)
        else:
            actual = matmul_toeplitz((c, r), b, check_finite)
        desired = toeplitz(c, r) @ b
        assert_allclose(actual, desired)
