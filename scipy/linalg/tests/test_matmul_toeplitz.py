"""Test functions for linalg.matmul_toeplitz function
"""

import numpy as np
from scipy.linalg import toeplitz, matmul_toeplitz

from pytest import raises as assert_raises
from numpy.testing import assert_allclose


class TestMatmulToeplitz:

    def setup_method(self):
        self.rng = np.random.RandomState(42)
        self.tolerance = 1e-10

    def test_real(self):
        cases = []

        n = 1
        c = self.rng.normal(size=n)
        r = self.rng.normal(size=n)
        b = self.rng.normal(size=(n, 1))
        cases.append((b, c, r, False))

        n = 2
        c = self.rng.normal(size=n)
        r = self.rng.normal(size=n)
        b = self.rng.normal(size=(n, 1))
        cases.append((b, c, r, False))

        n = 10
        c = self.rng.normal(size=n)
        r = self.rng.normal(size=n)
        b = self.rng.normal(size=(n, 1))
        cases.append((b, c, r, True))

        n = 1000
        c = self.rng.normal(size=n)
        r = self.rng.normal(size=n)
        b = self.rng.normal(size=(n, 1))
        cases.append((b, c, r, False))

        n = 100
        c = self.rng.normal(size=n)
        r = self.rng.normal(size=n)
        b = self.rng.normal(size=(n, self.rng.randint(1, 10)))
        cases.append((b, c, r, False))

        n = 100
        c = self.rng.normal(size=(n, 1))
        r = self.rng.normal(size=(n, 1))
        b = self.rng.normal(size=(n, self.rng.randint(1, 10)))
        cases.append((b, c, r, True))

        n = 100
        c = self.rng.normal(size=(n, 1))
        r = None
        b = self.rng.normal(size=(n, self.rng.randint(1, 10)))
        cases.append((b, c, r, True, -1))

        n = 100
        c = self.rng.normal(size=(n, 1))
        r = None
        b = self.rng.normal(size=n)
        cases.append((b, c, r, False))

        [self.do(*i) for i in cases]

    def test_complex(self):
        n = 127
        c = self.rng.normal(size=(n, 1)) + self.rng.normal(size=(n, 1))*1j
        r = self.rng.normal(size=(n, 1)) + self.rng.normal(size=(n, 1))*1j
        b = self.rng.normal(size=(n, 3)) + self.rng.normal(size=(n, 3))*1j
        self.do(b, c, r, False)

    def test_exceptions(self):

        n = 100
        c = self.rng.normal(size=2*n)
        r = self.rng.normal(size=n)
        b = self.rng.normal(size=n)
        assert_raises(ValueError, matmul_toeplitz, (c, r), b, True)

        n = 100
        c = self.rng.normal(size=n)
        r = self.rng.normal(size=2*n)
        b = self.rng.normal(size=n)
        assert_raises(ValueError, matmul_toeplitz, (c, r), b, True)

        n = 100
        c = self.rng.normal(size=n)
        r = self.rng.normal(size=n)
        b = self.rng.normal(size=n-1)
        assert_raises(ValueError, matmul_toeplitz, (c, r), b, True)

    # For toeplitz matrices, matmul_toeplitz() should be equivalent to @.
    def do(self, b, c, r=None, check_finite=False, workers=None):
        if r is None:
            actual = matmul_toeplitz(c, b, check_finite, workers)
        else:
            actual = matmul_toeplitz((c, r), b, check_finite)
        desired = toeplitz(c, r) @ b
        assert_allclose(actual, desired, rtol=self.tolerance)
