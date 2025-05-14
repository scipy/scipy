"""Test functions for the sparse.linalg._krylov_funm module."""
from functools import partial

import numpy as np
import pytest
from numpy.testing import (assert_allclose)
import scipy.sparse
import scipy.linalg
from scipy.sparse.linalg import aslinearoperator
from scipy.linalg import (expm, cosm, coshm, sinm, sinhm)

from scipy.sparse.linalg._krylov_funmv import krylov_funmv

REAL_DTYPES = (np.float32, np.float64)
COMPLEX_DTYPES = (np.complex64, np.complex128)
DTYPES = REAL_DTYPES + COMPLEX_DTYPES

# This the phi_1 function from exponential integrators: phi_1 (z) = (e^z - 1) / z
def custom(X):
    return scipy.linalg.solve(X, expm(X) - np.eye(X.shape[0]))

FUNCS = (expm, cosm, coshm, sinm, sinhm, custom)

class TestKrylovFunmv:

    @pytest.mark.parametrize("f", FUNCS)
    def test_krylov_funmv_nonhermitian_dense(self, f):
        rng = np.random.default_rng(1738151906092735)
        n = 60
        nsamples = 10

        for i in range(nsamples):
            A = rng.standard_normal((n, n))
            b = rng.standard_normal(n)
            t = rng.random()

            fA = f(t * A)
            expected = fA @ b
            observed = krylov_funmv(f, t, A, b)
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)
            observed = krylov_funmv(f, t, aslinearoperator(A), b)
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)

    @pytest.mark.parametrize("f", FUNCS)
    def test_krylov_funmv_nonhermitian_sparse(self, f):

        rng = np.random.default_rng(1738151906092735)
        n = 100
        nsamples = 10

        for i in range(nsamples):
            D = scipy.sparse.diags(rng.standard_normal(n))
            A = scipy.sparse.random_array((n, n), density = 0.01, rng = rng) + D
            denseA = A.todense()
            b = rng.standard_normal(n)
            t = rng.random()

            fA = f(t * denseA)
            expected = fA @ b
            observed = krylov_funmv(f, t, A, b)
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)
            observed = krylov_funmv(f, t, aslinearoperator(A), b)
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)

    @pytest.mark.parametrize("f", FUNCS)
    def test_krylov_funmv_hermitian_dense(self, f):

        rng = np.random.default_rng(1738151906092735)
        n = 60
        nsamples = 10

        for i in range(nsamples):
            R = np.triu(rng.standard_normal((n, n)))
            A = R.T + R
            b = rng.standard_normal(n)
            t = rng.random()

            fA = f(t * A)
            expected = fA @ b
            observed = krylov_funmv(f, t, A, b, ortho_method = 'lanczos')
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)
            observed = krylov_funmv(f, t, aslinearoperator(A), b,
                                    ortho_method = 'lanczos')
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)

    @pytest.mark.parametrize("f", FUNCS)
    def test_krylov_funmv_hermitian_sparse(self, f):

        rng = np.random.default_rng(1738151906092735)
        n = 100
        nsamples = 10

        for i in range(nsamples):
            D = scipy.sparse.diags(rng.standard_normal(n))
            A = scipy.sparse.random_array((n, n), density = 0.01, rng = rng)
            R = scipy.sparse.triu(A)
            A = R + R.T + D
            denseA = A.todense()
            b = rng.standard_normal(n)
            t = rng.random()

            fA = f(t * denseA)
            expected = fA @ b
            observed = krylov_funmv(f, t, A, b, ortho_method = 'lanczos')
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)

            observed = krylov_funmv(f, t, aslinearoperator(A), b,
                                    ortho_method='lanczos')
            assert_allclose(observed, expected, rtol = 1E-6, atol = 1E-8)


@pytest.mark.parametrize("dtype_a", DTYPES)
@pytest.mark.parametrize("dtype_b", DTYPES)
def test_krylov_funmv_types(dtype_a, dtype_b):
    assert_allclose_ = (partial(assert_allclose, rtol = 1.8e-3, atol = 1e-5)
                        if {dtype_a, dtype_b} else assert_allclose)

    rng = np.random.default_rng(1738151906092735)
    n = 50
    t = rng.random()

    if dtype_a in REAL_DTYPES:
        A = rng.random([n, n]).astype(dtype_a)
    else:
        A = (rng.random([n, n]) + 1j * rng.random([n, n])).astype(dtype_a)

    if dtype_b in REAL_DTYPES:
        b = (2 * rng.random(n)).astype(dtype_b)
    else:
        b = (rng.random(n) + 1j * rng.random(n)).astype(dtype_b)

        expA = expm(t * A)
        expected = expA @ b
        observed = krylov_funmv(expm, t, A, b)
        assert_allclose_(observed, expected)
        observed = krylov_funmv(expm, t, aslinearoperator(A), b)
        assert_allclose_(observed, expected)
