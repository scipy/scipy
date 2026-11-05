"""Test functions for linalg._solve_toeplitz module
"""
import numpy as np
from scipy.linalg._solve_toeplitz import levinson, levinson_inverse_columns
from scipy.linalg import solve, toeplitz, solve_toeplitz
from scipy.linalg._basic import _solve_toeplitz_multiple_rhs
from numpy.testing import assert_equal, assert_allclose
from scipy._lib._array_api import xp_assert_close

import pytest
from pytest import raises as assert_raises


def test_solve_equivalence():
    # For toeplitz matrices, solve_toeplitz() should be equivalent to solve().
    random = np.random.RandomState(1234)
    for n in (1, 2, 3, 10):
        c = random.randn(n)
        if random.rand() < 0.5:
            c = c + 1j * random.randn(n)
        r = random.randn(n)
        if random.rand() < 0.5:
            r = r + 1j * random.randn(n)
        y = random.randn(n)
        if random.rand() < 0.5:
            y = y + 1j * random.randn(n)

        # Check equivalence when both the column and row are provided.
        actual = solve_toeplitz((c,r), y)
        desired = solve(toeplitz(c, r=r), y)
        assert_allclose(actual, desired)

        # Check equivalence when the column is provided but not the row.
        actual = solve_toeplitz(c, b=y)
        desired = solve(toeplitz(c), y)
        assert_allclose(actual, desired)


def test_multiple_rhs():
    random = np.random.RandomState(1234)
    c = random.randn(4)
    r = random.randn(4)
    for offset in [0, 1j]:
        for yshape in ((4,), (4, 3)):
            y = random.randn(*yshape) + offset
            actual = solve_toeplitz((c,r), b=y)
            desired = solve(toeplitz(c, r=r), y)
            assert_equal(actual.shape, yshape)
            assert_equal(desired.shape, yshape)
            assert_allclose(actual, desired)


def test_native_list_arguments():
    c = [1,2,4,7]
    r = [1,3,9,12]
    y = [5,1,4,2]
    actual = solve_toeplitz((c,r), y)
    desired = solve(toeplitz(c, r=r), y)
    assert_allclose(actual, desired)


def test_zero_diag_error():
    # The Levinson-Durbin implementation fails when the diagonal is zero.
    random = np.random.RandomState(1234)
    n = 4
    c = random.randn(n)
    r = random.randn(n)
    y = random.randn(n)
    c[0] = 0
    assert_raises(np.linalg.LinAlgError,
        solve_toeplitz, (c, r), b=y)


def test_wikipedia_counterexample():
    # The Levinson-Durbin implementation also fails in other cases.
    # This example is from the talk page of the wikipedia article.
    random = np.random.RandomState(1234)
    c = [2, 2, 1]
    y = random.randn(3)
    assert_raises(np.linalg.LinAlgError, solve_toeplitz, c, b=y)


def test_reflection_coeffs():
    # check that the partial solutions are given by the reflection
    # coefficients

    random = np.random.RandomState(1234)
    y_d = random.randn(10)
    y_z = random.randn(10) + 1j
    reflection_coeffs_d = [1]
    reflection_coeffs_z = [1]
    for i in range(2, 10):
        reflection_coeffs_d.append(solve_toeplitz(y_d[:(i-1)], b=y_d[1:i])[-1])
        reflection_coeffs_z.append(solve_toeplitz(y_z[:(i-1)], b=y_z[1:i])[-1])

    y_d_concat = np.concatenate((y_d[-2:0:-1], y_d[:-1]))
    y_z_concat = np.concatenate((y_z[-2:0:-1].conj(), y_z[:-1]))
    _, ref_d = levinson(y_d_concat, b=y_d[1:])
    _, ref_z = levinson(y_z_concat, b=y_z[1:])

    assert_allclose(reflection_coeffs_d, ref_d[:-1])
    assert_allclose(reflection_coeffs_z, ref_z[:-1])


@pytest.mark.xfail(reason='Instability of Levinson iteration')
def test_unstable():
    # this is a "Gaussian Toeplitz matrix", as mentioned in Example 2 of
    # I. Gohbert, T. Kailath and V. Olshevsky "Fast Gaussian Elimination with
    # Partial Pivoting for Matrices with Displacement Structure"
    # Mathematics of Computation, 64, 212 (1995), pp 1557-1576
    # which can be unstable for levinson recursion.

    # other fast toeplitz solvers such as GKO or Burg should be better.
    random = np.random.RandomState(1234)
    n = 100
    c = 0.9 ** (np.arange(n)**2)
    y = random.randn(n)

    solution1 = solve_toeplitz(c, b=y)
    solution2 = solve(toeplitz(c), y)

    assert_allclose(solution1, solution2)


@pytest.mark.parametrize('dt_c', [int, float, np.float32, complex, np.complex64])
@pytest.mark.parametrize('dt_b', [int, float, np.float32, complex, np.complex64])
def test_empty(dt_c, dt_b):
    c = np.array([], dtype=dt_c)
    b = np.array([], dtype=dt_b)
    x = solve_toeplitz(c, b)
    assert x.shape == (0,)
    if not (dt_b in {np.float32, np.complex64} and dt_c in {np.float32, np.complex64}):
        # see https://github.com/scipy/scipy/issues/25964#issuecomment-5521651145
        assert x.dtype == solve_toeplitz(np.array([2, 1], dtype=dt_c),
                                        np.ones(2, dtype=dt_b)).dtype

    b = np.empty((0, 0), dtype=dt_b)
    x1 = solve_toeplitz(c, b)
    assert x1.shape == (0, 0)
    assert x1.dtype == x.dtype


_test_dtypes = [np.int64, np.float32, np.float64, np.complex64, np.complex128]
@pytest.mark.parametrize("dt_c", _test_dtypes)
@pytest.mark.parametrize("dt_r", _test_dtypes)
@pytest.mark.parametrize("dt_b", _test_dtypes)
@pytest.mark.parametrize("n", [0, 5])
@pytest.mark.parametrize("nrhs", [0, 1, 2])
def test_shape_dtype(dt_c, dt_r, dt_b, n, nrhs):
    rng = np.random.default_rng(seed=12345)

    ## Setup
    c = rng.normal(size=(n,))
    if np.isdtype(dt_c, "complex floating"):
        c = c + 1j * rng.normal(size=(n,))
    c = np.astype(c, dt_c)

    r = rng.normal(size=(n,))
    if np.isdtype(dt_r, "complex floating"):
        r = r + 1j * rng.normal(size=(n,))
    r = np.astype(r, dt_r)

    shape_b = (n,) if nrhs == 0 else (n, nrhs)
    b = rng.normal(size=shape_b)
    if np.isdtype(dt_b, "complex floating"):
        b = b + 1j * rng.normal(size=shape_b)
    b = np.astype(b, dt_b)

    ## Compute solutions
    a = toeplitz(c, r)
    ref = solve(a, b, assume_a="general")
    res = solve_toeplitz((c, r), b)

    xp_assert_close(res, ref, atol=1e-15)


@pytest.mark.parametrize('dtype', [np.float32, np.float64, np.complex64, np.complex128])
@pytest.mark.parametrize('n', [1, 2, 3, 10])
def test_inverse_columns(dtype, n):
    # levinson_inverse_columns returns the first and last columns of the
    # inverse of the Toeplitz matrix
    rng = np.random.default_rng(1234)
    c = rng.standard_normal(n)
    r = rng.standard_normal(n)
    if np.issubdtype(dtype, np.complexfloating):
        c = c + 1j * rng.standard_normal(n)
        r = r + 1j * rng.standard_normal(n)
    c[0] = r[0] = 2 * n
    vals = np.concatenate((r[-1:0:-1], c)).astype(dtype)
    x, y = levinson_inverse_columns(vals)
    assert x.dtype == y.dtype == dtype
    inv = np.linalg.inv(toeplitz(c, r))
    rtol = 1e-5 if dtype in (np.float32, np.complex64) else 1e-12
    assert_allclose(x, inv[:, 0], rtol=rtol)
    assert_allclose(y, inv[:, -1], rtol=rtol)


def test_inverse_columns_singular():
    # Same failures as Levinson recursion
    for c in ([0., 1, 2], [2., 2, 1]):
        c = np.array(c)
        vals = np.concatenate((c[-1:0:-1], c))
        assert_raises(np.linalg.LinAlgError, levinson_inverse_columns, vals)


@pytest.mark.parametrize('dtype', [np.float32, np.float64, np.complex64, np.complex128])
@pytest.mark.parametrize('n', [1, 2, 3, 300])
@pytest.mark.parametrize('nrhs', [4, 7])
@pytest.mark.parametrize('block_size', [1, 2, 256])
def test_multiple_rhs_inverse_rows(dtype, n, nrhs, block_size):
    # Several right-hand sides are solved with the rows of the inverse, block
    # by block, which must match solving each column with Levinson recursion
    rng = np.random.default_rng(1234)
    c = rng.standard_normal(n)
    r = rng.standard_normal(n)
    b = rng.standard_normal((n, nrhs))
    if np.issubdtype(dtype, np.complexfloating):
        c = c + 1j * rng.standard_normal(n)
        r = r + 1j * rng.standard_normal(n)
        b = b + 1j * rng.standard_normal((n, nrhs))
    c[0] = r[0] = 2 * np.sqrt(n)
    vals = np.concatenate((r[-1:0:-1], c)).astype(dtype)
    b = b.astype(dtype)
    x = _solve_toeplitz_multiple_rhs(vals, b, block_size=block_size)
    assert x.dtype == dtype
    assert x.shape == b.shape
    desired = np.column_stack([levinson(vals, np.ascontiguousarray(b[:, i]))[0]
                               for i in range(nrhs)])
    rtol = 1e-4 if dtype in (np.float32, np.complex64) else 1e-11
    assert_allclose(x, desired, rtol=rtol, atol=rtol * np.abs(desired).max())


def test_multiple_rhs_matches_single_rhs():
    # Many columns go through the inverse rows, a single one through Levinson
    # recursion; both must solve the system. The matrix is a positive definite
    # covariance, as in Gaussian process regression on a regular grid.
    rng = np.random.default_rng(1234)
    n = 500
    c = np.exp(-(np.arange(n) / 10.)**2)
    c[0] += 1e-6
    b = rng.standard_normal((n, 20))
    x = solve_toeplitz(c, b)
    # The matrix has a condition number of about 2e7, so the two methods
    # agree to that accuracy in norm rather than element by element.
    for i in (0, 7, 19):
        xi = solve_toeplitz(c, b[:, i])
        assert np.linalg.norm(x[:, i] - xi) <= 1e-8 * np.linalg.norm(xi)
    xd = solve(toeplitz(c), b, assume_a='pos')
    assert np.linalg.norm(x - xd) <= 1e-8 * np.linalg.norm(xd)


def test_multiple_rhs_batch():
    # Batches of matrices with enough right-hand sides each to use the rows
    # of the inverse
    rng = np.random.default_rng(1234)
    c = rng.standard_normal((3, 6))
    c[:, 0] = 10
    b = rng.standard_normal((3, 6, 30))
    x = solve_toeplitz(c, b)
    for k in range(3):
        assert_allclose(x[k], solve(toeplitz(c[k]), b[k]), rtol=1e-12)
