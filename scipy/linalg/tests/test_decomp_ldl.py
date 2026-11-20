import numpy as np
from numpy.testing import (assert_array_almost_equal, assert_array_equal,
                           assert_allclose, assert_)
from numpy import (array, eye, zeros, empty_like, empty, tril_indices_from,
                   tril, triu_indices_from, spacing, float32, float64,
                   complex64, complex128)
from numpy.exceptions import ComplexWarning
from scipy.linalg import ldl
import pytest


def test_args():
    A = eye(3)
    # Nonsquare array
    with pytest.raises(ValueError):
        ldl(A[:, :2])
    # Complex matrix with imaginary diagonal entries with "hermitian=True"
    with pytest.warns(ComplexWarning):
        ldl(A*1j)


def test_empty_array():
    a = empty((0, 0), dtype=complex)
    l, d, p = ldl(empty((0, 0)))
    assert_array_almost_equal(l, empty_like(a))
    assert_array_almost_equal(d, empty_like(a))
    assert_array_almost_equal(p, array([], dtype=int))


def test_simple():
    a = array([[-0.39-0.71j, 5.14-0.64j, -7.86-2.96j, 3.80+0.92j],
               [5.14-0.64j, 8.86+1.81j, -3.52+0.58j, 5.32-1.59j],
               [-7.86-2.96j, -3.52+0.58j, -2.83-0.03j, -1.54-2.86j],
               [3.80+0.92j, 5.32-1.59j, -1.54-2.86j, -0.56+0.12j]])
    b = array([[5., 10, 1, 18],
               [10., 2, 11, 1],
               [1., 11, 19, 9],
               [18., 1, 9, 0]])
    c = array([[52., 97, 112, 107, 50],
               [97., 114, 89, 98, 13],
               [112., 89, 64, 33, 6],
               [107., 98, 33, 60, 73],
               [50., 13, 6, 73, 77]])

    d = array([[2., 2, -4, 0, 4],
               [2., -2, -2, 10, -8],
               [-4., -2, 6, -8, -4],
               [0., 10, -8, 6, -6],
               [4., -8, -4, -6, 10]])
    e = array([[-1.36+0.00j, 0+0j, 0+0j, 0+0j],
               [1.58-0.90j, -8.87+0j, 0+0j, 0+0j],
               [2.21+0.21j, -1.84+0.03j, -4.63+0j, 0+0j],
               [3.91-1.50j, -1.78-1.18j, 0.11-0.11j, -1.84+0.00j]])
    for x in (b, c, d):
        l, d, p = ldl(x)
        assert_allclose(l.dot(d).dot(l.T), x, atol=spacing(1000.), rtol=0)

        u, d, p = ldl(x, lower=False)
        assert_allclose(u.dot(d).dot(u.T), x, atol=spacing(1000.), rtol=0)

    l, d, p = ldl(a, hermitian=False)
    assert_allclose(l.dot(d).dot(l.T), a, atol=spacing(1000.), rtol=0)

    u, d, p = ldl(a, lower=False, hermitian=False)
    assert_allclose(u.dot(d).dot(u.T), a, atol=spacing(1000.), rtol=0)

    # Use upper part for the computation and use the lower part for comparison
    l, d, p = ldl(e.conj().T, lower=0)
    assert_allclose(tril(l.dot(d).dot(l.conj().T)-e), zeros((4, 4)),
                    atol=spacing(1000.), rtol=0)


def test_permutations():
    rng = np.random.default_rng(1234)
    for _ in range(10):
        n = rng.integers(1, 100)
        # Random real/complex array
        x = rng.random((n, n)) + 0 if rng.integers(2) else rng.random((n, n))*1j
        x = x + x.conj().T
        x += eye(n)*rng.integers(5, 1e6)
        l_ind = tril_indices_from(x, k=-1)
        u_ind = triu_indices_from(x, k=1)

        # Test whether permutations lead to a triangular array
        u, d, p = ldl(x, lower=0)
        # lower part should be zero
        assert_(not any(u[p, :][l_ind]), f'Spin {_} failed')

        l, d, p = ldl(x, lower=1)
        # upper part should be zero
        assert_(not any(l[p, :][u_ind]), f'Spin {_} failed')


@pytest.mark.parametrize("dtype", [float32, float64])
@pytest.mark.parametrize("n", [30, 150])
def test_ldl_type_size_combinations_real(n, dtype):
    rng = np.random.default_rng(1234)
    msg = (f"Failed for size: {n}, dtype: {dtype}")

    x = rng.random((n, n)).astype(dtype)
    x = x + x.T
    x += eye(n, dtype=dtype)*dtype(rng.integers(5, 1e6))

    l, d1, p = ldl(x)
    u, d2, p = ldl(x, lower=0)
    rtol = 1e-4 if dtype is float32 else 1e-10
    assert_allclose(l.dot(d1).dot(l.T), x, rtol=rtol, err_msg=msg)
    assert_allclose(u.dot(d2).dot(u.T), x, rtol=rtol, err_msg=msg)


@pytest.mark.parametrize("dtype", [complex64, complex128])
@pytest.mark.parametrize("n", [30, 150])
def test_ldl_type_size_combinations_complex(n, dtype):
    rng = np.random.default_rng(1234)
    msg1 = (f"Her failed for size: {n}, dtype: {dtype}")
    msg2 = (f"Sym failed for size: {n}, dtype: {dtype}")

    # Complex hermitian upper/lower
    x = (rng.random((n, n))+1j*rng.random((n, n))).astype(dtype)
    x = x+x.conj().T
    x += eye(n, dtype=dtype)*dtype(rng.integers(5, 1e6))

    l, d1, p = ldl(x)
    u, d2, p = ldl(x, lower=0)
    rtol = 2e-4 if dtype is complex64 else 1e-10
    assert_allclose(l.dot(d1).dot(l.conj().T), x, rtol=rtol, err_msg=msg1)
    assert_allclose(u.dot(d2).dot(u.conj().T), x, rtol=rtol, err_msg=msg1)

    # Complex symmetric upper/lower
    x = (rng.random((n, n))+1j*rng.random((n, n))).astype(dtype)
    x = x+x.T
    x += eye(n, dtype=dtype)*dtype(rng.integers(5, 1e6))

    l, d1, p = ldl(x, hermitian=0)
    u, d2, p = ldl(x, lower=0, hermitian=0)
    assert_allclose(l.dot(d1).dot(l.T), x, rtol=rtol, err_msg=msg2)
    assert_allclose(u.dot(d2).dot(u.T), x, rtol=rtol, err_msg=msg2)


@pytest.mark.parametrize("dtype, hermitian", [
    (float32, True), (float64, True),
    (complex64, True), (complex128, True),
    (complex64, False), (complex128, False),
])
@pytest.mark.parametrize("n", [64, 65, 130])
@pytest.mark.parametrize("lower", [False, True])
@pytest.mark.parametrize("overwrite_a", [False, True])
@pytest.mark.parametrize("order", ["C", "F"])
def test_singular_factorization(dtype, hermitian, n, lower, overwrite_a, order):
    # gh-26330: older blocked LAPACK routines leave stale values in D when
    # a pivot becomes zero during elimination. Include two such pivots,
    # at opposite ends, as INFO only identifies the first one encountered.
    a = np.array(np.eye(n) * 0.04, dtype=dtype, order=order)
    v = np.array([1, 1j] if np.iscomplexobj(a) else [1, 1], dtype=dtype)
    block = np.outer(v, v.conj() if hermitian else v) * 0.0625
    a[:2, :2] = a[-2:, -2:] = block
    # Exercise a two-by-two pivot as well as the singular one-by-one pivots.
    a[3:5, 3:5] = [[0, 0.125], [0.125, 0]]
    expected = a.copy()

    lu, d, perm = ldl(a, lower=lower, hermitian=hermitian,
                      overwrite_a=overwrite_a)

    transpose = lu.conj().T if hermitian else lu.T
    assert_allclose(lu @ d @ transpose, expected, rtol=0, atol=1e-7)
    triangular = np.tril if lower else np.triu
    assert_array_equal(lu[perm], triangular(lu[perm]))
    assert lu.dtype == d.dtype == a.dtype
    if not overwrite_a:
        assert_array_equal(a, expected)


@pytest.mark.parametrize("lower", [False, True])
def test_singular_factorization_selected_triangle(lower):
    a = np.eye(130)
    a[:2, :2] = a[-2:, -2:] = 0.25
    # Scatter the singular blocks so triangularization requires permutation.
    perm = np.random.default_rng(1749).permutation(len(a))
    expected = a[np.ix_(perm, perm)]
    a = np.array(expected, order="F")
    unused = (np.triu_indices_from(a, 1) if lower
              else np.tril_indices_from(a, -1))
    a[unused] = np.nan

    lu, d, perm = ldl(a, lower=lower, overwrite_a=True, check_finite=False)

    assert_allclose(lu @ d @ lu.T, expected, rtol=0, atol=1e-14)
    triangular = np.tril if lower else np.triu
    assert_array_equal(lu[perm], triangular(lu[perm]))


def test_singular_factorization_batch():
    a = np.stack([np.eye(65), np.eye(65)])
    a[0, :2, :2] = 0.25
    a[1, -2:, -2:] = 0.5
    expected = a.copy()

    lu, d, perm = ldl(a)

    assert_allclose(lu @ d @ lu.swapaxes(-1, -2), expected, rtol=0, atol=1e-14)
    assert_array_equal(a, expected)
