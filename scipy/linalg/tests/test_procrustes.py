from itertools import product, permutations

import numpy as np
from numpy.testing import assert_array_less, assert_allclose, assert_raises

from scipy.linalg import inv, eigh, norm
from scipy.linalg import orthogonal_procrustes


def test_orthogonal_procrustes_ndim_too_large():
    np.random.seed(1234)
    A = np.random.randn(3, 4, 5)
    B = np.random.randn(3, 4, 5)
    assert_raises(ValueError, orthogonal_procrustes, A, B)


def test_orthogonal_procrustes_ndim_too_small():
    np.random.seed(1234)
    A = np.random.randn(3)
    B = np.random.randn(3)
    assert_raises(ValueError, orthogonal_procrustes, A, B)


def test_orthogonal_procrustes_shape_mismatch():
    np.random.seed(1234)
    shapes = ((3, 3), (3, 4), (4, 3), (4, 4))
    for a, b in permutations(shapes, 2):
        A = np.random.randn(*a)
        B = np.random.randn(*b)
        assert_raises(ValueError, orthogonal_procrustes, A, B)


def test_orthogonal_procrustes_checkfinite_exception():
    np.random.seed(1234)
    m, n = 2, 3
    A_good = np.random.randn(m, n)
    B_good = np.random.randn(m, n)
    for bad_value in np.inf, -np.inf, np.nan:
        A_bad = A_good.copy()
        A_bad[1, 2] = bad_value
        B_bad = B_good.copy()
        B_bad[1, 2] = bad_value
        for A, B in ((A_good, B_bad), (A_bad, B_good), (A_bad, B_bad)):
            assert_raises(ValueError, orthogonal_procrustes, A, B)


def test_orthogonal_procrustes_array_conversion():
    np.random.seed(1234)
    for m, n in ((6, 4), (4, 4), (4, 6)):
        A_arr = np.random.randn(m, n)
        B_arr = np.random.randn(m, n)
        As = (A_arr, A_arr.tolist(), np.matrix(A_arr))
        Bs = (B_arr, B_arr.tolist(), np.matrix(B_arr))
        kwarg_dicts = ({}, dict(check_finite=True), dict(check_finite=False))
        R_arr = orthogonal_procrustes(A_arr, B_arr)
        AR_arr = A_arr.dot(R_arr)
        for A, B, kwargs in product(As, Bs, kwarg_dicts):
            R = orthogonal_procrustes(A, B, **kwargs)
            AR = A_arr.dot(R)
            assert_allclose(AR, AR_arr)


def test_orthogonal_procrustes():
    np.random.seed(1234)
    for m, n in ((6, 4), (4, 4), (4, 6)):
        # Sample a random target matrix.
        B = np.random.randn(m, n)
        # Sample a random orthogonal matrix
        # by computing eigh of a sampled symmetric matrix.
        X = np.random.randn(n, n)
        w, V = eigh(X.T + X)
        assert_allclose(inv(V), V.T)
        # Compute a matrix with a known orthogonal transformation that gives B.
        A = np.dot(B, V.T)
        # Check that an orthogonal transformation from A to B can be recovered.
        R = orthogonal_procrustes(A, B)
        assert_allclose(inv(R), R.T)
        assert_allclose(A.dot(R), B)
        # Create a perturbed input matrix.
        A_perturbed = A + 1e-2 * np.random.randn(m, n)
        # Check that the orthogonal procrustes function can find an orthogonal
        # transformation that is better than the orthogonal transformation
        # computed from the original input matrix.
        R_prime = orthogonal_procrustes(A_perturbed, B)
        assert_allclose(inv(R_prime), R_prime.T)
        # Compute the naive and optimal transformations of the perturbed input.
        naive_approx = A_perturbed.dot(R)
        optim_approx = A_perturbed.dot(R_prime)
        # Compute the Frobenius norm errors of the matrix approximations.
        naive_approx_error = norm(naive_approx - B, ord='fro')
        optim_approx_error = norm(optim_approx - B, ord='fro')
        # Check that the orthogonal Procrustes approximation is better.
        assert_array_less(optim_approx_error, naive_approx_error)

