import numpy as np
from numpy.linalg import norm
from numpy.testing import assert_equal, assert_allclose, assert_
from scipy.sparse.linalg._isolve import minres

from pytest import raises as assert_raises


def get_sample_problem(make_complex=False):
    # A random 10 x 10 symmetric matrix
    rng = np.random.RandomState(1234)
    if make_complex:
        matrix = rng.rand(10, 10) + 1j * rng.rand(10, 10)
    else:
        matrix = rng.rand(10, 10)
    if make_complex:
        matrix = matrix + matrix.T.conj()
    else:
        matrix = matrix + matrix.T
    # A random vector of length 10
    vector = rng.rand(10)
    return matrix, vector

def test_singular():
    for val in [False, True]:
        A, b = get_sample_problem(make_complex=val)
        A[0, ] = 0
        b[0] = 0
        xp, info = minres(A, b)
        assert_equal(info, 0)
        assert norm(A @ xp - b) <= 1e-5 * norm(b)


def test_x0_is_used_by():
    for val in [False, True]:
        A, b = get_sample_problem(make_complex=val)
        # Random x0 to feed minres
        rng = np.random.RandomState(12345)
        x0 = rng.rand(10)
        trace = []

        def trace_iterates(xk):
            trace.append(xk)
        minres(A, b, x0=x0, callback=trace_iterates)
        trace_with_x0 = trace

        trace = []
        minres(A, b, callback=trace_iterates)
        assert_(not np.array_equal(trace_with_x0[0], trace[0]))


def test_shift():
    for val in [False, True]:
        A, b = get_sample_problem(make_complex=val)
        shift = 0.5
        shifted_A = A - shift * np.eye(10)
        x1, info1 = minres(A, b, shift=shift)
        x2, info2 = minres(shifted_A, b)
        assert_equal(info1, 0)
        assert_allclose(x1, x2, rtol=1e-5)


def test_asymmetric_fail():
    """Asymmetric matrix should raise `ValueError` when check=True"""
    for val in [False, True]:
        A, b = get_sample_problem(make_complex=val)
        A[1, 2] = 1
        A[2, 1] = 2
        with assert_raises(ValueError):
            xp, info = minres(A, b, check=True)


def test_minres_non_default_x0():
    rng = np.random.RandomState(1234)
    rtol = 1e-6
    for val in [False, True]:
        a = rng.randn(5, 5) if not val else rng.randn(5, 5) + 1j * rng.randn(5, 5)
        a = np.dot(a, a.T) if not val else np.dot(a, a.T.conj())
        b = rng.randn(5)
        c = rng.randn(5)
        x = minres(a, b, x0=c, rtol=rtol)[0]
        assert norm(a @ x - b) <= rtol * norm(b)


def test_minres_precond_non_default_x0():
    rng = np.random.RandomState(12345)
    rtol = 1e-6
    for val in [False, True]:
        a = rng.randn(5, 5) if not val else rng.randn(5, 5) + 1j * rng.randn(5, 5)
        a = np.dot(a, a.T) if not val else np.dot(a, a.T.conj())
        b = rng.randn(5)
        c = rng.randn(5)
        m = rng.randn(5, 5)
        m = np.dot(m, m.T)
        x = minres(a, b, M=m, x0=c, rtol=rtol)[0]
        assert norm(a @ x - b) <= rtol * norm(b)


def test_minres_precond_exact_x0():
    rng = np.random.RandomState(1234)
    rtol = 1e-6
    a = np.eye(10)
    b = np.ones(10)
    c = np.ones(10)
    m = rng.randn(10, 10)
    m = np.dot(m, m.T)
    x = minres(a, b, M=m, x0=c, rtol=rtol)[0]
    assert norm(a @ x - b) <= rtol * norm(b)