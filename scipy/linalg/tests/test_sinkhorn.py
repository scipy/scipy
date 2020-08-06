import pytest
from pytest import raises as assert_raises
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from scipy.linalg import sinkhorn


def test_square():
    with assert_raises(ValueError, match="square"):
        sinkhorn(np.ones(4))

    with assert_raises(ValueError, match="square"):
        sinkhorn(np.ones((4, 5)))


def test_positive():
    with assert_raises(ValueError, match="strictly positive"):
        sinkhorn(-np.ones((5, 5)))


def test_complex():
    with assert_raises(ValueError, match="real"):
        sinkhorn(np.ones((5, 5)) + 1j)


@pytest.mark.parametrize("value", [float("inf"), float("NaN")])
def test_inf_nan(value):
    with assert_raises(ValueError, match="must not contain infs or NaNs"):
        a = np.ones((5, 5))
        a[0, 0] = value
        sinkhorn(a)


@pytest.mark.parametrize("seed", range(4))
def test_row_columns_sums(seed):
    n = 10
    rng = np.random.RandomState(seed=seed)
    A = rng.uniform(0, 10, (n, n))
    D1, S, D2 = sinkhorn(A)

    # test correct product
    assert_allclose(A, D1 @ S @ D2)

    # test doubly stochastic
    assert_allclose(np.sum(S, axis=0), 1)
    assert_allclose(np.sum(S, axis=1), 1)

    # test diagonal matrices
    assert_equal(D1, np.diag(np.diag(D1)))
    assert_equal(D2, np.diag(np.diag(D2)))
