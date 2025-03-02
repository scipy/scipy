import os
import numpy as np
from numpy.testing import assert_allclose

import pytest

import scipy.special as sc


@pytest.mark.parametrize('x, expected', [
    (np.array([1000, 1]), np.array([0, -999])),

    # Expected value computed using mpmath (with mpmath.mp.dps = 200) and then
    # converted to float.
    (np.arange(4), np.array([-3.4401896985611953,
                             -2.4401896985611953,
                             -1.4401896985611953,
                             -0.44018969856119533]))
])
def test_log_softmax(x, expected):
    assert_allclose(sc.log_softmax(x), expected, rtol=1e-13)


@pytest.mark.parametrize('x, expected', [
    # we shouldn't return zero on the smallest subnormal input
    (np.array([-np.log(np.finfo(np.float32).smallest_subnormal), 0], dtype=np.float32),
     np.array([float.fromhex('-0x1.00000p-149'), float.fromhex('-0x1.9d1dap+6')],
              dtype=np.float32)),
    (np.array([-np.log(np.finfo(np.float64).smallest_subnormal), 0], dtype=np.float64),
     np.array([float.fromhex('-0x0.0000000000001p-1022'),
               float.fromhex('-0x1.74385446d71c3p+9')],
              dtype=np.float64)),
])
def test_log_softmax_use_atol_on_windows(x, expected):
    # Windows has a bugged implementation of floating point arithmetic, see
    # https://github.com/scipy/scipy/pull/19549#discussion_r1458551340
    if os.name == 'nt':  # Check if the operating system is Windows
        atol = np.finfo(x.dtype).smallest_subnormal
        assert_allclose(sc.log_softmax(x), expected, atol=atol)
    else:
        assert_allclose(sc.log_softmax(x), expected, rtol=1e-13)


@pytest.fixture
def log_softmax_x():
    x = np.arange(4)
    return x


@pytest.fixture
def log_softmax_expected():
    # Expected value computed using mpmath (with mpmath.mp.dps = 200) and then
    # converted to float.
    expected = np.array([-3.4401896985611953,
                         -2.4401896985611953,
                         -1.4401896985611953,
                         -0.44018969856119533])
    return expected


def test_log_softmax_translation(log_softmax_x, log_softmax_expected):
    # Translation property.  If all the values are changed by the same amount,
    # the softmax result does not change.
    x = log_softmax_x + 100
    expected = log_softmax_expected
    assert_allclose(sc.log_softmax(x), expected, rtol=1e-13)


def test_log_softmax_noneaxis(log_softmax_x, log_softmax_expected):
    # When axis=None, softmax operates on the entire array, and preserves
    # the shape.
    x = log_softmax_x.reshape(2, 2)
    expected = log_softmax_expected.reshape(2, 2)
    assert_allclose(sc.log_softmax(x), expected, rtol=1e-13)


@pytest.mark.parametrize('axis_2d, expected_2d', [
    (0, np.log(0.5) * np.ones((2, 2))),
    (1, np.array([[0, -999], [0, -999]]))
])
def test_axes(axis_2d, expected_2d):
    assert_allclose(
        sc.log_softmax([[1000, 1], [1000, 1]], axis=axis_2d),
        expected_2d,
        rtol=1e-13,
    )


@pytest.fixture
def log_softmax_2d_x():
    x = np.arange(8).reshape(2, 4)
    return x


@pytest.fixture
def log_softmax_2d_expected():
    # Expected value computed using mpmath (with mpmath.mp.dps = 200) and then
    # converted to float.
    expected = np.array([[-3.4401896985611953,
                         -2.4401896985611953,
                         -1.4401896985611953,
                         -0.44018969856119533],
                        [-3.4401896985611953,
                         -2.4401896985611953,
                         -1.4401896985611953,
                         -0.44018969856119533]])
    return expected


def test_log_softmax_2d_axis1(log_softmax_2d_x, log_softmax_2d_expected):
    x = log_softmax_2d_x
    expected = log_softmax_2d_expected
    assert_allclose(sc.log_softmax(x, axis=1), expected, rtol=1e-13)


def test_log_softmax_2d_axis0(log_softmax_2d_x, log_softmax_2d_expected):
    x = log_softmax_2d_x.T
    expected = log_softmax_2d_expected.T
    assert_allclose(sc.log_softmax(x, axis=0), expected, rtol=1e-13)


def test_log_softmax_3d(log_softmax_2d_x, log_softmax_2d_expected):
    # 3-d input, with a tuple for the axis.
    x_3d = log_softmax_2d_x.reshape(2, 2, 2)
    expected_3d = log_softmax_2d_expected.reshape(2, 2, 2)
    assert_allclose(sc.log_softmax(x_3d, axis=(1, 2)), expected_3d, rtol=1e-13)


def test_log_softmax_scalar():
    assert_allclose(sc.log_softmax(1.0), 0.0, rtol=1e-13)
