import numpy as np
from scipy.linalg import get_array_bandwidth
import pytest
from pytest import raises


def test_get_array_bandwidth_dtypes():
    n = 5
    for t in np.typecodes['All']:
        A = np.zeros([n, n], dtype=t)
        if t in 'eUVOMm':
            raises(TypeError, get_array_bandwidth, A)
        else:  # no-op
            _ = get_array_bandwidth(A)


def test_get_array_bandwidth_non2d_input():
    A = np.array([1, 2, 3])
    raises(ValueError, get_array_bandwidth, A)
    A = np.array([[[1, 2, 3], [4, 5, 6]]])
    raises(ValueError, get_array_bandwidth, A)


@pytest.mark.parametrize('T', [x for x in np.typecodes['All']
                               if x not in 'eUVOMm'])
def test_get_array_bandwidth_square_inputs(T):
    n = 20
    k = 4
    R = np.zeros([n, n], dtype=T, order='F')
    # form a banded matrix inplace
    R[[x for x in range(n)], [x for x in range(n)]] = 1
    R[[x for x in range(n-k)], [x for x in range(k, n)]] = 1
    R[[x for x in range(1, n)], [x for x in range(n-1)]] = 1
    R[[x for x in range(k, n)], [x for x in range(n-k)]] = 1
    assert get_array_bandwidth(R) == (k, k)


@pytest.mark.parametrize('T', [x for x in np.typecodes['All']
                               if x not in 'eUVOMm'])
def test_get_array_bandwidth_rect_inputs(T):
    n, m = 10, 20
    k = 5
    R = np.zeros([n, m], dtype=T, order='F')
    # form a banded matrix inplace
    R[[x for x in range(n)], [x for x in range(m)]] = 1
    R[[x for x in range(n-k)], [x for x in range(k, m)]] = 1
    R[[x for x in range(1, n)], [x for x in range(m-1)]] = 1
    R[[x for x in range(k, n)], [x for x in range(m-k)]] = 1
    assert get_array_bandwidth(R) == (k, k)

