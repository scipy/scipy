import pytest
import numpy as np
from numpy.testing import assert_almost_equal, assert_allclose

from scipy.stats import ortho_group, unitary_group
from scipy.linalg import cossin

REAL_DTYPES = [np.float32, np.float64]
COMPLEX_DTYPES = [np.complex64, np.complex128]
DTYPES = REAL_DTYPES + COMPLEX_DTYPES


@pytest.mark.parametrize('dtype_', DTYPES)
@pytest.mark.parametrize('m, p, q',
                         [
                             (2, 1, 1),
                             (3, 2, 1),
                             (3, 1, 2),
                             (4, 2, 2),
                             (4, 1, 2),
                             (4, 2, 1),
                             (4, 3, 1),
                             (4, 1, 3),
                             (100, 50, 1),
                             (100, 50, 50),
                         ])
def test_cossin(dtype_, m, p, q):
    if dtype_ in COMPLEX_DTYPES:
        x = np.array(unitary_group.rvs(m), dtype=dtype_)
    else:
        x = np.array(ortho_group.rvs(m), dtype=dtype_)

    u, cs, vh = cossin(x, (p, q))
    assert_allclose(x, u @ cs @ vh, rtol=0., atol=1e4 * np.finfo(dtype_).eps)
    assert u.dtype == dtype_
    assert cs.dtype == dtype_
    assert vh.dtype == dtype_

    u, cs, vh = cossin([x[:p, :q], x[:p, q:], x[p:, :q], x[p:, q:]])
    assert_allclose(x, u @ cs @ vh, rtol=0., atol=1e4 * np.finfo(dtype_).eps)
    assert u.dtype == dtype_
    assert cs.dtype == dtype_
    assert vh.dtype == dtype_

    _, cs2, vh2 = cossin(x, (p, q), compute_u=False)
    assert_allclose(cs, cs2, rtol=0., atol=0.)
    assert_allclose(vh, vh2, rtol=0., atol=0.)

    u2, cs2, _ = cossin(x, (p, q), compute_vh=False)
    assert_allclose(u, u2, rtol=0., atol=0.)
    assert_allclose(cs, cs2, rtol=0., atol=0.)

    _, cs2, _ = cossin(x, (p, q), compute_u=False, compute_vh=False)
    assert_allclose(cs, cs2, rtol=0., atol=0.)

def test_cossin_mixed_types():
    x = np.array(ortho_group.rvs(4), dtype=np.float)
    u, cs, vh = cossin([x[:2, :2],
                        np.array(x[:2, 2:], dtype=np.complex128),
                        x[2:, :2],
                        x[2:, 2:]])

    assert u.dtype == np.complex128
    assert cs.dtype == np.complex128
    assert vh.dtype == np.complex128
    assert_allclose(x, u @ cs @ vh, rtol=0.,
                    atol=1e4 * np.finfo(np.complex128).eps)

def test_cossin_error_incorrect_subblocks():
    with pytest.raises(ValueError, match="invalid submatrix dimensions.*"
                                           "[(1, 2), (1, 3), (1, 2), (1, 3)]"):
        cossin(([1, 2], [3, 4, 5], [6, 7], [8, 9, 10]))


def test_cossin_error_missing_partitioning():
    with pytest.raises(ValueError, match=".*exactly four submatrices.* got 2"):
        cossin(unitary_group.rvs(2))

    with pytest.raises(ValueError, match="you forgot.*partitioning"):
        cossin(unitary_group.rvs(4))


def test_cossin_error_non_iterable():
    with pytest.raises(ValueError, match="must be Iterable"):
        cossin(12j)


def test_cossin_error_non_square():
    with pytest.raises(ValueError, match="only supports square"):
        cossin([[1, 2]], (1, 1))

