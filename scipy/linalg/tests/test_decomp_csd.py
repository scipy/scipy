import pytest
import numpy as np
from numpy.linalg import eigh
from numpy.testing import assert_almost_equal
from scipy.stats import ortho_group

from scipy.linalg import decomp_csd


@pytest.mark.parametrize('m, p, q, dtype, decimals',
                         [
                             (3, 2, 1, np.complex128, 12),
                             (4, 2, 1, np.complex128, 12),
                             (10, 7, 3, np.complex128, 12),
                             (10, 6, 4, np.complex128, 12),
                             (9, 6, 3, np.complex128, 12),
                             (10, 3, 7, np.complex128, 12),
                             (10, 4, 6, np.complex128, 12),
                             (9, 3, 6, np.complex128, 12),
                             (9, 5, 4, np.complex128, 12),
                             (100, 50, 50, np.complex128, 12),
                             (100, 50, 50, np.complex64, 6)])
def test_cs_decomp_unitary(m, p, q, dtype, decimals):
    u = _random_unitary(m)

    u1, u2, v1h, v2h, theta = decomp_csd.cs_decomp(
        np.array(u, dtype=dtype), p, q)

    ud = _block_diag_square_matrices(u1, u2)
    vdh = _block_diag_square_matrices(v1h, v2h)
    cs = _construct_cs(theta, m, p, q)

    assert_almost_equal(u, ud @ cs @ vdh, decimals)


def _random_unitary(m):
    h = np.random.rand(m, m) + 1.9j * np.random.rand(m, m)
    h = h + h.conj().T
    _, u = eigh(h)
    return u


@pytest.mark.parametrize('m, p, q, dtype, decimals',
                         [
                             (3, 2, 1, np.float, 12),
                             (4, 2, 1, np.float, 12),
                             (10, 7, 3, np.float, 12),
                             (10, 6, 4, np.float, 12),
                             (9, 6, 3, np.float, 12),
                             (10, 3, 7, np.float, 12),
                             (10, 4, 6, np.float, 12),
                             (9, 3, 6, np.float, 12),
                             (9, 5, 4, np.float, 12),
                             (10, 5, 5, np.float, 6),
                             (10, 5, 5, np.double, 12),
                             (100, 50, 50, np.float, 6),
                             (100, 50, 50, np.double, 12)])
def test_cs_decomp_orthogonal(m, p, q, dtype, decimals):
    o = ortho_group.rvs(m)
    u1, u2, v1h, v2h, theta = decomp_csd.cs_decomp(
        np.array(o, dtype=dtype), p, q)
    ud = _block_diag_square_matrices(u1, u2)
    vdh = _block_diag_square_matrices(v1h, v2h)
    cs = _construct_cs(theta, m, p, q)

    assert_almost_equal(o, ud @ cs @ vdh, decimals)


def _block_diag_square_matrices(m1, m2):
    l1, l2 = len(m1), len(m2)
    z = lambda r, c: np.zeros((r, c))
    return np.block([
        [m1, z(l1, l2)],
        [z(l2, l1), m2]
    ])

# TODO(balintp): maybe it would be useful to actually return this matrix
#  from the cs_decomp routine itself or expose it as part of the package?
def _construct_cs(theta, m, p, q):
    r = min(p, m - p, q, m - q)
    assert len(theta) == r
    z = lambda r, c: np.zeros((r, c))
    i = lambda d: np.eye(d)
    COSINE = np.diag(np.cos(theta))
    SINE = np.diag(np.sin(theta))

    # identity matrix dimensions based on
    # https://www.nag.com/numeric/fl/nagdoc_fl26/pdf/f08/f08rnf.pdf
    i11 = min(p, q) - r
    i12 = min(p, m - q) - r
    i21 = min(m - p, q) - r
    i22 = min(m - p, m - q) - r
    CS = np.block([
        [i(i11), z(i11, r), z(i11, i21), z(i11, i22), z(i11, r), z(i11, i12)],
        [z(r, i11), COSINE, z(r, i21), z(r, i22), -SINE, z(r, i12)],
        [z(i12, i11), z(i12, r), z(i12, i21), z(i12, i22), z(i12, r), -i(i12)],
        [z(i22, i11), z(i22, r), z(i22, i21), i(i22), z(i22, r), z(i22, i12)],
        [z(r, i11), SINE, z(r, i21), z(r, i22), COSINE, z(r, i12)],
        [z(i21, i11), z(i21, r), i(i21), z(i21, i22), z(i21, r), z(i21, i12)],
    ]
    )
    return CS
