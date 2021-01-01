import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_raises

from scipy.sparse.linalg.propack import svdp
from scipy.sparse import csr_matrix, csc_matrix

TOLS = {
    np.float32: 1e-4,
    np.float64: 1e-8,
    np.complex64: 1e-4,
    np.complex128: 1e-8,
}


def is_complex_type(dtype):
    return np.dtype(dtype).char.isupper()


def generate_matrix(constructor, n, m, f,
                    dtype=float, rseed=0, **kwargs):
    """Generate a random sparse matrix"""
    rng = np.random.RandomState(rseed)
    if is_complex_type(dtype):
        M = (- 5 + 10 * rng.rand(n, m)
             - 5j + 10j * rng.rand(n, m)).astype(dtype)
    else:
        M = (-5 + 10 * rng.rand(n, m)).astype(dtype)
    M[M.real > 10 * f - 5] = 0
    return constructor(M, **kwargs)


def assert_orthogonal(u1, u2, rtol, atol):
    """Check that the first k rows of u1 and u2 are orthogonal"""
    I = abs(np.dot(u1.conj().T, u2))
    assert_allclose(I, np.eye(u1.shape[1], u2.shape[1]), rtol=rtol, atol=atol)


def check_svdp(n, m, constructor, dtype, k, irl_mode, which, f=0.8):
    tol = TOLS[dtype]

    M = generate_matrix(np.asarray, n, m, f, dtype)
    Msp = constructor(M)

    u1, sigma1, vt1 = np.linalg.svd(M, full_matrices=False)
    u2, sigma2, vt2 = svdp(Msp, k=k, which=which, irl_mode=irl_mode, tol=tol)

    # check the which
    if which.upper() == 'SM':
        u1 = np.roll(u1, k, 1)
        vt1 = np.roll(vt1, k, 0)
        sigma1 = np.roll(sigma1, k)
    elif which.upper() == 'LM':
        pass
    else:
        raise ValueError(f"which = '{which}' not recognized")

    # check that singular values agree
    assert_allclose(sigma1[:k], sigma2, rtol=tol, atol=tol)

    # check that singular vectors are orthogonal
    assert_orthogonal(u1, u2, rtol=tol, atol=tol)
    assert_orthogonal(vt1.T, vt2.T, rtol=tol, atol=tol)


@pytest.mark.parametrize('ctor', (np.array, csr_matrix, csc_matrix))
@pytest.mark.parametrize('dtype', (np.float32, np.float64, np.complex64, np.complex128))
@pytest.mark.parametrize('irl', (True, False))
@pytest.mark.parametrize('which', ('LM', 'SM'))
def test_svdp(ctor, dtype, irl, which):
    np.random.seed(0)
    n, m, k = 10, 20, 3
    if which == 'SM' and not irl:
        with assert_raises(ValueError):
            check_svdp(n, m, ctor, dtype, k, irl, which)
    else:
        check_svdp(n, m, ctor, dtype, k, irl, which)
