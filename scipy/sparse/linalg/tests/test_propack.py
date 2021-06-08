import os
import pytest

import numpy as np
from numpy.testing import assert_allclose, assert_raises

from scipy.sparse.linalg.propack import svdp
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix


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


def load_coord(folder, precision, file="illc1850.coord"):
    dtype = {"single": np.float32, "double": np.float64}[precision]
    path = os.path.join(folder, precision, file)
    with open(path) as f:
        m, n, nnz = (int(val) for val in f.readline().split())
        coord = np.array([[float(val) for val in line.split()] for line in f])
    i = coord[:, 0].astype(int) - 1
    j = coord[:, 1].astype(int) - 1
    data = coord[:, 2].astype(dtype)
    A = coo_matrix((data, (i, j)))
    return A


def load_sigma(folder, precision="double", irl=False):
    dtype = {"single": np.float32, "double": np.float64}[precision]
    s_name = "Sigma_IRL.ascii" if irl else "Sigma.ascii"
    path = os.path.join(folder, precision, s_name)
    with open(path) as f:
        data = np.array([float(line.split()[1]) for line in f], dtype=dtype)
    return data


def load_uv(folder, precision="double", uv="U", irl=False):
    dtype = {"single": np.float32, "double": np.float64}[precision]
    uv_name = (uv + "_IRL.ascii") if irl else (uv + ".ascii")
    path = os.path.join(folder, precision, uv_name)
    with open(path) as f:
        m, n = (int(val) for val in f.readline().split())
        data = np.array([float(val.strip()) for val in f], dtype=dtype)
    return data.reshape((n, m)).T


@pytest.mark.parametrize('precision', ('single', 'double'))
def test_examples(precision):
    atol = {'single': 1e-3, 'double': 1e-11}[precision]

    path_prefix = os.path.dirname(__file__)
    relative_path = "propack_examples"
    folder = os.path.join(path_prefix, relative_path)

    A = load_coord(folder, precision)
    s_expected = load_sigma(folder, precision)
    u_expected = load_uv(folder, precision, "U")
    vt_expected = load_uv(folder, precision, "V").T

    k = len(s_expected)
    u, s, vt = svdp(A, k)

    assert_allclose(s, s_expected, atol)
    assert_allclose(np.abs(u), np.abs(u_expected), atol=atol)
    assert_allclose(np.abs(vt), np.abs(vt_expected), atol=atol)
