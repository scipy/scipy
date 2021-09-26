import os
import pytest
import sys

import numpy as np
from numpy.testing import assert_allclose
from pytest import raises as assert_raises
from scipy.sparse.linalg._svdp import _svdp
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix

TOLS = {
    np.float32: 1e-4,
    np.float64: 1e-8,
    np.complex64: 1e-4,
    np.complex128: 1e-8,
}

_dtype_map = {
    "single": np.float32,
    "double": np.float64,
    "complex8": np.complex64,
    "complex16": np.complex128,
}


def is_complex_type(dtype):
    return np.dtype(dtype).char.isupper()


def is_32bit():
    return sys.maxsize <= 2**32  # (usually 2**31-1 on 32-bit)


_dtype_testing = []
for dtype in _dtype_map:
    if 'complex' in dtype and is_32bit():
        # PROPACK has issues w/ complex on 32-bit; see gh-14433
        marks = [pytest.mark.skip]
    elif 'complex' in dtype:
        marks = [pytest.mark.slow]  # type: ignore[list-item]
    else:
        marks = []
    _dtype_testing.append(pytest.param(dtype, marks=marks))
_dtype_testing = tuple(_dtype_testing)  # type: ignore[assignment]


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
    A = abs(np.dot(u1.conj().T, u2))
    assert_allclose(A, np.eye(u1.shape[1], u2.shape[1]), rtol=rtol, atol=atol)


def check_svdp(n, m, constructor, dtype, k, irl_mode, which, f=0.8):
    tol = TOLS[dtype]

    M = generate_matrix(np.asarray, n, m, f, dtype)
    Msp = constructor(M)

    u1, sigma1, vt1 = np.linalg.svd(M, full_matrices=False)
    u2, sigma2, vt2, _ = _svdp(Msp, k=k, which=which, irl_mode=irl_mode,
                               tol=tol)

    # check the which
    if which.upper() == 'SM':
        u1 = np.roll(u1, k, 1)
        vt1 = np.roll(vt1, k, 0)
        sigma1 = np.roll(sigma1, k)

    # check that singular values agree
    assert_allclose(sigma1[:k], sigma2, rtol=tol, atol=tol)

    # check that singular vectors are orthogonal
    assert_orthogonal(u1, u2, rtol=tol, atol=tol)
    assert_orthogonal(vt1.T, vt2.T, rtol=tol, atol=tol)


@pytest.mark.parametrize('ctor', (np.array, csr_matrix, csc_matrix))
@pytest.mark.parametrize('precision', _dtype_testing)
@pytest.mark.parametrize('irl', (True, False))
@pytest.mark.parametrize('which', ('LM', 'SM'))
def test_svdp(ctor, precision, irl, which):
    np.random.seed(0)
    dtype = _dtype_map[precision]
    n, m, k = 10, 20, 3
    if which == 'SM' and not irl:
        message = "`which`='SM' requires irl_mode=True"
        with assert_raises(ValueError, match=message):
            check_svdp(n, m, ctor, dtype, k, irl, which)
    else:
        if is_32bit() and is_complex_type(dtype):
            message = 'PROPACK complex-valued SVD methods not available '
            with assert_raises(TypeError, match=message):
                check_svdp(n, m, ctor, dtype, k, irl, which)
        else:
            check_svdp(n, m, ctor, dtype, k, irl, which)


def load_real(folder, precision, file='illc1850.coord'):
    dtype = _dtype_map[precision]
    path = os.path.join(folder, precision, 'Examples', file)
    # Coordinate Text File
    with open(path) as f:
        m, n, nnz = (int(val) for val in f.readline().split())
        coord = np.array([[float(val) for val in line.split()] for line in f])
    i = coord[:, 0].astype(int) - 1
    j = coord[:, 1].astype(int) - 1
    data = coord[:, 2].astype(dtype)
    A = coo_matrix((data, (i, j)))
    return A


def load_complex(folder, precision):
    file = 'mhd1280b.cua'

    dtype = _dtype_map[precision]
    path = os.path.join(folder, precision, 'Examples', file)
    with open(path, "r") as f:
        contents = f.readlines()
    file_metadata = contents[1].split()
    matrix_metadata = contents[2].split()
    datum_length = 15  # hard code rather than getting from contents[3]

    n_header = 4
    n_total, n_indptr, n_indices, n_data, _ = (int(n) for n in file_metadata)
    m, n, nnz, _ = (int(n) for n in matrix_metadata[1:])

    line_indptr = n_header
    line_indices = line_indptr + n_indptr
    line_data = line_indices + n_indices

    def _concatenate_lines(lines):
        return "".join([line.rstrip() for line in lines])

    indptr = _concatenate_lines(contents[line_indptr:line_indices])
    indptr = np.asarray([int(i) for i in indptr.split()])-1

    indices = _concatenate_lines(contents[line_indices:line_data])
    indices = np.asarray([int(i) for i in indices.split()])-1

    data = _concatenate_lines(contents[line_data:])
    data = np.asarray([float(data[i:i+datum_length])
                       for i in range(0, len(data), datum_length)])
    real, imag = data[::2], data[1::2]
    data = real + imag*1.0j
    data.astype(dtype)

    return csc_matrix((data, indices, indptr), (m, n))


def load_sigma(folder, precision="double", irl=False):
    dtype = _dtype_map[precision]
    s_name = "Sigma_IRL.ascii" if irl else "Sigma.ascii"
    path = os.path.join(folder, precision, 'Examples', 'Output', s_name)
    pydtype = {
        np.complex64: complex,
        np.complex128: complex,
        np.float32: float,
        np.float64: float,
    }[dtype]
    with open(path) as f:
        data = np.array([pydtype(line.split()[1]) for line in f]).astype(dtype)
    return data


def load_uv(folder, precision="double", uv="U", irl=False):
    dtype = _dtype_map[precision]
    uv_name = (uv + "_IRL.ascii") if irl else (uv + ".ascii")
    path = os.path.join(folder, precision, 'Examples', 'Output', uv_name)
    with open(path) as f:
        m, n = (int(val) for val in f.readline().split())
        if precision in {'single', 'double'}:
            data = np.array([dtype(val.strip()) for val in f], dtype=dtype)
        else:
            data = np.loadtxt(
                path, dtype=dtype, skiprows=1, delimiter='\n',
                converters={0: lambda s: complex(
                    *[float(n) for n in s.decode()[1:-1].split(',')])})
    return data.reshape((n, m)).T


@pytest.mark.parametrize('precision', _dtype_testing)
@pytest.mark.parametrize('irl', (False, True))
def test_examples(precision, irl):
    atol = {
        'single': 1e-3,
        'double': 1e-9,
        'complex8': 1e-4,
        'complex16': 1e-9,
    }[precision]

    path_prefix = os.path.dirname(__file__)
    relative_path = ['..', '_propack', 'PROPACK']
    folder = os.path.join(path_prefix, *relative_path)

    A = {
        'single': load_real,
        'double': load_real,
        'complex8': load_complex,
        'complex16': load_complex,
    }[precision](folder, precision)
    s_expected = load_sigma(folder, precision, irl=irl)
    u_expected = load_uv(folder, precision, "U", irl=irl)
    vh_expected = load_uv(folder, precision, "V", irl=irl).T

    k = len(s_expected)
    u, s, vh, _ = _svdp(A, k, irl_mode=irl, random_state=0)

    # Check singular values
    assert_allclose(s, s_expected.real, atol=atol)

    # complex example matrix has many repeated singular values, so check only
    # beginning non-repeated singular vectors to avoid permutations
    sv_check = 27 if precision in {'complex8', 'complex16'} else k
    u = u[:, :sv_check]
    u_expected = u_expected[:, :sv_check]
    vh = vh[:sv_check, :]
    vh_expected = vh_expected[:sv_check, :]
    s = s[:sv_check]

    assert_allclose(np.abs(u), np.abs(u_expected), atol=atol)
    assert_allclose(np.abs(vh), np.abs(vh_expected), atol=atol)

    # Check orthogonality of singular vectors
    assert_allclose(np.eye(u.shape[1]), u.conj().T @ u, atol=atol)
    assert_allclose(np.eye(vh.shape[0]), vh @ vh.conj().T, atol=atol)

    # Ensure the norm of the difference between the np.linalg.svd and
    # PROPACK reconstructed matrices is small
    u3, s3, vh3 = np.linalg.svd(A.todense())
    u3 = u3[:, :sv_check]
    s3 = s3[:sv_check]
    vh3 = vh3[:sv_check, :]
    A3 = u3 @ np.diag(s3) @ vh3
    recon = u @ np.diag(s) @ vh
    assert_allclose(np.linalg.norm(A3 - recon), 0, atol=atol)


@pytest.mark.parametrize('shifts', (None, -10, 0, 1, 10, 70))
@pytest.mark.parametrize('precision', ('single', 'double'))
def test_shifts(shifts, precision):
    np.random.seed(0)
    n, k = 70, 10
    A = np.random.random((n, n))
    if shifts is not None and ((shifts < 0) or (k > min(n-1-shifts, n))):
        with pytest.raises(ValueError):
            _svdp(A, k, shifts=shifts, kmax=5*k, irl_mode=True)
    else:
        _svdp(A, k, shifts=shifts, kmax=5*k, irl_mode=True)


@pytest.mark.slow
@pytest.mark.xfail()
def test_shifts_accuracy():
    np.random.seed(0)
    n, k = 70, 10
    A = np.random.random((n, n)).astype(np.double)
    u1, s1, vt1, _ = _svdp(A, k, shifts=None, which='SM', irl_mode=True)
    u2, s2, vt2, _ = _svdp(A, k, shifts=32, which='SM', irl_mode=True)
    # shifts <= 32 doesn't agree with shifts > 32
    # Does agree when which='LM' instead of 'SM'
    assert_allclose(s1, s2)
