import numpy as np
import scipy.linalg
from scipy.linalg import set_backend
from scipy.linalg.tests import mock_backend

from numpy.testing import assert_equal
import pytest


fnames = [
    # solvers
    'solve_sylvester',
    'solve_continuous_lyapunov', 'solve_discrete_lyapunov',
    'solve_lyapunov',
    'solve_continuous_are', 'solve_discrete_are',
    # decomp (eigen value problems)
    'eig', 'eigvals', 'eigh', 'eigvalsh',
    'eig_banded', 'eigvals_banded',
    'eigh_tridiagonal', 'eigvalsh_tridiagonal', 'hessenberg', 'cdf2rdf',
    # matrix functions
    'expm', 'cosm', 'sinm', 'tanm', 'coshm', 'sinhm',
    'tanhm', 'logm', 'funm', 'signm', 'sqrtm',
    'expm_frechet', 'expm_cond', 'fractional_matrix_power',
    'khatri_rao',
    # sketches
    'clarkson_woodruff_transform',
    # special matrices
    'tri', 'tril', 'triu', 'toeplitz', 'circulant', 'hankel',
    'hadamard', 'leslie', 'kron', 'companion',
    'fiedler', 'fiedler_companion', 'convolution_matrix',
    # decompositions
    'cholesky', 'cho_factor', 'cho_solve', 'cholesky_banded',
    'cho_solve_banded', 'ldl', 'lu', 'lu_solve', 'lu_factor',
    'polar', 'qr', 'qr_multiply', 'rq', 'qz', 'ordqz', 'schur', 'rsf2csf',
    'svd', 'svdvals', 'diagsvd', 'orth', 'subspace_angles', 'null_space',
    'qr_delete', 'qr_insert', 'qr_update',
    # basic
    'solve', 'solve_triangular', 'solveh_banded', 'solve_banded',
    'solve_toeplitz', 'solve_circulant', 'inv', 'det', 'lstsq',
    'pinv', 'pinv2', 'pinvh', 'matrix_balance', 'matmul_toeplitz'
]


funcs = [getattr(scipy.linalg, fname) for fname in fnames]
mocks = [getattr(mock_backend, fname) for fname in fnames]


@pytest.mark.parametrize("func, mock", zip(funcs, mocks))
def test_backend_call(func, mock):
    """
    Checks fake backend dispatch.
    """
    x = np.array([[0, 0], [1, 1]])

    with set_backend(mock_backend, only=True):
        mock.number_calls = 0
        y = func(x)
        assert_equal(y, mock.return_value)
        assert_equal(mock.number_calls, 1)
