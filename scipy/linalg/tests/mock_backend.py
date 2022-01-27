import numpy as np


class _MockFunction:
    def __init__(self, return_value=None):
        self.number_calls = 0
        self.return_value = return_value
        self.last_args = ([], {})

    def __call__(self, *args, **kwargs):
        self.number_calls += 1
        self.last_args = (args, kwargs)
        return self.return_value


method_names = [
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

for name in method_names:
    globals()[name] = _MockFunction(np.array([[0, 0], [1, 1]]))


__ua_domain__ = "numpy.scipy.linalg"


def __ua_function__(method, args, kwargs):
    fn = globals().get(method.__name__)
    return (fn(*args, **kwargs) if fn is not None else NotImplemented)
