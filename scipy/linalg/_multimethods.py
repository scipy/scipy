import functools
import numpy as np
from scipy._lib.uarray import Dispatchable, all_of_type, create_multimethod
from scipy.linalg import _api
from scipy.linalg._backend import scalar_tuple_callable_array

# they don't need to be dispatchabled
from ._api import (hilbert, helmert, invhilbert, pascal, invpascal, dft,
                   block_diag, cossin)



__all__ = [
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
    'hadamard', 'leslie', 'kron', 'block_diag', 'companion',
    'helmert', 'hilbert', 'invhilbert', 'pascal', 'invpascal', 'dft',
    'fiedler', 'fiedler_companion', 'convolution_matrix',
    # decompositions
    'cholesky', 'cho_factor', 'cho_solve', 'cholesky_banded',
    'cho_solve_banded', 'cossin', 'ldl', 'lu', 'lu_solve', 'lu_factor',
    'polar', 'qr', 'qr_multiply', 'rq', 'qz', 'ordqz', 'schur', 'rsf2csf',
    'svd', 'svdvals', 'diagsvd', 'orth', 'subspace_angles', 'null_space',
    'qr_delete', 'qr_insert', 'qr_update',
    # basic
    'solve', 'solve_triangular', 'solveh_banded', 'solve_banded',
    'solve_toeplitz', 'solve_circulant', 'inv', 'det', 'lstsq',
    'pinv', 'pinv2', 'pinvh', 'matrix_balance', 'matmul_toeplitz'
]


_create_linalg = functools.partial(
                    create_multimethod,
                    domain="numpy.scipy.linalg"
                 )


_mark_scalar_tuple_callable_array = functools.partial(
                            Dispatchable,
                            dispatch_type=scalar_tuple_callable_array,
                            coercible=True
                        )


def _get_docs(func):
    func.__doc__ = getattr(_api, func.__name__).__doc__
    return func


def _a_b_q_replacer(args, kwargs, dispatchables):
    def self_method(a, b, q, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_a_b_q_replacer)
@all_of_type(np.ndarray)
@_get_docs
def solve_sylvester(a, b, q):
    return a, b, q


def _a_q_replacer(args, kwargs, dispatchables):
    def self_method(a, q, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_a_q_replacer)
@all_of_type(np.ndarray)
@_get_docs
def solve_continuous_lyapunov(a, q):
    return a, q


@_create_linalg(_a_q_replacer)
@all_of_type(np.ndarray)
@_get_docs
def solve_lyapunov(a, q):
    return a, q


@_create_linalg(_a_q_replacer)
@all_of_type(np.ndarray)
@_get_docs
def solve_discrete_lyapunov(a, q, method=None):
    return a, q


def _a_b_q_r_e_s_replacer(args, kwargs, dispatchables):
    def self_method(a, b, q, r, e=None, s=None, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_a_b_q_r_e_s_replacer)
@all_of_type(np.ndarray)
@_get_docs
def solve_continuous_are(a, b, q, r, e=None, s=None, balanced=True):
    return a, b, q, r, e, s


@_create_linalg(_a_b_q_r_e_s_replacer)
@all_of_type(np.ndarray)
@_get_docs
def solve_discrete_are(a, b, q, r, e=None, s=None, balanced=True):
    return a, b, q, r, e, s


###################### decomp (eigenvalue problems) ############################


def _a_b_replacer(args, kwargs, dispatchables):
    def self_method(a, b=None, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_a_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def eig(a, b=None, left=False, right=True, overwrite_a=False,
        overwrite_b=False, check_finite=True, homogeneous_eigvals=False):
    return a, b


@_create_linalg(_a_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def eigh(a, b=None, lower=True, eigvals_only=False, overwrite_a=False,
         overwrite_b=False, turbo=True, eigvals=None, type=1,
         check_finite=True, subset_by_index=None, subset_by_value=None,
         driver=None):
    return a, b


def _aband_replacer(args, kwargs, dispatchables):
    def self_method(a_band, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_aband_replacer)
@all_of_type(np.ndarray)
@_get_docs
def eig_banded(a_band, lower=False, eigvals_only=False, overwrite_a_band=False,
               select='a', select_range=None, max_ev=0, check_finite=True):
    return (a_band, )


@_create_linalg(_a_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def eigvals(a, b=None, overwrite_a=False, check_finite=True,
            homogeneous_eigvals=False):
    return a, b


@_create_linalg(_a_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def eigvalsh(a, b=None, lower=True, overwrite_a=False,
             overwrite_b=False, turbo=True, eigvals=None, type=1,
             check_finite=True, subset_by_index=None, subset_by_value=None,
             driver=None):
    return a, b


@_create_linalg(_aband_replacer)
@all_of_type(np.ndarray)
@_get_docs
def eigvals_banded(a_band, lower=False, overwrite_a_band=False,
                   select='a', select_range=None, check_finite=True):
    return (a_band, )


def _d_e_replacer(args, kwargs, dispatchables):
    def self_method(d, e, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_d_e_replacer)
@all_of_type(np.ndarray)
@_get_docs
def eigvalsh_tridiagonal(d, e, select='a', select_range=None,
                         check_finite=True, tol=0., lapack_driver='auto'):
    return d, e


@_create_linalg(_d_e_replacer)
@all_of_type(np.ndarray)
@_get_docs
def eigh_tridiagonal(d, e, eigvals_only=False, select='a', select_range=None,
                     check_finite=True, tol=0., lapack_driver='auto'):
    return d, e


def _a_replacer(args, kwargs, dispatchables):
    def self_method(a, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def hessenberg(a, calc_q=False, overwrite_a=False, check_finite=True):
    return (a, )


def _w_v_replacer(args, kwargs, dispatchables):
    def self_method(w, v, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_w_v_replacer)
@all_of_type(np.ndarray)
@_get_docs
def cdf2rdf(w, v):
    return w, v


############################## Matrix functions ################################

def _A_replacer(args, kwargs, dispatchables):
    def self_method(A, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def fractional_matrix_power(A, t):
    return (A, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def logm(A, disp=True):
    return (A, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def expm(A):
    return (A, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def cosm(A):
    return (A, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def sinm(A):
    return (A, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def tanm(A):
    return (A, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def coshm(A):
    return (A, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def sinhm(A):
    return (A, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def tanhm(A):
    return (A, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def funm(A, func, disp=True):
    return (A, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def signm(A, disp=True):
    return (A, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def sqrtm(A, disp=True, blocksize=64):
    return (A, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def expm_cond(A, check_finite=True):
    return (A, )


def _A_E_replacer(args, kwargs, dispatchables):
    def self_method(A, E, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_A_E_replacer)
@all_of_type(np.ndarray)
@_get_docs
def expm_frechet(A, E, method=None, compute_expm=True, check_finite=True):
    return A, E


@_create_linalg(_a_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def khatri_rao(a, b):
    return a, b

############################### sketches #######################################


def _inputmatrix_replacer(args, kwargs, dispatchables):
    def self_method(input_matrix, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_inputmatrix_replacer)
@all_of_type(np.ndarray)
@_get_docs
def clarkson_woodruff_transform(input_matrix, sketch_size, seed=None):
    return (input_matrix, )


############################### special matrices ###############################

def _N_M_k_dtype_replacer(args, kwargs, dispatchables):
    def self_method(N, M=None, k=0, dtype=None, *args, **kwargs):
        return (N, M, k, dispatchables[0]) + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_N_M_k_dtype_replacer)
@all_of_type(np.dtype)
@_get_docs
def tri(N, M=None, k=0, dtype=None):
    return (dtype, )


def _m_replacer(args, kwargs, dispatchables):
    def self_method(m, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_m_replacer)
@all_of_type(np.ndarray)
@_get_docs
def tril(m, k=0):
    return (m, )


@_create_linalg(_m_replacer)
@all_of_type(np.ndarray)
@_get_docs
def triu(m, k=0):
    return (m, )


def _c_r_replacer(args, kwargs, dispatchables):
    def self_method(c, r=None, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_c_r_replacer)
@all_of_type(np.ndarray)
@_get_docs
def toeplitz(c, r=None):
    return c, r


@_create_linalg(_c_r_replacer)
@all_of_type(np.ndarray)
@_get_docs
def hankel(c, r=None):
    return c, r


def _c_replacer(args, kwargs, dispatchables):
    def self_method(c, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_c_replacer)
@all_of_type(np.ndarray)
@_get_docs
def circulant(c):
    return (c, )


def _n_dtype_replacer(args, kwargs, dispatchables):
    def self_method(n, dtype=int, *args, **kwargs):
        return (n, dispatchables[0]) + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_n_dtype_replacer)
@all_of_type(np.dtype)
@_get_docs
def hadamard(n, dtype=int):
    return (dtype, )


def _f_s_replacer(args, kwargs, dispatchables):
    def self_method(f, s, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_f_s_replacer)
@all_of_type(np.ndarray)
@_get_docs
def leslie(f, s):
    return f, s


@_create_linalg(_a_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def kron(a, b):
    return a, b


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def companion(a):
    return (a, )


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def fiedler(a):
    return (a, )


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def fiedler_companion(a):
    return (a, )


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def convolution_matrix(a, n, mode='full'):
    return (a, )


################################# decompositions ###############################


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def cholesky(a, lower=False, overwrite_a=False, check_finite=True):
    return (a, )


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def cho_factor(a, lower=False, overwrite_a=False, check_finite=True):
    return (a, )


def _clower_b_replacer(args, kwargs, dispatchables):
    def self_method(c_and_lower, b, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_clower_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def cho_solve(c_and_lower, b, overwrite_b=False, check_finite=True):
    return _mark_scalar_tuple_callable_array(c_and_lower), b


def _ab_replacer(args, kwargs, dispatchables):
    def self_method(ab, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_ab_replacer)
@all_of_type(np.ndarray)
@_get_docs
def cholesky_banded(ab, overwrite_ab=False, lower=False, check_finite=True):
    return (ab, )


def _cblower_b_replacer(args, kwargs, dispatchables):
    def self_method(cb_and_lower, b, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_cblower_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def cho_solve_banded(cb_and_lower, b, overwrite_b=False, check_finite=True):
    return _mark_scalar_tuple_callable_array(cb_and_lower), b


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def ldl(A, lower=True, hermitian=True, overwrite_a=False, check_finite=True):
    return (A, )


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def lu_factor(a, overwrite_a=False, check_finite=True):
    return (a, )


def _lupiv_b_replacer(args, kwargs, dispatchables):
    def self_method(lu_and_piv, b, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_lupiv_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def lu_solve(lu_and_piv, b, trans=0, overwrite_b=False, check_finite=True):
    return _mark_scalar_tuple_callable_array(lu_and_piv), b


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def lu(a, permute_l=False, overwrite_a=False, check_finite=True):
    return (a, )


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def polar(a, side="right"):
    return (a, )


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def qr(a, overwrite_a=False, lwork=None, mode='full', pivoting=False,
       check_finite=True):
    return (a, )


def _a_c_replacer(args, kwargs, dispatchables):
    def self_method(a, c, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_a_c_replacer)
@all_of_type(np.ndarray)
@_get_docs
def qr_multiply(a, c, mode='right', pivoting=False, conjugate=False,
                overwrite_a=False, overwrite_c=False):
    return a, c


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def rq(a, overwrite_a=False, lwork=None, mode='full', check_finite=True):
    return (a, )


def _A_B_replacer(args, kwargs, dispatchables):
    def self_method(A, B, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_A_B_replacer)
@all_of_type(np.ndarray)
@_get_docs
def qz(A, B, output='real', lwork=None, sort=None, overwrite_a=False,
       overwrite_b=False, check_finite=True):
    return A, B


@_create_linalg(_A_B_replacer)
@all_of_type(np.ndarray)
@_get_docs
def ordqz(A, B, sort='lhp', output='real', overwrite_a=False,
          overwrite_b=False, check_finite=True):
    return A, B


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def schur(a, output='real', lwork=None, overwrite_a=False, sort=None,
          check_finite=True):
    return (a, )


def _T_Z_replacer(args, kwargs, dispatchables):
    def self_method(T, Z, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_T_Z_replacer)
@all_of_type(np.ndarray)
@_get_docs
def rsf2csf(T, Z, check_finite=True):
    return T, Z


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def svd(a, full_matrices=True, compute_uv=True, overwrite_a=False,
        check_finite=True, lapack_driver='gesdd'):
    return (a, )


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def svdvals(a, overwrite_a=False, check_finite=True):
    return (a, )


def _s_replacer(args, kwargs, dispatchables):
    def self_method(s, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_s_replacer)
@all_of_type(np.ndarray)
@_get_docs
def diagsvd(s, M, N):
    return (s, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def orth(A, rcond=None):
    return (A, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def null_space(A, rcond=None):
    return (A, )


@_create_linalg(_A_B_replacer)
@all_of_type(np.ndarray)
@_get_docs
def subspace_angles(A, B):
    return A, B


def _Q_R_replacer(args, kwargs, dispatchables):
    def self_method(Q, R, k, *args, **kwargs):
        return dispatchables + (k, ) + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_Q_R_replacer)
@all_of_type(np.ndarray)
@_get_docs
def qr_delete(Q, R, k, p=1, which='row', overwrite_qr=False,
              check_finite=True):
    return Q, R


def _Q_R_u_replacer(args, kwargs, dispatchables):
    def self_method(Q, R, u, k, *args, **kwargs):
        return dispatchables + (k, ) + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_Q_R_u_replacer)
@all_of_type(np.ndarray)
@_get_docs
def qr_insert(Q, R, u, k, which='row', rcond=None, overwrite_qru=False,
              check_finite=True):
    return Q, R, u


def _Q_R_u_v_replacer(args, kwargs, dispatchables):
    def self_method(Q, R, u, v, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_Q_R_u_v_replacer)
@all_of_type(np.ndarray)
@_get_docs
def qr_update(Q, R, u, v, overwrite_qruv=False, check_finite=True):
    return Q, R, u, v


################################ basic #########################################


@_create_linalg(_a_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def solve(a, b, sym_pos=False, lower=False, overwrite_a=False,
          overwrite_b=False, debug=None, check_finite=True, assume_a='gen',
          transposed=False):
    return a, b


@_create_linalg(_a_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def solve_triangular(a, b, trans=0, lower=False, unit_diagonal=False,
                     overwrite_b=False, debug=None, check_finite=True):
    return a, b


def _lu_ab_b_replacer(args, kwargs, dispatchables):
    def self_method(l_and_u, ab, b, *args, **kwargs):
        return (l_and_u, ) + dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_lu_ab_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def solve_banded(l_and_u, ab, b, overwrite_ab=False, overwrite_b=False,
                 debug=None, check_finite=True):
    return ab, b


def _ab_b_replacer(args, kwargs, dispatchables):
    def self_method(ab, b, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_ab_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def solveh_banded(ab, b, overwrite_ab=False, overwrite_b=False, lower=False,
                  check_finite=True):
    return ab, b


def _corcr_b_replacer(args, kwargs, dispatchables):
    def self_method(c_or_cr, b, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_corcr_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def solve_toeplitz(c_or_cr, b, check_finite=True):
    return _mark_scalar_tuple_callable_array(c_or_cr), b


def _c_b_replacer(args, kwargs, dispatchables):
    def self_method(c, b, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_c_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def solve_circulant(c, b, singular='raise', tol=None,
                    caxis=-1, baxis=0, outaxis=0):
    return c, b


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def inv(a, overwrite_a=False, check_finite=True):
    return (a, )


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def det(a, overwrite_a=False, check_finite=True):
    return (a, )


@_create_linalg(_a_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def lstsq(a, b, cond=None, overwrite_a=False, overwrite_b=False,
          check_finite=True, lapack_driver=None):
    return a, b


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def pinv(a, atol=None, rtol=None, return_rank=False, check_finite=True,
         cond=None, rcond=None):
    return (a, )


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def pinv2(a, cond=None, rcond=None, return_rank=False, check_finite=True):
    return (a, )


@_create_linalg(_a_replacer)
@all_of_type(np.ndarray)
@_get_docs
def pinvh(a, atol=None, rtol=None, lower=True, return_rank=False,
          check_finite=True, cond=None, rcond=None):
    return (a, )


@_create_linalg(_A_replacer)
@all_of_type(np.ndarray)
@_get_docs
def matrix_balance(A, permute=True, scale=True, separate=False,
                   overwrite_a=False):
    return (A, )


def _corcr_x_replacer(args, kwargs, dispatchables):
    def self_method(c_or_cr, x, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_linalg(_corcr_x_replacer)
@all_of_type(np.ndarray)
@_get_docs
def matmul_toeplitz(c_or_cr, x, check_finite=False, workers=None):
    return _mark_scalar_tuple_callable_array(c_or_cr), x
