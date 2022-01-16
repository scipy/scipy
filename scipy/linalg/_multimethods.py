import functools
import numpy as np
from scipy._lib.uarray import Dispatchable, all_of_type, create_multimethod
from scipy.linalg import _api
from scipy.linalg._backend import scalar_tuple_callable_array

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
    'khatri_rao'
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
