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
    'solve_continuous_are', 'solve_discrete_are'
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
