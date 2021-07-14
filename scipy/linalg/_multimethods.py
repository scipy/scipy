import functools
import numpy as np
from scipy._lib.uarray import all_of_type, create_multimethod
from scipy.linalg import _api


__all__ = ['eig', 'eigvals', 'eigh', 'eigvalsh',
           'eig_banded', 'eigvals_banded',
           'eigh_tridiagonal', 'eigvalsh_tridiagonal']


_create_linalg = functools.partial(create_multimethod,
                                   domain="numpy.scipy.linalg")


def _get_docs(func):
    func.__doc__ = getattr(_api, func.__name__).__doc__
    return func


############################################
# decomp (eigenvalue problems) multimethods
############################################


def _a_b_replacer(args, kwargs, dispatchables):
    def self_method(a, b=None, *args, **kwargs):
        return (dispatchables[0], dispatchables[1]) + args, kwargs

    return self_method(*args, **kwargs)


def _a_band_replacer(args, kwargs, dispatchables):
    def self_method(a_band, *args, **kwargs):
        return (dispatchables[0],) + args, kwargs

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


@_create_linalg(_a_band_replacer)
@all_of_type(np.ndarray)
@_get_docs
def eig_banded(a_band, lower=False, eigvals_only=False, overwrite_a_band=False,
               select='a', select_range=None, max_ev=0, check_finite=True):
    return a_band


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


@_create_linalg(_a_band_replacer)
@all_of_type(np.ndarray)
@_get_docs
def eigvals_banded(a_band, lower=False, overwrite_a_band=False,
                   select='a', select_range=None, check_finite=True):
    return a_band


@_create_linalg(_a_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def eigvalsh_tridiagonal(d, e, select='a', select_range=None,
                         check_finite=True, tol=0., lapack_driver='auto'):
    return d, e


@_create_linalg(_a_b_replacer)
@all_of_type(np.ndarray)
@_get_docs
def eigh_tridiagonal(d, e, eigvals_only=False, select='a', select_range=None,
                     check_finite=True, tol=0., lapack_driver='auto'):
    return d, e