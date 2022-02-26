import functools
import numpy as np
from scipy._lib.uarray import Dispatchable, all_of_type, create_multimethod
from scipy.special import _api
from scipy.special._backend import scalar_tuple_callable_array


__all__ = [
    # Error function and Fresnel integrals (ufuncs)
    'erf', 'erfc', 'erfcx', 'erfi'
]


_create_special = functools.partial(
                    create_multimethod,
                    domain="numpy.scipy.special"
                 )


_mark_scalar_tuple_callable_array = functools.partial(
                            Dispatchable,
                            dispatch_type=scalar_tuple_callable_array,
                            coercible=True
                        )


def _get_docs(func):
    func.__doc__ = getattr(_api, func.__name__).__doc__
    return func


def _z_replacer(args, kwargs, dispatchables):
    def self_method(z, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


def _x_out_replacer(args, kwargs, dispatchables):
    def self_method(x, out=None, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


def _z_out_replacer(args, kwargs, dispatchables):
    def self_method(z, out=None, *args, **kwargs):
        return dispatchables + args, kwargs

    return self_method(*args, **kwargs)


@_create_special(_z_replacer)
@all_of_type(np.ndarray)
@_get_docs
def erf(z):
    return (_mark_scalar_tuple_callable_array(z), )


@_create_special(_x_out_replacer)
@all_of_type(np.ndarray)
@_get_docs
def erfc(x, out=None):
    return _mark_scalar_tuple_callable_array(x), out


@_create_special(_x_out_replacer)
@all_of_type(np.ndarray)
@_get_docs
def erfcx(x, out=None):
    return _mark_scalar_tuple_callable_array(x), out


@_create_special(_z_out_replacer)
@all_of_type(np.ndarray)
@_get_docs
def erfi(z, out=None):
    return _mark_scalar_tuple_callable_array(z), out
