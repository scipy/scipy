import functools
import operator
from dataclasses import dataclass
from types import ModuleType

import numpy as np
from scipy._lib._array_api import (
    array_namespace, scipy_namespace_for, is_numpy, is_dask, is_marray,
    xp_promote, SCIPY_ARRAY_API
)
import scipy._lib.array_api_extra as xpx
from . import _basic, _spfun_stats, _ufuncs
# These are imported so that type checkers know the signatures of the exported
# functions.
from ._basic import digamma, polygamma, zeta
from ._spfun_stats import multigammaln
from ._ufuncs import (
    beta, betainc, betaincc, betaln, binom, chdtr, chdtrc, entr, erf, erfc, erfinv,
    expi, expit, expn, gamma, gammainc, gammaincc, gammaln, i0, i0e, i1, i1e, log_ndtr,
    logit, ndtr, ndtri, rel_entr, stdtr, stdtrit, xlogy
)


__all__ = [
    'beta', 'betainc', 'betaincc', 'betaln', 'binom', 'chdtr',  # noqa: F822
    'chdtrc', 'digamma', 'entr', 'erf', 'erfc', 'erfinv', 'expi',  # noqa: F822
    'expit', 'expn', 'gamma', 'gammainc', 'gammaincc', 'gammaln',  # noqa: F822
    'i0', 'i0e', 'i1', 'i1e', 'log_ndtr', 'logit', 'multigammaln',  # noqa: F822
    'ndtr', 'ndtri', 'polygamma', 'rel_entr', 'stdtr', 'stdtrit',  # noqa: F822
    'xlogy', 'zeta']  # noqa: F822


array_api_compat_prefix = "scipy._lib.array_api_compat"


def get_array_special_func(f_name, f_module, xp):
    if is_numpy(xp):
        return getattr(f_module, f_name)

    spx = scipy_namespace_for(xp)
    if spx is not None:
        f = getattr(spx.special, f_name, None)
        if f is not None:
            return f

    # if generic array-API implementation is available, use that;
    # otherwise, fall back to NumPy/SciPy
    if f_name in _generic_implementations:
        f = _generic_implementations[f_name](xp=xp, spx=spx)
        if f is not None:
            return f

    def f(*args, **kwargs):
        if is_marray(xp):
            _f = globals()[f_name]  # Allow nested wrapping
            data_args = [arg.data for arg in args]
            out = _f(*data_args, **kwargs)
            mask = functools.reduce(operator.or_, (arg.mask for arg in args))
            return xp.asarray(out, mask=mask)

        elif is_dask(xp):
            # IMPORTANT: map_blocks works only because all ufuncs in this module
            # are elementwise. It would be a grave mistake to apply this to gufuncs
            # or any other function with reductions, as they would change their
            # output depending on chunking!

            _f = globals()[f_name]  # Allow nested wrapping
            # Hide dtype kwarg from map_blocks
            return xp.map_blocks(functools.partial(_f, **kwargs), *args)

        else:
            _f = getattr(f_module, f_name)
            args = [np.asarray(arg) for arg in args]
            out = _f(*args, **kwargs)
            return xp.asarray(out)

    return f


def _rel_entr(xp, spx):
    def __rel_entr(x, y, *, xp=xp):
        # https://github.com/data-apis/array-api-extra/issues/160
        mxp = array_namespace(x._meta, y._meta) if is_dask(xp) else xp
        x, y = xp_promote(x, y, broadcast=True, force_floating=True, xp=xp)
        xy_pos = (x > 0) & (y > 0)
        xy_inf = xp.isinf(x) & xp.isinf(y)
        res = xpx.apply_where(
            xy_pos & ~xy_inf,
            (x, y),
            # Note: for very large x, this can overflow.
            lambda x, y: x * (mxp.log(x) - mxp.log(y)),
            fill_value=xp.inf
        )
        res = xpx.at(res)[(x == 0) & (y >= 0)].set(0)
        res = xpx.at(res)[xp.isnan(x) | xp.isnan(y) | (xy_pos & xy_inf)].set(xp.nan)
        return res

    return __rel_entr


def _xlogy(xp, spx):
    def __xlogy(x, y, *, xp=xp):
        with np.errstate(divide='ignore', invalid='ignore'):
            temp = x * xp.log(y)
        return xp.where(x == 0., 0., temp)
    return __xlogy


def _get_native_func(xp, spx, f_name):
    f = getattr(spx.special, f_name, None) if spx else None
    if f is None and hasattr(xp, 'special'):
        f = getattr(xp.special, f_name, None)
    return f


def _chdtr(xp, spx):
    # The difference between this and just using `gammainc`
    # defined by `get_array_special_func` is that if `gammainc`
    # isn't found, we don't want to use the SciPy version; we'll
    # return None here and use the SciPy version of `chdtr`.
    gammainc = _get_native_func(xp, spx, 'gammainc')  # noqa: F811
    if gammainc is None:
        return None

    def __chdtr(v, x):
        res = gammainc(v / 2, x / 2)  # this is almost all we need
        # The rest can be removed when google/jax#20507 is resolved
        mask = (v == 0) & (x > 0)  # JAX returns NaN
        res = xp.where(mask, 1., res)
        mask = xp.isinf(v) & xp.isinf(x)  # JAX returns 1.0
        return xp.where(mask, xp.nan, res)
    return __chdtr


def _chdtrc(xp, spx):
    # The difference between this and just using `gammaincc`
    # defined by `get_array_special_func` is that if `gammaincc`
    # isn't found, we don't want to use the SciPy version; we'll
    # return None here and use the SciPy version of `chdtrc`.
    gammaincc = _get_native_func(xp, spx, 'gammaincc')  # noqa: F811
    if gammaincc is None:
        return None

    def __chdtrc(v, x):
        res = xp.where(x >= 0, gammaincc(v/2, x/2), 1)
        i_nan = ((x == 0) & (v == 0)) | xp.isnan(x) | xp.isnan(v) | (v <= 0)
        res = xp.where(i_nan, xp.nan, res)
        return res
    return __chdtrc


def _betaincc(xp, spx):
    betainc = _get_native_func(xp, spx, 'betainc')  # noqa: F811
    if betainc is None:
        return None

    def __betaincc(a, b, x):
        # not perfect; might want to just rely on SciPy
        return betainc(b, a, 1-x)
    return __betaincc


def _stdtr(xp, spx):
    betainc = _get_native_func(xp, spx, 'betainc')  # noqa: F811
    if betainc is None:
        return None

    def __stdtr(df, t):
        x = df / (t ** 2 + df)
        tail = betainc(df / 2, 0.5, x) / 2
        return xp.where(t < 0, tail, 1 - tail)

    return __stdtr


def _stdtrit(xp, spx):
    betainc = _get_native_func(xp, spx, 'betainc')  # noqa: F811
    # If betainc is not defined, the root-finding would be done with `xp`
    # despite `stdtr` being evaluated with SciPy/NumPy `stdtr`. Save the
    # conversions: in this case, just evaluate `stdtrit` with SciPy/NumPy.
    if betainc is None:
        return None

    from scipy.optimize.elementwise import bracket_root, find_root

    def __stdtrit(df, p):
        def fun(t, df, p):  return stdtr(df, t) - p
        res_bracket = bracket_root(fun, xp.zeros_like(p), args=(df, p))
        res_root = find_root(fun, res_bracket.bracket, args=(df, p))
        return res_root.x

    return __stdtrit


_generic_implementations = {'rel_entr': _rel_entr,
                            'xlogy': _xlogy,
                            'chdtr': _chdtr,
                            'chdtrc': _chdtrc,
                            'betaincc': _betaincc,
                            'stdtr': _stdtr,
                            'stdtrit': _stdtrit,
                            }


# functools.wraps doesn't work because:
# 'numpy.ufunc' object has no attribute '__module__'
def support_alternative_backends(f_name, f_module):
    func = getattr(f_module, f_name)

    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        xp = array_namespace(*args)
        f = get_array_special_func(f_name, f_module, xp)
        return f(*args, **kwargs)

    if hasattr(func, 'types'):
        # Some tests use the types attribute to generate test cases.
        wrapped.types = func.types  # type: ignore

    return wrapped


# function name: number of args (for testing purposes)
@dataclass
class FunctionInfo:
    n_args: int
    module: ModuleType = _ufuncs


array_special_func_map = {
    'logit': FunctionInfo(1),
    'expit': FunctionInfo(1),
    'log_ndtr': FunctionInfo(1),
    'ndtr': FunctionInfo(1),
    'ndtri': FunctionInfo(1),
    'digamma': FunctionInfo(1, _basic),
    'polygamma': FunctionInfo(2, _basic),
    'multigammaln': FunctionInfo(2, _spfun_stats),
    'gammaln': FunctionInfo(1),
    'gamma': FunctionInfo(1),
    'gammainc': FunctionInfo(2),
    'gammaincc': FunctionInfo(2),
    'betaln': FunctionInfo(2),
    'beta': FunctionInfo(2),
    'betainc': FunctionInfo(3),
    'erf': FunctionInfo(1),
    'erfc': FunctionInfo(1),
    'erfinv': FunctionInfo(1),
    'zeta': FunctionInfo(2, _basic),
    'binom': FunctionInfo(2),
    'expi': FunctionInfo(1),
    'expn': FunctionInfo(2),
    'i0': FunctionInfo(1),
    'i0e': FunctionInfo(1),
    'i1': FunctionInfo(1),
    'i1e': FunctionInfo(1),
    'entr': FunctionInfo(1),
    'rel_entr': FunctionInfo(2),
    'xlogy': FunctionInfo(2),
    'chdtr': FunctionInfo(2),
    'chdtrc': FunctionInfo(2),
    'betaincc': FunctionInfo(3),
    'stdtr': FunctionInfo(2),
    'stdtrit': FunctionInfo(2),
}

globals().update(
    {f_name: support_alternative_backends(f_name, info.module)
     if SCIPY_ARRAY_API
     else getattr(info.module, f_name)
     for f_name, info in array_special_func_map.items()}
)
