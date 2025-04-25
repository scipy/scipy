"""Partial replacements for numpy polynomial routines, with Array API compatibility.

This module contains both "old-style", np.poly1d, routines from the main numpy
namespace, and "new-style", np.polynomial.polynomial, routines.

To distinguish the two sets, the "new-style" routine names start with `npp_`
"""
import scipy._lib.array_api_extra as xpx


def polyroots(coef, *, xp):
    """numpy.roots, best-effor replacement
    """
    if coef.shape[0] < 2:
        return xp.asarray([], dtype=coef.dtype)

    root_func = getattr(xp, 'roots', None)
    if root_func:
        # OK on numpy and jax, is broken on cupy
        return root_func(coef)

    # companion matrix
    n = coef.shape[0]
    a = xp.eye(n - 1, n - 1, k=-1)
    a[:, -1] = -xp.flip(coef[1:]) / coef[0]

    # non-symmetric eigenvalue problem is not in the spec but is available on e.g. torch
    if hasattr(xp.linalg, 'eigvals'):
        return xp.linalg.eigvals(a)
    else:
        import numpy as np
        return xp.asarray(np.linalg.eigvals(np.asarray(a)))


# ### Old-style routines ###


# https://github.com/numpy/numpy/blob/v2.2.0/numpy/lib/_polynomial_impl.py#L702-L779
def polyval(p, x, *, xp):
    """ Old-style polynomial, `np.polyval`
    """
    p = xp.asarray(p)
    x = xp.asarray(x)
    y = xp.zeros_like(x)

    for pv in p:
        y = y * x + pv
    return y


# ### New-style routines ###


# https://github.com/numpy/numpy/blob/v2.2.0/numpy/polynomial/polynomial.py#L758-L842
def npp_polyvalfromroots(x, r, *, xp, tensor=True):
    r = xp.asarray(r, copy=None)
    r = xpx.atleast_nd(r, ndim=1, xp=xp)
    # if r.dtype.char in '?bBhHiIlLqQpP':
    #    r = r.astype(np.double)

    if isinstance(x, (tuple, list)):
        x = np.asarray(x)

    if tensor:
        r = xp.reshape(r, r.shape + (1,) * x.ndim)
    elif x.ndim >= r.ndim:
        raise ValueError("x.ndim must be < r.ndim when tensor == False")
    return xp.prod(x - r, axis=0)
