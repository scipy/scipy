__all__ = ['splrep', 'splprep', 'splev', 'splint', 'sproot', 'spalde',
           'bisplrep', 'bisplev', 'insert', 'splder', 'splantider']


import numpy as np

# These are in the API for fitpack even if not used in fitpack.py itself.
from ._fitpack_impl import bisplrep, bisplev, dblint  # noqa: F401
from . import _fitpack_impl as _impl
from ._bsplines import BSpline


def splprep(x, w=None, u=None, ub=None, ue=None, k=3, task=0, s=None, t=None,
            full_output=0, nest=None, per=0, quiet=1):
    res = _impl.splprep(x, w, u, ub, ue, k, task, s, t, full_output, nest, per,
                        quiet)
    return res
splprep.__doc__ = _impl.splprep.__doc__


def splrep(x, y, w=None, xb=None, xe=None, k=3, task=0, s=None, t=None,
           full_output=0, per=0, quiet=1):
    res = _impl.splrep(x, y, w, xb, xe, k, task, s, t, full_output, per, quiet)
    return res
splrep.__doc__ = _impl.splrep.__doc__


def splev(x, tck, der=0, ext=0):
    if isinstance(tck, BSpline):
        if tck.c.ndim > 1:
            mesg = ("Calling splev() with BSpline objects with c.ndim > 1 is "
                    "not allowed. Use BSpline.__call__(x) instead.")
            raise ValueError(mesg)

        # remap the out-of-bounds behavior
        try:
            extrapolate = {0: True, }[ext]
        except KeyError as e:
            raise ValueError("Extrapolation mode %s is not supported "
                             "by BSpline." % ext) from e

        return tck(x, der, extrapolate=extrapolate)
    else:
        return _impl.splev(x, tck, der, ext)
splev.__doc__ = _impl.splev.__doc__


def splint(a, b, tck, full_output=0):
    if isinstance(tck, BSpline):
        if tck.c.ndim > 1:
            mesg = ("Calling splint() with BSpline objects with c.ndim > 1 is "
                    "not allowed. Use BSpline.integrate() instead.")
            raise ValueError(mesg)

        if full_output != 0:
            mesg = ("full_output = %s is not supported. Proceeding as if "
                    "full_output = 0" % full_output)

        return tck.integrate(a, b, extrapolate=False)
    else:
        return _impl.splint(a, b, tck, full_output)
splint.__doc__ = _impl.splint.__doc__


def sproot(tck, mest=10):
    if isinstance(tck, BSpline):
        if tck.c.ndim > 1:
            mesg = ("Calling sproot() with BSpline objects with c.ndim > 1 is "
                    "not allowed.")
            raise ValueError(mesg)

        t, c, k = tck.tck

        # _impl.sproot expects the interpolation axis to be last, so roll it.
        # NB: This transpose is a no-op if c is 1D.
        sh = tuple(range(c.ndim))
        c = c.transpose(sh[1:] + (0,))
        return _impl.sproot((t, c, k), mest)
    else:
        return _impl.sproot(tck, mest)
sproot.__doc__ = _impl.sproot.__doc__


def spalde(x, tck):
    if isinstance(tck, BSpline):
        raise TypeError("spalde does not accept BSpline instances.")
    else:
        return _impl.spalde(x, tck)
spalde.__doc__ = _impl.spalde.__doc__


def insert(x, tck, m=1, per=0):
    if isinstance(tck, BSpline):

        t, c, k = tck.tck

        # FITPACK expects the interpolation axis to be last, so roll it over
        # NB: if c array is 1D, transposes are no-ops
        sh = tuple(range(c.ndim))
        c = c.transpose(sh[1:] + (0,))
        t_, c_, k_ = _impl.insert(x, (t, c, k), m, per)

        # and roll the last axis back
        c_ = np.asarray(c_)
        c_ = c_.transpose((sh[-1],) + sh[:-1])
        return BSpline(t_, c_, k_)
    else:
        return _impl.insert(x, tck, m, per)
insert.__doc__ = _impl.insert.__doc__


def splder(tck, n=1):
    if isinstance(tck, BSpline):
        return tck.derivative(n)
    else:
        return _impl.splder(tck, n)
splder.__doc__ = _impl.splder.__doc__


def splantider(tck, n=1):
    if isinstance(tck, BSpline):
        return tck.antiderivative(n)
    else:
        return _impl.splantider(tck, n)
splantider.__doc__ = _impl.splantider.__doc__
