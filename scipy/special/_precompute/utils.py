from __future__ import division, print_function, absolute_import
import os
import numpy as np
from numpy.testing import dec, assert_

try:
    import mpmath as mp
except ImportError:
    try:
        import sympy.mpmath as mp
    except ImportError:
        pass

try:
    from sympy.abc import x
except ImportError:
    pass


def skip():
    msg = "Set environment variable SCIPY_PRECOMPUTE=1 to run precompute tests."

    def deco(func):
        try:
            if bool(os.environ['SCIPY_PRECOMPUTE']):
                return func
        except (ValueError, KeyError):
            pass
        return dec.skipif(True, msg)(func)
    return deco


def mpf_assert_allclose(res, std, atol=0, rtol=1e-17):
    n = len(std)
    if len(res) != n:
        raise AssertionError("Lengths of inputs not equal.")

    failures = []
    for k in range(n):
        try:
            assert_(mp.fabs(res[k] - std[k]) <= atol + rtol*mp.fabs(std[k]))
        except AssertionError:
            failures.append(k)

    ndigits = int(abs(np.log10(rtol)))
    msg = [""]
    msg.append("Bad results ({} out of {}) for the following points:"
               .format(len(failures), n))
    for k in failures:
        resrep = mp.nstr(res[k], ndigits, min_fixed=0, max_fixed=0)
        stdrep = mp.nstr(std[k], ndigits, min_fixed=0, max_fixed=0)
        if std[k] == 0:
            rdiff = "inf"
        else:
            rdiff = mp.fabs((res[k] - std[k])/std[k])
            rdiff = mp.nstr(rdiff, 3)
        msg.append("{}: {} != {} (rdiff {})".format(k, resrep, stdrep, rdiff))
    if failures:
        assert_(False, "\n".join(msg))


def lagrange_inversion(a):
    """Given a series

    f(x) = a[1]*x + a[2]*x**2 + ... + a[n-1]*x**(n - 1),

    use the Lagrange inversion formula to compute a series

    g(x) = b[1]*x + b[2]*x**2 + ... + b[n-1]*x**(n - 1)

    so that f(g(x)) = g(f(x)) = x mod x**n. We must have a[0] = 0, so
    necessarily b[0] = 0 too.

    The algorithm is naive and could be improved, but speed isn't an
    issue here and it's easy to read.

    """
    n = len(a)
    f = sum(a[i]*x**i for i in range(len(a)))
    h = (x/f).series(x, 0, n).removeO()
    hpower = [h**0]
    for k in range(n):
        hpower.append((hpower[-1]*h).expand())
    b = [mp.mpf(0)]
    for k in range(1, n):
        b.append(hpower[k].coeff(x, k - 1)/k)
    b = map(lambda x: mp.mpf(x), b)
    return b
