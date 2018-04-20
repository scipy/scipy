from __future__ import division, print_function, absolute_import

import warnings

from . import _zeros
from numpy import finfo, isfinite, abs

_iter = 100
_xtol = 2e-12
_rtol = 4 * finfo(float).eps

__all__ = ['newton', 'sidi', 'bisect', 'ridder', 'brentq', 'brenth']

CONVERGED = 'converged'
SIGNERR = 'sign error'
CONVERR = 'convergence error'
_ECONVERGED = 0
_ESIGNERR = -1
_ECONVERR = -2

flag_map = {_ECONVERGED: CONVERGED, _ESIGNERR: SIGNERR, _ECONVERR: CONVERR}


class RootResults(object):
    """ Represents the root finding result.
    Attributes
    ----------
    root : float
        Estimated root location.
    iterations : int
        Number of iterations needed to find the root.
    function_calls : int
        Number of times the function was called.
    converged : bool
        True if the routine converged.
    flag : str
        Description of the cause of termination.
    """

    def __init__(self, root, iterations, function_calls, flag):
        self.root = root
        self.iterations = iterations
        self.function_calls = function_calls
        self.converged = flag == 0
        try:
            self.flag = flag_map[flag]
        except KeyError:
            self.flag = 'unknown error %d' % (flag,)

    def __repr__(self):
        attrs = ['converged', 'flag', 'function_calls',
                 'iterations', 'root']
        m = max(map(len, attrs)) + 1
        return '\n'.join([a.rjust(m) + ': ' + repr(getattr(self, a))
                          for a in attrs])


def results_c(full_output, r):
    if full_output:
        x, funcalls, iterations, flag = r
        results = RootResults(root=x,
                              iterations=iterations,
                              function_calls=funcalls,
                              flag=flag)
        return x, results
    else:
        return r


def _results_select(full_output, r):
    r"""Select from a tuple of (root, funccalls, iterations, flag)"""
    x, funcalls, iterations, flag = r
    if full_output:
        results = RootResults(root=x,
                              iterations=iterations,
                              function_calls=funcalls,
                              flag=flag)
        return x, results
    return x


def _withinTolerance(x, y, atol, rtol):
    diff = abs(x - y)
    z = abs(y)
    result = (diff <= (atol + rtol * z))
    return result


# Newton-Raphson method
def newton(func, x0, fprime=None, args=(), tol=1.48e-8, maxiter=50,
           fprime2=None,
           full_output=False, disp=True):
    """
    Find a zero using the Newton-Raphson or secant method.

    Find a zero of the function `func` given a nearby starting point `x0`.
    The Newton-Raphson method is used if the derivative `fprime` of `func`
    is provided, otherwise the secant method is used.  If the second order
    derivative `fprime2` of `func` is provided, then Halley's method is used.

    Parameters
    ----------
    func : function
        The function whose zero is wanted. It must be a function of a
        single variable of the form f(x,a,b,c...), where a,b,c... are extra
        arguments that can be passed in the `args` parameter.
    x0 : float
        An initial estimate of the zero that should be somewhere near the
        actual zero.
    fprime : function, optional
        The derivative of the function when available and convenient. If it
        is None (default), then the secant method is used.
    args : tuple, optional
        Extra arguments to be used in the function call.
    tol : float, optional
        The allowable error of the zero value.
    maxiter : int, optional
        Maximum number of iterations.
    fprime2 : function, optional
        The second order derivative of the function when available and
        convenient. If it is None (default), then the normal Newton-Raphson
        or the secant method is used. If it is not None, then Halley's method
        is used.

    Returns
    -------
    zero : float
        Estimated location where function is zero.

    See Also
    --------
    brentq, brenth, ridder, bisect
    fsolve : find zeroes in n dimensions.

    Notes
    -----
    The convergence rate of the Newton-Raphson method is quadratic,
    the Halley method is cubic, and the secant method is
    sub-quadratic.  This means that if the function is well behaved
    the actual error in the estimated zero is approximately the square
    (cube for Halley) of the requested tolerance up to roundoff
    error. However, the stopping criterion used here is the step size
    and there is no guarantee that a zero has been found. Consequently
    the result should be verified. Safer algorithms are brentq,
    brenth, ridder, and bisect, but they all require that the root
    first be bracketed in an interval where the function changes
    sign. The brentq algorithm is recommended for general use in one
    dimensional problems when such an interval has been found.

    Examples
    --------

    >>> def f(x):
    ...     return (x**3 - 1)  # only one real root at x = 1
    
    >>> from scipy import optimize

    ``fprime`` not provided, use secant method
    
    >>> root = optimize.newton(f, 1.5)
    >>> root
    1.0000000000000016
    >>> root = optimize.newton(f, 1.5, fprime2=lambda x: 6 * x)
    >>> root
    1.0000000000000016

    Only ``fprime`` provided, use Newton Raphson method
    
    >>> root = optimize.newton(f, 1.5, fprime=lambda x: 3 * x**2)
    >>> root
    1.0
    
    Both ``fprime2`` and ``fprime`` provided, use Halley's method

    >>> root = optimize.newton(f, 1.5, fprime=lambda x: 3 * x**2,
    ...                        fprime2=lambda x: 6 * x)
    >>> root
    1.0

    """
    if tol <= 0:
        raise ValueError("tol too small (%g <= 0)" % tol)
    if maxiter < 1:
        raise ValueError("maxiter must be greater than 0")
    # Multiply by 1.0 to convert to floating point.  We don't use float(x0)
    # so it still works if x0 is complex.
    p0 = 1.0 * x0
    funcalls = 0
    if fprime is not None:
        # Newton-Raphson method
        for itr in range(maxiter):
            fder = fprime(p0, *args)
            funcalls += 1
            if fder == 0:
                if disp:
                    msg = "derivative was zero."
                    warnings.warn(msg, RuntimeWarning)
                return _results_select(full_output, (p0, funcalls, itr + 1, _ECONVERR))
            fval = func(p0, *args)
            funcalls += 1
            newton_step = fval / fder
            if fprime2 is None:
                # Newton step
                p = p0 - newton_step
            else:
                fder2 = fprime2(p0, *args)
                # Halley's method
                p = p0 - newton_step / (1.0 - 0.5 * newton_step * fder2 / fder)
            if abs(p - p0) < tol:
                return _results_select(full_output, (p, funcalls, itr + 1, _ECONVERGED))
            p0 = p
    else:
        # Secant method
        if x0 >= 0:
            p1 = x0 * (1 + 1e-4) + 1e-4
        else:
            p1 = x0 * (1 + 1e-4) - 1e-4
        q0 = func(p0, *args)
        funcalls += 1
        q1 = func(p1, *args)
        funcalls += 1
        for itr in range(maxiter):
            if q1 == q0:
                if p1 != p0:
                    if disp:
                        msg = "Tolerance of %s reached" % (p1 - p0)
                        warnings.warn(msg, RuntimeWarning)
                    p = (p1 + p0) / 2.0
                    return _results_select(full_output, (p, funcalls, itr + 1, _ECONVERGED))
            else:
                p = p1 - q1 * (p1 - p0) / (q1 - q0)
            if abs(p - p1) < tol:
                return _results_select(full_output, (p, funcalls, itr + 1, _ECONVERGED))
            p0 = p1
            q0 = q1
            p1 = p
            q1 = func(p1, *args)
            funcalls += 1
    if disp:
        msg = "Failed to converge after %d iterations, value is %s" % (itr + 1, p)
        raise RuntimeError(msg)

    return _results_select(full_output, (p, funcalls, itr + 1, _ECONVERR))


# Helpers for Sidi's method
def _updateDividedDiffs(xx, divdiffs, x, fx, k):
    r'''Update the divided differences "matrix"'''
    # Instead of keeping a full triangular array of f values and
    # divided differences, just keep the trailing diagonal
    # I.e. [f7, f67, f567]  ->
    #    [f8, f78=(f7-f8)/(x7-x8), f678=(f67-f78)/(x6-x8)]
    # See diagram on [Sidi2008, p118]
    xx.insert(0, x)
    divdiffs.insert(0, fx)
    # for i in range(1, min(k + 2, len(xx))):
    for i in range(1, len(xx)):
        divdiffs[i] = (divdiffs[i] - divdiffs[i - 1]) / (xx[i] - xx[0])
    if len(xx) > k + 1:
        xx.pop(-1)
        divdiffs.pop(-1)


def _updateBracket(brkt, fbrkt, x, fx):
    r"""Update [x0, x1], [f0, f1] with the knowledge that func(x) = fx"""
    if not brkt:
        return
    if brkt[0] < x < brkt[1]:
        repidx = (0 if (fx < 0) == (fbrkt[0] < 0) else 1)
        brkt[repidx] = x
        fbrkt[repidx] = fx


# Sidi's method
def sidi(f, a, b, args=(), k=2,
         xtol=_xtol, rtol=_rtol, maxiter=_iter, safe=True,
         full_output=False, disp=True):
    """
    Find a zero using Sidi's method.

    Uses Sidi's method [Sidi2008] to find a zero of the function `f` on
    the sign changing interval [a, b].  The method is a generalization
    of the secant method.  Newton's interpolation formula constructs
    a polynomial p_k(x) through k+1 points.  The derivative of p_k(x0)
    is used as a replacement for f'(x0) in the Newton-Raphson formula.
    Divided differences enable an efficient computation of p_k'(x).

    Parameters
    ----------
    f : function
        Python function returning a number.  The function :math:`f`
        must be continuous, and :math:`f(a)` and :math:`f(b)`
        have opposite signs.
    a : number
        One end of the bracketing interval :math:`[a, b]`.
    b : number
        The other end of the bracketing interval :math:`[a, b]`.
    k : number, optional
        The degree of the interpolating polynomial, k>=1.  The default is 2.
        k=1 is the Secant method.  Higher degrees use more divided
        difference terms.
    xtol : number, optional
        The computed root ``x0`` will satisfy ``np.allclose(x, x0,
        atol=xtol, rtol=rtol)``, where ``x`` is the exact root. The
        parameter must be nonnegative.
    rtol : number, optional
        The computed root ``x0`` will satisfy ``np.allclose(x, x0,
        atol=xtol, rtol=rtol)``, where ``x`` is the exact root.
    maxiter : number, optional
        if convergence is not achieved in maxiter iterations, an error is
        raised.  Must be >= 0.
    args : tuple, optional
        containing extra arguments for the function `f`.
        `f` is called by ``apply(f, (x)+args)``.
    full_output : bool, optional
        If `full_output` is False, the root is returned.  If `full_output` is
        True, the return value is ``(x, r)``, where `x` is the root, and `r` is
        a RootResults object.
    disp : bool, optional
        If True, raise RuntimeError if the algorithm didn't converge.
    safe : bool, optional
        If True, whenever an intermediate derivative is 0, a bisection
        step is performsd.
        If False and [a, b] is a bracketing interval, iterates outside
        the interval report a convergence error.
        If False and [a, b] is not a bracketing interval, a and b are
        just the first two iterates, and iterations proceed until
        either convergence, or the maximum number of iterations is exceeded.

    Returns
    -------
    x0 : float
        Zero of `f` between `a` and `b`.
    r : RootResults (present if ``full_output = True``)
        Object containing information about the convergence.  In particular,
        ``r.converged`` is True if the routine converged.

    See Also
    --------
    brentq, brenth, ridder, bisect, newton
    fsolve : find zeroes in n dimensions.

    Notes
    -----
    `f` must be continuous.
    If safe is True, f(a) and f(b) must have opposite signs.
    Otherwise [a, b] is treated as a bracketing interval iff f(a)
    and f(b) have opposite signs.
    If f has at least (k+1) continuous derivatives, then Sidi's method
    of degree k has approximate order of convergence s_k, where
    s_1, s_2, s_3, ... = 1.618, 1.839, 1.928, ...
    k=1 is just the secant method.

    Examples
    --------

    >>> def f(x):
    ...     return (x**3 - 1)  # only one real root at x = 1

    >>> from scipy import optimize
    >>> root, results = optimize.sidi(f, 0, 2, full_output=True)
    >>> root
    1.0
    >>> results
          converged: True
               flag: 'converged'
     function_calls: 11
         iterations: 9
               root: 1.0


    References
    ----------
    .. [Sidi2008]
       Sidi, Avram,
       *Generalization of the secant method for nonlinear equations.*.
       Applied Mathematics E-Notes,
       http://eudml.org/doc/55912,

    """
    if xtol <= 0:
        raise ValueError("xtol too small (%g <= 0)" % xtol)
    if rtol < _rtol:
        raise ValueError("rtol too small (%g < %g)" % (rtol, _rtol))
    if maxiter < 1:
        raise ValueError("maxiter must be greater than 0")
    if not isfinite(a):
        raise ValueError("a is not finite" % a)
    if not isfinite(b):
        raise ValueError("b is not finite" % b)
    if k < 1:
        raise ValueError("k too small (%d < %d" % (k, 1))
    if b <= a and safe:
        raise ValueError("Invalid bracket %s" % [a, b])

    if not isinstance(args, tuple):
        args = (args,)

    funcalls = 0

    fa = f(a, *args)
    funcalls += 1
    if fa == 0:
        return _results_select(full_output, (a, funcalls, 0, _ECONVERGED))
    fb = f(b, *args)
    funcalls += 1
    if fb == 0:
        return _results_select(full_output, (b, funcalls, 0, _ECONVERGED))

    # Ensure starting from the endpoint with the smaller of the two absolute fvalues
    if abs(fa) > abs(fb):
        xvals, divdiffs = [a], [fa]
        _updateDividedDiffs(xvals, divdiffs, b, fb, k)
    else:
        xvals, divdiffs = [b], [fb]
        _updateDividedDiffs(xvals, divdiffs, a, fa, k)

    # Is this a bracketing interval?
    if fa * fb > 0:
        if safe:
            return _results_select(full_output, (a, funcalls, 0, _ESIGNERR))
        brkt, fbrk= [], []
    else:
        brkt, fbrk = [a, b], [fa, fb]

    for itr in range(maxiter):
        denom = divdiffs[1]
        m = 1.0
        for i in range(1, min(k + 1, len(xvals) - 1)):
            m *= xvals[0] - xvals[i]
            denom += divdiffs[i + 1] * m
        if denom == 0:
            if not safe:
                if disp:
                    msg = "Denominator is zero during iteration %d, value is %s (%s)" % (
                        itr, xvals[0], divdiffs)
                    raise RuntimeError(msg)
                return _results_select(full_output, (xvals[0], funcalls, itr, _ECONVERR))
            xn = sum(brkt) / 2.0
        else:
            delta = divdiffs[0] / denom
            assert isfinite(delta)
            xn = xvals[0] - delta
            # If the movement was too small to register...
            if delta != 0 and xn == xvals[0]:
                return _results_select(full_output, (xn, funcalls, itr + 1, _ECONVERGED))

        # If outside bracket
        if brkt and not (brkt[0] <= xn <= brkt[1]):
            if safe:
                # Disallow iterates outside of [a, b]
                # Disallow two consecutive iterates to be outside the brkt
                if not (a <= xn <= b) or not (brkt[0] <= xvals[0] <= brkt[1]):
                    # reset
                    xvals, divdiffs = [brkt[0]], [fbrk[0]]
                    _updateDividedDiffs(xvals, divdiffs, brkt[1], fbrk[1], k)
                    #  use mid-point of brkt
                    xn = sum(brkt[:2]) / 2.0
            elif not (a <= xn <= b):
                # [a, b] is a containing bracket, we've gone outside
                if disp:
                    msg = ("Iterate %f outside bracket during iteration %d, "
                           "last value is %s") % (xn, itr, xvals[0])
                    raise RuntimeError(msg)
                break

        # Ensure the next value is not already in the list
        if safe and xn in xvals:
            xn = sum(brkt) / 2.0

        fxn = f(xn, *args)
        funcalls += 1
        if fxn == 0:
            return _results_select(full_output, (xn, funcalls, itr + 1, _ECONVERGED))

        _updateDividedDiffs(xvals, divdiffs, xn, fxn, k)
        _updateBracket(brkt, fbrk, xn, fxn)

        if brkt:
            if (brkt[0] <= xn <= brkt[1]) and _withinTolerance(brkt[0], brkt[1], xtol, rtol):
                return _results_select(full_output, (xn, funcalls, itr + 1, _ECONVERGED))
        else:
            if _withinTolerance(xn, xvals[1], xtol, rtol):
                return _results_select(full_output, (xn, funcalls, itr + 1, _ECONVERGED))

    if disp:
        msg = "Failed to converge after %d iterations, value is %s" % (
            itr + 1, xvals[0])
        raise RuntimeError(msg)

    return _results_select(full_output, (xvals[0], funcalls, itr + 1, _ECONVERR))


def bisect(f, a, b, args=(),
           xtol=_xtol, rtol=_rtol, maxiter=_iter,
           full_output=False, disp=True):
    """
    Find root of a function within an interval.

    Basic bisection routine to find a zero of the function `f` between the
    arguments `a` and `b`. `f(a)` and `f(b)` cannot have the same signs.
    Slow but sure.

    Parameters
    ----------
    f : function
        Python function returning a number.  `f` must be continuous, and
        f(a) and f(b) must have opposite signs.
    a : number
        One end of the bracketing interval [a,b].
    b : number
        The other end of the bracketing interval [a,b].
    xtol : number, optional
        The computed root ``x0`` will satisfy ``np.allclose(x, x0,
        atol=xtol, rtol=rtol)``, where ``x`` is the exact root. The
        parameter must be nonnegative.
    rtol : number, optional
        The computed root ``x0`` will satisfy ``np.allclose(x, x0,
        atol=xtol, rtol=rtol)``, where ``x`` is the exact root. The
        parameter cannot be smaller than its default value of
        ``4*np.finfo(float).eps``.
    maxiter : number, optional
        if convergence is not achieved in `maxiter` iterations, an error is
        raised.  Must be >= 0.
    args : tuple, optional
        containing extra arguments for the function `f`.
        `f` is called by ``apply(f, (x)+args)``.
    full_output : bool, optional
        If `full_output` is False, the root is returned.  If `full_output` is
        True, the return value is ``(x, r)``, where x is the root, and r is
        a `RootResults` object.
    disp : bool, optional
        If True, raise RuntimeError if the algorithm didn't converge.

    Returns
    -------
    x0 : float
        Zero of `f` between `a` and `b`.
    r : RootResults (present if ``full_output = True``)
        Object containing information about the convergence.  In particular,
        ``r.converged`` is True if the routine converged.

    Examples
    --------

    >>> def f(x):
    ...     return (x**2 - 1)

    >>> from scipy import optimize

    >>> root = optimize.bisect(f, 0, 2)
    >>> root
    1.0

    >>> root = optimize.bisect(f, -2, 0)
    >>> root
    -1.0

    See Also
    --------
    brentq, brenth, bisect, newton
    fixed_point : scalar fixed-point finder
    fsolve : n-dimensional root-finding

    """
    if not isinstance(args, tuple):
        args = (args,)
    if xtol <= 0:
        raise ValueError("xtol too small (%g <= 0)" % xtol)
    if rtol < _rtol:
        raise ValueError("rtol too small (%g < %g)" % (rtol, _rtol))
    r = _zeros._bisect(f, a, b, xtol, rtol, maxiter, args, full_output, disp)
    return results_c(full_output, r)


def ridder(f, a, b, args=(),
           xtol=_xtol, rtol=_rtol, maxiter=_iter,
           full_output=False, disp=True):
    """
    Find a root of a function in an interval.

    Parameters
    ----------
    f : function
        Python function returning a number.  f must be continuous, and f(a) and
        f(b) must have opposite signs.
    a : number
        One end of the bracketing interval [a,b].
    b : number
        The other end of the bracketing interval [a,b].
    xtol : number, optional
        The computed root ``x0`` will satisfy ``np.allclose(x, x0,
        atol=xtol, rtol=rtol)``, where ``x`` is the exact root. The
        parameter must be nonnegative.
    rtol : number, optional
        The computed root ``x0`` will satisfy ``np.allclose(x, x0,
        atol=xtol, rtol=rtol)``, where ``x`` is the exact root. The
        parameter cannot be smaller than its default value of
        ``4*np.finfo(float).eps``.
    maxiter : number, optional
        if convergence is not achieved in maxiter iterations, an error is
        raised.  Must be >= 0.
    args : tuple, optional
        containing extra arguments for the function `f`.
        `f` is called by ``apply(f, (x)+args)``.
    full_output : bool, optional
        If `full_output` is False, the root is returned.  If `full_output` is
        True, the return value is ``(x, r)``, where `x` is the root, and `r` is
        a RootResults object.
    disp : bool, optional
        If True, raise RuntimeError if the algorithm didn't converge.

    Returns
    -------
    x0 : float
        Zero of `f` between `a` and `b`.
    r : RootResults (present if ``full_output = True``)
        Object containing information about the convergence.
        In particular, ``r.converged`` is True if the routine converged.

    See Also
    --------
    brentq, brenth, bisect, newton : one-dimensional root-finding
    fixed_point : scalar fixed-point finder

    Notes
    -----
    Uses [Ridders1979]_ method to find a zero of the function `f` between the
    arguments `a` and `b`. Ridders' method is faster than bisection, but not
    generally as fast as the Brent routines. [Ridders1979]_ provides the
    classic description and source of the algorithm. A description can also be
    found in any recent edition of Numerical Recipes.

    The routine used here diverges slightly from standard presentations in
    order to be a bit more careful of tolerance.

    Examples
    --------

    >>> def f(x):
    ...     return (x**2 - 1)

    >>> from scipy import optimize

    >>> root = optimize.ridder(f, 0, 2)
    >>> root
    1.0

    >>> root = optimize.ridder(f, -2, 0)
    >>> root
    -1.0

    References
    ----------
    .. [Ridders1979]
       Ridders, C. F. J. "A New Algorithm for Computing a
       Single Root of a Real Continuous Function."
       IEEE Trans. Circuits Systems 26, 979-980, 1979.

    """
    if not isinstance(args, tuple):
        args = (args,)
    if xtol <= 0:
        raise ValueError("xtol too small (%g <= 0)" % xtol)
    if rtol < _rtol:
        raise ValueError("rtol too small (%g < %g)" % (rtol, _rtol))
    r = _zeros._ridder(f, a, b, xtol, rtol, maxiter, args, full_output, disp)
    return results_c(full_output, r)


def brentq(f, a, b, args=(),
           xtol=_xtol, rtol=_rtol, maxiter=_iter,
           full_output=False, disp=True):
    """
    Find a root of a function in a bracketing interval using Brent's method.

    Uses the classic Brent's method to find a zero of the function `f` on
    the sign changing interval [a , b].  Generally considered the best of the
    rootfinding routines here.  It is a safe version of the secant method that
    uses inverse quadratic extrapolation.  Brent's method combines root
    bracketing, interval bisection, and inverse quadratic interpolation.  It is
    sometimes known as the van Wijngaarden-Dekker-Brent method.  Brent (1973)
    claims convergence is guaranteed for functions computable within [a,b].

    [Brent1973]_ provides the classic description of the algorithm.  Another
    description can be found in a recent edition of Numerical Recipes, including
    [PressEtal1992]_.  Another description is at
    http://mathworld.wolfram.com/BrentsMethod.html.  It should be easy to
    understand the algorithm just by reading our code.  Our code diverges a bit
    from standard presentations: we choose a different formula for the
    extrapolation step.

    Parameters
    ----------
    f : function
        Python function returning a number.  The function :math:`f`
        must be continuous, and :math:`f(a)` and :math:`f(b)` must
        have opposite signs.
    a : number
        One end of the bracketing interval :math:`[a, b]`.
    b : number
        The other end of the bracketing interval :math:`[a, b]`.
    xtol : number, optional
        The computed root ``x0`` will satisfy ``np.allclose(x, x0,
        atol=xtol, rtol=rtol)``, where ``x`` is the exact root. The
        parameter must be nonnegative. For nice functions, Brent's
        method will often satisfy the above condition with ``xtol/2``
        and ``rtol/2``. [Brent1973]_
    rtol : number, optional
        The computed root ``x0`` will satisfy ``np.allclose(x, x0,
        atol=xtol, rtol=rtol)``, where ``x`` is the exact root. The
        parameter cannot be smaller than its default value of
        ``4*np.finfo(float).eps``. For nice functions, Brent's
        method will often satisfy the above condition with ``xtol/2``
        and ``rtol/2``. [Brent1973]_
    maxiter : number, optional
        if convergence is not achieved in maxiter iterations, an error is
        raised.  Must be >= 0.
    args : tuple, optional
        containing extra arguments for the function `f`.
        `f` is called by ``apply(f, (x)+args)``.
    full_output : bool, optional
        If `full_output` is False, the root is returned.  If `full_output` is
        True, the return value is ``(x, r)``, where `x` is the root, and `r` is
        a RootResults object.
    disp : bool, optional
        If True, raise RuntimeError if the algorithm didn't converge.

    Returns
    -------
    x0 : float
        Zero of `f` between `a` and `b`.
    r : RootResults (present if ``full_output = True``)
        Object containing information about the convergence.  In particular,
        ``r.converged`` is True if the routine converged.

    See Also
    --------
    multivariate local optimizers
      `fmin`, `fmin_powell`, `fmin_cg`, `fmin_bfgs`, `fmin_ncg`
    nonlinear least squares minimizer
      `leastsq`
    constrained multivariate optimizers
      `fmin_l_bfgs_b`, `fmin_tnc`, `fmin_cobyla`
    global optimizers
      `basinhopping`, `brute`, `differential_evolution`
    local scalar minimizers
      `fminbound`, `brent`, `golden`, `bracket`
    n-dimensional root-finding
      `fsolve`
    one-dimensional root-finding
      `brenth`, `ridder`, `bisect`, `newton`
    scalar fixed-point finder
      `fixed_point`

    Notes
    -----
    `f` must be continuous.  f(a) and f(b) must have opposite signs.

    Examples
    --------
    >>> def f(x):
    ...     return (x**2 - 1)

    >>> from scipy import optimize

    >>> root = optimize.brentq(f, -2, 0)
    >>> root
    -1.0

    >>> root = optimize.brentq(f, 0, 2)
    >>> root
    1.0

    References
    ----------
    .. [Brent1973]
       Brent, R. P.,
       *Algorithms for Minimization Without Derivatives*.
       Englewood Cliffs, NJ: Prentice-Hall, 1973. Ch. 3-4.

    .. [PressEtal1992]
       Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T.
       *Numerical Recipes in FORTRAN: The Art of Scientific Computing*, 2nd ed.
       Cambridge, England: Cambridge University Press, pp. 352-355, 1992.
       Section 9.3:  "Van Wijngaarden-Dekker-Brent Method."

    """
    if not isinstance(args, tuple):
        args = (args,)
    if xtol <= 0:
        raise ValueError("xtol too small (%g <= 0)" % xtol)
    if rtol < _rtol:
        raise ValueError("rtol too small (%g < %g)" % (rtol, _rtol))
    r = _zeros._brentq(f, a, b, xtol, rtol, maxiter, args, full_output, disp)
    return results_c(full_output, r)


def brenth(f, a, b, args=(),
           xtol=_xtol, rtol=_rtol, maxiter=_iter,
           full_output=False, disp=True):
    """Find root of f in [a,b].

    A variation on the classic Brent routine to find a zero of the function f
    between the arguments a and b that uses hyperbolic extrapolation instead of
    inverse quadratic extrapolation. There was a paper back in the 1980's ...
    f(a) and f(b) cannot have the same signs. Generally on a par with the
    brent routine, but not as heavily tested.  It is a safe version of the
    secant method that uses hyperbolic extrapolation. The version here is by
    Chuck Harris.

    Parameters
    ----------
    f : function
        Python function returning a number.  f must be continuous, and f(a) and
        f(b) must have opposite signs.
    a : number
        One end of the bracketing interval [a,b].
    b : number
        The other end of the bracketing interval [a,b].
    xtol : number, optional
        The computed root ``x0`` will satisfy ``np.allclose(x, x0,
        atol=xtol, rtol=rtol)``, where ``x`` is the exact root. The
        parameter must be nonnegative. As with `brentq`, for nice
        functions the method will often satisfy the above condition
        with ``xtol/2`` and ``rtol/2``.
    rtol : number, optional
        The computed root ``x0`` will satisfy ``np.allclose(x, x0,
        atol=xtol, rtol=rtol)``, where ``x`` is the exact root. The
        parameter cannot be smaller than its default value of
        ``4*np.finfo(float).eps``. As with `brentq`, for nice functions
        the method will often satisfy the above condition with
        ``xtol/2`` and ``rtol/2``.
    maxiter : number, optional
        if convergence is not achieved in maxiter iterations, an error is
        raised.  Must be >= 0.
    args : tuple, optional
        containing extra arguments for the function `f`.
        `f` is called by ``apply(f, (x)+args)``.
    full_output : bool, optional
        If `full_output` is False, the root is returned.  If `full_output` is
        True, the return value is ``(x, r)``, where `x` is the root, and `r` is
        a RootResults object.
    disp : bool, optional
        If True, raise RuntimeError if the algorithm didn't converge.

    Returns
    -------
    x0 : float
        Zero of `f` between `a` and `b`.
    r : RootResults (present if ``full_output = True``)
        Object containing information about the convergence.  In particular,
        ``r.converged`` is True if the routine converged.

    Examples
    --------
    >>> def f(x):
    ...     return (x**2 - 1)

    >>> from scipy import optimize

    >>> root = optimize.brenth(f, -2, 0)
    >>> root
    -1.0

    >>> root = optimize.brenth(f, 0, 2)
    >>> root
    1.0

    See Also
    --------
    fmin, fmin_powell, fmin_cg,
           fmin_bfgs, fmin_ncg : multivariate local optimizers

    leastsq : nonlinear least squares minimizer

    fmin_l_bfgs_b, fmin_tnc, fmin_cobyla : constrained multivariate optimizers

    basinhopping, differential_evolution, brute : global optimizers

    fminbound, brent, golden, bracket : local scalar minimizers

    fsolve : n-dimensional root-finding

    brentq, brenth, ridder, bisect, newton : one-dimensional root-finding

    fixed_point : scalar fixed-point finder

    """
    if not isinstance(args, tuple):
        args = (args,)
    if xtol <= 0:
        raise ValueError("xtol too small (%g <= 0)" % xtol)
    if rtol < _rtol:
        raise ValueError("rtol too small (%g < %g)" % (rtol, _rtol))
    r = _zeros._brenth(f, a, b, xtol, rtol, maxiter, args, full_output, disp)
    return results_c(full_output, r)
