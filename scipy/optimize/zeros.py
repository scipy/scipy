## Automatically adapted for scipy Oct 07, 2005 by convertcode.py

import _zeros
from numpy import finfo

_iter = 100
_xtol = 1e-12
# not actually used at the moment
_rtol = finfo(float).eps * 2

__all__ = ['bisect','ridder','brentq','brenth']

CONVERGED = 'converged'
SIGNERR = 'sign error'
CONVERR = 'convergence error'
flag_map = {0 : CONVERGED, -1 : SIGNERR, -2 : CONVERR}

class RootResults(object):
    def __init__(self, root, iterations, function_calls, flag):
        self.root = root
        self.iterations = iterations
        self.function_calls = function_calls
        self.converged = flag == 0
        try:
            self.flag = flag_map[flag]
        except KeyError:
            self.flag = 'unknown error %d' % (flag,)

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

def bisect(f, a, b, args=(),
           xtol=_xtol, rtol=_rtol, maxiter=_iter,
           full_output=False, disp=True):
    """Find root of f in [a,b].

    Basic bisection routine to find a zero of the function f between the
    arguments a and b. f(a) and f(b) can not have the same signs. Slow but
    sure.

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
        The routine converges when a root is known to lie within xtol of the
        value return. Should be >= 0.  The routine modifies this to take into
        account the relative precision of doubles.
    maxiter : number, optional
        if convergence is not achieved in maxiter iterations, and error is
        raised.  Must be >= 0.
    args : tuple, optional
        containing extra arguments for the function `f`.
        `f` is called by ``apply(f, (x)+args)``.
    full_output : bool, optional
        If `full_output` is False, the root is returned.  If `full_output` is
        True, the return value is ``(x, r)``, where `x` is the root, and `r` is
        a RootResults object.
    disp : {True, bool} optional
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
        brentq, brenth, bisect, newton : one-dimensional root-finding
        fixed_point : scalar fixed-point finder
        fsolve -- n-dimenstional root-finding

    """
    if type(args) != type(()) :
        args = (args,)
    r = _zeros._bisect(f,a,b,xtol,maxiter,args,full_output,disp)
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
        The routine converges when a root is known to lie within xtol of the
        value return. Should be >= 0.  The routine modifies this to take into
        account the relative precision of doubles.
    maxiter : number, optional
        if convergence is not achieved in maxiter iterations, and error is
        raised.  Must be >= 0.
    args : tuple, optional
        containing extra arguments for the function `f`.
        `f` is called by ``apply(f, (x)+args)``.
    full_output : bool, optional
        If `full_output` is False, the root is returned.  If `full_output` is
        True, the return value is ``(x, r)``, where `x` is the root, and `r` is
        a RootResults object.
    disp : {True, bool} optional
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
    generally as fast as the Brent rountines. [Ridders1979]_ provides the
    classic description and source of the algorithm. A description can also be
    found in any recent edition of Numerical Recipes.

    The routine used here diverges slightly from standard presentations in
    order to be a bit more careful of tolerance.

    References
    ----------
    .. [Ridders1979]
       Ridders, C. F. J. "A New Algorithm for Computing a
       Single Root of a Real Continuous Function."
       IEEE Trans. Circuits Systems 26, 979-980, 1979.

    """
    if type(args) != type(()) :
        args = (args,)
    r = _zeros._ridder(f,a,b,xtol,maxiter,args,full_output,disp)
    return results_c(full_output, r)

def brentq(f, a, b, args=(),
           xtol=_xtol, rtol=_rtol, maxiter=_iter,
           full_output=False, disp=True):
    """
    Find a root of a function in given interval.

    Return float, a zero of `f` between `a` and `b`.  `f` must be a continuous
    function, and [a,b] must be a sign changing interval.

    Description:
    Uses the classic Brent (1973) method to find a zero of the function `f` on
    the sign changing interval [a , b].  Generally considered the best of the
    rootfinding routines here.  It is a safe version of the secant method that
    uses inverse quadratic extrapolation.  Brent's method combines root
    bracketing, interval bisection, and inverse quadratic interpolation.  It is
    sometimes known as the van Wijngaarden-Deker-Brent method.  Brent (1973)
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
        Python function returning a number.  f must be continuous, and f(a) and
        f(b) must have opposite signs.
    a : number
        One end of the bracketing interval [a,b].
    b : number
        The other end of the bracketing interval [a,b].
    xtol : number, optional
        The routine converges when a root is known to lie within xtol of the
        value return. Should be >= 0.  The routine modifies this to take into
        account the relative precision of doubles.
    maxiter : number, optional
        if convergence is not achieved in maxiter iterations, and error is
        raised.  Must be >= 0.
    args : tuple, optional
        containing extra arguments for the function `f`.
        `f` is called by ``apply(f, (x)+args)``.
    full_output : bool, optional
        If `full_output` is False, the root is returned.  If `full_output` is
        True, the return value is ``(x, r)``, where `x` is the root, and `r` is
        a RootResults object.
    disp : {True, bool} optional
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
      `anneal`, `brute`
    local scalar minimizers
      `fminbound`, `brent`, `golden`, `bracket`
    n-dimenstional root-finding
      `fsolve`
    one-dimensional root-finding
      `brentq`, `brenth`, `ridder`, `bisect`, `newton`
    scalar fixed-point finder
      `fixed_point`

    Notes
    -----

    f must be continuous.  f(a) and f(b) must have opposite signs.


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
    if type(args) != type(()) :
        args = (args,)
    r = _zeros._brentq(f,a,b,xtol,maxiter,args,full_output,disp)
    return results_c(full_output, r)

def brenth(f, a, b, args=(),
           xtol=_xtol, rtol=_rtol, maxiter=_iter,
           full_output=False, disp=True):
    """Find root of f in [a,b].

    A variation on the classic Brent routine to find a zero of the function f
    between the arguments a and b that uses hyperbolic extrapolation instead of
    inverse quadratic extrapolation. There was a paper back in the 1980's ...
    f(a) and f(b) can not have the same signs. Generally on a par with the
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
        The routine converges when a root is known to lie within xtol of the
        value return. Should be >= 0.  The routine modifies this to take into
        account the relative precision of doubles.
    maxiter : number, optional
        if convergence is not achieved in maxiter iterations, and error is
        raised.  Must be >= 0.
    args : tuple, optional
        containing extra arguments for the function `f`.
        `f` is called by ``apply(f, (x)+args)``.
    full_output : bool, optional
        If `full_output` is False, the root is returned.  If `full_output` is
        True, the return value is ``(x, r)``, where `x` is the root, and `r` is
        a RootResults object.
    disp : {True, bool} optional
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
      fmin, fmin_powell, fmin_cg,
             fmin_bfgs, fmin_ncg -- multivariate local optimizers
      leastsq -- nonlinear least squares minimizer

      fmin_l_bfgs_b, fmin_tnc,
             fmin_cobyla -- constrained multivariate optimizers

      anneal, brute -- global optimizers

      fminbound, brent, golden, bracket -- local scalar minimizers

      fsolve -- n-dimenstional root-finding

      brentq, brenth, ridder, bisect, newton -- one-dimensional root-finding

      fixed_point -- scalar fixed-point finder

    """
    if type(args) != type(()) :
        args = (args,)
    r = _zeros._brenth(f,a, b, xtol, maxiter, args, full_output, disp)
    return results_c(full_output, r)
