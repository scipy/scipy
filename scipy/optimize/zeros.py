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
           full_output=False, disp=False):
    """Find root of f in [a,b]

    Basic bisection routine to find a zero of the function
    f between the arguments a and b. f(a) and f(b) can not
    have the same signs. Slow but sure.

    f : Python function returning a number.

    a : Number, one end of the bracketing interval.

    b : Number, the other end of the bracketing interval.

    xtol : Number, the routine converges when a root is known
        to lie within xtol of the value return. Should be >= 0.
        The routine modifies this to take into account the relative
        precision of doubles.

    maxiter : Number, if convergence is not achieved in
        maxiter iterations, and error is raised. Must be
        >= 0.

    args : tuple containing extra arguments for the function f.
        f is called by apply(f,(x)+args).

    If full_output is False, the root is returned.

    If full_output is True, the return value is (x, r), where x
    is the root, and r is a RootResults object containing information
    about the convergence. In particular, r.converged is True if the
    the routine converged.

    See also:

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
    r = _zeros._bisect(f,a,b,xtol,maxiter,args,full_output,disp)
    return results_c(full_output, r)

def ridder(f, a, b, args=(),
           xtol=_xtol, rtol=_rtol, maxiter=_iter,
           full_output=False, disp=False):
    """Find root of f in [a,b]

    Ridder routine to find a zero of the function
    f between the arguments a and b. f(a) and f(b) can not
    have the same signs. Faster than bisection, but not
    generaly as fast as the brent rountines. A description
    may be found in a recent edition of Numerical Recipes.
    The routine here is modified a bit to be more careful
    of tolerance.

    f : Python function returning a number.

    a : Number, one end of the bracketing interval.

    b : Number, the other end of the bracketing interval.

    xtol : Number, the routine converges when a root is known
        to lie within xtol of the value return. Should be >= 0.
        The routine modifies this to take into account the relative
        precision of doubles.

    maxiter : Number, if convergence is not achieved in
        maxiter iterations, and error is raised. Must be
        >= 0.

    args : tuple containing extra arguments for the function f.
        f is called by apply(f,(x)+args).

    If full_output is False, the root is returned.

    If full_output is True, the return value is (x, r), where x
    is the root, and r is a RootResults object containing information
    about the convergence. In particular, r.converged is True if the
    the routine converged.

    See also:

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
    r = _zeros._ridder(f,a,b,xtol,maxiter,args,full_output,disp)
    return results_c(full_output, r)

def brentq(f, a, b, args=(),
           xtol=_xtol, rtol=_rtol, maxiter=_iter,
           full_output=False, disp=False):
    """Find root of f in [a,b]

    The classic Brent routine to find a zero of the function
    f between the arguments a and b. f(a) and f(b) can not
    have the same signs. Generally the best of the routines here.
    It is a safe version of the secant method that uses inverse
    quadratic extrapolation. The version here is a slight
    modification that uses a different formula in the extrapolation
    step. A description may be found in Numerical Recipes, but the
    code here is probably easier to understand.

    f : Python function returning a number.

    a : Number, one end of the bracketing interval.

    b : Number, the other end of the bracketing interval.

    xtol : Number, the routine converges when a root is known
        to lie within xtol of the value return. Should be >= 0.
        The routine modifies this to take into account the relative
        precision of doubles.

    maxiter : Number, if convergence is not achieved in
        maxiter iterations, and error is raised. Must be
        >= 0.

    args : tuple containing extra arguments for the function f.
        f is called by apply(f,(x)+args).

    If full_output is False, the root is returned.

    If full_output is True, the return value is (x, r), where x
    is the root, and r is a RootResults object containing information
    about the convergence. In particular, r.converged is True if the
    the routine converged.

    See also:

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
    r = _zeros._brentq(f,a,b,xtol,maxiter,args,full_output,disp)
    return results_c(full_output, r)

def brenth(f, a, b, args=(),
           xtol=_xtol, rtol=_rtol, maxiter=_iter,
           full_output=False, disp=False):
    """Find root of f in [a,b]

    A variation on the classic Brent routine to find a zero
    of the function f between the arguments a and b that uses
    hyperbolic extrapolation instead of inverse quadratic
    extrapolation. There was a paper back in the 1980's ...
    f(a) and f(b) can not have the same signs. Generally on a
    par with the brent routine, but not as heavily tested.
    It is a safe version of the secant method that uses hyperbolic
    extrapolation. The version here is by Chuck Harris.

    f : Python function returning a number.

    a : Number, one end of the bracketing interval.

    b : Number, the other end of the bracketing interval.

    xtol : Number, the routine converges when a root is known
        to lie within xtol of the value return. Should be >= 0.
        The routine modifies this to take into account the relative
        precision of doubles.

    maxiter : Number, if convergence is not achieved in
        maxiter iterations, and error is raised. Must be
        >= 0.

    args : tuple containing extra arguments for the function f.
        f is called by apply(f,(x)+args).

    If full_output is False, the root is returned.

    If full_output is True, the return value is (x, r), where x
    is the root, and r is a RootResults object containing information
    about the convergence. In particular, r.converged is True if the
    the routine converged.

    See also:

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
    r = _zeros._brenth(f,a,b,xtol,maxiter,args,full_output,disp)
    return results_c(full_output, r)
