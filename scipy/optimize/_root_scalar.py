"""
Unified interfaces to root finding algorithms for real or complex
scalar functions.

Functions
---------
- root : find a root of a scalar function.
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy._lib.six import callable

from . import zeros as optzeros

__all__ = ['root_scalar']


class MemoizeDer(object):
    """Decorator that caches the value and derivative(s) of function each
    time it is called.

    This is a simplistic memoizer that calls and caches a single value
    of `f(x, *args)`.
    It assumes that `args` does not change between invocations.
    It supports the use case of a root-finder where `args` is fixed,
    `x` changes, and only rarely, if at all, does x assume the same value
    more than once."""
    def __init__(self, fun):
        self.fun = fun
        self.vals = None
        self.x = None
        self.n_calls = 0

    def __call__(self, x, *args):
        r"""Calculate f or used cached value if available"""
        # Derivative may be requested before the function itself, always check
        if self.vals is None or x != self.x:
            fg = self.fun(x, *args)
            self.x = x
            self.n_calls += 1
            self.vals = fg[:]
        return self.vals[0]

    def fprime(self, x, *args):
        r"""Calculate f' or used a cached value if available"""
        if self.vals is None or x != self.x:
            self(x, *args)
        return self.vals[1]

    def fprime2(self, x, *args):
        r"""Calculate f'' or used a cached value if available"""
        if self.vals is None or x != self.x:
            self(x, *args)
        return self.vals[2]

    def ncalls(self):
        return self.n_calls


def root_scalar(f, args=(), method=None, bracket=None,
                fprime=None, fprime2=None,
                x0=None, x1=None,
                xtol=None, rtol=None, maxiter=None,
                options=None):
    """
    Find a root of a scalar function.

    Parameters
    ----------
    f : callable
        A function to find a root of.
    args : tuple, optional
        Extra arguments passed to the objective function and its derivative(s).
    method : str, optional
        Type of solver.  Should be one of

            - 'bisect'    :ref:`(see here) <optimize.root_scalar-bisect>`
            - 'brentq'    :ref:`(see here) <optimize.root_scalar-brentq>`
            - 'brenth'    :ref:`(see here) <optimize.root_scalar-brenth>`
            - 'ridder'    :ref:`(see here) <optimize.root_scalar-ridder>`
            - 'newton'    :ref:`(see here) <optimize.root_scalar-newton>`

    bracket: A sequence of 2 floats, optional
        An interval bracketing a root.  `f(x, *args)` must have different
        signs at the two endpoints.
    x0 : float, optional
        Initial guess.
    x1 : float, optional
        A second guess.
    fprime : bool or callable, optional
        If `fprime` is a boolean and is True, `f` is assumed to return the
        value of derivative along with the objective function.
        `fprime` can also be a callable returning the derivative of `f`. In
        this case, it must accept the same arguments as `f`.
    fprime2 : bool or callable, optional
        If `fprime2` is a boolean and is True, `f` is assumed to return the
        value of 1st and 2nd derivatives along with the objective function.
        `fprime2` can also be a callable returning the 2nd derivative of `f`.
        In this case, it must accept the same arguments as `f`.
    xtol : float, optional
        Tolerance (absolute) for termination.
    rtol : float, optional
        Tolerance (relative) for termination.
    maxiter : int, optional
        Maximum number of iterations.
    options : dict, optional
        A dictionary of solver options. E.g. `k` or `maxiter`, see
        :obj:`show_options()` for details.

    Returns
    -------
    sol : RootResults
        The solution represented as a ``RootResults`` object.
        Important attributes are: ``root`` the solution , ``converged`` a
        boolean flag indicating if the algorithm exited successfully and
        ``flag`` which describes the cause of the termination. See
        `RootResults` for a description of other attributes.

    See also
    --------
    show_options : Additional options accepted by the solvers

    Examples
    --------

    Find the root of a simple cubic

    >>> from scipy import optimize
    >>> def f(x):
    ...     return (x**3 - 1)  # only one real root at x = 1

    >>> def fprime(x):
    ...     return 3*x**2

    brentq takes as input a bracket

    >>> sol = optimize.root_scalar(f, bracket=[0, 3], method='brentq')
    >>> sol.root, sol.iterations, sol.function_calls
    (1.0, 10, 11)

    newton starts takes as input a single point and uses the derivative(s)

    >>> sol = optimize.root_scalar(f, x0=0.2, fprime=fprime, method='newton')
    >>> sol.root, sol.iterations, sol.function_calls
    (1.0, 11, 22)

    The function can provide the value and derivative in a single call.

    >>> def f_p_pp(x):
    ...     return (x**3 - 1), 3*x**2, 6*x

    >>> sol = optimize.root_scalar(f_p_pp, x0=0.2, fprime=True, method='newton')
    >>> sol.root, sol.iterations, sol.function_calls
    (1.0, 11, 11)

    >>> sol = optimize.root_scalar(f_p_pp, x0=0.2, fprime=True, fprime2=True, method='newton')
    >>> sol.root, sol.iterations, sol.function_calls
    (1.0, 6, 6)


    """
    if not isinstance(args, tuple):
        args = (args,)

    if options is None:
        options = {}

    # fun also returns the derivative(s)
    is_memoized = False
    if fprime2 is not None and not callable(fprime2):
        if bool(fprime2):
            f = MemoizeDer(f)
            is_memoized = True
            fprime2 = f.fprime2
            fprime = f.fprime
        else:
            fprime2 = None
    if fprime is not None and not callable(fprime):
        if bool(fprime):
            f = MemoizeDer(f)
            is_memoized = True
            fprime = f.fprime
        else:
            fprime = None

    # respect solver-specific default tolerances - only pass in if actually set
    kwargs = {}
    for k in ['xtol', 'rtol', 'maxiter']:
        v = locals().get(k)
        if v is not None:
            kwargs[k] = v

    # Set any solver-specific options
    if options:
        kwargs.update(options)
    # Always want full_output
    kwargs.update(full_output=True, disp=False)

    # Pick a method if not specified
    if not method:
        if bracket:
            method = 'brentq'
        elif x0 is not None:
            method = 'newton'
        else:
            raise ValueError('Unable to select a solver as neither bracket '
                             'nor starting point provided.')
    meth = method.lower()

    try:
        methodc = getattr(optzeros, meth)
    except AttributeError:
        raise ValueError('Unknown solver %s' % meth)

    if meth in ['bisect', 'ridder', 'brentq', 'brenth']:
        if not isinstance(bracket, (list, tuple, np.ndarray)):
            raise ValueError('Bracket needed for %s' % method)

        a, b = bracket[:2]
        r, sol = methodc(f, a, b, args=args, **kwargs)
    elif meth == 'newton':
        if x0 is None:
            raise ValueError('x0 musxt not be None for %s' % method)
        if 'xtol' in kwargs:
            kwargs['tol'] = kwargs.pop('xtol')
        r, sol = methodc(f, x0, args=args, fprime=fprime, fprime2=fprime2,
                         x1=x1, **kwargs)
    else:
        raise ValueError('Unknown solver %s' % method)

    if is_memoized:
        # Replace the function_calls count with the memoized count.
        n_calls = f.n_calls
        sol.function_calls = n_calls

    return sol
