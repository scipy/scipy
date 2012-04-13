"""
Unified interfaces to root finding algorithms.

Functions
---------
- root : find a root of a vector function.
"""

__all__ = ['root']

from optimize import MemoizeJac
from minpack import _root_hybr, leastsq

def root(fun, x0, args=(), method='hybr', jac=None, options=None,
         full_output=False):
    """
    Find a root of a vector function.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    fun : callable
        A vector function to find a root of.
    x0 : ndarray
        Initial guess.
    args : tuple, optional
        Extra arguments passed to the objective function and its Jacobian.
    method : str, optional
        Type of solver.  Should be one of: {'hybr', 'lm'}
    jac : bool or callable, optional
        If `jac` is a Boolean and is True, `fun` is assumed to return the
        value of Jacobian along with the objective function. If False, the
        Jacobian will be estimated numerically.
        `jac` can also be a callable returning the Jacobian of the
        objective. In this case, it must accept the same arguments as
        `fun`.
    options : dict, optional
        A dictionary of solver options. All methods accept the following
        generic options:
            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.
        For method-specific options, see `show_minimize_options`.
    full_output : bool, optional
        If True, return optional outputs.  Default is False.

    Returns
    -------
    x : ndarray
        The solution.
    info : dict
        A dictionary of optional outputs (depending on the chosen method)
        with the keys:
            solution : ndarray
                The solution (same as `x`).
            success : bool
                Boolean flag indicating if a solution was found.
            status : int
                An integer flag indicating the type of termination.  Its
                value depends on the underlying solver.  Refer to `message`
                for more information.
            message : str
                A string message giving information about the cause of the
                termination.
            fun, jac : ndarray
                Values of objective function and Jacobian (if available).
            nfev, njev: int
                Number of evaluations of the objective functions and of its
                jacobian.
            nit: int
                Number of iterations.

    Notes
    -----
    This section describes the available solvers that can be selected by the
    'method' parameter. The default method is *hybr*.

    Method *hybr* uses a modification of the Powell hybrid method as
    implemented in MINPACK [1]_.

    Method *lm* solves the system of nonlinear equations in a least squares
    sense using a modification of the Levenberg-Marquardt algorithm as
    implemented in MINPACK [1]_.

    References
    ----------
    .. [1] More, Jorge J., Burton S. Garbow, and Kenneth E. Hillstrom.
       1980. User Guide for MINPACK-1.

    Examples
    --------
    The following functions define a system of nonlinear equations and its
    jacobian.
    >>> def fun(x):
    ...     return [x[0]  + 0.5 * (x[0] - x[1])**3 - 1.0,
    ...             0.5 * (x[1] - x[0])**3 + x[1]]
    >>> def jac(x):
    ...     return np.array([[1 + 1.5 * (x[0] - x[1])**2,
    ...                       -1.5 * (x[0] - x[1])**2],
    ...                      [-1.5 * (x[1] - x[0])**2,
    ...                       1 + 1.5 * (x[1] - x[0])**2]])

    A solution can be obtained as follows.
    >>> x, info = root(fun, [0, 0], jac=jac, method='hybr')
    >>> x
    array([ 0.8411639,  0.1588361])
    """
    meth = method.lower()
    if options is None:
        options = {}

    # fun also returns the jacobian
    if not callable(jac):
        if bool(jac):
            fun = MemoizeJac(fun)
            jac = fun.derivative
        else:
            jac = None

    if meth == 'hybr':
        out = _root_hybr(fun, x0, args=args, jac=jac, options=options,
                          full_output=full_output)
        if full_output:
            x, info = out
        else:
            x = out
    elif meth == 'lm':
        col_deriv = options.get('col_deriv', 0)
        xtol      = options.get('xtol', 1.49012e-08)
        ftol      = options.get('ftol', 1.49012e-08)
        gtol      = options.get('gtol', 0.0)
        maxfev    = options.get('maxfev', 0)
        epsfcn    = options.get('epsfcn', 0.0)
        factor    = options.get('factor', 100)
        diag      = options.get('diag', None)
        out = leastsq(fun, x0, args=args, Dfun=jac,
                      full_output=full_output, col_deriv=col_deriv,
                      xtol=xtol, ftol=ftol, gtol=gtol, maxfev=maxfev,
                      epsfcn=epsfcn, factor=factor, diag=diag)
        if full_output:
            x, cov_x, infodict, mesg, ier = out
            info = infodict
            info['message'] = mesg
            info['status'] = ier
            info['success'] = ier in (1, 2, 3, 4)
            info['cov_x'] = cov_x
        else:
            x = out[0]
    else:
        raise ValueError('Unknown solver %s' % method)

    if full_output:
        return x, info
    else:
        return x

