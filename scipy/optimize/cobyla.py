"""Interface to Constrained Optimization By Linear Approximation

Functions:
fmin_coblya(func, x0, cons, args=(), consargs=None, rhobeg=1.0, rhoend=1e-4,
            iprint=1, maxfun=1000)
    Minimize a function using the Constrained Optimization BY Linear
    Approximation (COBYLA) method

"""

import _cobyla
from numpy import copy
def fmin_cobyla(func, x0, cons, args=(), consargs=None, rhobeg=1.0, rhoend=1e-4,
                iprint=1, maxfun=1000, disp=None):
    """
    Minimize a function using the Constrained Optimization BY Linear
    Approximation (COBYLA) method.

    Parameters
    ----------
    func : callable
        Function to minimize. In the form func(x, \\*args).
    x0 : ndarray
        Initial guess.
    cons : sequence
        Constraint functions; must all be ``>=0`` (a single function
        if only 1 constraint). Each function takes the parameters `x`
        as its first argument.
    args : tuple
        Extra arguments to pass to function.
    consargs : tuple
        Extra arguments to pass to constraint functions (default of None means
        use same extra arguments as those passed to func).
        Use ``()`` for no extra arguments.
    rhobeg :
        Reasonable initial changes to the variables.
    rhoend :
        Final accuracy in the optimization (not precisely guaranteed).
    iprint : {0, 1, 2, 3}
        Controls the frequency of output; 0 implies no output.  Deprecated.
    disp : {0, 1, 2, 3}
        Over-rides the iprint interface.  Preferred.
    maxfun : int
        Maximum number of function evaluations.

    Returns
    -------
    x : ndarray
        The argument that minimises `f`.

    Examples
    --------
    Minimize the objective function f(x,y) = x*y subject
    to the constraints x**2 + y**2 < 1 and y > 0::

        >>> def objective(x):
        ...     return x[0]*x[1]
        ...
        >>> def constr1(x):
        ...     return 1 - (x[0]**2 + x[1]**2)
        ...
        >>> def constr2(x):
        ...     return x[1]
        ...
        >>> fmin_cobyla(objective, [0.0, 0.1], [constr1, constr2], rhoend=1e-7)

           Normal return from subroutine COBYLA

           NFVALS =   64   F =-5.000000E-01    MAXCV = 1.998401E-14
           X =-7.071069E-01   7.071067E-01
        array([-0.70710685,  0.70710671])

    The exact solution is (-sqrt(2)/2, sqrt(2)/2).

    """
    err = "cons must be a sequence of callable functions or a single"\
          " callable function."
    try:
        m = len(cons)
    except TypeError:
        if callable(cons):
            m = 1
            cons = [cons]
        else:
            raise TypeError(err)
    else:
        for thisfunc in cons:
            if not callable(thisfunc):
                raise TypeError(err)

    if consargs is None:
        consargs = args

    if disp is not None:
        iprint = disp

    def calcfc(x, con):
        f = func(x, *args)
        k = 0
        for constraints in cons:
            con[k] = constraints(x, *consargs)
            k += 1
        return f

    xopt = _cobyla.minimize(calcfc, m=m, x=copy(x0), rhobeg=rhobeg,
                            rhoend=rhoend, iprint=iprint, maxfun=maxfun)

    return xopt
