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
                iprint=1, maxfun=1000):
    """
    Minimize a function using the Constrained Optimization BY Linear
    Approximation (COBYLA) method.

    Parameters
    ----------
    func : callable f(x, *args)
        Function to minimize.
    x0 : ndarray
        Initial guess.
    cons : sequence
        Constraint functions; must all be ``>=0`` (a single function
        if only 1 constraint).
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
        Controls the frequency of output; 0 implies no output.
    maxfun : int
        Maximum number of function evaluations.

    Returns
    -------
    x : ndarray
        The argument that minimises `f`.

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
