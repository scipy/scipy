"""Interface to Constrained Optimization By Linear Approximation

Functions:
fmin_coblya(func, x0, cons, args=(), consargs=None, rhobeg=1.0, rhoend=1e-4,
            iprint=1, maxfun=1000)
    Minimize a function using the Constrained Optimization BY Linear
    Approximation (COBYLA) method

"""

import _cobyla
from numpy import copy

__all__ = ['fmin_cobyla']


def fmin_cobyla(func, x0, cons, args=(), consargs=None, rhobeg=1.0, rhoend=1e-4,
                iprint=1, maxfun=1000, disp=None):
    """
    Minimize a function using the Constrained Optimization BY Linear
    Approximation (COBYLA) method. This method wraps a FORTRAN
    implentation of the algorithm.

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
        Final accuracy in the optimization (not precisely guaranteed). This
        is a lower bound on the size of the trust region.
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

    Notes
    -----
    This algorithm is based on linear approximations to the objective
    function and each constraint. We briefly describe the algorithm.

    Suppose the function is being minimized over k variables. At the
    jth iteration the algorithm has k+1 points v_1, ..., v_(k+1),
    an approximate solution x_j, and a radius RHO_j.
    (i.e. linear plus a constant) approximations to the objective
    function and constraint functions such that their function values
    agree with the linear approximation on the k+1 points v_1,.., v_(k+1).
    This gives a linear program to solve (where the linear approximations
    of the constraint functions are constrained to be non-negative).

    However the linear approximations are likely only good
    approximations near the current simplex, so the linear program is
    given the further requirement that the solution, which
    will become x_(j+1), must be within RHO_j from x_j. RHO_j only
    decreases, never increases. The initial RHO_j is rhobeg and the
    final RHO_j is rhoend. In this way COBYLA's iterations behave
    like a trust region algorithm.

    Additionally, the linear program may be inconsistent, or the
    approximation may give poor improvement. For details about
    how these issues are resolved, as well as how the points v_i are
    updated, refer to the source code or the references below.


    References
    ----------
    Powell M.J.D. (1994), "A direct search optimization method that models
    the objective and constraint functions by linear interpolation.", in
    Advances in Optimization and Numerical Analysis, eds. S. Gomez and
    J-P Hennart, Kluwer Academic (Dordrecht), pp. 51-67

    Powell M.J.D. (1998), "Direct search algorithms for optimization
    calculations", Acta Numerica 7, 287-336

    Powell M.J.D. (2007), "A view of algorithms for optimization without
    derivatives", Cambridge University Technical Report DAMTP 2007/NA03


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
