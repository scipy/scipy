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
    Approximation (COBYLA) method

    Arguments:

    func     -- function to minimize. Called as func(x, *args)

    x0       -- initial guess to minimum

    cons     -- a sequence of functions that all must be >=0 (a single function
                if only 1 constraint)

    args     -- extra arguments to pass to function

    consargs -- extra arguments to pass to constraints (default of None means
                use same extra arguments as those passed to func).
                Use () for no extra arguments.

    rhobeg --  reasonable initial changes to the variables

    rhoend --  final accuracy in the optimization (not precisely guaranteed)

    iprint  -- controls the frequency of output: 0 (no output),1,2,3

    maxfun  -- maximum number of function evaluations.


    Returns:

    x -- the minimum

    See also:

        scikits.openopt, which offers a unified syntax to call this and other solvers

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

    xopt = _cobyla.minimize(calcfc, m=m, x=copy(x0), rhobeg=rhobeg, rhoend=rhoend,
                            iprint=iprint, maxfun=maxfun)

    return xopt
