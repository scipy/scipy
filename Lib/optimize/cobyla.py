# Interface to Constrained Optimization By Linear Approximation

from __future__ import nested_scopes
import _cobyla

def fmin_cobyla(func, x0, cons, args=(), consargs=None, rhobeg=1.0, rhoend=1e-4,
                iprint=1, maxfun=1000):
    """
    Minimize a function using the Contrained Optimization BY Linear Approximation
    (COBYLA) method

    Arguments:

    func     -- function to minimize. Called as func(x, *args)

    x0       -- initial guess to minimum

    cons     -- a list of functions that all must be >=0 (a single function
                if only 1 constraint)

    args     -- extra arguments to pass to function

    consargs -- extra arguments to pass to constraints (default of None means
                use same extra arguments as those passed to func).  Use () for no
                extra arguments.

    rhobeg --  reasonable initial changes to the variables

    rhoend --  final accuracy in the optimization (not precisely guaranteed)

    iprint  -- controls the frequency of output: 0 (no output),1,2,3

    maxfun  -- maximum number of function evaluations.


    Returns:

    x -- the minimum
    
    """
    err = "cons must be a list of callable functions or a single"\
              " callable function."
    n = len(x0)
    if isinstance(cons, list):
        m = len(cons)
        for thisfunc in cons:
            if not callable(thisfunc):
                raise TypeError, err
    elif callable(cons):
        m = 1
        cons = [cons]
    else:
        raise TypeError, "cons must be a list of callable functions or a single"\
              " callable function."

    if consargs is None:
        consargs = args
        
    def calcfc(x, con):
        f = func(x, *args)
        k = 0
        for constraints in cons:
            con[k] = constraints(x,*consargs)
            k += 1
        return f

    xopt = _cobyla.minimize(calcfc, m=m, x=x0,rhobeg=rhobeg,rhoend=rhoend,iprint=iprint,
                            maxfun=maxfun)
    
    return xopt
