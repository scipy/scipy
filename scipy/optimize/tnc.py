# TNC Python interface
# @(#) $Jeannot: tnc.py,v 1.11 2005/01/28 18:27:31 js Exp $

# Copyright (c) 2004-2005, Jean-Sebastien Roy (js@jeannot.org)

# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:

# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
TNC: A python interface to the TNC non-linear optimizer

TNC is a non-linear optimizer. To use it, you must provide a function to
minimize. The function must take one argument: the list of coordinates where to
evaluate the function; and it must return either a tuple, whose first element is the
value of the function, and whose second argument is the gradient of the function
(as a list of values); or None, to abort the minimization.
"""
from scipy.optimize import moduleTNC
from numpy import asarray, inf, array

__all__ = ['fmin_tnc']


MSG_NONE = 0 # No messages
MSG_ITER = 1 # One line per iteration
MSG_INFO = 2 # Informational messages
MSG_VERS = 4 # Version info
MSG_EXIT = 8 # Exit reasons
MSG_ALL = MSG_ITER + MSG_INFO + MSG_VERS + MSG_EXIT

MSGS = {
        MSG_NONE : "No messages",
        MSG_ITER : "One line per iteration",
        MSG_INFO : "Informational messages",
        MSG_VERS : "Version info",
        MSG_EXIT : "Exit reasons",
        MSG_ALL  : "All messages"
}

INFEASIBLE   = -1 # Infeasible (low > up)
LOCALMINIMUM =  0 # Local minima reach (|pg| ~= 0)
FCONVERGED   =  1 # Converged (|f_n-f_(n-1)| ~= 0)
XCONVERGED   =  2 # Converged (|x_n-x_(n-1)| ~= 0)
MAXFUN       =  3 # Max. number of function evaluations reach
LSFAIL       =  4 # Linear search failed
CONSTANT     =  5 # All lower bounds are equal to the upper bounds
NOPROGRESS   =  6 # Unable to progress
USERABORT    =  7 # User requested end of minimization

RCSTRINGS = {
        INFEASIBLE   : "Infeasible (low > up)",
        LOCALMINIMUM : "Local minima reach (|pg| ~= 0)",
        FCONVERGED   : "Converged (|f_n-f_(n-1)| ~= 0)",
        XCONVERGED   : "Converged (|x_n-x_(n-1)| ~= 0)",
        MAXFUN       : "Max. number of function evaluations reach",
        LSFAIL       : "Linear search failed",
        CONSTANT     : "All lower bounds are equal to the upper bounds",
        NOPROGRESS   : "Unable to progress",
        USERABORT    : "User requested end of minimization"
}

# Changes to interface made by Travis Oliphant, Apr. 2004 for inclusion in
#  SciPy

import optimize
approx_fprime = optimize.approx_fprime

def fmin_tnc(func, x0, fprime=None, args=(), approx_grad=0,
             bounds=None, epsilon=1e-8, scale=None, offset=None,
             messages=MSG_ALL, maxCGit=-1, maxfun=None, eta=-1,
             stepmx=0, accuracy=0, fmin=0, ftol=-1, xtol=-1, pgtol=-1,
             rescale=-1, disp=None):
    """
    Minimize a function with variables subject to bounds, using
    gradient information in a truncated Newton algorithm. This
    method wraps a C implementation of the algorithm.

    Parameters
    ----------
    func : callable ``func(x, *args)``
        Function to minimize.  Must do one of
        1. Return f and g, where f is
        the value of the function and g its gradient (a list of
        floats).
        2. Return the function value but supply gradient function
        seperately as fprime
        3. Return the function value and set approx_grad=True.
        If the function returns None, the minimization
        is aborted.
    x0 : list of floats
        Initial estimate of minimum.
    fprime : callable ``fprime(x, *args)``
        Gradient of func. If None, then either func must return the
        function value and the gradient (``f,g = func(x, *args)``)
        or approx_grad must be True.
    args : tuple
        Arguments to pass to function.
    approx_grad : bool
        If true, approximate the gradient numerically.
    bounds : list
        (min, max) pairs for each element in x0, defining the
        bounds on that parameter. Use None or +/-inf for one of
        min or max when there is no bound in that direction.
    epsilon: float
        Used if approx_grad is True. The stepsize in a finite
        difference approximation for fprime.
    scale : list of floats
        Scaling factors to apply to each variable.  If None, the
        factors are up-low for interval bounded variables and
        1+|x] fo the others.  Defaults to None
    offset : float
        Value to substract from each variable.  If None, the
        offsets are (up+low)/2 for interval bounded variables
        and x for the others.
    messages :
        Bit mask used to select messages display during
        minimization values defined in the MSGS dict.  Defaults to
        MGS_ALL.
    disp : int
        Integer interface to messages.  0 = no message, 5 = all messages
    maxCGit : int
        Maximum number of hessian*vector evaluations per main
        iteration.  If maxCGit == 0, the direction chosen is
        -gradient if maxCGit < 0, maxCGit is set to
        max(1,min(50,n/2)).  Defaults to -1.
    maxfun : int
        Maximum number of function evaluation.  if None, maxfun is
        set to max(100, 10*len(x0)).  Defaults to None.
    eta : float
        Severity of the line search. if < 0 or > 1, set to 0.25.
        Defaults to -1.
    stepmx : float
        Maximum step for the line search.  May be increased during
        call.  If too small, it will be set to 10.0.  Defaults to 0.
    accuracy : float
        Relative precision for finite difference calculations.  If
        <= machine_precision, set to sqrt(machine_precision).
        Defaults to 0.
    fmin : float
        Minimum function value estimate.  Defaults to 0.
    ftol : float
        Precision goal for the value of f in the stoping criterion.
        If ftol < 0.0, ftol is set to 0.0 defaults to -1.
    xtol : float
        Precision goal for the value of x in the stopping
        criterion (after applying x scaling factors).  If xtol <
        0.0, xtol is set to sqrt(machine_precision).  Defaults to
        -1.
    pgtol : float
        Precision goal for the value of the projected gradient in
        the stopping criterion (after applying x scaling factors).
        If pgtol < 0.0, pgtol is set to 1e-2 * sqrt(accuracy).
        Setting it to 0.0 is not recommended.  Defaults to -1.
    rescale : float
        Scaling factor (in log10) used to trigger f value
        rescaling.  If 0, rescale at each iteration.  If a large
        value, never rescale.  If < 0, rescale is set to 1.3.

    Returns
    -------
    x : list of floats
        The solution.
    nfeval : int
        The number of function evaluations.
    rc : int
        Return code as defined in the RCSTRINGS dict.


    Notes
    -----
    The underlying algorithm is truncated Newton, also called
    Newton Conjugate-Gradient. This method differs from
    scipy.optimize.fmin_ncg in that

    1. It wraps a C implementation of the algorithm
    2. It allows each variable to be given an upper and lower bound.


    The algorithm incoporates the bound constraints by determining
    the descent direction as in an unconstrained truncated Newton,
    but never taking a step-size large enough to leave the space
    of feasible x's. The algorithm keeps track of a set of
    currently active constraints, and ignores them when computing
    the minimum allowable step size. (The x's associated with the
    active constraint are kept fixed.) If the maximum allowable
    step size is zero then a new constraint is added. At the end
    of each iteration one of the constraints may be deemed no
    longer active and removed. A constraint is considered
    no longer active is if it is currently active
    but the gradient for that variable points inward from the
    constraint. The specific constraint removed is the one
    associated with the variable of largest index whose
    constraint is no longer active.


    References
    ----------
    Wright S., Nocedal J. (2006), 'Numerical Optimization'

    Nash S.G. (1984), "Newton-Type Minimization Via the Lanczos Method",
    SIAM Journal of Numerical Analysis 21, pp. 770-778


    """
    x0 = asarray(x0, dtype=float).tolist()
    n = len(x0)

    if bounds is None:
        bounds = [(None,None)] * n
    if len(bounds) != n:
        raise ValueError('length of x0 != length of bounds')

    if disp is not None:
        messages = {0:MSG_NONE, 1:MSG_ITER, 2:MSG_INFO, 3:MSG_VERS,
                    4:MSG_EXIT, 5:MSG_ALL}.get(disp, MSG_ALL)

    if approx_grad:
        def func_and_grad(x):
            x = asarray(x)
            f = func(x, *args)
            g = approx_fprime(x, func, epsilon, *args)
            return f, list(g)
    elif fprime is None:
        def func_and_grad(x):
            x = asarray(x)
            f, g = func(x, *args)
            return f, list(g)
    else:
        def func_and_grad(x):
            x = asarray(x)
            f = func(x, *args)
            g = fprime(x, *args)
            return f, list(g)

    """
    low, up   : the bounds (lists of floats)
                if low is None, the lower bounds are removed.
                if up is None, the upper bounds are removed.
                low and up defaults to None
    """
    low = [0]*n
    up = [0]*n
    for i in range(n):
        if bounds[i] is None: l, u = -inf, inf
        else:
            l,u = bounds[i]
            if l is None:
                low[i] = -inf
            else:
                low[i] = l
            if u is None:
                up[i] = inf
            else:
                up[i] = u

    if scale is None:
        scale = []

    if offset is None:
        offset = []

    if maxfun is None:
        maxfun = max(100, 10*len(x0))

    rc, nf, x = moduleTNC.minimize(func_and_grad, x0, low, up, scale, offset,
            messages, maxCGit, maxfun, eta, stepmx, accuracy,
            fmin, ftol, xtol, pgtol, rescale)
    return array(x), nf, rc

if __name__ == '__main__':
    # Examples for TNC

    def example():
        print "Example"
        # A function to minimize
        def function(x):
            f = pow(x[0],2.0)+pow(abs(x[1]),3.0)
            g = [0,0]
            g[0] = 2.0*x[0]
            g[1] = 3.0*pow(abs(x[1]),2.0)
            if x[1]<0:
                g[1] = -g[1]
            return f, g

        # Optimizer call
        x, nf, rc = fmin_tnc(function, [-7, 3], bounds=([-10, 1], [10, 10]))

        print "After", nf, "function evaluations, TNC returned:", RCSTRINGS[rc]
        print "x =", x
        print "exact value = [0, 1]"
        print

    example()
