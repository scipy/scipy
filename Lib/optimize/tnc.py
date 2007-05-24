## Automatically adapted for scipy Oct 07, 2005 by convertcode.py

# TNC Python interface
# @(#) $Jeannot: tnc.py,v 1.7 2004/04/02 20:40:21 js Exp $

# Copyright (c) 2004, Jean-Sebastien Roy (js@jeannot.org)

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
from numpy import asarray, inf

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

EINVAL       = -2 # Invalid parameters (n<1)
INFEASIBLE   = -1 # Infeasible (low > up)
LOCALMINIMUM =  0 # Local minima reach (|pg| ~= 0)
CONVERGED    =  1 # Converged (|f_n-f_(n-1)| ~= 0)
MAXFUN       =  2 # Max. number of function evaluations reach
LSFAIL       =  3 # Linear search failed
CONSTANT     =  4 # All lower bounds are equal to the upper bounds
NOPROGRESS   =  5 # Unable to progress
USERABORT    =  6 # User requested end of minimization

RCSTRINGS = {
        EINVAL       : "Invalid parameters (n<1)",
        INFEASIBLE   : "Infeasible (low > up)",
        LOCALMINIMUM : "Local minima reach (|pg| ~= 0)",
        CONVERGED    : "Converged (|f_n-f_(n-1)| ~= 0)",
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


def fmin_tnc(func, x0, fprime=None, args=(), approx_grad=0, bounds=None, epsilon=1e-8,
             scale=None, messages=MSG_ALL, maxCGit=-1, maxfun=None, eta=-1,
             stepmx=0, accuracy=0, fmin=0, ftol=0, rescale=-1):
    """Minimize a function with variables subject to bounds, using gradient
    information.

    :Parameters:

        func : callable func(x, *args)
            Function to minimize.
        x0 : float
            Initial guess to minimum.
        fprime : callable fprime(x, *args)
            Gradient of func. If None, then func must return the
            function value and the gradient, e.g.
            f,g = func(x,*args).
        args : tuple
            Arguments to pass to function.
        approx_grad : bool
            If true, approximate the gradient numerically.
        bounds : list
               (min, max) pairs for each element in x, defining the
               bounds on that parameter. Use None for one of min or
               max when there is no bound in that direction
        scale : list of floats
            Scaling factors to apply to each variable.  If None, the
            factors are up-low for interval bounded variables and
            1+|x] fo the others.  Defaults to None.
        messages :
            Bit mask used to select messages display during
            minimization values defined in the optimize.tnc.MSGS dict.
            defaults to optimize.tnc.MGS_ALL.
        maxCGit : int
            Maximum number of hessian*vector evaluation per main
            iteration.  If maxCGit == 0, the direction chosen is
            -gradient.  If maxCGit < 0, maxCGit is set to
            max(1,min(50,n/2)).  Defaults to -1.
        maxfun : int
            Maximum number of function evaluation.  If None, maxfun
            is set to max(1000, 100*len(x0)).  Defaults to None.
        eta : float
            Severity of the line search. if < 0 or > 1, set to 0.25.
            Defaults to -1.
        stepmx : float
            Maximum step for the line search.  May be increased during
            call.  If too small, will be set to 10.0.  Defaults to 0.
        accuracy : float
            Relative precision for finite difference calculations.  If
            <= machine_precision, set to sqrt(machine_precision).
            Defaults to 0.
        fmin : float
            Minimum function value estimate.  Defaults to 0.
        ftol : float
            Precision goal for the value of f in the stoping criterion
            relative to the machine precision and the value of f.  If
            ftol < 0.0, ftol is set to 0.0.  Defaults to 0.
        rescale : float
            Scaling factor (in log10) used to trigger rescaling.  If
            0, rescale at each iteration.  If a large value, never
            rescale.  If < 0, rescale is set to 1.3.

    :Returns:

        x : list of floats
            The solution.
        nfeval : int
            The number of function evaluations.
        rc :
            Return code (corresponding message in optimize.tnc.RCSTRINGS).

    :SeeAlso:

      - fmin, fmin_powell, fmin_cg, fmin_bfgs, fmin_ncg : multivariate
        local optimizers

      - leastsq : nonlinear least squares minimizer

      - fmin_l_bfgs_b, fmin_tnc, fmin_cobyla : constrained
        multivariate optimizers

      - anneal, brute : global optimizers

      - fminbound, brent, golden, bracket : local scalar minimizers

      - fsolve : n-dimenstional root-finding

      - brentq, brenth, ridder, bisect, newton : one-dimensional root-finding

      - fixed_point : scalar fixed-point finder

    """

    n = len(x0)

    if bounds is None:
        bounds = [(None,None)] * n
    if len(bounds) != n:
        raise ValueError('length of x0 != length of bounds')

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

    low = [0]*n
    up = [0]*n
    for i in range(n):
        l,u = bounds[i]
        if l is None:
            low[i] = -inf
        else:
            low[i] = l
        if u is None:
            up[i] = inf
        else:
            up[i] = u

    if scale == None:
        scale = []

    if maxfun == None:
        maxfun = max(1000, 100*len(x0))

    return moduleTNC.minimize(func_and_grad, x0, low, up, scale, messages,
                              maxCGit, maxfun, eta, stepmx, accuracy,
                              fmin, ftol, rescale)

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
        x, nf, rc = fmin_tnc(function, [-7, 3], bounds=([-10, 10], [1, 10]))

        print "After", nf, "function evaluations, TNC returned:", RCSTRINGS[rc]
        print "x =", x
        print "exact value = [0, 1]"
        print

    example()
