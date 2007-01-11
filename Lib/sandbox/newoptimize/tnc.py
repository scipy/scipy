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

import moduleTNC

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

HUGE_VAL=1e200*1e200 # No standard representation of Infinity in Python 2.3.3

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

def minimize(function, x, low = None, up = None, scale = None, offset = None,
        messages = MSG_ALL, maxCGit = -1, maxnfeval = None, eta = -1, stepmx = 0,
        accuracy = 0, fmin = 0, ftol = -1, xtol = -1, pgtol = -1, rescale = -1):
    """Minimize a function with variables subject to bounds, using gradient
    information.

    returns (rc, nfeval, x).

    Inputs:
    x         : initial estimate (a list of floats)
    function  : the function to minimize. Must take one argument, x and return
                f and g, where f is the value of the function and g its
                gradient (a list of floats).
                if the function returns None, the minimization is aborted.
    low, up   : the bounds (lists of floats)
                set low[i] to -HUGE_VAL to remove the lower bound
                set up[i] to HUGE_VAL to remove the upper bound
                if low == None, the lower bounds are removed.
                if up == None, the upper bounds are removed.
                low and up defaults to None
    scale     : scaling factors to apply to each variable (a list of floats)
                if None, the factors are up-low for interval bounded variables
                and 1+|x] fo the others.
                defaults to None
    offset    : constant to substract to each variable
                if None, the constant are (up+low)/2 for interval bounded
                variables and x for the others.
    messages  : bit mask used to select messages display during minimization
                values defined in the MSGS dict.
                defaults to MGS_ALL
    maxCGit   : max. number of hessian*vector evaluation per main iteration
                if maxCGit == 0, the direction chosen is -gradient
                if maxCGit < 0, maxCGit is set to max(1,min(50,n/2))
                defaults to -1
    maxnfeval : max. number of function evaluation
                if None, maxnfeval is set to max(100, 10*len(x0))
                defaults to None
    eta       : severity of the line search. if < 0 or > 1, set to 0.25
                defaults to -1
    stepmx    : maximum step for the line search. may be increased during call
                if too small, will be set to 10.0
                defaults to 0
    accuracy  : relative precision for finite difference calculations
                if <= machine_precision, set to sqrt(machine_precision)
                defaults to 0
    fmin      : minimum function value estimate
                defaults to 0
    ftol      : precision goal for the value of f in the stoping criterion
                if ftol < 0.0, ftol is set to 0.0
                defaults to -1
    xtol      : precision goal for the value of x in the stopping criterion
                (after applying x scaling factors)
                if xtol < 0.0, xtol is set to sqrt(machine_precision)
                defaults to -1
    pgtol     : precision goal for the value of the projected gradient in the
                stopping criterion (after applying x scaling factors)
                if pgtol < 0.0, pgtol is set to 1e-2 * sqrt(accuracy)
                setting it to 0.0 is not recommended.
                defaults to -1
    rescale   : f scaling factor (in log10) used to trigger f value rescaling
                if 0, rescale at each iteration
                if a large value, never rescale
                if < 0, rescale is set to 1.3

    Outputs:
    x         : the solution (a list of floats)
    nfeval    : the number of function evaluations
    rc        : return code as defined in the RCSTRINGS dict"""

    if low == None:
        low = [-HUGE_VAL]*len(x)

    if up == None:
        up = [HUGE_VAL]*len(x)

    if scale == None:
        scale = []

    if offset == None:
        offset = []

    if maxnfeval == None:
        maxnfeval = max(100, 10*len(x))

    return moduleTNC.minimize(function, x, low, up, scale, offset,
            messages, maxCGit, maxnfeval, eta, stepmx, accuracy,
            fmin, ftol, xtol, pgtol, rescale)

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
        rc, nf, x = minimize(function, [-7, 3], [-10, 1], [10, 10])

        print "After", nf, "function evaluations, TNC returned:", RCSTRINGS[rc]
        print "x =", x
        print "exact value = [0, 1]"
        print

    example()

    # Tests
    # These tests are taken from Prof. K. Schittkowski test examples for
    # constrained nonlinear programming.
    # http://www.uni-bayreuth.de/departments/math/~kschittkowski/home.htm
    tests = []
    def test1fg(x):
        f = 100.0*pow((x[1]-pow(x[0],2)),2)+pow(1.0-x[0],2)
        dif = [0,0]
        dif[1] = 200.0*(x[1]-pow(x[0],2))
        dif[0] = -2.0*(x[0]*(dif[1]-1.0)+1.0)
        return f, dif
    tests.append ((test1fg, [-2,1], [-HUGE_VAL, -1.5], None, [1,1]))

    def test2fg(x):
        f = 100.0*pow((x[1]-pow(x[0],2)),2)+pow(1.0-x[0],2)
        dif = [0,0]
        dif[1] = 200.0*(x[1]-pow(x[0],2))
        dif[0] = -2.0*(x[0]*(dif[1]-1.0)+1.0)
        return f, dif
    tests.append ((test2fg, [-2,1], [-HUGE_VAL, 1.5], None, [-1.2210262419616387,1.5]))

    def test3fg(x):
        f = x[1]+pow(x[1]-x[0],2)*1.0e-5
        dif = [0,0]
        dif[0] = -2.0*(x[1]-x[0])*1.0e-5
        dif[1] = 1.0-dif[0]
        return f, dif
    tests.append ((test3fg, [10,1], [-HUGE_VAL, 0.0], None, [0,0]))

    def test4fg(x):
        f = pow(x[0]+1.0,3)/3.0+x[1]
        dif = [0,0]
        dif[0] = pow(x[0]+1.0,2)
        dif[1] = 1.0
        return f, dif
    tests.append ((test4fg, [1.125,0.125], [1, 0], None, [1,0]))

    from math import *

    def test5fg(x):
        f = sin(x[0]+x[1])+pow(x[0]-x[1],2)-1.5*x[0]+2.5*x[1]+1.0
        dif = [0,0]
        v1 = cos(x[0]+x[1]);
        v2 = 2.0*(x[0]-x[1]);

        dif[0] = v1+v2-1.5;
        dif[1] = v1-v2+2.5;
        return f, dif
    tests.append ((test5fg, [0,0], [-1.5, -3], [4,3], [-0.54719755119659763, -1.5471975511965976]))

    def test38fg(x):
        f = (100.0*pow(x[1]-pow(x[0],2),2)+pow(1.0-x[0],2)+90.0*pow(x[3]-pow(x[2],2),2) \
                +pow(1.0-x[2],2)+10.1*(pow(x[1]-1.0,2)+pow(x[3]-1.0,2)) \
                +19.8*(x[1]-1.0)*(x[3]-1.0))*1.0e-5
        dif = [0,0,0,0]
        dif[0] = (-400.0*x[0]*(x[1]-pow(x[0],2))-2.0*(1.0-x[0]))*1.0e-5
        dif[1] = (200.0*(x[1]-pow(x[0],2))+20.2*(x[1]-1.0)+19.8*(x[3]-1.0))*1.0e-5
        dif[2] = (-360.0*x[2]*(x[3]-pow(x[2],2))-2.0*(1.0-x[2]))*1.0e-5
        dif[3] = (180.0*(x[3]-pow(x[2],2))+20.2*(x[3]-1.0)+19.8*(x[1]-1.0))*1.0e-5
        return f, dif
    tests.append ((test38fg, [-3,-1,-3,-1], [-10]*4, [10]*4, [1]*4))

    def test45fg(x):
        f = 2.0-x[0]*x[1]*x[2]*x[3]*x[4]/120.0
        dif = [0]*5
        dif[0] = -x[1]*x[2]*x[3]*x[4]/120.0
        dif[1] = -x[0]*x[2]*x[3]*x[4]/120.0
        dif[2] = -x[0]*x[1]*x[3]*x[4]/120.0
        dif[3] = -x[0]*x[1]*x[2]*x[4]/120.0
        dif[4] = -x[0]*x[1]*x[2]*x[3]/120.0
        return f, dif
    tests.append ((test45fg, [2]*5, [0]*5, [1,2,3,4,5], [1,2,3,4,5]))

    def test(fg, x, low, up, xopt):
        print "** Test", fg.__name__
        rc, nf, x = minimize(fg, x, low, up, messages = MSG_NONE, maxnfeval = 200)
        print "After", nf, "function evaluations, TNC returned:", RCSTRINGS[rc]
        print "x =", x
        print "exact value =", xopt
        enorm = 0.0
        norm = 1.0
        for y,yo in zip(x, xopt):
            enorm += (y-yo)*(y-yo)
            norm += yo*yo
        ex = pow(enorm/norm, 0.5)
        print "X Error =", ex
        ef = abs(fg(xopt)[0] - fg(x)[0])
        print "F Error =", ef
        if ef > 1e-8:
            raise "Test "+fg.__name__+" failed"

    for fg, x, low, up, xopt in tests:
        test(fg, x, low, up, xopt)

    print
    print "** All TNC tests passed."
