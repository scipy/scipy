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

import moduleTNC
from scipy_base import asarray

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

HUGE_VAL=1e500 # No standard representation of Infinity in Python 2.3.3

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


def fmin_tnc(func, x0, fprime=None, args=(), approx_grad=False, bounds=None, epsilon=1e-8,
             scale=None, messages=MSG_ALL, maxCGit=-1, maxfun=None, eta=-1,
             stepmx=0, accuracy=0, fmin=0, ftol=0, rescale=-1):
    """Minimize a function with variables subject to bounds, using gradient 
    information.
    
    returns (rc, nfeval, x).
    
    Inputs:

    func    -- function to minimize. Called as func(x, *args)
    
    x0      -- initial guess to minimum
    
    fprime  -- gradient of func. If None, then func returns the function
               value and the gradient ( f, g = func(x, *args) ).
               Called as fprime(x, *args)
                   
    args    -- arguments to pass to function

    approx_grad -- if true, approximate the gradient numerically    

    bounds  -- a list of (min, max) pairs for each element in x, defining
               the bounds on that parameter. Use None for one of min or max
               when there is no bound in that direction

        scale     : scaling factors to apply to each variable (a list of floats)
                    if None, the factors are up-low for interval bounded variables
                    and 1+|x] fo the others.
                    defaults to None
        messages  : bit mask used to select messages display during minimization
                    values defined in the optimize.tnc.MSGS dict.
                    defaults to optimize.tnc.MGS_ALL
        maxCGit   : max. number of hessian*vector evaluation per main iteration
                    if maxCGit == 0, the direction chosen is -gradient
                    if maxCGit < 0, maxCGit is set to max(1,min(50,n/2))
                    defaults to -1
        maxnfeval : max. number of function evaluation
                    if None, maxnfeval is set to max(1000, 100*len(x0))
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
                    relative to the machine precision and the value of f.
                    if ftol < 0.0, ftol is set to 0.0
                    defaults to 0
        rescale   : Scaling factor (in log10) used to trigger rescaling
                    if 0, rescale at each iteration
                    if a large value, never rescale
                    if < 0, rescale is set to 1.3

    Outputs:

        x         : the solution (a list of floats)
        nfeval    : the number of function evaluations
        rc        : return code (corresponding message in optimize.tnc.RCSTRINGS)
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
            low[i] = -HUGE_VAL
        else:
            low[i] = l
        if u is None:
            up[i] = HUGE_VAL
        else:
            up[i] = l
        
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
                rc, nf, x = minimize(function, [-7, 3], bounds=([-10, 10], [1, 10]))

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
        tests.append ((test1fg, [-2,1], ([-HUGE_VAL,None],[-1.5,None]), [1,1]))

        def test2fg(x):
                f = 100.0*pow((x[1]-pow(x[0],2)),2)+pow(1.0-x[0],2)
                dif = [0,0]
                dif[1] = 200.0*(x[1]-pow(x[0],2))
                dif[0] = -2.0*(x[0]*(dif[1]-1.0)+1.0)
                return f, dif
        tests.append ((test2fg, [-2,1], [(-HUGE_VAL,None),(1.5,None)], [-1.2210262419616387,1.5]))

        def test3fg(x):
                f = x[1]+pow(x[1]-x[0],2)*1.0e-5
                dif = [0,0]
                dif[0] = -2.0*(x[1]-x[0])*1.0e-5
                dif[1] = 1.0-dif[0]
                return f, dif
        tests.append ((test3fg, [10,1], [(-HUGE_VAL,None),(0.0, None)], [0,0]))

        def test4fg(x):
                f = pow(x[0]+1.0,3)/3.0+x[1]
                dif = [0,0]
                dif[0] = pow(x[0]+1.0,2)
                dif[1] = 1.0
                return f, dif
        tests.append ((test4fg, [1.125,0.125], [(1, None),(0, None)], [1,0]))

        from math import *

        def test5fg(x):
                f = sin(x[0]+x[1])+pow(x[0]-x[1],2)-1.5*x[0]+2.5*x[1]+1.0
                dif = [0,0]
                v1 = cos(x[0]+x[1]);
                v2 = 2.0*(x[0]-x[1]);

                dif[0] = v1+v2-1.5;
                dif[1] = v1-v2+2.5;
                return f, dif
        tests.append ((test5fg, [0,0], [(-1.5, 4),(-3,3)], [-0.54719755119659763, -1.5471975511965976]))

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
        tests.append ((test38fg, [-3,-1,-3,-1], [(-10,10)]*4, [1]*4))

        def test45fg(x):
                f = 2.0-x[0]*x[1]*x[2]*x[3]*x[4]/120.0
                dif = [0]*5
                dif[0] = -x[1]*x[2]*x[3]*x[4]/120.0
                dif[1] = -x[0]*x[2]*x[3]*x[4]/120.0
                dif[2] = -x[0]*x[1]*x[3]*x[4]/120.0
                dif[3] = -x[0]*x[1]*x[2]*x[4]/120.0
                dif[4] = -x[0]*x[1]*x[2]*x[3]/120.0
                return f, dif
        tests.append ((test45fg, [2]*5, [(0,1),(0,2),(0,3),(0,4),(0,5)], [1,2,3,4,5]))

        def test(fg, x, bounds, xopt):
                print "** Test", fg.__name__
                rc, nf, x = minimize(fg, x, bounds=bounds,  messages = MSG_NONE, maxnfeval = 200)
                print "After", nf, "function evaluations, TNC returned:", RCSTRINGS[rc]
                print "x =", x
                print "exact value =", xopt
                enorm = 0.0
                norm = 1.0
                for y,yo in zip(x, xopt):
                        enorm += (y-yo)*(y-yo)
                        norm += yo*yo
                e = pow(enorm/norm, 0.5)
                print "Error =", e
                if e > 1e-8:
                        raise "Test "+fg.__name__+" failed"

        for fg, x, bounds, xopt in tests:
                test(fg, x, bounds, xopt)
        
        print
        print "** All TNC tests passed."
