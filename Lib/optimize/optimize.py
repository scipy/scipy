from __future__ import nested_scopes
# ******NOTICE***************
# optimize.py module by Travis E. Oliphant
#
# You may copy and use this module as you see fit with no
# guarantee implied provided you keep this notice in all copies.
# *****END NOTICE************

# A collection of optimization algorithms.  Version 0.5
# CHANGES
#  Added fminbound (July 2001)
#  Added brute (Aug. 2002)
#  Finished line search satisfying strong Wolfe conditions (Mar. 2004)
#  Updated strong Wolfe conditions line search to use cubic-interpolation (Mar. 2004)

# Minimization routines
"""optimize.py

A collection of general-purpose optimization routines using Numeric


N-D Algorithms
===============
fmin        ---      Nelder-Mead Simplex algorithm (uses only function calls).
fmin_powell ---      Powell (direction set) method (uses only function calls).
fmin_cg     ---      Non-linear (Polak-Rubiere) conjugate gradient algorithm
                       (can use function and gradient).
fmin_bfgs   ---      Quasi-Newton method (Broyden-Fletcher-Goldfarb-Shanno
                       (can use function and gradient).
fmin_ncg    ---      Line-search Newton Conjugate Gradient (can use function, 
                       gradient and hessian).
                     

brute       ---      Perform a brute force search for the minimum
                       with final optimization if desired.

1-D Algorithms
===============
brent       ---      Use Brent's method (does not need inital guess)
fminbound   ---      Bounded minimization for scalar functions on an
                     interval using Brent's parabolic/golden_mean method.
golden      --       Use Golden Section method (does not need initial guess)

bracket     ---      Find a bracket containing the minimum.



"""


__all__ = ['fmin', 'fmin_powell','fmin_bfgs', 'fmin_ncg', 'fmin_cg',
           'fminbound','brent', 'golden','bracket','rosen','rosen_der',
           'rosen_hess', 'rosen_hess_prod', 'brute', 'approx_fprime',
           'line_search', 'check_grad']

import Numeric
import MLab
from scipy_base import atleast_1d, eye, mgrid, argmin, zeros, shape, \
     squeeze, isscalar, vectorize, asarray, absolute, sqrt, Inf
import scipy_base
import linesearch
Num = Numeric
max = MLab.max
min = MLab.min
abs = absolute
import __builtin__
pymin = __builtin__.min
pymax = __builtin__.max
__version__="0.7"
_epsilon = sqrt(scipy_base.limits.double_epsilon)

def vecnorm(x, ord=2):
    if ord == Inf:
        return scipy_base.amax(abs(x))
    elif ord == -Inf:
        return scipy_base.amin(abs(x))
    else:
        return scipy_base.sum(abs(x)**ord)**(1.0/ord)
        
def rosen(x):  # The Rosenbrock function
    x = asarray(x)
    return MLab.sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)

def rosen_der(x):
    x = asarray(x)
    xm = x[1:-1]
    xm_m1 = x[:-2]
    xm_p1 = x[2:]
    der = MLab.zeros(x.shape,x.typecode())
    der[1:-1] = 200*(xm-xm_m1**2) - 400*(xm_p1 - xm**2)*xm - 2*(1-xm)
    der[0] = -400*x[0]*(x[1]-x[0]**2) - 2*(1-x[0])
    der[-1] = 200*(x[-1]-x[-2]**2)
    return der

def rosen_hess(x):
    x = atleast_1d(x)
    H = MLab.diag(-400*x[:-1],1) - MLab.diag(400*x[:-1],-1)
    diagonal = Num.zeros(len(x),x.typecode())
    diagonal[0] = 1200*x[0]-400*x[1]+2
    diagonal[-1] = 200
    diagonal[1:-1] = 202 + 1200*x[1:-1]**2 - 400*x[2:]
    H = H + MLab.diag(diagonal)
    return H

def rosen_hess_prod(x,p):
    x = atleast_1d(x)
    Hp = Num.zeros(len(x),x.typecode())
    Hp[0] = (1200*x[0]**2 - 400*x[1] + 2)*p[0] - 400*x[0]*p[1]
    Hp[1:-1] = -400*x[:-2]*p[:-2]+(202+1200*x[1:-1]**2-400*x[2:])*p[1:-1] \
               -400*x[1:-1]*p[2:]
    Hp[-1] = -400*x[-2]*p[-2] + 200*p[-1]
    return Hp
        
def fmin(func, x0, args=(), xtol=1e-4, ftol=1e-4, maxiter=None, maxfun=None, 
         full_output=0, disp=1, retall=0):
    """Minimize a function using the downhill simplex algorithm.

    Description:
    
      Uses a Nelder-Mead simplex algorithm to find the minimum of function
      of one or more variables.

    Inputs:

      func -- the Python function or method to be minimized.
      x0 -- the initial guess.
      args -- extra arguments for func.

    Outputs: (xopt, {fopt, iter, funcalls, warnflag})

      xopt -- minimizer of function

      fopt -- value of function at minimum: fopt = func(xopt)
      iter -- number of iterations
      funcalls -- number of function calls
      warnflag -- Integer warning flag:
                  1 : 'Maximum number of function evaluations.'
                  2 : 'Maximum number of iterations.'
      allvecs  -- a list of solutions at each iteration

    Additional Inputs:

      xtol -- acceptable relative error in xopt for convergence.
      ftol -- acceptable relative error in func(xopt) for convergence.
      maxiter -- the maximum number of iterations to perform.
      maxfun -- the maximum number of function evaluations.
      full_output -- non-zero if fval and warnflag outputs are desired.
      disp -- non-zero to print convergence messages.
      retall -- non-zero to return list of solutions at each iteration
      
      """
    x0 = asarray(x0)
    N = len(x0)
    rank = len(x0.shape)
    if not -1 < rank < 2:
        raise ValueError, "Initial guess must be a scalar or rank-1 sequence."
    if maxiter is None:
        maxiter = N * 200
    if maxfun is None:
        maxfun = N * 200

    rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
    one2np1 = range(1,N+1)

    if rank == 0:
        sim = Num.zeros((N+1,),x0.typecode())
    else:        
        sim = Num.zeros((N+1,N),x0.typecode())
    fsim = Num.zeros((N+1,),'d')
    sim[0] = x0
    if retall:
        allvecs = [sim[0]]
    fsim[0] = apply(func,(x0,)+args)
    nonzdelt = 0.05
    zdelt = 0.00025
    for k in range(0,N):
        y = Num.array(x0,copy=1)
        if y[k] != 0:
            y[k] = (1+nonzdelt)*y[k]
        else:
            y[k] = zdelt

        sim[k+1] = y
        f = apply(func,(y,)+args)
        fsim[k+1] = f

    ind = Num.argsort(fsim)
    fsim = Num.take(fsim,ind)  # sort so sim[0,:] has the lowest function value
    sim = Num.take(sim,ind,0)
    
    iterations = 1
    funcalls = N+1
    
    while (funcalls < maxfun and iterations < maxiter):
        if (max(Num.ravel(abs(sim[1:]-sim[0]))) <= xtol \
            and max(abs(fsim[0]-fsim[1:])) <= ftol):
            break

        xbar = Num.add.reduce(sim[:-1],0) / N
        xr = (1+rho)*xbar - rho*sim[-1]
        fxr = apply(func,(xr,)+args)
        funcalls = funcalls + 1
        doshrink = 0

        if fxr < fsim[0]:
            xe = (1+rho*chi)*xbar - rho*chi*sim[-1]
            fxe = apply(func,(xe,)+args)
            funcalls = funcalls + 1

            if fxe < fxr:
                sim[-1] = xe
                fsim[-1] = fxe
            else:
                sim[-1] = xr
                fsim[-1] = fxr
        else: # fsim[0] <= fxr
            if fxr < fsim[-2]:
                sim[-1] = xr
                fsim[-1] = fxr
            else: # fxr >= fsim[-2]
                # Perform contraction
                if fxr < fsim[-1]:
                    xc = (1+psi*rho)*xbar - psi*rho*sim[-1]
                    fxc = apply(func,(xc,)+args)
                    funcalls = funcalls + 1

                    if fxc <= fxr:
                        sim[-1] = xc
                        fsim[-1] = fxc
                    else:
                        doshrink=1
                else:
                    # Perform an inside contraction
                    xcc = (1-psi)*xbar + psi*sim[-1]
                    fxcc = apply(func,(xcc,)+args)
                    funcalls = funcalls + 1

                    if fxcc < fsim[-1]:
                        sim[-1] = xcc
                        fsim[-1] = fxcc
                    else:
                        doshrink = 1

                if doshrink:
                    for j in one2np1:
                        sim[j] = sim[0] + sigma*(sim[j] - sim[0])
                        fsim[j] = apply(func,(sim[j],)+args)
                    funcalls = funcalls + N

        ind = Num.argsort(fsim)
        sim = Num.take(sim,ind,0)
        fsim = Num.take(fsim,ind)
        iterations = iterations + 1
        if retall:
            allvecs.append(sim[0])        

    x = sim[0]
    fval = min(fsim)
    warnflag = 0

    if funcalls >= maxfun:
        warnflag = 1
        if disp:
            print "Warning: Maximum number of function evaluations has "\
                  "been exceeded."
    elif iterations >= maxiter:
        warnflag = 2
        if disp:
            print "Warning: Maximum number of iterations has been exceeded"
    else:
        if disp:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % iterations
            print "         Function evaluations: %d" % funcalls


    if full_output:
        retlist = x, fval, iterations, funcalls, warnflag
        if retall:
            retlist += (allvecs,)
    else: 
        retlist = x
        if retall:
            retlist = (x, allvecs)

    return retlist



def _cubicmin(a,fa,fpa,b,fb,c,fc):
    # finds the minimizer for a cubic polynomial that goes through the
    #  points (a,fa), (b,fb), and (c,fc) with derivative at a of fpa.
    #
    # if no minimizer can be found return None
    #
    # f(x) = A *(x-a)^3 + B*(x-a)^2 + C*(x-a) + D

    C = fpa
    D = fa
    db = b-a
    dc = c-a
    if (db == 0) or (dc == 0) or (b==c): return None
    denom = (db*dc)**2 * (db-dc)
    [A,B] = Num.dot([[dc**2, -db**2],[-dc**3, db**3]],[fb-fa-C*db,fc-fa-C*dc])
    A /= denom
    B /= denom
    radical = B*B-3*A*C
    if radical < 0:  return None
    if (A == 0): return None
    xmin = a + (-B + sqrt(radical))/(3*A)
    return xmin

    
def _quadmin(a,fa,fpa,b,fb):
    # finds the minimizer for a quadratic polynomial that goes through
    #  the points (a,fa), (b,fb) with derivative at a of fpa
    # f(x) = B*(x-a)^2 + C*(x-a) + D
    D = fa
    C = fpa
    db = b-a*1.0
    if (db==0): return None
    B = (fb-D-C*db)/(db*db)
    if (B <= 0): return None
    xmin = a  - C / (2.0*B)
    return xmin

def zoom(a_lo, a_hi, phi_lo, phi_hi, derphi_lo,
         phi, derphi, phi0, derphi0, c1, c2):
    maxiter = 10
    i = 0
    delta1 = 0.2  # cubic interpolant check
    delta2 = 0.1  # quadratic interpolant check
    phi_rec = phi0
    a_rec = 0
    while 1:
        # interpolate to find a trial step length between a_lo and a_hi
        # Need to choose interpolation here.  Use cubic interpolation and then if the
        #  result is within delta * dalpha or outside of the interval bounded by a_lo or a_hi 
        #  then use quadratic interpolation, if the result is still too close, then use bisection

        dalpha = a_hi-a_lo;
        if dalpha < 0: a,b = a_hi,a_lo
        else: a,b = a_lo, a_hi

        # minimizer of cubic interpolant
        #    (uses phi_lo, derphi_lo, phi_hi, and the most recent value of phi)
        #      if the result is too close to the end points (or out of the interval)
        #         then use quadratic interpolation with phi_lo, derphi_lo and phi_hi
        #      if the result is stil too close to the end points (or out of the interval)
        #         then use bisection

        if (i > 0):
            cchk = delta1*dalpha
            a_j = _cubicmin(a_lo, phi_lo, derphi_lo, a_hi, phi_hi, a_rec, phi_rec)
        if (i==0) or (a_j is None) or (a_j > b-cchk) or (a_j < a+cchk):
            qchk = delta2*dalpha
            a_j = _quadmin(a_lo, phi_lo, derphi_lo, a_hi, phi_hi)
            if (a_j is None) or (a_j > b-qchk) or (a_j < a+qchk):
                a_j = a_lo + 0.5*dalpha
#                print "Using bisection."
#            else: print "Using quadratic."
#        else: print "Using cubic."

        # Check new value of a_j 

        phi_aj = phi(a_j)
        if (phi_aj > phi0 + c1*a_j*derphi0) or (phi_aj >= phi_lo):
            phi_rec = phi_hi
            a_rec = a_hi
            a_hi = a_j
            phi_hi = phi_aj
        else:
            derphi_aj = derphi(a_j)            
            if abs(derphi_aj) <= -c2*derphi0:
                a_star = a_j
                val_star = phi_aj
                valprime_star = derphi_aj
                break
            if derphi_aj*(a_hi - a_lo) >= 0:
                phi_rec = phi_hi
                a_rec = a_hi
                a_hi = a_lo
                phi_hi = phi_lo
            else:                
                phi_rec = phi_lo
                a_rec = a_lo
            a_lo = a_j
            phi_lo = phi_aj
            derphi_lo = derphi_aj
        i += 1
        if (i > maxiter):
            a_star = a_j
            val_star = phi_aj
            valprime_star = None
            break    
    return a_star, val_star, valprime_star

def line_search(f, myfprime, xk, pk, gfk, old_fval, old_old_fval,
                args=(), c1=1e-4, c2=0.9, amax=50):
    """Find alpha that satisfies strong Wolfe conditions. 
    
    Uses the line search algorithm to enforce strong Wolfe conditions 
    Wright and Nocedal, 'Numerical Optimization', 1999, pg. 59-60

    For the zoom phase it uses an algorithm by 
    Outputs: (alpha0, gc, fc)
    """

    global _ls_fc, _ls_gc, _ls_ingfk
    _ls_fc = 0
    _ls_gc = 0
    _ls_ingfk = None
    def phi(alpha):
        global _ls_fc
        _ls_fc += 1
        return f(xk+alpha*pk,*args)

    if isinstance(myfprime,type(())):
        def phiprime(alpha):
            global _ls_fc, _ls_ingfk
            _ls_fc += len(xk)+1
            eps = myfprime[1]
            fprime = myfprime[0]
            newargs = (f,eps) + args
            _ls_ingfk = fprime(xk+alpha*pk,*newargs)  # store for later use
            return Num.dot(_ls_ingfk,pk)
    else:
        fprime = myfprime
        def phiprime(alpha):
            global _ls_gc, _ls_ingfk
            _ls_gc += 1
            _ls_ingfk = fprime(xk+alpha*pk,*args)  # store for later use
            return Num.dot(_ls_ingfk,pk)

    alpha0 = 0
    phi0 = old_fval
    derphi0 = Num.dot(gfk,pk)

    alpha1 = pymin(1.0,1.01*2*(phi0-old_old_fval)/derphi0)
    phi_a1 = phi(alpha1)
    #derphi_a1 = phiprime(alpha1)  evaluated below

    phi_a0 = phi0
    derphi_a0 = derphi0

    i = 1
    maxiter = 10
    while 1:         # bracketing phase 
        if (phi_a1 > phi0 + c1*alpha1*derphi0) or \
           ((phi_a1 >= phi_a0) and (i > 1)):
            alpha_star, fval_star, fprime_star = \
                        zoom(alpha0, alpha1, phi_a0,
                             phi_a1, derphi_a0, phi, phiprime,
                             phi0, derphi0, c1, c2)
            break

        derphi_a1 = phiprime(alpha1)
        if (abs(derphi_a1) <= -c2*derphi0):
            alpha_star = alpha1
            fval_star = phi_a1
            fprime_star = derphi_a1
            break

        if (derphi_a1 >= 0):
            alpha_star, fval_star, fprime_star = \
                        zoom(alpha1, alpha0, phi_a1,
                             phi_a0, derphi_a1, phi, phiprime,
                             phi0, derphi0, c1, c2)
            break

        alpha2 = 2 * alpha1   # increase by factor of two on each iteration
        i = i + 1
        alpha0 = alpha1
        alpha1 = alpha2
        phi_a0 = phi_a1
        phi_a1 = phi(alpha1)
        derphi_a0 = derphi_a1

        # stopping test if lower function not found
        if (i > maxiter):          
            alpha_star = alpha1
            fval_star = phi_a1
            fprime_star = None
            break

    if fprime_star is not None:
        # fprime_star is a number (derphi) -- so use the most recently calculated gradient
        #                                     used in computing it derphi = gfk*pk
        #                                     this is the gradient at the next step
        #                                     no need to compute it again in the outer loop.
        fprime_star = _ls_ingfk
        
    return alpha_star, _ls_fc, _ls_gc, fval_star, old_fval, fprime_star
    

def line_search_BFGS(f, xk, pk, gfk, old_fval, args=(), c1=1e-4, alpha0=1):
    """Minimize over alpha, the function f(xk+alpha pk)

    Uses the interpolation algorithm (Armiijo backtracking) as suggested by
    Wright and Nocedal in 'Numerical Optimization', 1999, pg. 56-57
    
    Outputs: (alpha, fc, gc)
    """

    fc = 0
    phi0 = old_fval                            # compute f(xk) -- done in past loop
    phi_a0 = apply(f,(xk+alpha0*pk,)+args)     # compute f
    fc = fc + 1
    derphi0 = Num.dot(gfk,pk)

    if (phi_a0 <= phi0 + c1*alpha0*derphi0):
        return alpha0, fc, 0, phi_a0

    # Otherwise compute the minimizer of a quadratic interpolant:

    alpha1 = -(derphi0) * alpha0**2 / 2.0 / (phi_a0 - phi0 - derphi0 * alpha0)
    phi_a1 = apply(f,(xk+alpha1*pk,)+args)
    fc = fc + 1

    if (phi_a1 <= phi0 + c1*alpha1*derphi0):
        return alpha1, fc, 0, phi_a1

    # Otherwise loop with cubic interpolation until we find an alpha which 
    # satifies the first Wolfe condition (since we are backtracking, we will
    # assume that the value of alpha is not too small and satisfies the second
    # condition.

    while 1:       # we are assuming pk is a descent direction
        factor = alpha0**2 * alpha1**2 * (alpha1-alpha0)
        a = alpha0**2 * (phi_a1 - phi0 - derphi0*alpha1) - \
            alpha1**2 * (phi_a0 - phi0 - derphi0*alpha0)
        a = a / factor
        b = -alpha0**3 * (phi_a1 - phi0 - derphi0*alpha1) + \
            alpha1**3 * (phi_a0 - phi0 - derphi0*alpha0)
        b = b / factor

        alpha2 = (-b + Num.sqrt(abs(b**2 - 3 * a * derphi0))) / (3.0*a)
        phi_a2 = apply(f,(xk+alpha2*pk,)+args)
        fc = fc + 1

        if (phi_a2 <= phi0 + c1*alpha2*derphi0):
            return alpha2, fc, 0, phi_a2

        if (alpha1 - alpha2) > alpha1 / 2.0 or (1 - alpha2/alpha1) < 0.96:
            alpha2 = alpha1 / 2.0

        alpha0 = alpha1
        alpha1 = alpha2
        phi_a0 = phi_a1
        phi_a1 = phi_a2


def approx_fprime(xk,f,epsilon,*args):
    f0 = apply(f,(xk,)+args)
    grad = Num.zeros((len(xk),),'d')
    ei = Num.zeros((len(xk),),'d')
    for k in range(len(xk)):
        ei[k] = epsilon
        grad[k] = (apply(f,(xk+ei,)+args) - f0)/epsilon
        ei[k] = 0.0
    return grad

def check_grad(func, grad, x0, *args):
    return sqrt(sum((grad(x0,*args)-approx_fprime(x0,func,_epsilon,*args))**2))

def approx_fhess_p(x0,p,fprime,epsilon,*args):
    f2 = apply(fprime,(x0+epsilon*p,)+args)
    f1 = apply(fprime,(x0,)+args)
    return (f2 - f1)/epsilon


def fmin_bfgs(f, x0, fprime=None, args=(), gtol=1e-5, norm=Inf,
              epsilon=_epsilon, maxiter=None, full_output=0, disp=1,
              retall=0):
    """Minimize a function using the BFGS algorithm.

    Description:

      Optimize the function, f, whose gradient is given by fprime using the
      quasi-Newton method of Broyden, Fletcher, Goldfarb, and Shanno (BFGS)
      See Wright, and Nocedal 'Numerical Optimization', 1999, pg. 198.

    Inputs:

      f -- the Python function or method to be minimized.
      x0 -- the initial guess for the minimizer.

      fprime -- a function to compute the gradient of f.
      args -- extra arguments to f and fprime.
      gtol -- gradient norm must be less than gtol before succesful termination
      norm -- order of norm (Inf is max, -Inf is min)
      epsilon -- if fprime is approximated use this value for
                 the step size (can be scalar or vector)

    Outputs: (xopt, {fopt, gopt, Hopt, func_calls, grad_calls, warnflag}, <allvecs>)

      xopt -- the minimizer of f.

      fopt -- the value of f(xopt).
      gopt -- the value of f'(xopt).  (Should be near 0)
      Bopt -- the value of 1/f''(xopt).  (inverse hessian matrix)
      func_calls -- the number of function_calls.
      grad_calls -- the number of gradient calls.
      warnflag -- an integer warning flag:
                  1 : 'Maximum number of iterations exceeded.'
                  2 : 'Gradient and/or function calls not changing'
      allvecs  --  a list of all iterates  (only returned if retall==1)

    Additional Inputs:

      maxiter -- the maximum number of iterations.
      full_output -- if non-zero then return fopt, func_calls, grad_calls,
                     and warnflag in addition to xopt.
      disp -- print convergence message if non-zero.
      retall -- return a list of results at each iteration if non-zero
      """
    app_fprime = 0
    if fprime is None:
        app_fprime = 1

    x0 = asarray(x0)
    if maxiter is None:
        maxiter = len(x0)*200
    func_calls = 0
    grad_calls = 0
    k = 0
    N = len(x0)
    I = MLab.eye(N)
    Hk = I
    old_fval = f(x0,*args)
    old_old_fval = old_fval + 5000
    func_calls += 1
    if app_fprime:
        gfk = apply(approx_fprime,(x0,f,epsilon)+args)
        myfprime = (approx_fprime,epsilon)
        func_calls = func_calls + len(x0) + 1
    else:
        gfk = apply(fprime,(x0,)+args)
        myfprime = fprime
        grad_calls = grad_calls + 1
    xk = x0
    if retall:
        allvecs = [x0]
    sk = [2*gtol]
    warnflag = 0
    gnorm = vecnorm(gfk,ord=norm)
    while (gnorm > gtol) and (k < maxiter):
        pk = -Num.dot(Hk,gfk)
        alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = \
           linesearch.line_search(f,myfprime,xk,pk,gfk,
                                  old_fval,old_old_fval,args=args)
        if alpha_k is None:  # line search failed try different one.
            func_calls = func_calls + fc
            grad_calls = grad_calls + gc            
            alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = \
                     line_search(f,myfprime,xk,pk,gfk,
                                 old_fval,old_old_fval,args=args)
            
        func_calls = func_calls + fc
        grad_calls = grad_calls + gc
        xkp1 = xk + alpha_k * pk
        if retall:
            allvecs.append(xkp1)
        sk = xkp1 - xk
        xk = xkp1
        if gfkp1 is None:
            if app_fprime:
                gfkp1 = apply(approx_fprime,(xkp1,f,epsilon)+args)
                func_calls = func_calls + len(x0) + 1
            else:
                gfkp1 = apply(fprime,(xkp1,)+args)
                grad_calls = grad_calls + 1

        yk = gfkp1 - gfk
        gfk = gfkp1        
        k = k + 1
        gnorm = vecnorm(gfk,ord=norm)
        if (gnorm <= gtol):
            break

        try:
            rhok = 1 / Num.dot(yk,sk)
        except ZeroDivisionError:
            warnflag = 2
            break
            #print "Divide by zero encountered:  Hessian calculation reset."
            #Hk = I
        else:
            A1 = I - sk[:,Num.NewAxis] * yk[Num.NewAxis,:] * rhok
            A2 = I - yk[:,Num.NewAxis] * sk[Num.NewAxis,:] * rhok
            Hk = Num.dot(A1,Num.dot(Hk,A2)) + rhok * sk[:,Num.NewAxis] \
                 * sk[Num.NewAxis,:]

    if disp or full_output:
        fval = old_fval
    if warnflag == 2:
        if disp:
            print "Warning: Desired error not necessarily achieved due to precision loss"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls
        
    elif k >= maxiter:
        warnflag = 1
        if disp:
            print "Warning: Maximum number of iterations has been exceeded"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls
    else:
        if disp:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls

    if full_output:
        retlist = xk, fval, gfk, Hk, func_calls, grad_calls, warnflag
        if retall:
            retlist += (allvecs,)
    else: 
        retlist = xk
        if retall:
            retlist = (xk, allvecs)

    return retlist


def fmin_cg(f, x0, fprime=None, args=(), gtol=1e-5, norm=Inf, epsilon=_epsilon,
              maxiter=None, full_output=0, disp=1, retall=0):
    """Minimize a function with nonlinear conjugate gradient algorithm.

    Description:

      Optimize the function, f, whose gradient is given by fprime using the
      nonlinear conjugate gradient algorithm of Polak and Ribiere
      See Wright, and Nocedal 'Numerical Optimization', 1999, pg. 120-122.

    Inputs:

      f -- the Python function or method to be minimized.
      x0 -- the initial guess for the minimizer.

      fprime -- a function to compute the gradient of f.
      args -- extra arguments to f and fprime.
      gtol -- stop when norm of gradient is less than gtol
      norm -- order of vector norm to use
      epsilon -- if fprime is approximated use this value for
                 the step size (can be scalar or vector)

    Outputs: (xopt, {fopt, func_calls, grad_calls, warnflag}, {allvecs})

      xopt -- the minimizer of f.

      fopt -- the value of f(xopt).
      func_calls -- the number of function_calls.
      grad_calls -- the number of gradient calls.
      warnflag -- an integer warning flag:
                  1 : 'Maximum number of iterations exceeded.'
                  2 : 'Gradient and/or function calls not changing'
      allvecs  --  if retall then this vector of the iterates is returned

    Additional Inputs:

      maxiter -- the maximum number of iterations.
      full_output -- if non-zero then return fopt, func_calls, grad_calls,
                     and warnflag in addition to xopt.
      disp -- print convergence message if non-zero.
      retall -- return a list of results at each iteration if True
      """
    app_fprime = 0
    if fprime is None:
        app_fprime = 1

    x0 = asarray(x0)
    if maxiter is None:
        maxiter = len(x0)*200
    func_calls = 0
    grad_calls = 0
    k = 0
    N = len(x0)
    xk = x0
    old_fval = f(xk,*args)
    old_old_fval = old_fval + 5000
    func_calls +=1 

    if app_fprime:
        gfk = apply(approx_fprime,(x0,f,epsilon)+args)
        myfprime = (approx_fprime,epsilon)
        func_calls = func_calls + len(x0) + 1
    else:
        gfk = apply(fprime,(x0,)+args)
        myfprime = fprime
        grad_calls = grad_calls + 1
    if retall:
        allvecs = [xk]
    sk = [2*gtol]
    warnflag = 0
    pk = -gfk
    gnorm = vecnorm(gfk,ord=norm)
    while (gnorm > gtol) and (k < maxiter):
        deltak = Num.dot(gfk,gfk)
        alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = \
           linesearch.line_search(f,myfprime,xk,pk,gfk,old_fval,
                                  old_old_fval,args=args,c2=0.4)
        if alpha_k is None:  # line search failed -- use different one.
            func_calls += fc
            grad_calls += gc                  
            alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = \
                     line_search(f,myfprime,xk,pk,gfk,
                                 old_fval,old_old_fval,args=args)
        func_calls += fc
        grad_calls += gc        
        xk = xk + alpha_k*pk
        if retall:
            allvecs.append(xk)
        if gfkp1 is None:
            if app_fprime:
                gfkp1 = apply(approx_fprime,(xk,f,epsilon)+args)
                func_calls = func_calls + len(x0) + 1
            else:
                gfkp1 = apply(fprime,(xk,)+args)
                grad_calls = grad_calls + 1

        yk = gfkp1 - gfk
        beta_k = pymax(0,Num.dot(yk,gfkp1)/deltak)
        pk = -gfkp1 + beta_k * pk
        gfk = gfkp1
        gnorm = vecnorm(gfk,ord=norm)
        k = k + 1
        
        
    if disp or full_output:
        fval = old_fval
    if warnflag == 2:
        if disp:
            print "Warning: Desired error not necessarily achieved due to precision loss"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls
        
    elif k >= maxiter:
        warnflag = 1
        if disp:
            print "Warning: Maximum number of iterations has been exceeded"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls
    else:
        if disp:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls


    if full_output:
        retlist = xk, fval, func_calls, grad_calls, warnflag
        if retall:
            retlist += (allvecs,)
    else: 
        retlist = xk
        if retall:
            retlist = (xk, allvecs)

    return retlist

def fmin_ncg(f, x0, fprime, fhess_p=None, fhess=None, args=(), avextol=1e-5,
             epsilon=_epsilon, maxiter=None, full_output=0, disp=1, retall=0):
    """Description:

    Minimize the function, f, whose gradient is given by fprime using the
    Newton-CG method.  fhess_p must compute the hessian times an arbitrary
    vector. If it is not given, finite-differences on fprime are used to
    compute it. See Wright, and Nocedal 'Numerical Optimization', 1999,
    pg. 140.

  Inputs:

    f -- the Python function or method to be minimized.
    x0 -- the initial guess for the minimizer.
    fprime -- a function to compute the gradient of f: fprime(x, *args)
    fhess_p -- a function to compute the Hessian of f times an
               arbitrary vector: fhess_p (x, p, *args)
    fhess -- a function to compute the Hessian matrix of f.
    args -- extra arguments for f, fprime, fhess_p, and fhess (the same
            set of extra arguments is supplied to all of these functions).

    epsilon -- if fhess is approximated use this value for
                 the step size (can be scalar or vector)

  Outputs: (xopt, {fopt, fcalls, gcalls, hcalls, warnflag},{allvecs})

    xopt -- the minimizer of f
    
    fopt -- the value of the function at xopt: fopt = f(xopt)
    fcalls -- the number of function calls.
    gcalls -- the number of gradient calls.
    hcalls -- the number of hessian calls.
    warnflag -- algorithm warnings:
                1 : 'Maximum number of iterations exceeded.'
    allvecs -- a list of all tried iterates

  Additional Inputs:

    avextol -- Convergence is assumed when the average relative error in
               the minimizer falls below this amount.  
    maxiter -- Maximum number of iterations to allow.
    full_output -- If non-zero return the optional outputs.
    disp -- If non-zero print convergence message.
    retall -- return a list of results at each iteration if True    
    
  Remarks:

    Only one of fhess_p or fhess need be given.  If fhess is provided,
    then fhess_p will be ignored.  If neither fhess nor fhess_p is
    provided, then the hessian product will be approximated using finite
    differences on fprime.

    """

    x0 = asarray(x0)
    fcalls = 0
    gcalls = 0
    hcalls = 0
    if maxiter is None:
        maxiter = len(x0)*200
    
    xtol = len(x0)*avextol
    update = [2*xtol]
    xk = x0
    if retall:
        allvecs = [xk]
    k = 0
    old_fval = f(x0,*args)
    fcalls += 1
    while (Num.add.reduce(abs(update)) > xtol) and (k < maxiter):
        # Compute a search direction pk by applying the CG method to
        #  del2 f(xk) p = - grad f(xk) starting from 0.
        b = -apply(fprime,(xk,)+args)
        gcalls = gcalls + 1
        maggrad = Num.add.reduce(abs(b))
        eta = min([0.5,Num.sqrt(maggrad)])
        termcond = eta * maggrad
        xsupi = 0
        ri = -b
        psupi = -ri
        i = 0
        dri0 = Num.dot(ri,ri)

        if fhess is not None:             # you want to compute hessian once.
            A = apply(fhess,(xk,)+args)
            hcalls = hcalls + 1

        while Num.add.reduce(abs(ri)) > termcond:
            if fhess is None:
                if fhess_p is None:
                    Ap = apply(approx_fhess_p,(xk,psupi,fprime,epsilon)+args)
                    gcalls = gcalls + 2
                else:
                    Ap = apply(fhess_p,(xk,psupi)+args)
                    hcalls = hcalls + 1
            else:
                Ap = Num.dot(A,psupi)
            # check curvature
            curv = Num.dot(psupi,Ap)
            if (curv <= 0):
                if (i > 0):
                    break
                else:
                    xsupi = xsupi + dri0/curv * psupi
                    break
            alphai = dri0 / curv
            xsupi = xsupi + alphai * psupi
            ri = ri + alphai * Ap
            dri1 = Num.dot(ri,ri)
            betai = dri1 / dri0
            psupi = -ri + betai * psupi
            i = i + 1
            dri0 = dri1          # update Num.dot(ri,ri) for next time.
    
        pk = xsupi  # search direction is solution to system.
        gfk = -b    # gradient at xk
        alphak, fc, gc, old_fval = line_search_BFGS(f,xk,pk,gfk,old_fval,args)
        fcalls = fcalls + fc
        gcalls = gcalls + gc

        update = alphak * pk
        xk = xk + update
        if retall:
            allvecs.append(xk)
        k = k + 1

    if disp or full_output:
        fval = old_fval
    if k >= maxiter:
        warnflag = 1
        if disp:
            print "Warning: Maximum number of iterations has been exceeded"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % fcalls
            print "         Gradient evaluations: %d" % gcalls
            print "         Hessian evaluations: %d" % hcalls
    else:
        warnflag = 0
        if disp:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % fcalls
            print "         Gradient evaluations: %d" % gcalls
            print "         Hessian evaluations: %d" % hcalls
            
    if full_output:
        retlist = xk, fval, fcalls, gcalls, hcalls, warnflag
        if retall:
            retlist += (allvecs,)
    else: 
        retlist = xk
        if retall:
            retlist = (xk, allvecs)

    return retlist
    

def fminbound(func, x1, x2, args=(), xtol=1e-5, maxfun=500, 
              full_output=0, disp=1):
    """Bounded minimization for scalar functions.

    Description:

      Finds a local minimizer of the scalar function func in the interval
      x1 < xopt < x2 using Brent's method.  (See brent for auto-bracketing).

    Inputs:

      func -- the function to be minimized (must accept scalar input and return
              scalar output).
      x1, x2 -- the optimization bounds.
      args -- extra arguments to pass to function.
      xtol -- the convergence tolerance.
      maxfun -- maximum function evaluations.
      full_output -- Non-zero to return optional outputs. 
      disp -- Non-zero to print messages.
              0 : no message printing.
              1 : non-convergence notification messages only.
              2 : print a message on convergence too.
              3 : print iteration results. 


    Outputs: (xopt, {fval, ierr, numfunc})

      xopt -- The minimizer of the function over the interval.
      fval -- The function value at the minimum point.
      ierr -- An error flag (0 if converged, 1 if maximum number of
              function calls reached).
      numfunc -- The number of function calls.
    """

    if x1 > x2:
        raise ValueError, "The lower bound exceeds the upper bound."

    flag = 0
    header = ' Func-count     x          f(x)          Procedure'
    step='       initial'

    sqrt_eps = sqrt(2.2e-16)
    golden_mean = 0.5*(3.0-sqrt(5.0))
    a, b = x1, x2
    fulc = a + golden_mean*(b-a)
    nfc, xf = fulc, fulc
    rat = e = 0.0
    x = xf
    fx = func(x,*args)
    num = 1
    fmin_data = (1, xf, fx)

    ffulc = fnfc = fx
    xm = 0.5*(a+b)
    tol1 = sqrt_eps*abs(xf) + xtol / 3.0
    tol2 = 2.0*tol1

    if disp > 2:
        print (" ")
        print (header)
        print "%5.0f   %12.6g %12.6g %s" % (fmin_data + (step,))
    

    while ( abs(xf-xm) > (tol2 - 0.5*(b-a)) ):
        golden = 1
        # Check for parabolic fit
        if abs(e) > tol1:
            golden = 0
            r = (xf-nfc)*(fx-ffulc)
            q = (xf-fulc)*(fx-fnfc)
            p = (xf-fulc)*q - (xf-nfc)*r
            q = 2.0*(q-r)
            if q > 0.0: p = -p
            q = abs(q)
            r = e
            e = rat

            # Check for acceptability of parabola
            if ( (abs(p) < abs(0.5*q*r)) and (p > q*(a-xf)) and (p < q*(b-xf))):
                rat = (p+0.0) / q;
                x = xf + rat
                step = '       parabolic'

                if ((x-a) < tol2) or ((b-x) < tol2):
                    si = Numeric.sign(xm-xf) + ((xm-xf)==0)
                    rat = tol1*si
            else:      # do a golden section step
                golden = 1

        if golden:  # Do a golden-section step
            if xf >= xm:
                e=a-xf
            else:
                e=b-xf
            rat = golden_mean*e
            step = '       golden'

        si = Numeric.sign(rat) + (rat == 0)
        x = xf + si*max([abs(rat), tol1])
        fu = func(x,*args)
        num += 1
        fmin_data = (num, x, fu)
        if disp > 2:
            print "%5.0f   %12.6g %12.6g %s" % (fmin_data + (step,))
                
        if fu <= fx:
            if x >= xf:
                a = xf
            else:
                b = xf
            fulc, ffulc = nfc, fnfc
            nfc, fnfc = xf, fx
            xf, fx = x, fu
        else:
            if x < xf:
                a = x
            else:
                b = x
            if (fu <= fnfc) or (nfc == xf):
                fulc, ffulc = nfc, fnfc
                nfc, fnfc = x, fu
            elif (fu <= ffulc) or (fulc == xf) or (fulc == nfc):
                fulc, ffulc = x, fu

        xm = 0.5*(a+b)
        tol1 = sqrt_eps*abs(xf) + xtol/3.0
        tol2 = 2.0*tol1

        if num >= maxfun:
            flag = 1
            fval = fx
            if disp > 0:
                _endprint(x, flag, fval, maxfun, tol, disp)
            if full_output:
                return xf, fval, flag, num
            else:
                return xf

    fval = fx
    if disp > 0:
        _endprint(x, flag, fval, maxfun, xtol, disp)
    
    if full_output:
        return xf, fval, flag, num
    else:
        return xf

def brent(func, args=(), brack=None, tol=1.48e-8, full_output=0, maxiter=500):
    """ Given a function of one-variable and a possible bracketing interval,
    return the minimum of the function isolated to a fractional precision of
    tol. A bracketing interval is a triple (a,b,c) where (a<b<c) and
    func(b) < func(a),func(c).  If bracket is two numbers then they are
    assumed to be a starting interval for a downhill bracket search
    (see bracket)

    Uses inverse parabolic interpolation when possible to speed up convergence
    of golden section method.

    """
    _mintol = 1.0e-11
    _cg = 0.3819660
    if brack is None:
        xa,xb,xc,fa,fb,fc,funcalls = bracket(func, args=args)
    elif len(brack) == 2:
        xa,xb,xc,fa,fb,fc,funcalls = bracket(func, xa=brack[0], xb=brack[1], args=args)
    elif len(brack) == 3:
        xa,xb,xc = brack
        if (xa > xc):  # swap so xa < xc can be assumed
            dum = xa; xa=xc; xc=dum
        assert ((xa < xb) and (xb < xc)), "Not a bracketing interval."
        fa = apply(func, (xa,)+args)
        fb = apply(func, (xb,)+args)
        fc = apply(func, (xc,)+args)
        assert ((fb<fa) and (fb < fc)), "Not a bracketing interval."
        funcalls = 3
    else:
        raise ValuError, "Bracketing interval must be length 2 or 3 sequence."

    x=w=v=xb
    fw=fv=fx=apply(func, (x,)+args)
    if (xa < xc):
        a = xa; b = xc
    else:
        a = xc; b = xa
    deltax= 0.0
    funcalls = 1
    iter = 0
    while (iter < maxiter):
        tol1 = tol*abs(x) + _mintol
        tol2 = 2.0*tol1
        xmid = 0.5*(a+b)
        if abs(x-xmid) < (tol2-0.5*(b-a)):  # check for convergence
            xmin=x; fval=fx
            break
        if (abs(deltax) <= tol1):           
            if (x>=xmid): deltax=a-x       # do a golden section step
            else: deltax=b-x
            rat = _cg*deltax
        else:                              # do a parabolic step
            tmp1 = (x-w)*(fx-fv)
            tmp2 = (x-v)*(fx-fw)
            p = (x-v)*tmp2 - (x-w)*tmp1;
            tmp2 = 2.0*(tmp2-tmp1)
            if (tmp2 > 0.0): p = -p
            tmp2 = abs(tmp2)
            dx_temp = deltax
            deltax= rat
            # check parabolic fit
            if ((p > tmp2*(a-x)) and (p < tmp2*(b-x)) and (abs(p) < abs(0.5*tmp2*dx_temp))):
                rat = p*1.0/tmp2        # if parabolic step is useful.
                u = x + rat
                if ((u-a) < tol2 or (b-u) < tol2):
                    if xmid-x >= 0: rat = tol1
                    else: rat = -tol1
            else:
                if (x>=xmid): deltax=a-x # if it's not do a golden section step
                else: deltax=b-x
                rat = _cg*deltax

        if (abs(rat) < tol1):            # update by at least tol1
            if rat >= 0: u = x + tol1
            else: u = x - tol1
        else:
            u = x + rat
        fu = apply(func, (u,)+args)      # calculate new output value
        funcalls += 1

        if (fu > fx):                 # if it's bigger than current
            if (u<x): a=u
            else: b=u
            if (fu<=fw) or (w==x):
                v=w; w=u; fv=fw; fw=fu
            elif (fu<=fv) or (v==x) or (v==w):
                v=u; fv=fu
        else: 
            if (u >= x): a = x
            else: b = x
            v=w; w=x; x=u
            fv=fw; fw=fx; fx=fu
        
    xmin = x
    fval = fx
    if full_output:
        return xmin, fval, iter, funcalls
    else:
        return xmin
    
def golden(func, args=(), brack=None, tol=_epsilon, full_output=0):
    """ Given a function of one-variable and a possible bracketing interval,
    return the minimum of the function isolated to a fractional precision of
    tol. A bracketing interval is a triple (a,b,c) where (a<b<c) and
    func(b) < func(a),func(c).  If bracket is two numbers then they are
    assumed to be a starting interval for a downhill bracket search
    (see bracket)

    Uses analog of bisection method to decrease the bracketed interval.
    """
    if brack is None:
        xa,xb,xc,fa,fb,fc,funcalls = bracket(func, args=args)
    elif len(brack) == 2:
        xa,xb,xc,fa,fb,fc,funcalls = bracket(func, xa=brack[0], xb=brack[1], args=args)
    elif len(brack) == 3:
        xa,xb,xc = brack
        if (xa > xc):  # swap so xa < xc can be assumed
            dum = xa; xa=xc; xc=dum
        assert ((xa < xb) and (xb < xc)), "Not a bracketing interval."
        fa = apply(func, (xa,)+args)
        fb = apply(func, (xb,)+args)
        fc = apply(func, (xc,)+args)
        assert ((fb<fa) and (fb < fc)), "Not a bracketing interval."
        funcalls = 3
    else:
        raise ValuError, "Bracketing interval must be length 2 or 3 sequence."

    _gR = 0.61803399
    _gC = 1.0-_gR
    x3 = xc
    x0 = xa
    if (abs(xc-xb) > abs(xb-xa)):
        x1 = xb
        x2 = xb + _gC*(xc-xb)
    else:
        x2 = xb
        x1 = xb - _gC*(xb-xa)
    f1 = apply(func, (x1,)+args)
    f2 = apply(func, (x2,)+args)
    funcalls += 2
    while (abs(x3-x0) > tol*(abs(x1)+abs(x2))):
        if (f2 < f1):
            x0 = x1; x1 = x2; x2 = _gR*x1 + _gC*x3
            f1 = f2; f2 = apply(func, (x2,)+args)
        else:
            x3 = x2; x2 = x1; x1 = _gR*x2 + _gC*x0
            f2 = f1; f1 = apply(func, (x1,)+args)
        funcalls += 1
    if (f1 < f2):
        xmin = x1
        fval = f1
    else:
        xmin = x2
        fval = f2
    if full_output:
        return xmin, fval, funcalls
    else:
        return xmin    


def bracket(func, xa=0.0, xb=1.0, args=(), grow_limit=110.0):
    """Given a function and distinct initial points, search in the downhill
    direction (as defined by the initital points) and return new points
    xa, xb, xc that bracket the minimum of the function:
    f(xa) > f(xb) < f(xc)
    """
    _gold = 1.618034
    _verysmall_num = 1e-21
    fa = apply(func, (xa,)+args)
    fb = apply(func, (xb,)+args)
    if (fa < fb):                      # Switch so fa > fb 
        dum = xa; xa = xb; xb = dum
        dum = fa; fa = fb; fb = dum
    xc = xb + _gold*(xb-xa)
    fc = apply(func, (xc,)+args)
    funcalls = 3
    iter = 0
    while (fc < fb):
        tmp1 = (xb - xa)*(fb-fc)
        tmp2 = (xb - xc)*(fb-fa)
        val = tmp2-tmp1
        if abs(val) < _verysmall_num:
            denom = 2.0*_verysmall_num
        else:
            denom = 2.0*val
        w = xb - ((xb-xc)*tmp2-(xb-xa)*tmp1)/denom
        wlim = xb + grow_limit*(xc-xb)
        if iter > 1000:
            raise RunTimeError, "Too many iterations."
        if (w-xc)*(xb-w) > 0.0:
            fw = apply(func, (w,)+args)
            funcalls += 1
            if (fw < fc):
                xa = xb; xb=w; fa=fb; fb=fw
                return xa, xb, xc, fa, fb, fc, funcalls
            elif (fw > fb):
                xc = w; fc=fw
                return xa, xb, xc, fa, fb, fc, funcalls
            w = xc + _gold*(xc-xb)
            fw = apply(func, (w,)+args)
            funcalls += 1
        elif (w-wlim)*(wlim-xc) >= 0.0:
            w = wlim
            fw = apply(func, (w,)+args)
            funcalls += 1
        elif (w-wlim)*(xc-w) > 0.0:
            fw = apply(func, (w,)+args)
            funcalls += 1
            if (fw < fc):
                xb=xc; xc=w; w=xc+_gold*(xc-xb)
                fb=fc; fc=fw; fw=apply(func, (w,)+args)
                funcalls += 1
        else:
            w = xc + _gold*(xc-xb)
            fw = apply(func, (w,)+args)
            funcalls += 1
        xa=xb; xb=xc; xc=w
        fa=fb; fb=fc; fc=fw
    return xa, xb, xc, fa, fb, fc, funcalls
            
            
global _powell_funcalls

def _myfunc(alpha, func, x0, direc, args=()):
    funcargs = (x0 + alpha * direc,)+args
    return func(*funcargs)
    
def _linesearch_powell(func, p, xi, args=(), tol=1e-3):
    # line-search algorithm using fminbound
    #  find the minimium of the function
    #  func(x0+ alpha*direc)
    global _powell_funcalls
    extra_args = (func, p, xi) + args
    alpha_min, fret, iter, num = brent(_myfunc, args=extra_args,
                                           full_output=1, tol=tol)
    xi = alpha_min*xi
    _powell_funcalls += num
    return squeeze(fret), p+xi, xi
    

def fmin_powell(func, x0, args=(), xtol=1e-4, ftol=1e-4, maxiter=None,
                maxfun=None, full_output=0, disp=1, retall=0):
    """Minimize a function using modified Powell's method.

    Description:
    
      Uses a modification of Powell's method to find the minimum of a function
      of N variables

    Inputs:

      func -- the Python function or method to be minimized.
      x0 -- the initial guess.
      args -- extra arguments for func.

    Outputs: (xopt, {fopt, xi, direc, iter, funcalls, warnflag}, {allvecs})

      xopt -- minimizer of function

      fopt  -- value of function at minimum: fopt = func(xopt)
      direc -- current direction set
      iter -- number of iterations
      funcalls -- number of function calls
      warnflag -- Integer warning flag:
                  1 : 'Maximum number of function evaluations.'
                  2 : 'Maximum number of iterations.'
      allvecs -- a list of solutions at each iteration

    Additional Inputs:

      xtol -- line-search error tolerance.
      ftol -- acceptable relative error in func(xopt) for convergence.
      maxiter -- the maximum number of iterations to perform.
      maxfun -- the maximum number of function evaluations.
      full_output -- non-zero if fval and warnflag outputs are desired.
      disp -- non-zero to print convergence messages.
      retall -- non-zero to return a list of the solution at each iteration

      """
    global _powell_funcalls
    x = asarray(x0)
    if retall:
        allvecs = [x]
    N = len(x)
    rank = len(x.shape)
    if not -1 < rank < 2:
        raise ValueError, "Initial guess must be a scalar or rank-1 sequence."
    if maxiter is None:
        maxiter = N * 1000
    if maxfun is None:
        maxfun = N * 1000

    direc = eye(N,typecode='d')
    fval = squeeze(apply(func, (x,)+args))
    _powell_funcalls = 1
    x1 = x.copy()
    iter = 0;
    ilist = range(N)
    while 1:
        fx = fval
        bigind = 0
        delta = 0.0
        for i in ilist:
            direc1 = direc[i]
            fx2 = fval
            fval, x, direc1 = _linesearch_powell(func, x, direc1, args=args, tol=xtol*100)
            if (fx2 - fval) > delta:
                delta = fx2 - fval
                bigind = i
        iter += 1
        if retall:
            allvecs.append(x)
        if (2.0*(fx - fval) <= ftol*(abs(fx)+abs(fval))+1e-20): break
        if _powell_funcalls >= maxfun: break
        if iter >= maxiter: break

        # Construct the extrapolated point
        direc1 = x - x1
        x2 = 2*x - x1
        x1 = x.copy()
        fx2 = squeeze(apply(func, (x2,)+args))
        _powell_funcalls +=1

        if (fx > fx2):
            t = 2.0*(fx+fx2-2.0*fval)
            temp = (fx-fval-delta)
            t *= temp*temp
            temp = fx-fx2
            t -= delta*temp*temp
            if t < 0.0:
                fval, x, direc1 = _linesearch_powell(func, x, direc1, args=args, tol=xtol*100)
                direc[bigind] = direc[-1]
                direc[-1] = direc1            

    warnflag = 0
    if _powell_funcalls >= maxfun:
        warnflag = 1
        if disp:
            print "Warning: Maximum number of function evaluations has "\
                  "been exceeded."
    elif iter >= maxiter:
        warnflag = 2
        if disp:
            print "Warning: Maximum number of iterations has been exceeded"
    else:
        if disp:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % iter
            print "         Function evaluations: %d" % _powell_funcalls

    x = squeeze(x)

    if full_output:
        retlist = x, fval, direc, iter, _powell_funcalls, warnflag
        if retall:
            retlist += (allvecs,)
    else: 
        retlist = x
        if retall:
            retlist = (x, allvecs)

    return retlist
    
    


def _endprint(x, flag, fval, maxfun, xtol, disp):
    if flag == 0:
        if disp > 1:
            print "\nOptimization terminated successfully;\n the returned value" + \
                  " satisfies the termination criteria\n (using xtol = ", xtol, ")"
    if flag == 1:
        print "\nMaximum number of function evaluations exceeded --- increase maxfun argument.\n"
    return


def brute(func, ranges, args=(), Ns=20, full_output=0, finish=fmin):
    """Minimize a function over a given range by brute force.

    That is find the minimum of a function evaluated on a grid
    given by the tuple ranges.

    Inputs:

    func        -- Function to be optimized
    ranges       -- Tuple where each element is a tuple of parameters
                      or a slice object to be handed to scipy.mgrid

    args        -- Extra arguments to function.
    Ns          -- Default number of samples if not given
    full_output -- Nonzero to return evaluation grid.

    Outputs: (x0, fval, {grid, Jout})

    x0          -- Value of arguments giving minimum over the grird
    fval        -- Function value at minimum
    grid        -- tuple with same length as x0 representing the
                    evaluation grid
    Jout        -- Function values over grid:  Jout = func(*grid)
    """
    N = len(ranges)
    if N > 40:
        raise ValueError, "Brute Force not possible with more than 40 variables."
    lrange = list(ranges)
    for k in range(N):
        if type(lrange[k]) is not type(slice(None)):
            if len(lrange[k]) < 3:
                lrange[k] = tuple(lrange[k]) + (complex(Ns),)
            lrange[k] = slice(*lrange[k])
    if (N==1):
        lrange = lrange[0]
        
    def _scalarfunc(*params):
        params = squeeze(asarray(params))
        return func(params,*args)
        
    vecfunc = vectorize(_scalarfunc)
    grid = mgrid[lrange]
    if (N==1):
        grid = (grid,)
    Jout = vecfunc(*grid)
    Nshape = shape(Jout)
    indx = argmin(Jout.flat)
    Nindx = zeros(N)
    xmin = zeros(N,'d')
    for k in range(N-1,-1,-1):
        thisN = Nshape[k]
        Nindx[k] = indx % Nshape[k]
        indx = indx / thisN
    for k in range(N):
        xmin[k] = grid[k][tuple(Nindx)]

    Jmin = Jout[tuple(Nindx)]
    if (N==1):
        grid = grid[0]
        xmin = xmin[0]
    if callable(finish):
        vals = finish(func,xmin,args=args,full_output=1, disp=0)
        xmin = vals[0]
        Jmin = vals[1]
        if vals[-1] > 0:
            print "Warning: Final optimization did not succeed"        
    if full_output:
        return xmin, Jmin, grid, Jout
    else:
        return xmin

            
if __name__ == "__main__":
    import string
    import time

    
    times = []
    algor = []
    x0 = [0.8,1.2,0.7]
    print "Nelder-Mead Simplex"
    print "==================="
    start = time.time()
    x = fmin(rosen,x0)
    print x
    times.append(time.time() - start)
    algor.append('Nelder-Mead Simplex\t')

    print
    print "Powell Direction Set Method"
    print "==========================="
    start = time.time()
    x = fmin_powell(rosen,x0)
    print x
    times.append(time.time() - start)
    algor.append('Powell Direction Set Method.')

    print
    print "Nonlinear CG"
    print "============"
    start = time.time()
    x = fmin_cg(rosen, x0, fprime=rosen_der, maxiter=200)
    print x
    times.append(time.time() - start)
    algor.append('Nonlinear CG     \t')

    print
    print "BFGS Quasi-Newton"
    print "================="
    start = time.time()
    x = fmin_bfgs(rosen, x0, fprime=rosen_der, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('BFGS Quasi-Newton\t')

    print
    print "BFGS approximate gradient"
    print "========================="
    start = time.time()
    x = fmin_bfgs(rosen, x0, gtol=1e-4, maxiter=100)
    print x
    times.append(time.time() - start)
    algor.append('BFGS without gradient\t')


    print
    print "Newton-CG with Hessian product"
    print "=============================="
    start = time.time()
    x = fmin_ncg(rosen, x0, rosen_der, fhess_p=rosen_hess_prod, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('Newton-CG with hessian product')
    

    print
    print "Newton-CG with full Hessian"
    print "==========================="
    start = time.time()
    x = fmin_ncg(rosen, x0, rosen_der, fhess=rosen_hess, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('Newton-CG with full hessian')

    print
    print "\nMinimizing the Rosenbrock function of order 3\n"
    print " Algorithm \t\t\t       Seconds"
    print "===========\t\t\t      ========="
    for k in range(len(algor)):
        print algor[k], "\t -- ", times[k]
        














