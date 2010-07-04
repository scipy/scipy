#__docformat__ = "restructuredtext en"
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

__all__ = ['fmin', 'fmin_powell','fmin_bfgs', 'fmin_ncg', 'fmin_cg',
           'fminbound','brent', 'golden','bracket','rosen','rosen_der',
           'rosen_hess', 'rosen_hess_prod', 'brute', 'approx_fprime',
           'line_search', 'check_grad']

__docformat__ = "restructuredtext en"

import numpy
from numpy import atleast_1d, eye, mgrid, argmin, zeros, shape, empty, \
     squeeze, vectorize, asarray, absolute, sqrt, Inf, asfarray, isinf
import linesearch

# These have been copied from Numeric's MLab.py
# I don't think they made the transition to scipy_core
def max(m,axis=0):
    """max(m,axis=0) returns the maximum of m along dimension axis.
    """
    m = asarray(m)
    return numpy.maximum.reduce(m,axis)

def min(m,axis=0):
    """min(m,axis=0) returns the minimum of m along dimension axis.
    """
    m = asarray(m)
    return numpy.minimum.reduce(m,axis)

def is_array_scalar(x):
    """Test whether `x` is either a scalar or an array scalar.

    """
    return len(atleast_1d(x) == 1)

abs = absolute
import __builtin__
pymin = __builtin__.min
pymax = __builtin__.max
__version__="0.7"
_epsilon = sqrt(numpy.finfo(float).eps)

def vecnorm(x, ord=2):
    if ord == Inf:
        return numpy.amax(abs(x))
    elif ord == -Inf:
        return numpy.amin(abs(x))
    else:
        return numpy.sum(abs(x)**ord,axis=0)**(1.0/ord)

def rosen(x):  # The Rosenbrock function
    x = asarray(x)
    return numpy.sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0,axis=0)

def rosen_der(x):
    x = asarray(x)
    xm = x[1:-1]
    xm_m1 = x[:-2]
    xm_p1 = x[2:]
    der = numpy.zeros_like(x)
    der[1:-1] = 200*(xm-xm_m1**2) - 400*(xm_p1 - xm**2)*xm - 2*(1-xm)
    der[0] = -400*x[0]*(x[1]-x[0]**2) - 2*(1-x[0])
    der[-1] = 200*(x[-1]-x[-2]**2)
    return der

def rosen_hess(x):
    x = atleast_1d(x)
    H = numpy.diag(-400*x[:-1],1) - numpy.diag(400*x[:-1],-1)
    diagonal = numpy.zeros(len(x), dtype=x.dtype)
    diagonal[0] = 1200*x[0]-400*x[1]+2
    diagonal[-1] = 200
    diagonal[1:-1] = 202 + 1200*x[1:-1]**2 - 400*x[2:]
    H = H + numpy.diag(diagonal)
    return H

def rosen_hess_prod(x,p):
    x = atleast_1d(x)
    Hp = numpy.zeros(len(x), dtype=x.dtype)
    Hp[0] = (1200*x[0]**2 - 400*x[1] + 2)*p[0] - 400*x[0]*p[1]
    Hp[1:-1] = -400*x[:-2]*p[:-2]+(202+1200*x[1:-1]**2-400*x[2:])*p[1:-1] \
               -400*x[1:-1]*p[2:]
    Hp[-1] = -400*x[-2]*p[-2] + 200*p[-1]
    return Hp

def wrap_function(function, args):
    ncalls = [0]
    def function_wrapper(x):
        ncalls[0] += 1
        return function(x, *args)
    return ncalls, function_wrapper

def fmin(func, x0, args=(), xtol=1e-4, ftol=1e-4, maxiter=None, maxfun=None,
         full_output=0, disp=1, retall=0, callback=None):
    """Minimize a function using the downhill simplex algorithm.

    Parameters
    ----------
    func : callable func(x,*args)
        The objective function to be minimized.
    x0 : ndarray
        Initial guess.
    args : tuple
        Extra arguments passed to func, i.e. ``f(x,*args)``.
    callback : callable
        Called after each iteration, as callback(xk), where xk is the
        current parameter vector.

    Returns
    -------
    xopt : ndarray
        Parameter that minimizes function.
    fopt : float
        Value of function at minimum: ``fopt = func(xopt)``.
    iter : int
        Number of iterations performed.
    funcalls : int
        Number of function calls made.
    warnflag : int
        1 : Maximum number of function evaluations made.
        2 : Maximum number of iterations reached.
    allvecs : list
        Solution at each iteration.

    Other Parameters
    ----------------
    xtol : float
        Relative error in xopt acceptable for convergence.
    ftol : number
        Relative error in func(xopt) acceptable for convergence.
    maxiter : int
        Maximum number of iterations to perform.
    maxfun : number
        Maximum number of function evaluations to make.
    full_output : bool
        Set to True if fval and warnflag outputs are desired.
    disp : bool
        Set to True to print convergence messages.
    retall : bool
        Set to True to return list of solutions at each iteration.

    Notes
    -----
    Uses a Nelder-Mead simplex algorithm to find the minimum of
    a function of one or more variables.

    """
    fcalls, func = wrap_function(func, args)
    x0 = asfarray(x0).flatten()
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
        sim = numpy.zeros((N+1,), dtype=x0.dtype)
    else:
        sim = numpy.zeros((N+1,N), dtype=x0.dtype)
    fsim = numpy.zeros((N+1,), float)
    sim[0] = x0
    if retall:
        allvecs = [sim[0]]
    fsim[0] = func(x0)
    nonzdelt = 0.05
    zdelt = 0.00025
    for k in range(0,N):
        y = numpy.array(x0,copy=True)
        if y[k] != 0:
            y[k] = (1+nonzdelt)*y[k]
        else:
            y[k] = zdelt

        sim[k+1] = y
        f = func(y)
        fsim[k+1] = f

    ind = numpy.argsort(fsim)
    fsim = numpy.take(fsim,ind,0)
    # sort so sim[0,:] has the lowest function value
    sim = numpy.take(sim,ind,0)

    iterations = 1

    while (fcalls[0] < maxfun and iterations < maxiter):
        if (max(numpy.ravel(abs(sim[1:]-sim[0]))) <= xtol \
            and max(abs(fsim[0]-fsim[1:])) <= ftol):
            break

        xbar = numpy.add.reduce(sim[:-1],0) / N
        xr = (1+rho)*xbar - rho*sim[-1]
        fxr = func(xr)
        doshrink = 0

        if fxr < fsim[0]:
            xe = (1+rho*chi)*xbar - rho*chi*sim[-1]
            fxe = func(xe)

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
                    fxc = func(xc)

                    if fxc <= fxr:
                        sim[-1] = xc
                        fsim[-1] = fxc
                    else:
                        doshrink=1
                else:
                    # Perform an inside contraction
                    xcc = (1-psi)*xbar + psi*sim[-1]
                    fxcc = func(xcc)

                    if fxcc < fsim[-1]:
                        sim[-1] = xcc
                        fsim[-1] = fxcc
                    else:
                        doshrink = 1

                if doshrink:
                    for j in one2np1:
                        sim[j] = sim[0] + sigma*(sim[j] - sim[0])
                        fsim[j] = func(sim[j])

        ind = numpy.argsort(fsim)
        sim = numpy.take(sim,ind,0)
        fsim = numpy.take(fsim,ind,0)
        if callback is not None:
            callback(sim[0])
        iterations += 1
        if retall:
            allvecs.append(sim[0])

    x = sim[0]
    fval = min(fsim)
    warnflag = 0

    if fcalls[0] >= maxfun:
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
            print "         Function evaluations: %d" % fcalls[0]


    if full_output:
        retlist = x, fval, iterations, fcalls[0], warnflag
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
    d1 = empty((2,2))
    d1[0,0] = dc**2
    d1[0,1] = -db**2
    d1[1,0] = -dc**3
    d1[1,1] = db**3
    [A,B] = numpy.dot(d1,asarray([fb-fa-C*db,fc-fa-C*dc]).flatten())
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
        # interpolate to find a trial step length between a_lo and
        # a_hi Need to choose interpolation here.  Use cubic
        # interpolation and then if the result is within delta *
        # dalpha or outside of the interval bounded by a_lo or a_hi
        # then use quadratic interpolation, if the result is still too
        # close, then use bisection

        dalpha = a_hi-a_lo;
        if dalpha < 0: a,b = a_hi,a_lo
        else: a,b = a_lo, a_hi

        # minimizer of cubic interpolant
        #    (uses phi_lo, derphi_lo, phi_hi, and the most recent value of phi)
        #      if the result is too close to the end points (or out of
        #        the interval) then use quadratic interpolation with
        #        phi_lo, derphi_lo and phi_hi
        #      if the result is stil too close to the end points (or
        #        out of the interval) then use bisection

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

    Parameters
    ----------
    f : callable f(x,*args)
        Objective function.
    myfprime : callable f'(x,*args)
        Objective function gradient (can be None).
    xk : ndarray
        Starting point.
    pk : ndarray
        Search direction.
    gfk : ndarray
        Gradient value for x=xk (xk being the current parameter
        estimate).
    args : tuple
        Additional arguments passed to objective function.
    c1 : float
        Parameter for Armijo condition rule.
    c2 : float
        Parameter for curvature condition rule.

    Returns
    -------
    alpha0 : float
        Alpha for which ``x_new = x0 + alpha * pk``.
    fc : int
        Number of function evaluations made.
    gc : int
        Number of gradient evaluations made.

    Notes
    -----
    Uses the line search algorithm to enforce strong Wolfe
    conditions.  See Wright and Nocedal, 'Numerical Optimization',
    1999, pg. 59-60.

    For the zoom phase it uses an algorithm by [...].

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
            return numpy.dot(_ls_ingfk,pk)
    else:
        fprime = myfprime
        def phiprime(alpha):
            global _ls_gc, _ls_ingfk
            _ls_gc += 1
            _ls_ingfk = fprime(xk+alpha*pk,*args)  # store for later use
            return numpy.dot(_ls_ingfk,pk)

    alpha0 = 0
    phi0 = old_fval
    derphi0 = numpy.dot(gfk,pk)

    alpha1 = pymin(1.0,1.01*2*(phi0-old_old_fval)/derphi0)

    if alpha1 == 0:
        # This shouldn't happen. Perhaps the increment has slipped below
        # machine precision?  For now, set the return variables skip the
        # useless while loop, and raise warnflag=2 due to possible imprecision.
        alpha_star = None
        fval_star = old_fval
        old_fval = old_old_fval
        fprime_star = None

    phi_a1 = phi(alpha1)
    #derphi_a1 = phiprime(alpha1)  evaluated below

    phi_a0 = phi0
    derphi_a0 = derphi0

    i = 1
    maxiter = 10
    while 1:         # bracketing phase
        if alpha1 == 0:
            break
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
        # fprime_star is a number (derphi) -- so use the most recently
        # calculated gradient used in computing it derphi = gfk*pk
        # this is the gradient at the next step no need to compute it
        # again in the outer loop.
        fprime_star = _ls_ingfk

    return alpha_star, _ls_fc, _ls_gc, fval_star, old_fval, fprime_star


def line_search_BFGS(f, xk, pk, gfk, old_fval, args=(), c1=1e-4, alpha0=1):
    """Minimize over alpha, the function ``f(xk+alpha pk)``.

    Uses the interpolation algorithm (Armiijo backtracking) as suggested by
    Wright and Nocedal in 'Numerical Optimization', 1999, pg. 56-57

    :Returns: (alpha, fc, gc)

    """

    xk = atleast_1d(xk)
    fc = 0
    phi0 = old_fval # compute f(xk) -- done in past loop
    phi_a0 = f(*((xk+alpha0*pk,)+args))
    fc = fc + 1
    derphi0 = numpy.dot(gfk,pk)

    if (phi_a0 <= phi0 + c1*alpha0*derphi0):
        return alpha0, fc, 0, phi_a0

    # Otherwise compute the minimizer of a quadratic interpolant:

    alpha1 = -(derphi0) * alpha0**2 / 2.0 / (phi_a0 - phi0 - derphi0 * alpha0)
    phi_a1 = f(*((xk+alpha1*pk,)+args))
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

        alpha2 = (-b + numpy.sqrt(abs(b**2 - 3 * a * derphi0))) / (3.0*a)
        phi_a2 = f(*((xk+alpha2*pk,)+args))
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
    f0 = f(*((xk,)+args))
    grad = numpy.zeros((len(xk),), float)
    ei = numpy.zeros((len(xk),), float)
    for k in range(len(xk)):
        ei[k] = epsilon
        grad[k] = (f(*((xk+ei,)+args)) - f0)/epsilon
        ei[k] = 0.0
    return grad

def check_grad(func, grad, x0, *args):
    return sqrt(sum((grad(x0,*args)-approx_fprime(x0,func,_epsilon,*args))**2))

def approx_fhess_p(x0,p,fprime,epsilon,*args):
    f2 = fprime(*((x0+epsilon*p,)+args))
    f1 = fprime(*((x0,)+args))
    return (f2 - f1)/epsilon


def fmin_bfgs(f, x0, fprime=None, args=(), gtol=1e-5, norm=Inf,
              epsilon=_epsilon, maxiter=None, full_output=0, disp=1,
              retall=0, callback=None):
    """Minimize a function using the BFGS algorithm.

    Parameters
    ----------
    f : callable f(x,*args)
        Objective function to be minimized.
    x0 : ndarray
        Initial guess.
    fprime : callable f'(x,*args)
        Gradient of f.
    args : tuple
        Extra arguments passed to f and fprime.
    gtol : float
        Gradient norm must be less than gtol before succesful termination.
    norm : float
        Order of norm (Inf is max, -Inf is min)
    epsilon : int or ndarray
        If fprime is approximated, use this value for the step size.
    callback : callable
        An optional user-supplied function to call after each
        iteration.  Called as callback(xk), where xk is the
        current parameter vector.

    Returns
    -------
    xopt : ndarray
        Parameters which minimize f, i.e. f(xopt) == fopt.
    fopt : float
        Minimum value.
    gopt : ndarray
        Value of gradient at minimum, f'(xopt), which should be near 0.
    Bopt : ndarray
        Value of 1/f''(xopt), i.e. the inverse hessian matrix.
    func_calls : int
        Number of function_calls made.
    grad_calls : int
        Number of gradient calls made.
    warnflag : integer
        1 : Maximum number of iterations exceeded.
        2 : Gradient and/or function calls not changing.
    allvecs  :  list
        Results at each iteration.  Only returned if retall is True.

    Other Parameters
    ----------------
    maxiter : int
        Maximum number of iterations to perform.
    full_output : bool
        If True,return fopt, func_calls, grad_calls, and warnflag
        in addition to xopt.
    disp : bool
        Print convergence message if True.
    retall : bool
        Return a list of results at each iteration if True.

    Notes
    -----
    Optimize the function, f, whose gradient is given by fprime
    using the quasi-Newton method of Broyden, Fletcher, Goldfarb,
    and Shanno (BFGS) See Wright, and Nocedal 'Numerical
    Optimization', 1999, pg. 198.

    """
    x0 = asarray(x0).squeeze()
    if x0.ndim == 0:
        x0.shape = (1,)
    if maxiter is None:
        maxiter = len(x0)*200
    func_calls, f = wrap_function(f, args)
    if fprime is None:
        grad_calls, myfprime = wrap_function(approx_fprime, (f, epsilon))
    else:
        grad_calls, myfprime = wrap_function(fprime, args)
    gfk = myfprime(x0)
    k = 0
    N = len(x0)
    I = numpy.eye(N,dtype=int)
    Hk = I
    old_fval = f(x0)
    old_old_fval = old_fval + 5000
    xk = x0
    if retall:
        allvecs = [x0]
    sk = [2*gtol]
    warnflag = 0
    gnorm = vecnorm(gfk,ord=norm)
    while (gnorm > gtol) and (k < maxiter):
        pk = -numpy.dot(Hk,gfk)
        alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = \
           linesearch.line_search(f,myfprime,xk,pk,gfk,
                                  old_fval,old_old_fval)
        if alpha_k is None:  # line search failed try different one.
            alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = \
                     line_search(f,myfprime,xk,pk,gfk,
                                 old_fval,old_old_fval)
            if alpha_k is None:
                # This line search also failed to find a better solution.
                warnflag = 2
                break
        xkp1 = xk + alpha_k * pk
        if retall:
            allvecs.append(xkp1)
        sk = xkp1 - xk
        xk = xkp1
        if gfkp1 is None:
            gfkp1 = myfprime(xkp1)

        yk = gfkp1 - gfk
        gfk = gfkp1
        if callback is not None:
            callback(xk)
        k += 1
        gnorm = vecnorm(gfk,ord=norm)
        if (gnorm <= gtol):
            break

        try: # this was handled in numeric, let it remaines for more safety
            rhok = 1.0 / (numpy.dot(yk,sk))
        except ZeroDivisionError:
            rhok = 1000.0
            print "Divide-by-zero encountered: rhok assumed large"
        if isinf(rhok): # this is patch for numpy
            rhok = 1000.0
            print "Divide-by-zero encountered: rhok assumed large"
        A1 = I - sk[:,numpy.newaxis] * yk[numpy.newaxis,:] * rhok
        A2 = I - yk[:,numpy.newaxis] * sk[numpy.newaxis,:] * rhok
        Hk = numpy.dot(A1,numpy.dot(Hk,A2)) + rhok * sk[:,numpy.newaxis] \
                 * sk[numpy.newaxis,:]

    if disp or full_output:
        fval = old_fval
    if warnflag == 2:
        if disp:
            print "Warning: Desired error not necessarily achieved" \
                  "due to precision loss"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls[0]
            print "         Gradient evaluations: %d" % grad_calls[0]

    elif k >= maxiter:
        warnflag = 1
        if disp:
            print "Warning: Maximum number of iterations has been exceeded"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls[0]
            print "         Gradient evaluations: %d" % grad_calls[0]
    else:
        if disp:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls[0]
            print "         Gradient evaluations: %d" % grad_calls[0]

    if full_output:
        retlist = xk, fval, gfk, Hk, func_calls[0], grad_calls[0], warnflag
        if retall:
            retlist += (allvecs,)
    else:
        retlist = xk
        if retall:
            retlist = (xk, allvecs)

    return retlist


def fmin_cg(f, x0, fprime=None, args=(), gtol=1e-5, norm=Inf, epsilon=_epsilon,
              maxiter=None, full_output=0, disp=1, retall=0, callback=None):
    """Minimize a function using a nonlinear conjugate gradient algorithm.

    Parameters
    ----------
    f : callable f(x,*args)
        Objective function to be minimized.
    x0 : ndarray
        Initial guess.
    fprime : callable f'(x,*args)
        Function which computes the gradient of f.
    args : tuple
        Extra arguments passed to f and fprime.
    gtol : float
        Stop when norm of gradient is less than gtol.
    norm : float
        Order of vector norm to use.  -Inf is min, Inf is max.
    epsilon : float or ndarray
        If fprime is approximated, use this value for the step
        size (can be scalar or vector).
    callback : callable
        An optional user-supplied function, called after each
        iteration.  Called as callback(xk), where xk is the
        current parameter vector.

    Returns
    -------
    xopt : ndarray
        Parameters which minimize f, i.e. f(xopt) == fopt.
    fopt : float
        Minimum value found, f(xopt).
    func_calls : int
        The number of function_calls made.
    grad_calls : int
        The number of gradient calls made.
    warnflag : int
        1 : Maximum number of iterations exceeded.
        2 : Gradient and/or function calls not changing.
    allvecs : ndarray
        If retall is True (see other parameters below), then this
        vector containing the result at each iteration is returned.

    Other Parameters
    ----------------
    maxiter : int
        Maximum number of iterations to perform.
    full_output : bool
        If True then return fopt, func_calls, grad_calls, and
        warnflag in addition to xopt.
    disp : bool
        Print convergence message if True.
    retall : bool
        Return a list of results at each iteration if True.

    Notes
    -----
    Optimize the function, f, whose gradient is given by fprime
    using the nonlinear conjugate gradient algorithm of Polak and
    Ribiere. See Wright & Nocedal, 'Numerical Optimization',
    1999, pg. 120-122.

    """
    x0 = asarray(x0).flatten()
    if maxiter is None:
        maxiter = len(x0)*200
    func_calls, f = wrap_function(f, args)
    if fprime is None:
        grad_calls, myfprime = wrap_function(approx_fprime, (f, epsilon))
    else:
        grad_calls, myfprime = wrap_function(fprime, args)
    gfk = myfprime(x0)
    k = 0
    N = len(x0)
    xk = x0
    old_fval = f(xk)
    old_old_fval = old_fval + 5000

    if retall:
        allvecs = [xk]
    sk = [2*gtol]
    warnflag = 0
    pk = -gfk
    gnorm = vecnorm(gfk,ord=norm)
    while (gnorm > gtol) and (k < maxiter):
        deltak = numpy.dot(gfk,gfk)

        # These values are modified by the line search, even if it fails
        old_fval_backup = old_fval
        old_old_fval_backup = old_old_fval

        alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = \
           linesearch.line_search(f,myfprime,xk,pk,gfk,old_fval,
                                  old_old_fval,c2=0.4)
        if alpha_k is None:  # line search failed -- use different one.
            alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = \
                     line_search(f,myfprime,xk,pk,gfk,
                                 old_fval_backup,old_old_fval_backup)
            if alpha_k is None or alpha_k == 0:
                # This line search also failed to find a better solution.
                warnflag = 2
                break
        xk = xk + alpha_k*pk
        if retall:
            allvecs.append(xk)
        if gfkp1 is None:
            gfkp1 = myfprime(xk)
        yk = gfkp1 - gfk
        beta_k = pymax(0,numpy.dot(yk,gfkp1)/deltak)
        pk = -gfkp1 + beta_k * pk
        gfk = gfkp1
        gnorm = vecnorm(gfk,ord=norm)
        if callback is not None:
            callback(xk)
        k += 1


    if disp or full_output:
        fval = old_fval
    if warnflag == 2:
        if disp:
            print "Warning: Desired error not necessarily achieved due to precision loss"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls[0]
            print "         Gradient evaluations: %d" % grad_calls[0]

    elif k >= maxiter:
        warnflag = 1
        if disp:
            print "Warning: Maximum number of iterations has been exceeded"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls[0]
            print "         Gradient evaluations: %d" % grad_calls[0]
    else:
        if disp:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls[0]
            print "         Gradient evaluations: %d" % grad_calls[0]


    if full_output:
        retlist = xk, fval, func_calls[0], grad_calls[0], warnflag
        if retall:
            retlist += (allvecs,)
    else:
        retlist = xk
        if retall:
            retlist = (xk, allvecs)

    return retlist

def fmin_ncg(f, x0, fprime, fhess_p=None, fhess=None, args=(), avextol=1e-5,
             epsilon=_epsilon, maxiter=None, full_output=0, disp=1, retall=0,
             callback=None):
    """Minimize a function using the Newton-CG method.

    Parameters
    ----------
    f : callable f(x,*args)
        Objective function to be minimized.
    x0 : ndarray
        Initial guess.
    fprime : callable f'(x,*args)
        Gradient of f.
    fhess_p : callable fhess_p(x,p,*args)
        Function which computes the Hessian of f times an
        arbitrary vector, p.
    fhess : callable fhess(x,*args)
        Function to compute the Hessian matrix of f.
    args : tuple
        Extra arguments passed to f, fprime, fhess_p, and fhess
        (the same set of extra arguments is supplied to all of
        these functions).
    epsilon : float or ndarray
        If fhess is approximated, use this value for the step size.
    callback : callable
        An optional user-supplied function which is called after
        each iteration.  Called as callback(xk), where xk is the
        current parameter vector.

    Returns
    -------
    xopt : ndarray
        Parameters which minimizer f, i.e. ``f(xopt) == fopt``.
    fopt : float
        Value of the function at xopt, i.e. ``fopt = f(xopt)``.
    fcalls : int
        Number of function calls made.
    gcalls : int
        Number of gradient calls made.
    hcalls : int
        Number of hessian calls made.
    warnflag : int
        Warnings generated by the algorithm.
        1 : Maximum number of iterations exceeded.
    allvecs : list
        The result at each iteration, if retall is True (see below).

    Other Parameters
    ----------------
    avextol : float
        Convergence is assumed when the average relative error in
        the minimizer falls below this amount.
    maxiter : int
        Maximum number of iterations to perform.
    full_output : bool
        If True, return the optional outputs.
    disp : bool
        If True, print convergence message.
    retall : bool
        If True, return a list of results at each iteration.

    Notes
    -----
    Only one of `fhess_p` or `fhess` need to be given.  If `fhess`
    is provided, then `fhess_p` will be ignored.  If neither `fhess`
    nor `fhess_p` is provided, then the hessian product will be
    approximated using finite differences on `fprime`. `fhess_p`
    must compute the hessian times an arbitrary vector. If it is not
    given, finite-differences on `fprime` are used to compute
    it. See Wright & Nocedal, 'Numerical Optimization', 1999,
    pg. 140.

    """
    x0 = asarray(x0).flatten()
    fcalls, f = wrap_function(f, args)
    gcalls, fprime = wrap_function(fprime, args)
    hcalls = 0
    if maxiter is None:
        maxiter = len(x0)*200

    xtol = len(x0)*avextol
    update = [2*xtol]
    xk = x0
    if retall:
        allvecs = [xk]
    k = 0
    old_fval = f(x0)
    while (numpy.add.reduce(abs(update)) > xtol) and (k < maxiter):
        # Compute a search direction pk by applying the CG method to
        #  del2 f(xk) p = - grad f(xk) starting from 0.
        b = -fprime(xk)
        maggrad = numpy.add.reduce(abs(b))
        eta = min([0.5,numpy.sqrt(maggrad)])
        termcond = eta * maggrad
        xsupi = zeros(len(x0), dtype=x0.dtype)
        ri = -b
        psupi = -ri
        i = 0
        dri0 = numpy.dot(ri,ri)

        if fhess is not None:             # you want to compute hessian once.
            A = fhess(*(xk,)+args)
            hcalls = hcalls + 1

        while numpy.add.reduce(abs(ri)) > termcond:
            if fhess is None:
                if fhess_p is None:
                    Ap = approx_fhess_p(xk,psupi,fprime,epsilon)
                else:
                    Ap = fhess_p(xk,psupi, *args)
                    hcalls = hcalls + 1
            else:
                Ap = numpy.dot(A,psupi)
            # check curvature
            Ap = asarray(Ap).squeeze() # get rid of matrices...
            curv = numpy.dot(psupi,Ap)
            if 0 <= curv <= 3*numpy.finfo(numpy.float64).eps:
                break
            elif curv < 0:
                if (i > 0):
                    break
                else:
                    xsupi = xsupi + dri0/curv * psupi
                    break
            alphai = dri0 / curv
            xsupi = xsupi + alphai * psupi
            ri = ri + alphai * Ap
            dri1 = numpy.dot(ri,ri)
            betai = dri1 / dri0
            psupi = -ri + betai * psupi
            i = i + 1
            dri0 = dri1          # update numpy.dot(ri,ri) for next time.

        pk = xsupi  # search direction is solution to system.
        gfk = -b    # gradient at xk
        alphak, fc, gc, old_fval = line_search_BFGS(f,xk,pk,gfk,old_fval)

        update = alphak * pk
        xk = xk + update        # upcast if necessary
        if callback is not None:
            callback(xk)
        if retall:
            allvecs.append(xk)
        k += 1

    if disp or full_output:
        fval = old_fval
    if k >= maxiter:
        warnflag = 1
        if disp:
            print "Warning: Maximum number of iterations has been exceeded"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % fcalls[0]
            print "         Gradient evaluations: %d" % gcalls[0]
            print "         Hessian evaluations: %d" % hcalls
    else:
        warnflag = 0
        if disp:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % fcalls[0]
            print "         Gradient evaluations: %d" % gcalls[0]
            print "         Hessian evaluations: %d" % hcalls

    if full_output:
        retlist = xk, fval, fcalls[0], gcalls[0], hcalls, warnflag
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

    Parameters
    ----------
    func : callable f(x,*args)
        Objective function to be minimized (must accept and return scalars).
    x1, x2 : float or array scalar
        The optimization bounds.
    args : tuple
        Extra arguments passed to function.
    xtol : float
        The convergence tolerance.
    maxfun : int
        Maximum number of function evaluations allowed.
    full_output : bool
        If True, return optional outputs.
    disp : int
        If non-zero, print messages.
            0 : no message printing.
            1 : non-convergence notification messages only.
            2 : print a message on convergence too.
            3 : print iteration results.


    Returns
    -------
    xopt : ndarray
        Parameters (over given interval) which minimize the
        objective function.
    fval : number
        The function value at the minimum point.
    ierr : int
        An error flag (0 if converged, 1 if maximum number of
        function calls reached).
    numfunc : int
      The number of function calls made.

    Notes
    -----
    Finds a local minimizer of the scalar function `func` in the
    interval x1 < xopt < x2 using Brent's method.  (See `brent`
    for auto-bracketing).

    """
    # Test bounds are of correct form

    if not (is_array_scalar(x1) and is_array_scalar(x2)):
        raise ValueError("Optimisation bounds must be scalars"
                         " or array scalars.")
    if x1 > x2:
        raise ValueError("The lower bound exceeds the upper bound.")

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
            if ( (abs(p) < abs(0.5*q*r)) and (p > q*(a-xf)) and \
                 (p < q*(b-xf))):
                rat = (p+0.0) / q;
                x = xf + rat
                step = '       parabolic'

                if ((x-a) < tol2) or ((b-x) < tol2):
                    si = numpy.sign(xm-xf) + ((xm-xf)==0)
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

        si = numpy.sign(rat) + (rat == 0)
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
                _endprint(x, flag, fval, maxfun, xtol, disp)
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

class Brent:
    #need to rethink design of __init__
    def __init__(self, func, args=(), tol=1.48e-8, maxiter=500,
                 full_output=0):
        self.func = func
        self.args = args
        self.tol = tol
        self.maxiter = maxiter
        self._mintol = 1.0e-11
        self._cg = 0.3819660
        self.xmin = None
        self.fval = None
        self.iter = 0
        self.funcalls = 0

    #need to rethink design of set_bracket (new options, etc)
    def set_bracket(self, brack = None):
        self.brack = brack
    def get_bracket_info(self):
        #set up
        func = self.func
        args = self.args
        brack = self.brack
        ### BEGIN core bracket_info code ###
        ### carefully DOCUMENT any CHANGES in core ##
        if brack is None:
            xa,xb,xc,fa,fb,fc,funcalls = bracket(func, args=args)
        elif len(brack) == 2:
            xa,xb,xc,fa,fb,fc,funcalls = bracket(func, xa=brack[0],
                                                 xb=brack[1], args=args)
        elif len(brack) == 3:
            xa,xb,xc = brack
            if (xa > xc):  # swap so xa < xc can be assumed
                dum = xa; xa=xc; xc=dum
            assert ((xa < xb) and (xb < xc)), "Not a bracketing interval."
            fa = func(*((xa,)+args))
            fb = func(*((xb,)+args))
            fc = func(*((xc,)+args))
            assert ((fb<fa) and (fb < fc)), "Not a bracketing interval."
            funcalls = 3
        else:
            raise ValueError("Bracketing interval must be " \
                             "length 2 or 3 sequence.")
        ### END core bracket_info code ###

        return xa,xb,xc,fa,fb,fc,funcalls

    def optimize(self):
        #set up for optimization
        func = self.func
        xa,xb,xc,fa,fb,fc,funcalls = self.get_bracket_info()
        _mintol = self._mintol
        _cg = self._cg
        #################################
        #BEGIN CORE ALGORITHM
        #we are making NO CHANGES in this
        #################################
        x=w=v=xb
        fw=fv=fx=func(*((x,)+self.args))
        if (xa < xc):
            a = xa; b = xc
        else:
            a = xc; b = xa
        deltax= 0.0
        funcalls = 1
        iter = 0
        while (iter < self.maxiter):
            tol1 = self.tol*abs(x) + _mintol
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
            fu = func(*((u,)+self.args))      # calculate new output value
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

            iter += 1
        #################################
        #END CORE ALGORITHM
        #################################

        self.xmin = x
        self.fval = fx
        self.iter = iter
        self.funcalls = funcalls

    def get_result(self, full_output=False):
        if full_output:
            return self.xmin, self.fval, self.iter, self.funcalls
        else:
            return self.xmin


def brent(func, args=(), brack=None, tol=1.48e-8, full_output=0, maxiter=500):
    """Given a function of one-variable and a possible bracketing interval,
    return the minimum of the function isolated to a fractional precision of
    tol.

    Parameters
    ----------
    func : callable f(x,*args)
        Objective function.
    args
        Additional arguments (if present).
    brack : tuple
        Triple (a,b,c) where (a<b<c) and func(b) <
        func(a),func(c).  If bracket consists of two numbers (a,c)
        then they are assumed to be a starting interval for a
        downhill bracket search (see `bracket`); it doesn't always
        mean that the obtained solution will satisfy a<=x<=c.
    full_output : bool
        If True, return all output args (xmin, fval, iter,
        funcalls).

    Returns
    -------
    xmin : ndarray
        Optimum point.
    fval : float
        Optimum value.
    iter : int
        Number of iterations.
    funcalls : int
        Number of objective function evaluations made.

    Notes
    -----
    Uses inverse parabolic interpolation when possible to speed up
    convergence of golden section method.

    """

    brent = Brent(func=func, args=args, tol=tol,
                  full_output=full_output, maxiter=maxiter)
    brent.set_bracket(brack)
    brent.optimize()
    return brent.get_result(full_output=full_output)


def golden(func, args=(), brack=None, tol=_epsilon, full_output=0):
    """ Given a function of one-variable and a possible bracketing interval,
    return the minimum of the function isolated to a fractional precision of
    tol.

    Parameters
    ----------
    func : callable func(x,*args)
        Objective function to minimize.
    args : tuple
        Additional arguments (if present), passed to func.
    brack : tuple
        Triple (a,b,c), where (a<b<c) and func(b) <
        func(a),func(c).  If bracket consists of two numbers (a,
        c), then they are assumed to be a starting interval for a
        downhill bracket search (see `bracket`); it doesn't always
        mean that obtained solution will satisfy a<=x<=c.
    tol : float
        x tolerance stop criterion
    full_output : bool
        If True, return optional outputs.

    Notes
    -----
    Uses analog of bisection method to decrease the bracketed
    interval.

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
        fa = func(*((xa,)+args))
        fb = func(*((xb,)+args))
        fc = func(*((xc,)+args))
        assert ((fb<fa) and (fb < fc)), "Not a bracketing interval."
        funcalls = 3
    else:
        raise ValueError, "Bracketing interval must be length 2 or 3 sequence."

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
    f1 = func(*((x1,)+args))
    f2 = func(*((x2,)+args))
    funcalls += 2
    while (abs(x3-x0) > tol*(abs(x1)+abs(x2))):
        if (f2 < f1):
            x0 = x1; x1 = x2; x2 = _gR*x1 + _gC*x3
            f1 = f2; f2 = func(*((x2,)+args))
        else:
            x3 = x2; x2 = x1; x1 = _gR*x2 + _gC*x0
            f2 = f1; f1 = func(*((x1,)+args))
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


def bracket(func, xa=0.0, xb=1.0, args=(), grow_limit=110.0, maxiter=1000):
    """Given a function and distinct initial points, search in the
    downhill direction (as defined by the initital points) and return
    new points xa, xb, xc that bracket the minimum of the function
    f(xa) > f(xb) < f(xc). It doesn't always mean that obtained
    solution will satisfy xa<=x<=xb

    Parameters
    ----------
    func : callable f(x,*args)
        Objective function to minimize.
    xa, xb : float
        Bracketing interval.
    args : tuple
        Additional arguments (if present), passed to `func`.
    grow_limit : float
        Maximum grow limit.
    maxiter : int
        Maximum number of iterations to perform.

    Returns
    -------
    xa, xb, xc : float
        Bracket.
    fa, fb, fc : float
        Objective function values in bracket.
    funcalls : int
        Number of function evaluations made.

    """
    _gold = 1.618034
    _verysmall_num = 1e-21
    fa = func(*(xa,)+args)
    fb = func(*(xb,)+args)
    if (fa < fb):                      # Switch so fa > fb
        dum = xa; xa = xb; xb = dum
        dum = fa; fa = fb; fb = dum
    xc = xb + _gold*(xb-xa)
    fc = func(*((xc,)+args))
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
        if iter > maxiter:
            raise RuntimeError, "Too many iterations."
        iter += 1
        if (w-xc)*(xb-w) > 0.0:
            fw = func(*((w,)+args))
            funcalls += 1
            if (fw < fc):
                xa = xb; xb=w; fa=fb; fb=fw
                return xa, xb, xc, fa, fb, fc, funcalls
            elif (fw > fb):
                xc = w; fc=fw
                return xa, xb, xc, fa, fb, fc, funcalls
            w = xc + _gold*(xc-xb)
            fw = func(*((w,)+args))
            funcalls += 1
        elif (w-wlim)*(wlim-xc) >= 0.0:
            w = wlim
            fw = func(*((w,)+args))
            funcalls += 1
        elif (w-wlim)*(xc-w) > 0.0:
            fw = func(*((w,)+args))
            funcalls += 1
            if (fw < fc):
                xb=xc; xc=w; w=xc+_gold*(xc-xb)
                fb=fc; fc=fw; fw=func(*((w,)+args))
                funcalls += 1
        else:
            w = xc + _gold*(xc-xb)
            fw = func(*((w,)+args))
            funcalls += 1
        xa=xb; xb=xc; xc=w
        fa=fb; fb=fc; fc=fw
    return xa, xb, xc, fa, fb, fc, funcalls



def _linesearch_powell(func, p, xi, tol=1e-3):
    """Line-search algorithm using fminbound.

    Find the minimium of the function ``func(x0+ alpha*direc)``.

    """
    def myfunc(alpha):
        return func(p + alpha * xi)
    alpha_min, fret, iter, num = brent(myfunc, full_output=1, tol=tol)
    xi = alpha_min*xi
    return squeeze(fret), p+xi, xi


def fmin_powell(func, x0, args=(), xtol=1e-4, ftol=1e-4, maxiter=None,
                maxfun=None, full_output=0, disp=1, retall=0, callback=None,
                direc=None):
    """Minimize a function using modified Powell's method.

    Parameters
    ----------
    func : callable f(x,*args)
        Objective function to be minimized.
    x0 : ndarray
        Initial guess.
    args : tuple
        Eextra arguments passed to func.
    callback : callable
        An optional user-supplied function, called after each
        iteration.  Called as ``callback(xk)``, where ``xk`` is the
        current parameter vector.
    direc : ndarray
        Initial direction set.

    Returns
    -------
    xopt : ndarray
        Parameter which minimizes `func`.
    fopt : number
        Value of function at minimum: ``fopt = func(xopt)``.
    direc : ndarray
        Current direction set.
    iter : int
        Number of iterations.
    funcalls : int
        Number of function calls made.
    warnflag : int
        Integer warning flag:
            1 : Maximum number of function evaluations.
            2 : Maximum number of iterations.
    allvecs : list
        List of solutions at each iteration.

    Other Parameters
    ----------------
    xtol : float
        Line-search error tolerance.
    ftol : float
        Relative error in ``func(xopt)`` acceptable for convergence.
    maxiter : int
        Maximum number of iterations to perform.
    maxfun : int
        Maximum number of function evaluations to make.
    full_output : bool
        If True, fopt, xi, direc, iter, funcalls, and
        warnflag are returned.
    disp : bool
        If True, print convergence messages.
    retall : bool
        If True, return a list of the solution at each iteration.

    Notes
    -----
    Uses a modification of Powell's method to find the minimum of
    a function of N variables.

    """
    # we need to use a mutable object here that we can update in the
    # wrapper function
    fcalls, func = wrap_function(func, args)
    x = asarray(x0).flatten()
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


    if direc is None:
        direc = eye(N, dtype=float)
    else:
        direc = asarray(direc, dtype=float)

    fval = squeeze(func(x))
    x1 = x.copy()
    iter = 0;
    ilist = range(N)
    while True:
        fx = fval
        bigind = 0
        delta = 0.0
        for i in ilist:
            direc1 = direc[i]
            fx2 = fval
            fval, x, direc1 = _linesearch_powell(func, x, direc1, tol=xtol*100)
            if (fx2 - fval) > delta:
                delta = fx2 - fval
                bigind = i
        iter += 1
        if callback is not None:
            callback(x)
        if retall:
            allvecs.append(x)
        if (2.0*(fx - fval) <= ftol*(abs(fx)+abs(fval))+1e-20): break
        if fcalls[0] >= maxfun: break
        if iter >= maxiter: break

        # Construct the extrapolated point
        direc1 = x - x1
        x2 = 2*x - x1
        x1 = x.copy()
        fx2 = squeeze(func(x2))

        if (fx > fx2):
            t = 2.0*(fx+fx2-2.0*fval)
            temp = (fx-fval-delta)
            t *= temp*temp
            temp = fx-fx2
            t -= delta*temp*temp
            if t < 0.0:
                fval, x, direc1 = _linesearch_powell(func, x, direc1,
                                                     tol=xtol*100)
                direc[bigind] = direc[-1]
                direc[-1] = direc1

    warnflag = 0
    if fcalls[0] >= maxfun:
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
            print "         Function evaluations: %d" % fcalls[0]

    x = squeeze(x)

    if full_output:
        retlist = x, fval, direc, iter, fcalls[0], warnflag
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
            print "\nOptimization terminated successfully;\n" \
                  "The returned value satisfies the termination criteria\n" \
                  "(using xtol = ", xtol, ")"
    if flag == 1:
        print "\nMaximum number of function evaluations exceeded --- " \
              "increase maxfun argument.\n"
    return


def brute(func, ranges, args=(), Ns=20, full_output=0, finish=fmin):
    """Minimize a function over a given range by brute force.

    Parameters
    ----------
    func : callable ``f(x,*args)``
        Objective function to be minimized.
    ranges : tuple
        Each element is a tuple of parameters or a slice object to
        be handed to ``numpy.mgrid``.
    args : tuple
        Extra arguments passed to function.
    Ns : int
        Default number of samples, if those are not provided.
    full_output : bool
        If True, return the evaluation grid.

    Returns
    -------
    x0 : ndarray
        Value of arguments to `func`, giving minimum over the grid.
    fval : int
        Function value at minimum.
    grid : tuple
        Representation of the evaluation grid.  It has the same
        length as x0.
    Jout : ndarray
        Function values over grid:  ``Jout = func(*grid)``.

    Notes
    -----
    Find the minimum of a function evaluated on a grid given by
    the tuple ranges.

    """
    N = len(ranges)
    if N > 40:
        raise ValueError, "Brute Force not possible with more " \
              "than 40 variables."
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
    indx = argmin(Jout.ravel(),axis=-1)
    Nindx = zeros(N,int)
    xmin = zeros(N,float)
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


def main():
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

if __name__ == "__main__":
    main()
