"""
Functions
---------
.. autosummary::
   :toctree: generated/

    line_search_armijo
    line_search_wolfe1
    line_search_wolfe2
    scalar_search_wolfe1
    scalar_search_wolfe2

"""
from scipy.optimize import minpack2
import numpy as np
from numpy.compat import asbytes

__all__ = ['line_search_wolfe1', 'line_search_wolfe2',
           'scalar_search_wolfe1', 'scalar_search_wolfe2',
           'line_search_armijo']

#------------------------------------------------------------------------------
# Minpack's Wolfe line and scalar searches
#------------------------------------------------------------------------------

def line_search_wolfe1(f, fprime, xk, pk, gfk=None,
                       old_fval=None, old_old_fval=None,
                       args=(), c1=1e-4, c2=0.9, amax=50, amin=1e-8,
                       xtol=1e-14):
    """
    As `scalar_search_wolfe1` but do a line search to direction `pk`

    Parameters
    ----------
    f : callable
        Function `f(x)`
    fprime : callable
        Gradient of `f`
    xk : array_like
        Current point
    pk : array_like
        Search direction

    gfk : array_like, optional
        Gradient of `f` at point `xk`
    old_fval : float, optional
        Value of `f` at point `xk`
    old_old_fval : float, optional
        Value of `f` at point preceding `xk`

    The rest of the parameters are the same as for `scalar_search_wolfe1`.

    Returns
    -------
    stp, f_count, g_count, fval, old_fval
        As in `line_search_wolfe1`
    gval : array
        Gradient of `f` at the final point

    """
    if gfk is None:
        gfk = fprime(xk)

    if isinstance(fprime, tuple):
        eps = fprime[1]
        fprime = fprime[0]
        newargs = (f, eps) + args
        gradient = False
    else:
        newargs = args
        gradient = True

    gval = [gfk]
    gc = [0]
    fc = [0]

    def phi(s):
        fc[0] += 1
        return f(xk + s*pk, *args)

    def derphi(s):
        gval[0] = fprime(xk + s*pk, *newargs)
        if gradient:
            gc[0] += 1
        else:
            fc[0] += len(xk) + 1
        return np.dot(gval[0], pk)

    derphi0 = np.dot(gfk, pk)

    stp, fval, old_fval = scalar_search_wolfe1(
            phi, derphi, old_fval, old_old_fval, derphi0,
            c1=c1, c2=c2, amax=amax, amin=amin, xtol=xtol)

    return stp, fc[0], gc[0], fval, old_fval, gval[0]


def scalar_search_wolfe1(phi, derphi, phi0=None, old_phi0=None, derphi0=None,
                         c1=1e-4, c2=0.9,
                         amax=50, amin=1e-8, xtol=1e-14):
    """
    Scalar function search for alpha that satisfies strong Wolfe conditions

    alpha > 0 is assumed to be a descent direction.

    Parameters
    ----------
    phi : callable phi(alpha)
        Function at point `alpha`
    derphi : callable dphi(alpha)
        Derivative `d phi(alpha)/ds`. Returns a scalar.

    phi0 : float, optional
        Value of `f` at 0
    old_phi0 : float, optional
        Value of `f` at the previous point
    derphi0 : float, optional
        Value `derphi` at 0
    amax : float, optional
        Maximum step size
    c1, c2 : float, optional
        Wolfe parameters

    Returns
    -------
    alpha : float
        Step size, or None if no suitable step was found
    phi : float
        Value of `phi` at the new point `alpha`
    phi0 : float
        Value of `phi` at `alpha=0`

    Notes
    -----
    Uses routine DCSRCH from MINPACK.

    """

    if phi0 is None:
        phi0 = phi(0.)
    if derphi0 is None:
        derphi0 = derphi(0.)

    if old_phi0 is not None:
        alpha1 = min(1.0, 1.01*2*(phi0 - old_phi0)/derphi0)
        if alpha1 < 0:
            alpha1 = 1.0
    else:
        alpha1 = 1.0

    phi1 = phi0
    derphi1 = derphi0
    isave = np.zeros((2,), np.intc)
    dsave = np.zeros((13,), float)
    task = asbytes('START')

    maxiter=30
    for i in xrange(maxiter):
        stp, phi1, derphi1, task = minpack2.dcsrch(alpha1, phi1, derphi1,
                                                   c1, c2, xtol, task,
                                                   amin, amax, isave, dsave)
        if task[:2] == asbytes('FG'):
            alpha1 = stp
            phi1 = phi(stp)
            derphi1 = derphi(stp)
        else:
            break
    else:
        # maxiter reached, the line search did not converge
        stp=None

    if task[:5] == asbytes('ERROR') or task[:4] == asbytes('WARN'):
        stp = None  # failed

    return stp, phi1, phi0

line_search = line_search_wolfe1

#------------------------------------------------------------------------------
# Pure-Python Wolfe line and scalar searches
#------------------------------------------------------------------------------

def line_search_wolfe2(f, myfprime, xk, pk, gfk=None, old_fval=None,
                       old_old_fval=None, args=(), c1=1e-4, c2=0.9, amax=50):
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
    gfk : ndarray, optional
        Gradient value for x=xk (xk being the current parameter
        estimate). Will be recomputed if omitted.
    old_fval : float, optional
        Function value for x=xk. Will be recomputed if omitted.
    old_old_fval : float, optional
        Function value for the point preceding x=xk
    args : tuple, optional
        Additional arguments passed to objective function.
    c1 : float, optional
        Parameter for Armijo condition rule.
    c2 : float, optional
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
    fc = [0]
    gc = [0]
    gval = [None]

    def phi(alpha):
        fc[0] += 1
        return f(xk + alpha * pk, *args)

    if isinstance(myfprime, tuple):
        def derphi(alpha):
            fc[0] += len(xk)+1
            eps = myfprime[1]
            fprime = myfprime[0]
            newargs = (f,eps) + args
            gval[0] = fprime(xk+alpha*pk, *newargs)  # store for later use
            return np.dot(gval[0], pk)
    else:
        fprime = myfprime
        def derphi(alpha):
            gc[0] += 1
            gval[0] = fprime(xk+alpha*pk, *args)  # store for later use
            return np.dot(gval[0], pk)

    derphi0 = np.dot(gfk, pk)

    alpha_star, phi_star, old_fval, derphi_star = \
                scalar_search_wolfe2(phi, derphi, old_fval, old_old_fval,
                                     derphi0, c1, c2, amax)

    if derphi_star is not None:
        # derphi_star is a number (derphi) -- so use the most recently
        # calculated gradient used in computing it derphi = gfk*pk
        # this is the gradient at the next step no need to compute it
        # again in the outer loop.
        derphi_star = gval[0]

    return alpha_star, fc[0], gc[0], phi_star, old_fval, derphi_star


def scalar_search_wolfe2(phi, derphi=None, phi0=None,
                         old_phi0=None, derphi0=None,
                         c1=1e-4, c2=0.9, amax=50):
    """Find alpha that satisfies strong Wolfe conditions.

    alpha > 0 is assumed to be a descent direction.

    Parameters
    ----------
    phi : callable f(x,*args)
        Objective scalar function.

    derphi : callable f'(x,*args), optional
        Objective function derivative (can be None)
    phi0 : float, optional
        Value of phi at s=0
    old_phi0 : float, optional
        Value of phi at previous point
    derphi0 : float, optional
        Value of derphi at s=0
    args : tuple
        Additional arguments passed to objective function.
    c1 : float
        Parameter for Armijo condition rule.
    c2 : float
        Parameter for curvature condition rule.

    Returns
    -------
    alpha_star : float
        Best alpha
    phi_star
        phi at alpha_star
    phi0
        phi at 0
    derphi_star
        derphi at alpha_star

    Notes
    -----
    Uses the line search algorithm to enforce strong Wolfe
    conditions.  See Wright and Nocedal, 'Numerical Optimization',
    1999, pg. 59-60.

    For the zoom phase it uses an algorithm by [...].

    """

    if phi0 is None:
        phi0 = phi(0.)

    if derphi0 is None and derphi is not None:
        derphi0 = derphi(0.)

    alpha0 = 0
    if old_phi0 is not None:
        alpha1 = min(1.0, 1.01*2*(phi0 - old_phi0)/derphi0)
    else:
        alpha1 = 1.0

    if alpha1 < 0:
        alpha1 = 1.0

    if alpha1 == 0:
        # This shouldn't happen. Perhaps the increment has slipped below
        # machine precision?  For now, set the return variables skip the
        # useless while loop, and raise warnflag=2 due to possible imprecision.
        alpha_star = None
        phi_star = phi0
        phi0 = old_phi0
        derphi_star = None

    phi_a1 = phi(alpha1)
    #derphi_a1 = derphi(alpha1)  evaluated below

    phi_a0 = phi0
    derphi_a0 = derphi0

    i = 1
    maxiter = 10
    for i in xrange(maxiter):
        if alpha1 == 0:
            break
        if (phi_a1 > phi0 + c1*alpha1*derphi0) or \
           ((phi_a1 >= phi_a0) and (i > 1)):
            alpha_star, phi_star, derphi_star = \
                        _zoom(alpha0, alpha1, phi_a0,
                              phi_a1, derphi_a0, phi, derphi,
                              phi0, derphi0, c1, c2)
            break

        derphi_a1 = derphi(alpha1)
        if (abs(derphi_a1) <= -c2*derphi0):
            alpha_star = alpha1
            phi_star = phi_a1
            derphi_star = derphi_a1
            break

        if (derphi_a1 >= 0):
            alpha_star, phi_star, derphi_star = \
                        _zoom(alpha1, alpha0, phi_a1,
                              phi_a0, derphi_a1, phi, derphi,
                              phi0, derphi0, c1, c2)
            break

        alpha2 = 2 * alpha1   # increase by factor of two on each iteration
        i = i + 1
        alpha0 = alpha1
        alpha1 = alpha2
        phi_a0 = phi_a1
        phi_a1 = phi(alpha1)
        derphi_a0 = derphi_a1

    else:
        # stopping test maxiter reached
        alpha_star = alpha1
        phi_star = phi_a1
        derphi_star = None

    return alpha_star, phi_star, phi0, derphi_star


def _cubicmin(a,fa,fpa,b,fb,c,fc):
    """
    Finds the minimizer for a cubic polynomial that goes through the
    points (a,fa), (b,fb), and (c,fc) with derivative at a of fpa.

    If no minimizer can be found return None

    """
    # f(x) = A *(x-a)^3 + B*(x-a)^2 + C*(x-a) + D

    C = fpa
    D = fa
    db = b-a
    dc = c-a
    if (db == 0) or (dc == 0) or (b==c): return None
    denom = (db*dc)**2 * (db-dc)
    d1 = np.empty((2,2))
    d1[0,0] = dc**2
    d1[0,1] = -db**2
    d1[1,0] = -dc**3
    d1[1,1] = db**3
    [A,B] = np.dot(d1, np.asarray([fb-fa-C*db,fc-fa-C*dc]).flatten())
    A /= denom
    B /= denom
    radical = B*B-3*A*C
    if radical < 0:  return None
    if (A == 0): return None
    xmin = a + (-B + np.sqrt(radical))/(3*A)
    return xmin


def _quadmin(a,fa,fpa,b,fb):
    """
    Finds the minimizer for a quadratic polynomial that goes through
    the points (a,fa), (b,fb) with derivative at a of fpa,

    """
    # f(x) = B*(x-a)^2 + C*(x-a) + D
    D = fa
    C = fpa
    db = b-a*1.0
    if (db==0): return None
    B = (fb-D-C*db)/(db*db)
    if (B <= 0): return None
    xmin = a  - C / (2.0*B)
    return xmin

def _zoom(a_lo, a_hi, phi_lo, phi_hi, derphi_lo,
          phi, derphi, phi0, derphi0, c1, c2):
    """
    Part of the optimization algorithm in `scalar_search_wolfe2`.
    """

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
        # (uses phi_lo, derphi_lo, phi_hi, and the most recent value of phi)
        #
        # if the result is too close to the end points (or out of the
        # interval) then use quadratic interpolation with phi_lo,
        # derphi_lo and phi_hi if the result is stil too close to the
        # end points (or out of the interval) then use bisection

        if (i > 0):
            cchk = delta1*dalpha
            a_j = _cubicmin(a_lo, phi_lo, derphi_lo, a_hi, phi_hi, a_rec, phi_rec)
        if (i==0) or (a_j is None) or (a_j > b-cchk) or (a_j < a+cchk):
            qchk = delta2*dalpha
            a_j = _quadmin(a_lo, phi_lo, derphi_lo, a_hi, phi_hi)
            if (a_j is None) or (a_j > b-qchk) or (a_j < a+qchk):
                a_j = a_lo + 0.5*dalpha

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


#------------------------------------------------------------------------------
# Armijo line and scalar searches
#------------------------------------------------------------------------------

def line_search_armijo(f, xk, pk, gfk, old_fval, args=(), c1=1e-4, alpha0=1):
    """Minimize over alpha, the function ``f(xk+alpha pk)``.

    Parameters
    ----------
    f : callable
        Function to be minimized.
    xk : array_like
        Current point.
    pk : array_like
        Search direction.
    gfk : array_like, optional
        Gradient of `f` at point `xk`.
    old_fval : float
        Value of `f` at point `xk`.
    args : tuple, optional
        Optional arguments.
    c1 : float, optional
        Value to control stopping criterion.
    alpha0 : scalar, optional
        Value of `alpha` at start of the optimization.

    Returns
    -------
    alpha
    f_count
    f_val_at_alpha

    Notes
    -----
    Uses the interpolation algorithm (Armijo backtracking) as suggested by
    Wright and Nocedal in 'Numerical Optimization', 1999, pg. 56-57

    """
    xk = np.atleast_1d(xk)
    fc = [0]

    def phi(alpha1):
        fc[0] += 1
        return f(xk + alpha1*pk, *args)

    if old_fval is None:
        phi0 = phi(0.)
    else:
        phi0 = old_fval # compute f(xk) -- done in past loop

    derphi0 = np.dot(gfk, pk)
    alpha, phi1 = scalar_search_armijo(phi, phi0, derphi0, c1=c1, alpha0=alpha0)
    return alpha, fc[0], phi1

def line_search_BFGS(f, xk, pk, gfk, old_fval, args=(), c1=1e-4, alpha0=1):
    """
    Compatibility wrapper for `line_search_armijo`
    """
    r = line_search_armijo(f, xk, pk, gfk, old_fval, args=args, c1=c1,
                           alpha0=alpha0)
    return r[0], r[1], 0, r[2]

def scalar_search_armijo(phi, phi0, derphi0, c1=1e-4, alpha0=1, amin=0):
    """Minimize over alpha, the function ``phi(alpha)``.

    Uses the interpolation algorithm (Armijo backtracking) as suggested by
    Wright and Nocedal in 'Numerical Optimization', 1999, pg. 56-57

    alpha > 0 is assumed to be a descent direction.

    Returns
    -------
    alpha
    phi1

    """
    phi_a0 = phi(alpha0)
    if phi_a0 <= phi0 + c1*alpha0*derphi0:
        return alpha0, phi_a0

    # Otherwise compute the minimizer of a quadratic interpolant:

    alpha1 = -(derphi0) * alpha0**2 / 2.0 / (phi_a0 - phi0 - derphi0 * alpha0)
    phi_a1 = phi(alpha1)

    if (phi_a1 <= phi0 + c1*alpha1*derphi0):
        return alpha1, phi_a1

    # Otherwise loop with cubic interpolation until we find an alpha which
    # satifies the first Wolfe condition (since we are backtracking, we will
    # assume that the value of alpha is not too small and satisfies the second
    # condition.

    while alpha1 > amin:       # we are assuming alpha>0 is a descent direction
        factor = alpha0**2 * alpha1**2 * (alpha1-alpha0)
        a = alpha0**2 * (phi_a1 - phi0 - derphi0*alpha1) - \
            alpha1**2 * (phi_a0 - phi0 - derphi0*alpha0)
        a = a / factor
        b = -alpha0**3 * (phi_a1 - phi0 - derphi0*alpha1) + \
            alpha1**3 * (phi_a0 - phi0 - derphi0*alpha0)
        b = b / factor

        alpha2 = (-b + np.sqrt(abs(b**2 - 3 * a * derphi0))) / (3.0*a)
        phi_a2 = phi(alpha2)

        if (phi_a2 <= phi0 + c1*alpha2*derphi0):
            return alpha2, phi_a2

        if (alpha1 - alpha2) > alpha1 / 2.0 or (1 - alpha2/alpha1) < 0.96:
            alpha2 = alpha1 / 2.0

        alpha0 = alpha1
        alpha1 = alpha2
        phi_a0 = phi_a1
        phi_a1 = phi_a2

    # Failed to find a suitable step length
    return None, phi_a1
