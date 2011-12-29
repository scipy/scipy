import warnings
import _minpack

from numpy import atleast_1d, dot, take, triu, shape, eye, \
                  transpose, zeros, product, greater, array, \
                  all, where, isscalar, asarray, inf, abs

error = _minpack.error

__all__ = ['fsolve', 'leastsq', 'fixed_point', 'curve_fit']

def _check_func(checker, argname, thefunc, x0, args, numinputs, output_shape=None):
    res = atleast_1d(thefunc(*((x0[:numinputs],) + args)))
    if (output_shape is not None) and (shape(res) != output_shape):
        if (output_shape[0] != 1):
            if len(output_shape) > 1:
                if output_shape[1] == 1:
                    return shape(res)
            msg = "%s: there is a mismatch between the input and output " \
                  "shape of the '%s' argument" % (checker, argname)
            func_name = getattr(thefunc, 'func_name', None)
            if func_name:
                msg += " '%s'." % func_name
            else:
                msg += "."
            raise TypeError(msg)
    return shape(res)


def fsolve(func, x0, args=(), fprime=None, full_output=0,
           col_deriv=0, xtol=1.49012e-8, maxfev=0, band=None,
           epsfcn=0.0, factor=100, diag=None):
    """
    Find the roots of a function.

    Return the roots of the (non-linear) equations defined by
    ``func(x) = 0`` given a starting estimate.

    Parameters
    ----------
    func : callable f(x, *args)
        A function that takes at least one (possibly vector) argument.
    x0 : ndarray
        The starting estimate for the roots of ``func(x) = 0``.
    args : tuple
        Any extra arguments to `func`.
    fprime : callable(x)
        A function to compute the Jacobian of `func` with derivatives
        across the rows. By default, the Jacobian will be estimated.
    full_output : bool
        If True, return optional outputs.
    col_deriv : bool
        Specify whether the Jacobian function computes derivatives down
        the columns (faster, because there is no transpose operation).

    Returns
    -------
    x : ndarray
        The solution (or the result of the last iteration for
        an unsuccessful call).
    infodict : dict
        A dictionary of optional outputs with the keys::

          * 'nfev': number of function calls
          * 'njev': number of Jacobian calls
          * 'fvec': function evaluated at the output
          * 'fjac': the orthogonal matrix, q, produced by the QR
                    factorization of the final approximate Jacobian
                    matrix, stored column wise
          * 'r': upper triangular matrix produced by QR factorization of same
                 matrix
          * 'qtf': the vector (transpose(q) * fvec)

    ier : int
        An integer flag.  Set to 1 if a solution was found, otherwise refer
        to `mesg` for more information.
    mesg : str
        If no solution is found, `mesg` details the cause of failure.

    Other Parameters
    ----------------
    xtol : float
        The calculation will terminate if the relative error between two
        consecutive iterates is at most `xtol`.
    maxfev : int
        The maximum number of calls to the function. If zero, then
        ``100*(N+1)`` is the maximum where N is the number of elements
        in `x0`.
    band : tuple
        If set to a two-sequence containing the number of sub- and
        super-diagonals within the band of the Jacobi matrix, the
        Jacobi matrix is considered banded (only for ``fprime=None``).
    epsfcn : float
        A suitable step length for the forward-difference
        approximation of the Jacobian (for ``fprime=None``). If
        `epsfcn` is less than the machine precision, it is assumed
        that the relative errors in the functions are of the order of
        the machine precision.
    factor : float
        A parameter determining the initial step bound
        (``factor * || diag * x||``).  Should be in the interval
        ``(0.1, 100)``.
    diag : sequence
        N positive entries that serve as a scale factors for the
        variables.

    Notes
    -----
    ``fsolve`` is a wrapper around MINPACK's hybrd and hybrj algorithms.

    """
    x0 = array(x0, ndmin=1)
    n = len(x0)
    if type(args) != type(()): args = (args,)
    _check_func('fsolve', 'func', func, x0, args, n, (n,))
    Dfun = fprime
    if Dfun is None:
        if band is None:
            ml, mu = -10,-10
        else:
            ml, mu = band[:2]
        if (maxfev == 0):
            maxfev = 200*(n + 1)
        retval = _minpack._hybrd(func, x0, args, full_output, xtol,
                maxfev, ml, mu, epsfcn, factor, diag)
    else:
        _check_func('fsolve', 'fprime', Dfun, x0, args, n, (n,n))
        if (maxfev == 0):
            maxfev = 100*(n + 1)
        retval = _minpack._hybrj(func, Dfun, x0, args, full_output,
                col_deriv, xtol, maxfev, factor,diag)

    errors = {0:["Improper input parameters were entered.",TypeError],
              1:["The solution converged.", None],
              2:["The number of calls to function has "
                 "reached maxfev = %d." % maxfev, ValueError],
              3:["xtol=%f is too small, no further improvement "
                 "in the approximate\n  solution "
                 "is possible." % xtol, ValueError],
              4:["The iteration is not making good progress, as measured "
                 "by the \n  improvement from the last five "
                 "Jacobian evaluations.", ValueError],
              5:["The iteration is not making good progress, "
                 "as measured by the \n  improvement from the last "
                 "ten iterations.", ValueError],
              'unknown': ["An error occurred.", TypeError]}

    info = retval[-1]    # The FORTRAN return value
    if (info != 1 and not full_output):
        if info in [2,3,4,5]:
            msg = errors[info][0]
            warnings.warn(msg, RuntimeWarning)
        else:
            try:
                raise errors[info][1](errors[info][0])
            except KeyError:
                raise errors['unknown'][1](errors['unknown'][0])

    if full_output:
        try:
            return retval + (errors[info][0],)  # Return all + the message
        except KeyError:
            return retval + (errors['unknown'][0],)
    else:
        return retval[0]


def leastsq(func, x0, args=(), Dfun=None, full_output=0,
            col_deriv=0, ftol=1.49012e-8, xtol=1.49012e-8,
            gtol=0.0, maxfev=0, epsfcn=0.0, factor=100, diag=None):
    """
    Minimize the sum of squares of a set of equations.

    ::

        x = arg min(sum(func(y)**2,axis=0))
                 y

    Parameters
    ----------
    func : callable
        should take at least one (possibly length N vector) argument and
        returns M floating point numbers.
    x0 : ndarray
        The starting estimate for the minimization.
    args : tuple
        Any extra arguments to func are placed in this tuple.
    Dfun : callable
        A function or method to compute the Jacobian of func with derivatives
        across the rows. If this is None, the Jacobian will be estimated.
    full_output : bool
        non-zero to return all optional outputs.
    col_deriv : bool
        non-zero to specify that the Jacobian function computes derivatives
        down the columns (faster, because there is no transpose operation).
    ftol : float
        Relative error desired in the sum of squares.
    xtol : float
        Relative error desired in the approximate solution.
    gtol : float
        Orthogonality desired between the function vector and the columns of
        the Jacobian.
    maxfev : int
        The maximum number of calls to the function. If zero, then 100*(N+1) is
        the maximum where N is the number of elements in x0.
    epsfcn : float
        A suitable step length for the forward-difference approximation of the
        Jacobian (for Dfun=None). If epsfcn is less than the machine precision,
        it is assumed that the relative errors in the functions are of the
        order of the machine precision.
    factor : float
        A parameter determining the initial step bound
        (``factor * || diag * x||``). Should be in interval ``(0.1, 100)``.
    diag : sequence
        N positive entries that serve as a scale factors for the variables.

    Returns
    -------
    x : ndarray
        The solution (or the result of the last iteration for an unsuccessful
        call).
    cov_x : ndarray
        Uses the fjac and ipvt optional outputs to construct an
        estimate of the jacobian around the solution.  ``None`` if a
        singular matrix encountered (indicates very flat curvature in
        some direction).  This matrix must be multiplied by the
        residual standard deviation to get the covariance of the
        parameter estimates -- see curve_fit.
    infodict : dict
        a dictionary of optional outputs with the key s::

            - 'nfev' : the number of function calls
            - 'fvec' : the function evaluated at the output
            - 'fjac' : A permutation of the R matrix of a QR
                     factorization of the final approximate
                     Jacobian matrix, stored column wise.
                     Together with ipvt, the covariance of the
                     estimate can be approximated.
            - 'ipvt' : an integer array of length N which defines
                     a permutation matrix, p, such that
                     fjac*p = q*r, where r is upper triangular
                     with diagonal elements of nonincreasing
                     magnitude. Column j of p is column ipvt(j)
                     of the identity matrix.
            - 'qtf'  : the vector (transpose(q) * fvec).

    mesg : str
        A string message giving information about the cause of failure.
    ier : int
        An integer flag.  If it is equal to 1, 2, 3 or 4, the solution was
        found.  Otherwise, the solution was not found. In either case, the
        optional output variable 'mesg' gives more information.

    Notes
    -----
    "leastsq" is a wrapper around MINPACK's lmdif and lmder algorithms.

    cov_x is a Jacobian approximation to the Hessian of the least squares
    objective function.
    This approximation assumes that the objective function is based on the
    difference between some observed target data (ydata) and a (non-linear)
    function of the parameters `f(xdata, params)` ::

           func(params) = ydata - f(xdata, params)

    so that the objective function is ::

           min   sum((ydata - f(xdata, params))**2, axis=0)
         params

    """
    x0 = array(x0, ndmin=1)
    n = len(x0)
    if type(args) != type(()):
        args = (args,)
    m = _check_func('leastsq', 'func', func, x0, args, n)[0]
    if n > m:
        raise TypeError('Improper input: N=%s must not exceed M=%s' % (n,m))
    if Dfun is None:
        if (maxfev == 0):
            maxfev = 200*(n + 1)
        retval = _minpack._lmdif(func, x0, args, full_output, ftol, xtol,
                gtol, maxfev, epsfcn, factor, diag)
    else:
        if col_deriv:
            _check_func('leastsq', 'Dfun', Dfun, x0, args, n, (n,m))
        else:
            _check_func('leastsq', 'Dfun', Dfun, x0, args, n, (m,n))
        if (maxfev == 0):
            maxfev = 100*(n + 1)
        retval = _minpack._lmder(func, Dfun, x0, args, full_output, col_deriv,
                ftol, xtol, gtol, maxfev, factor, diag)

    errors = {0:["Improper input parameters.", TypeError],
              1:["Both actual and predicted relative reductions "
                 "in the sum of squares\n  are at most %f" % ftol, None],
              2:["The relative error between two consecutive "
                 "iterates is at most %f" % xtol, None],
              3:["Both actual and predicted relative reductions in "
                 "the sum of squares\n  are at most %f and the "
                 "relative error between two consecutive "
                 "iterates is at \n  most %f" % (ftol,xtol), None],
              4:["The cosine of the angle between func(x) and any "
                 "column of the\n  Jacobian is at most %f in "
                 "absolute value" % gtol, None],
              5:["Number of calls to function has reached "
                 "maxfev = %d." % maxfev, ValueError],
              6:["ftol=%f is too small, no further reduction "
                 "in the sum of squares\n  is possible.""" % ftol, ValueError],
              7:["xtol=%f is too small, no further improvement in "
                 "the approximate\n  solution is possible." % xtol, ValueError],
              8:["gtol=%f is too small, func(x) is orthogonal to the "
                 "columns of\n  the Jacobian to machine "
                 "precision." % gtol, ValueError],
              'unknown':["Unknown error.", TypeError]}

    info = retval[-1]    # The FORTRAN return value

    if (info not in [1,2,3,4] and not full_output):
        if info in [5,6,7,8]:
            warnings.warn(errors[info][0], RuntimeWarning)
        else:
            try:
                raise errors[info][1](errors[info][0])
            except KeyError:
                raise errors['unknown'][1](errors['unknown'][0])

    mesg = errors[info][0]
    if full_output:
        cov_x = None
        if info in [1,2,3,4]:
            from numpy.dual import inv
            from numpy.linalg import LinAlgError
            perm = take(eye(n),retval[1]['ipvt']-1,0)
            r = triu(transpose(retval[1]['fjac'])[:n,:])
            R = dot(r, perm)
            try:
                cov_x = inv(dot(transpose(R),R))
            except LinAlgError:
                pass
        return (retval[0], cov_x) + retval[1:-1] + (mesg, info)
    else:
        return (retval[0], info)

def _general_function(params, xdata, ydata, function):
    return function(xdata, *params) - ydata

def _weighted_general_function(params, xdata, ydata, function, weights):
    return weights * (function(xdata, *params) - ydata)

def curve_fit(f, xdata, ydata, p0=None, sigma=None, **kw):
    """
    Use non-linear least squares to fit a function, f, to data.

    Assumes ``ydata = f(xdata, *params) + eps``

    Parameters
    ----------
    f : callable
        The model function, f(x, ...).  It must take the independent
        variable as the first argument and the parameters to fit as
        separate remaining arguments.
    xdata : An N-length sequence or an (k,N)-shaped array
        for functions with k predictors.
        The independent variable where the data is measured.
    ydata : N-length sequence
        The dependent data --- nominally f(xdata, ...)
    p0 : None, scalar, or M-length sequence
        Initial guess for the parameters.  If None, then the initial
        values will all be 1 (if the number of parameters for the function
        can be determined using introspection, otherwise a ValueError
        is raised).
    sigma : None or N-length sequence
        If not None, it represents the standard-deviation of ydata.
        This vector, if given, will be used as weights in the
        least-squares problem.

    Returns
    -------
    popt : array
        Optimal values for the parameters so that the sum of the squared error
        of ``f(xdata, *popt) - ydata`` is minimized
    pcov : 2d array
        The estimated covariance of popt.  The diagonals provide the variance
        of the parameter estimate.

    See Also
    --------
    leastsq

    Notes
    -----
    The algorithm uses the Levenburg-Marquardt algorithm through `leastsq`.
    Additional keyword arguments are passed directly to that algorithm.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.optimize import curve_fit
    >>> def func(x, a, b, c):
    ...     return a*np.exp(-b*x) + c

    >>> x = np.linspace(0,4,50)
    >>> y = func(x, 2.5, 1.3, 0.5)
    >>> yn = y + 0.2*np.random.normal(size=len(x))

    >>> popt, pcov = curve_fit(func, x, yn)

    """
    if p0 is None:
        # determine number of parameters by inspecting the function
        import inspect
        args, varargs, varkw, defaults = inspect.getargspec(f)
        if len(args) < 2:
            msg = "Unable to determine number of fit parameters."
            raise ValueError(msg)
        if 'self' in args:
            p0 = [1.0] * (len(args)-2)
        else:
            p0 = [1.0] * (len(args)-1)

    if isscalar(p0):
        p0 = array([p0])

    args = (xdata, ydata, f)
    if sigma is None:
        func = _general_function
    else:
        func = _weighted_general_function
        args += (1.0/asarray(sigma),)

    # Remove full_output from kw, otherwise we're passing it in twice.
    return_full = kw.pop('full_output', False)
    res = leastsq(func, p0, args=args, full_output=1, **kw)
    (popt, pcov, infodict, errmsg, ier) = res

    if ier not in [1,2,3,4]:
        msg = "Optimal parameters not found: " + errmsg
        raise RuntimeError(msg)

    if (len(ydata) > len(p0)) and pcov is not None:
        s_sq = (func(popt, *args)**2).sum()/(len(ydata)-len(p0))
        pcov = pcov * s_sq
    else:
        pcov = inf

    if return_full:
        return popt, pcov, infodict, errmsg, ier
    else:
        return popt, pcov

def check_gradient(fcn, Dfcn, x0, args=(), col_deriv=0):
    """Perform a simple check on the gradient for correctness.

    """

    x = atleast_1d(x0)
    n = len(x)
    x = x.reshape((n,))
    fvec = atleast_1d(fcn(x,*args))
    m = len(fvec)
    fvec = fvec.reshape((m,))
    ldfjac = m
    fjac = atleast_1d(Dfcn(x,*args))
    fjac = fjac.reshape((m,n))
    if col_deriv == 0:
        fjac = transpose(fjac)

    xp = zeros((n,), float)
    err = zeros((m,), float)
    fvecp = None
    _minpack._chkder(m, n, x, fvec, fjac, ldfjac, xp, fvecp, 1, err)

    fvecp = atleast_1d(fcn(xp,*args))
    fvecp = fvecp.reshape((m,))
    _minpack._chkder(m, n, x, fvec, fjac, ldfjac, xp, fvecp, 2, err)

    good = (product(greater(err, 0.5), axis=0))

    return (good, err)


# Steffensen's Method using Aitken's Del^2 convergence acceleration.
def fixed_point(func, x0, args=(), xtol=1e-8, maxiter=500):
    """Find the point where func(x) == x

    Given a function of one or more variables and a starting point, find a
    fixed-point of the function: i.e. where func(x)=x.

    Uses Steffensen's Method using Aitken's Del^2 convergence acceleration.
    See Burden, Faires, "Numerical Analysis", 5th edition, pg. 80

    Examples
    --------
    >>> from numpy import sqrt, array
    >>> from scipy.optimize import fixed_point
    >>> def func(x, c1, c2):
            return sqrt(c1/(x+c2))
    >>> c1 = array([10,12.])
    >>> c2 = array([3, 5.])
    >>> fixed_point(func, [1.2, 1.3], args=(c1,c2))
    array([ 1.4920333 ,  1.37228132])

    """
    if not isscalar(x0):
        x0 = asarray(x0)
        p0 = x0
        for iter in range(maxiter):
            p1 = func(p0, *args)
            p2 = func(p1, *args)
            d = p2 - 2.0 * p1 + p0
            p = where(d == 0, p2, p0 - (p1 - p0)*(p1 - p0) / d)
            relerr = where(p0 == 0, p, (p-p0)/p0)
            if all(abs(relerr) < xtol):
                return p
            p0 = p
    else:
        p0 = x0
        for iter in range(maxiter):
            p1 = func(p0, *args)
            p2 = func(p1, *args)
            d = p2 - 2.0 * p1 + p0
            if d == 0.0:
                return p2
            else:
                p = p0 - (p1 - p0)*(p1 - p0) / d
            if p0 == 0:
                relerr = p
            else:
                relerr = (p - p0)/p0
            if abs(relerr) < xtol:
                return p
            p0 = p
    msg = "Failed to converge after %d iterations, value is %s" % (maxiter, p)
    raise RuntimeError(msg)
