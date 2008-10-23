import _minpack

from numpy import atleast_1d, dot, take, triu, shape, eye, \
                  transpose, zeros, product, greater, array, \
                  all, where, isscalar, asarray

error = _minpack.error

__all__ = ['fsolve', 'leastsq', 'newton', 'fixed_point','bisection']

def check_func(thefunc, x0, args, numinputs, output_shape=None):
    res = atleast_1d(thefunc(*((x0[:numinputs],)+args)))
    if (output_shape is not None) and (shape(res) != output_shape):
        if (output_shape[0] != 1):
            if len(output_shape) > 1:
                if output_shape[1] == 1:
                    return shape(res)
            raise TypeError, "There is a mismatch between the input and output shape of %s." % thefunc.func_name
    return shape(res)


def fsolve(func,x0,args=(),fprime=None,full_output=0,col_deriv=0,xtol=1.49012e-8,maxfev=0,band=None,epsfcn=0.0,factor=100,diag=None, warning=True):
    """Find the roots of a function.

  Description:

    Return the roots of the (non-linear) equations defined by
    func(x)=0 given a starting estimate.

  Inputs:

    func -- A Python function or method which takes at least one
            (possibly vector) argument.
    x0 -- The starting estimate for the roots of func(x)=0.
    args -- Any extra arguments to func are placed in this tuple.
    fprime -- A function or method to compute the Jacobian of func with
            derivatives across the rows. If this is None, the
            Jacobian will be estimated.
    full_output -- non-zero to return the optional outputs.
    col_deriv -- non-zero to specify that the Jacobian function
                 computes derivatives down the columns (faster, because
                 there is no transpose operation).
    warning -- True to print a warning message when the call is
                unsuccessful; False to suppress the warning message.
  Outputs: (x, {infodict, ier, mesg})

    x -- the solution (or the result of the last iteration for an
         unsuccessful call.

    infodict -- a dictionary of optional outputs with the keys:
                'nfev' : the number of function calls
                'njev' : the number of jacobian calls
                'fvec' : the function evaluated at the output
                'fjac' : the orthogonal matrix, q, produced by the
                         QR factorization of the final approximate
                         Jacobian matrix, stored column wise.
                'r'    : upper triangular matrix produced by QR
                         factorization of same matrix.
                'qtf'  : the vector (transpose(q) * fvec).
    ier -- an integer flag.  If it is equal to 1 the solution was
           found.  If it is not equal to 1, the solution was not
           found and the following message gives more information.
    mesg -- a string message giving information about the cause of
            failure.

  Extended Inputs:

   xtol -- The calculation will terminate if the relative error
           between two consecutive iterates is at most xtol.
   maxfev -- The maximum number of calls to the function. If zero,
             then 100*(N+1) is the maximum where N is the number
             of elements in x0.
   band -- If set to a two-sequence containing the number of sub-
           and superdiagonals within the band of the Jacobi matrix,
           the Jacobi matrix is considered banded (only for fprime=None).
   epsfcn -- A suitable step length for the forward-difference
             approximation of the Jacobian (for fprime=None). If
             epsfcn is less than the machine precision, it is assumed
             that the relative errors in the functions are of
             the order of the machine precision.
   factor -- A parameter determining the initial step bound
             (factor * || diag * x||). Should be in interval (0.1,100).
   diag -- A sequency of N positive entries that serve as a
           scale factors for the variables.

  Remarks:

    "fsolve" is a wrapper around MINPACK's hybrd and hybrj algorithms.

  See also:

      scikits.openopt, which offers a unified syntax to call this and other solvers

      fmin, fmin_powell, fmin_cg,
             fmin_bfgs, fmin_ncg -- multivariate local optimizers
      leastsq -- nonlinear least squares minimizer

      fmin_l_bfgs_b, fmin_tnc,
             fmin_cobyla -- constrained multivariate optimizers

      anneal, brute -- global optimizers

      fminbound, brent, golden, bracket -- local scalar minimizers

      brentq, brenth, ridder, bisect, newton -- one-dimensional root-finding

      fixed_point -- scalar and vector fixed-point finder

    """
    x0 = array(x0,ndmin=1)
    n = len(x0)
    if type(args) != type(()): args = (args,)
    check_func(func,x0,args,n,(n,))
    Dfun = fprime
    if Dfun is None:
        if band is None:
            ml,mu = -10,-10
        else:
            ml,mu = band[:2]
        if (maxfev == 0):
            maxfev = 200*(n+1)
        retval = _minpack._hybrd(func,x0,args,full_output,xtol,maxfev,ml,mu,epsfcn,factor,diag)
    else:
        check_func(Dfun,x0,args,n,(n,n))
        if (maxfev == 0):
            maxfev = 100*(n+1)
        retval = _minpack._hybrj(func,Dfun,x0,args,full_output,col_deriv,xtol,maxfev,factor,diag)

    errors = {0:["Improper input parameters were entered.",TypeError],
              1:["The solution converged.",None],
              2:["The number of calls to function has reached maxfev = %d." % maxfev, ValueError],
              3:["xtol=%f is too small, no further improvement in the approximate\n  solution is possible." % xtol, ValueError],
              4:["The iteration is not making good progress, as measured by the \n  improvement from the last five Jacobian evaluations.", ValueError],
              5:["The iteration is not making good progress, as measured by the \n  improvement from the last ten iterations.", ValueError],
              'unknown': ["An error occurred.", TypeError]}

    info = retval[-1]    # The FORTRAN return value
    if (info != 1 and not full_output):
        if info in [2,3,4,5]:
            if warning:  print "Warning: " + errors[info][0]
        else:
            try:
                raise errors[info][1], errors[info][0]
            except KeyError:
                raise errors['unknown'][1], errors['unknown'][0]

    if n == 1:
        retval = (retval[0][0],) + retval[1:]

    if full_output:
        try:
            return retval + (errors[info][0],)  # Return all + the message
        except KeyError:
            return retval + (errors['unknown'][0],)
    else:
        return retval[0]


def leastsq(func,x0,args=(),Dfun=None,full_output=0,col_deriv=0,ftol=1.49012e-8,xtol=1.49012e-8,gtol=0.0,maxfev=0,epsfcn=0.0,factor=100,diag=None,warning=True):
    """Minimize the sum of squares of a set of equations.

  Description:

    Return the point which minimizes the sum of squares of M
    (non-linear) equations in N unknowns given a starting estimate, x0,
    using a modification of the Levenberg-Marquardt algorithm.

                    x = arg min(sum(func(y)**2,axis=0))
                             y

  Inputs:

    func -- A Python function or method which takes at least one
            (possibly length N vector) argument and returns M
            floating point numbers.
    x0 -- The starting estimate for the minimization.
    args -- Any extra arguments to func are placed in this tuple.
    Dfun -- A function or method to compute the Jacobian of func with
            derivatives across the rows. If this is None, the
            Jacobian will be estimated.
    full_output -- non-zero to return all optional outputs.
    col_deriv -- non-zero to specify that the Jacobian function
                 computes derivatives down the columns (faster, because
                 there is no transpose operation).
        warning -- True to print a warning message when the call is
             unsuccessful; False to suppress the warning message.

  Outputs: (x, {cov_x, infodict, mesg}, ier)

    x -- the solution (or the result of the last iteration for an
         unsuccessful call.

    cov_x -- uses the fjac and ipvt optional outputs to construct an
             estimate of the covariance matrix of the solution.
             None if a singular matrix encountered (indicates
             infinite covariance in some direction).
    infodict -- a dictionary of optional outputs with the keys:
                'nfev' : the number of function calls
                'fvec' : the function evaluated at the output
                'fjac' : A permutation of the R matrix of a QR
                         factorization of the final approximate
                         Jacobian matrix, stored column wise.
                         Together with ipvt, the covariance of the
                         estimate can be approximated.
                'ipvt' : an integer array of length N which defines
                         a permutation matrix, p, such that
                         fjac*p = q*r, where r is upper triangular
                         with diagonal elements of nonincreasing
                         magnitude. Column j of p is column ipvt(j)
                         of the identity matrix.
                'qtf'  : the vector (transpose(q) * fvec).
    mesg -- a string message giving information about the cause of failure.
    ier -- an integer flag.  If it is equal to 1, 2, 3 or 4, the
           solution was found.  Otherwise, the solution was not
           found. In either case, the optional output variable 'mesg'
           gives more information.


  Extended Inputs:

   ftol -- Relative error desired in the sum of squares.
   xtol -- Relative error desired in the approximate solution.
   gtol -- Orthogonality desired between the function vector
           and the columns of the Jacobian.
   maxfev -- The maximum number of calls to the function. If zero,
             then 100*(N+1) is the maximum where N is the number
             of elements in x0.
   epsfcn -- A suitable step length for the forward-difference
             approximation of the Jacobian (for Dfun=None). If
             epsfcn is less than the machine precision, it is assumed
             that the relative errors in the functions are of
             the order of the machine precision.
   factor -- A parameter determining the initial step bound
             (factor * || diag * x||). Should be in interval (0.1,100).
   diag -- A sequency of N positive entries that serve as a
           scale factors for the variables.

  Remarks:

    "leastsq" is a wrapper around MINPACK's lmdif and lmder algorithms.

  See also:

      scikits.openopt, which offers a unified syntax to call this and other solvers

      fmin, fmin_powell, fmin_cg,
             fmin_bfgs, fmin_ncg -- multivariate local optimizers

      fmin_l_bfgs_b, fmin_tnc,
             fmin_cobyla -- constrained multivariate optimizers

      anneal, brute -- global optimizers

      fminbound, brent, golden, bracket -- local scalar minimizers

      fsolve -- n-dimenstional root-finding

      brentq, brenth, ridder, bisect, newton -- one-dimensional root-finding

      fixed_point -- scalar and vector fixed-point finder

    """
    x0 = array(x0,ndmin=1)
    n = len(x0)
    if type(args) != type(()): args = (args,)
    m = check_func(func,x0,args,n)[0]
    if Dfun is None:
        if (maxfev == 0):
            maxfev = 200*(n+1)
        retval = _minpack._lmdif(func,x0,args,full_output,ftol,xtol,gtol,maxfev,epsfcn,factor,diag)
    else:
        if col_deriv:
            check_func(Dfun,x0,args,n,(n,m))
        else:
            check_func(Dfun,x0,args,n,(m,n))
        if (maxfev == 0):
            maxfev = 100*(n+1)
        retval = _minpack._lmder(func,Dfun,x0,args,full_output,col_deriv,ftol,xtol,gtol,maxfev,factor,diag)

    errors = {0:["Improper input parameters.", TypeError],
              1:["Both actual and predicted relative reductions in the sum of squares\n  are at most %f" % ftol, None],
              2:["The relative error between two consecutive iterates is at most %f" % xtol, None],
              3:["Both actual and predicted relative reductions in the sum of squares\n  are at most %f and the relative error between two consecutive iterates is at \n  most %f" % (ftol,xtol), None],
              4:["The cosine of the angle between func(x) and any column of the\n  Jacobian is at most %f in absolute value" % gtol, None],
              5:["Number of calls to function has reached maxfev = %d." % maxfev, ValueError],
              6:["ftol=%f is too small, no further reduction in the sum of squares\n  is possible.""" % ftol, ValueError],
              7:["xtol=%f is too small, no further improvement in the approximate\n  solution is possible." % xtol, ValueError],
              8:["gtol=%f is too small, func(x) is orthogonal to the columns of\n  the Jacobian to machine precision." % gtol, ValueError],
              'unknown':["Unknown error.", TypeError]}

    info = retval[-1]    # The FORTRAN return value

    if (info not in [1,2,3,4] and not full_output):
        if info in [5,6,7,8]:
            if warning:  print "Warning: " + errors[info][0]
        else:
            try:
                raise errors[info][1], errors[info][0]
            except KeyError:
                raise errors['unknown'][1], errors['unknown'][0]

    if n == 1:
        retval = (retval[0][0],) + retval[1:]

    mesg = errors[info][0]
    if full_output:
        from numpy.dual import inv
        from numpy.linalg import LinAlgError
        perm = take(eye(n),retval[1]['ipvt']-1,0)
        r = triu(transpose(retval[1]['fjac'])[:n,:])
        R = dot(r, perm)
        try:
            cov_x = inv(dot(transpose(R),R))
        except LinAlgError:
            cov_x = None
        return (retval[0], cov_x) + retval[1:-1] + (mesg,info)
    else:
        return (retval[0], info)


def check_gradient(fcn,Dfcn,x0,args=(),col_deriv=0):
    """Perform a simple check on the gradient for correctness.
    """

    x = atleast_1d(x0)
    n = len(x)
    x=x.reshape((n,))
    fvec = atleast_1d(fcn(x,*args))
    m = len(fvec)
    fvec=fvec.reshape((m,))
    ldfjac = m
    fjac = atleast_1d(Dfcn(x,*args))
    fjac=fjac.reshape((m,n))
    if col_deriv == 0:
        fjac = transpose(fjac)

    xp = zeros((n,), float)
    err = zeros((m,), float)
    fvecp = None
    _minpack._chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,1,err)

    fvecp = atleast_1d(fcn(xp,*args))
    fvecp=fvecp.reshape((m,))
    _minpack._chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,2,err)

    good = (product(greater(err,0.5),axis=0))

    return (good,err)


# Netwon-Raphson method
def newton(func, x0, fprime=None, args=(), tol=1.48e-8, maxiter=50):
    """Given a function of a single variable and a starting point,
    find a nearby zero using Newton-Raphson.

    fprime is the derivative of the function.  If not given, the
    Secant method is used.

    See also:

      fmin, fmin_powell, fmin_cg,
             fmin_bfgs, fmin_ncg -- multivariate local optimizers
      leastsq -- nonlinear least squares minimizer

      fmin_l_bfgs_b, fmin_tnc,
             fmin_cobyla -- constrained multivariate optimizers

      anneal, brute -- global optimizers

      fminbound, brent, golden, bracket -- local scalar minimizers

      fsolve -- n-dimenstional root-finding

      brentq, brenth, ridder, bisect, newton -- one-dimensional root-finding

      fixed_point -- scalar and vector fixed-point finder

    """

    if fprime is not None:
        p0 = x0
        for iter in range(maxiter):
            myargs = (p0,)+args
            fval = func(*myargs)
            fpval = fprime(*myargs)
            if fpval == 0:
                print "Warning: zero-derivative encountered."
                return p0
            p = p0 - func(*myargs)/fprime(*myargs)
            if abs(p-p0) < tol:
                return p
            p0 = p
    else: # Secant method
        p0 = x0
        p1 = x0*(1+1e-4)
        q0 = func(*((p0,)+args))
        q1 = func(*((p1,)+args))
        for iter in range(maxiter):
            if q1 == q0:
                if p1 != p0:
                    print "Tolerance of %s reached" % (p1-p0)
                return (p1+p0)/2.0
            else:
                p = p1 - q1*(p1-p0)/(q1-q0)
            if abs(p-p1) < tol:
                return p
            p0 = p1
            q0 = q1
            p1 = p
            q1 = func(*((p1,)+args))
    raise RuntimeError, "Failed to converge after %d iterations, value is %s" % (maxiter,p)


# Steffensen's Method using Aitken's Del^2 convergence acceleration.
def fixed_point(func, x0, args=(), xtol=1e-8, maxiter=500):
    """Find the point where func(x) == x

    Given a function of one or more variables and a starting point, find a
    fixed-point of the function: i.e. where func(x)=x.

    Uses Steffensen's Method using Aitken's Del^2 convergence acceleration.
    See Burden, Faires, "Numerical Analysis", 5th edition, pg. 80

    Example
    -------
    >>> from numpy import sqrt, array
    >>> from scipy.optimize import fixed_point
    >>> def func(x, c1, c2):
            return sqrt(c1/(x+c2))
    >>> c1 = array([10,12.])
    >>> c2 = array([3, 5.])
    >>> fixed_point(func, [1.2, 1.3], args=(c1,c2))
    array([ 1.4920333 ,  1.37228132])

    See also:

      fmin, fmin_powell, fmin_cg,
             fmin_bfgs, fmin_ncg -- multivariate local optimizers
      leastsq -- nonlinear least squares minimizer

      fmin_l_bfgs_b, fmin_tnc,
             fmin_cobyla -- constrained multivariate optimizers

      anneal, brute -- global optimizers

      fminbound, brent, golden, bracket -- local scalar minimizers

      fsolve -- n-dimenstional root-finding

      brentq, brenth, ridder, bisect, newton -- one-dimensional root-finding

    """
    if not isscalar(x0):
        x0 = asarray(x0)
        p0 = x0
        for iter in range(maxiter):
            p1 = func(p0, *args)
            p2 = func(p1, *args)
            d = p2 - 2.0 * p1 + p0
            p = where(d == 0, p2, p0 - (p1 - p0)*(p1-p0) / d)
            relerr = where(p0 == 0, p, (p-p0)/p0)
            if all(relerr < xtol):
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
                p = p0 - (p1 - p0)*(p1-p0) / d
            if p0 == 0:
                relerr = p
            else:
                relerr = (p-p0)/p0
            if relerr < xtol:
                return p
            p0 = p
    raise RuntimeError, "Failed to converge after %d iterations, value is %s" % (maxiter,p)


def bisection(func, a, b, args=(), xtol=1e-10, maxiter=400):
    """Bisection root-finding method.  Given a function and an interval with
    func(a) * func(b) < 0, find the root between a and b.

    See also:

      fmin, fmin_powell, fmin_cg,
             fmin_bfgs, fmin_ncg -- multivariate local optimizers
      leastsq -- nonlinear least squares minimizer

      fmin_l_bfgs_b, fmin_tnc,
             fmin_cobyla -- constrained multivariate optimizers

      anneal, brute -- global optimizers

      fminbound, brent, golden, bracket -- local scalar minimizers

      fsolve -- n-dimenstional root-finding

      brentq, brenth, ridder, bisect, newton -- one-dimensional root-finding

      fixed_point -- scalar and vector fixed-point finder

    """
    i = 1
    eva = func(a,*args)
    evb = func(b,*args)
    assert (eva*evb < 0), "Must start with interval with func(a) * func(b) <0"
    while i<=maxiter:
        dist = (b-a)/2.0
        p = a + dist
        if dist < xtol:
            return p
        ev = func(p,*args)
        if ev == 0:
            return p
        i += 1
        if ev*eva > 0:
            a = p
            eva = ev
        else:
            b = p
    print "Warning: Method failed after %d iterations." % maxiter
    return p
