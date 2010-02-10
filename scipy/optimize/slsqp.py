"""This module implements the Sequential Least SQuares Programming optimization
algorithm (SLSQP), orginally developed by Dieter Kraft.

See http://www.netlib.org/toms/733

"""

__all__ = ['approx_jacobian','fmin_slsqp']

from _slsqp import slsqp
from numpy import zeros, array, linalg, append, asfarray, concatenate, finfo, \
                  sqrt, vstack
from optimize import approx_fprime, wrap_function

__docformat__ = "restructuredtext en"

_epsilon = sqrt(finfo(float).eps)

def approx_jacobian(x,func,epsilon,*args):
    """Approximate the Jacobian matrix of a callable function.

    Parameters
    ----------
    x : array_like
        The state vector at which to compute the Jacobian matrix.
    func : callable f(x, *args)
        The vector-valued function.
    epsilon : float\
        The peturbation used to determine the partial derivatives.
    *args : tuple
        Additional arguments passed to func.

    Returns
    -------
    An array of dimensions ``(lenf, lenx)`` where ``lenf`` is the length
    of the outputs of `func`, and ``lenx`` is the number of elements in
    `x`.

    Notes
    -----
    The approximation is done using forward differences.

    """
    x0 = asfarray(x)
    f0 = func(*((x0,)+args))
    jac = zeros([len(x0),len(f0)])
    dx = zeros(len(x0))
    for i in range(len(x0)):
        dx[i] = epsilon
        jac[i] = (func(*((x0+dx,)+args)) - f0)/epsilon
        dx[i] = 0.0
    return jac.transpose()


def fmin_slsqp( func, x0 , eqcons=[], f_eqcons=None, ieqcons=[], f_ieqcons=None,
                bounds = [], fprime = None, fprime_eqcons=None,
                fprime_ieqcons=None, args = (), iter = 100, acc = 1.0E-6,
                iprint = 1, full_output = 0, epsilon = _epsilon ):
    """
    Minimize a function using Sequential Least SQuares Programming

    Python interface function for the SLSQP Optimization subroutine
    originally implemented by Dieter Kraft.

    Parameters
    ----------
    func : callable f(x,*args)
        Objective function.
    x0 : ndarray of float
        Initial guess for the independent variable(s).
    eqcons : list
        A list of functions of length n such that
        eqcons[j](x0,*args) == 0.0 in a successfully optimized
        problem.
    f_eqcons : callable f(x,*args)
        Returns an array in which each element must equal 0.0 in a
        successfully optimized problem.  If f_eqcons is specified,
        eqcons is ignored.
    ieqcons : list
        A list of functions of length n such that
        ieqcons[j](x0,*args) >= 0.0 in a successfully optimized
        problem.
    f_ieqcons : callable f(x0,*args)
        Returns an array in which each element must be greater or
        equal to 0.0 in a successfully optimized problem.  If
        f_ieqcons is specified, ieqcons is ignored.
    bounds : list
        A list of tuples specifying the lower and upper bound
        for each independent variable [(xl0, xu0),(xl1, xu1),...]
    fprime : callable `f(x,*args)`
        A function that evaluates the partial derivatives of func.
    fprime_eqcons : callable `f(x,*args)`
        A function of the form `f(x, *args)` that returns the m by n
        array of equality constraint normals.  If not provided,
        the normals will be approximated. The array returned by
        fprime_eqcons should be sized as ( len(eqcons), len(x0) ).
    fprime_ieqcons : callable `f(x,*args)`
        A function of the form `f(x, *args)` that returns the m by n
        array of inequality constraint normals.  If not provided,
        the normals will be approximated. The array returned by
        fprime_ieqcons should be sized as ( len(ieqcons), len(x0) ).
    args : sequence
        Additional arguments passed to func and fprime.
    iter : int
        The maximum number of iterations.
    acc : float
        Requested accuracy.
    iprint : int
        The verbosity of fmin_slsqp :

        * iprint <= 0 : Silent operation
        * iprint == 1 : Print summary upon completion (default)
        * iprint >= 2 : Print status of each iterate and summary
    full_output : bool
        If False, return only the minimizer of func (default).
        Otherwise, output final objective function and summary
        information.
    epsilon : float
        The step size for finite-difference derivative estimates.

    Returns
    -------
    x : ndarray of float
        The final minimizer of func.
    fx : ndarray of float, if full_output is true
        The final value of the objective function.
    its : int, if full_output is true
        The number of iterations.
    imode : int, if full_output is true
        The exit mode from the optimizer (see below).
    smode : string, if full_output is true
        Message describing the exit mode from the optimizer.

    Notes
    -----
    Exit modes are defined as follows ::

        -1 : Gradient evaluation required (g & a)
         0 : Optimization terminated successfully.
         1 : Function evaluation required (f & c)
         2 : More equality constraints than independent variables
         3 : More than 3*n iterations in LSQ subproblem
         4 : Inequality constraints incompatible
         5 : Singular matrix E in LSQ subproblem
         6 : Singular matrix C in LSQ subproblem
         7 : Rank-deficient equality constraint subproblem HFTI
         8 : Positive directional derivative for linesearch
         9 : Iteration limit exceeded

    Examples
    --------
    Examples are given :ref:`in the tutorial <tutorial-sqlsp>`.

    """

    exit_modes = { -1 : "Gradient evaluation required (g & a)",
                    0 : "Optimization terminated successfully.",
                    1 : "Function evaluation required (f & c)",
                    2 : "More equality constraints than independent variables",
                    3 : "More than 3*n iterations in LSQ subproblem",
                    4 : "Inequality constraints incompatible",
                    5 : "Singular matrix E in LSQ subproblem",
                    6 : "Singular matrix C in LSQ subproblem",
                    7 : "Rank-deficient equality constraint subproblem HFTI",
                    8 : "Positive directional derivative for linesearch",
                    9 : "Iteration limit exceeded" }

    # Now do a lot of function wrapping

    # Wrap func
    feval, func = wrap_function(func, args)
    # Wrap fprime, if provided, or approx_fprime if not
    if fprime:
        geval, fprime = wrap_function(fprime,args)
    else:
        geval, fprime = wrap_function(approx_fprime,(func,epsilon))

    if f_eqcons:
        # Equality constraints provided via f_eqcons
        ceval, f_eqcons = wrap_function(f_eqcons,args)
        if fprime_eqcons:
            # Wrap fprime_eqcons
            geval, fprime_eqcons = wrap_function(fprime_eqcons,args)
        else:
            # Wrap approx_jacobian
            geval, fprime_eqcons = wrap_function(approx_jacobian,
                                                 (f_eqcons,epsilon))
    else:
        # Equality constraints provided via eqcons[]
        eqcons_prime = []
        for i in range(len(eqcons)):
            eqcons_prime.append(None)
            if eqcons[i]:
                # Wrap eqcons and eqcons_prime
                ceval, eqcons[i] = wrap_function(eqcons[i],args)
                geval, eqcons_prime[i] = wrap_function(approx_fprime,
                                                       (eqcons[i],epsilon))

    if f_ieqcons:
        # Inequality constraints provided via f_ieqcons
        ceval, f_ieqcons = wrap_function(f_ieqcons,args)
        if fprime_ieqcons:
            # Wrap fprime_ieqcons
            geval, fprime_ieqcons = wrap_function(fprime_ieqcons,args)
        else:
            # Wrap approx_jacobian
            geval, fprime_ieqcons = wrap_function(approx_jacobian,
                                                  (f_ieqcons,epsilon))
    else:
        # Inequality constraints provided via ieqcons[]
        ieqcons_prime = []
        for i in range(len(ieqcons)):
            ieqcons_prime.append(None)
            if ieqcons[i]:
                # Wrap ieqcons and ieqcons_prime
                ceval, ieqcons[i] = wrap_function(ieqcons[i],args)
                geval, ieqcons_prime[i] = wrap_function(approx_fprime,
                                                        (ieqcons[i],epsilon))


    # Transform x0 into an array.
    x = asfarray(x0).flatten()

    # Set the parameters that SLSQP will need
    # meq = The number of equality constraints
    if f_eqcons:
        meq = len(f_eqcons(x))
    else:
        meq = len(eqcons)
    if f_ieqcons:
        mieq = len(f_ieqcons(x))
    else:
        mieq = len(ieqcons)
    # m = The total number of constraints
    m = meq + mieq
    # la = The number of constraints, or 1 if there are no constraints
    la = array([1,m]).max()
    # n = The number of independent variables
    n = len(x)

    # Define the workspaces for SLSQP
    n1 = n+1
    mineq = m - meq + n1 + n1
    len_w = (3*n1+m)*(n1+1)+(n1-meq+1)*(mineq+2) + 2*mineq+(n1+mineq)*(n1-meq) \
            + 2*meq + n1 +(n+1)*n/2 + 2*m + 3*n + 3*n1 + 1
    len_jw = mineq
    w = zeros(len_w)
    jw = zeros(len_jw)

    # Decompose bounds into xl and xu
    if len(bounds) == 0:
        bounds = [(-1.0E12, 1.0E12) for i in range(n)]
    elif len(bounds) != n:
        raise IndexError, \
        'SLSQP Error:  If bounds is specified, len(bounds) == len(x0)'
    else:
        for i in range(len(bounds)):
            if bounds[i][0] > bounds[i][1]:
                raise ValueError, \
                'SLSQP Error: lb > ub in bounds[' + str(i) +']  ' + str(bounds[4])

    xl = array( [ b[0] for b in bounds ] )
    xu = array( [ b[1] for b in bounds ] )



    # Initialize the iteration counter and the mode value
    mode = array(0,int)
    acc = array(acc,float)
    majiter = array(iter,int)
    majiter_prev = 0

    # Print the header if iprint >= 2
    if iprint >= 2:
        print "%5s %5s %16s %16s" % ("NIT","FC","OBJFUN","GNORM")

    while 1:

        if mode == 0 or mode == 1: # objective and constraint evaluation requird

            # Compute objective function
            fx = func(x)
            # Compute the constraints
            if f_eqcons:
                c_eq = f_eqcons(x)
            else:
                c_eq = array([ eqcons[i](x) for i in range(meq) ])
            if f_ieqcons:
                c_ieq = f_ieqcons(x)
            else:
                c_ieq = array([ ieqcons[i](x) for i in range(len(ieqcons)) ])

            # Now combine c_eq and c_ieq into a single matrix
            if m == 0:
                # no constraints
                c = zeros([la])
            else:
                # constraints exist
                if meq > 0 and mieq == 0:
                    # only equality constraints
                    c = c_eq
                if meq == 0 and mieq > 0:
                    # only inequality constraints
                    c = c_ieq
                if meq > 0 and mieq > 0:
                    # both equality and inequality constraints exist
                    c = append(c_eq, c_ieq)

        if mode == 0 or mode == -1: # gradient evaluation required

            # Compute the derivatives of the objective function
            # For some reason SLSQP wants g dimensioned to n+1
            g = append(fprime(x),0.0)

            # Compute the normals of the constraints
            if fprime_eqcons:
                a_eq = fprime_eqcons(x)
            else:
                a_eq = zeros([meq,n])
                for i in range(meq):
                    a_eq[i] = eqcons_prime[i](x)

            if fprime_ieqcons:
                a_ieq = fprime_ieqcons(x)
            else:
                a_ieq = zeros([mieq,n])
                for i in range(mieq):
                    a_ieq[i] = ieqcons_prime[i](x)

            # Now combine a_eq and a_ieq into a single a matrix
            if m == 0:
                # no constraints
                a = zeros([la,n])
            elif meq > 0 and mieq == 0:
                # only equality constraints
                a = a_eq
            elif meq == 0 and mieq > 0:
                # only inequality constraints
                a = a_ieq
            elif meq > 0 and mieq > 0:
                # both equality and inequality constraints exist
                a = vstack((a_eq,a_ieq))
            a = concatenate((a,zeros([la,1])),1)

        # Call SLSQP
        slsqp(m, meq, x, xl, xu, fx, c, g, a, acc, majiter, mode, w, jw)

        # Print the status of the current iterate if iprint > 2 and the
        # major iteration has incremented
        if iprint >= 2 and majiter > majiter_prev:
            print "%5i %5i % 16.6E % 16.6E" % (majiter,feval[0],
                                               fx,linalg.norm(g))

        # If exit mode is not -1 or 1, slsqp has completed
        if abs(mode) != 1:
            break

        majiter_prev = int(majiter)

    # Optimization loop complete.  Print status if requested
    if iprint >= 1:
        print exit_modes[int(mode)] + "    (Exit mode " + str(mode) + ')'
        print "            Current function value:", fx
        print "            Iterations:", majiter
        print "            Function evaluations:", feval[0]
        print "            Gradient evaluations:", geval[0]

    if not full_output:
        return x
    else:
        return [list(x),
                float(fx),
                int(majiter),
                int(mode),
                exit_modes[int(mode)] ]
