"""
This module implements the Sequential Least SQuares Programming optimization
algorithm (SLSQP), orginally developed by Dieter Kraft.
See http://www.netlib.org/toms/733

Functions
---------
.. autosummary::
   :toctree: generated/

    approx_jacobian
    fmin_slsqp

"""

from __future__ import division, print_function, absolute_import

__all__ = ['approx_jacobian','fmin_slsqp']

from scipy.optimize._slsqp import slsqp
from numpy import zeros, array, linalg, append, asfarray, concatenate, finfo, \
                  sqrt, vstack, exp, inf, where, isinf, atleast_1d
from .optimize import wrap_function, Result, _check_unknown_options

__docformat__ = "restructuredtext en"

_epsilon = sqrt(finfo(float).eps)

def approx_jacobian(x,func,epsilon,*args):
    """
    Approximate the Jacobian matrix of a callable function.

    Parameters
    ----------
    x : array_like
        The state vector at which to compute the Jacobian matrix.
    func : callable f(x,*args)
        The vector-valued function.
    epsilon : float
        The perturbation used to determine the partial derivatives.
    args : sequence
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
    f0 = atleast_1d(func(*((x0,)+args)))
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
                iprint = 1, disp = None, full_output = 0, epsilon = _epsilon ):
    """
    Minimize a function using Sequential Least SQuares Programming

    Python interface function for the SLSQP Optimization subroutine
    originally implemented by Dieter Kraft.

    Parameters
    ----------
    func : callable f(x,*args)
        Objective function.
    x0 : 1-D ndarray of float
        Initial guess for the independent variable(s).
    eqcons : list
        A list of functions of length n such that
        eqcons[j](x,*args) == 0.0 in a successfully optimized
        problem.
    f_eqcons : callable f(x,*args)
        Returns a 1-D array in which each element must equal 0.0 in a
        successfully optimized problem.  If f_eqcons is specified,
        eqcons is ignored.
    ieqcons : list
        A list of functions of length n such that
        ieqcons[j](x,*args) >= 0.0 in a successfully optimized
        problem.
    f_ieqcons : callable f(x,*args)
        Returns a 1-D ndarray in which each element must be greater or
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
    disp : int
        Over-rides the iprint interface (preferred).
    full_output : bool
        If False, return only the minimizer of func (default).
        Otherwise, output final objective function and summary
        information.
    epsilon : float
        The step size for finite-difference derivative estimates.

    Returns
    -------
    out : ndarray of float
        The final minimizer of func.
    fx : ndarray of float, if full_output is true
        The final value of the objective function.
    its : int, if full_output is true
        The number of iterations.
    imode : int, if full_output is true
        The exit mode from the optimizer (see below).
    smode : string, if full_output is true
        Message describing the exit mode from the optimizer.

    See also
    --------
    minimize: Interface to minimization algorithms for multivariate
        functions. See the 'SLSQP' `method` in particular.

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
    if disp is not None:
        iprint = disp
    opts = {'maxiter': iter,
            'ftol'   : acc,
            'iprint' : iprint,
            'disp'   : iprint != 0,
            'eps'    : epsilon}

    # Build the constraints as a tuple of dictionaries
    cons = ()
    # 1. constraints of the 1st kind (eqcons, ieqcons); no jacobian; take
    #    the same extra arguments as the objective function.
    cons += tuple({'type': 'eq', 'fun' : c, 'args': args} for c in eqcons)
    cons += tuple({'type': 'ineq', 'fun' : c, 'args': args} for c in ieqcons)
    # 2. constraints of the 2nd kind (f_eqcons, f_ieqcons) and their jacobian
    #    (fprime_eqcons, fprime_ieqcons); also take the same extra arguments
    #    as the objective function.
    if f_eqcons:
        cons += ({'type': 'eq', 'fun': f_eqcons, 'jac': fprime_eqcons,
                  'args': args}, )
    if f_ieqcons:
        cons += ({'type': 'ineq', 'fun': f_ieqcons, 'jac': fprime_ieqcons,
                  'args': args}, )

    res = _minimize_slsqp(func, x0, args, jac=fprime, bounds=bounds,
                          constraints=cons, **opts)
    if full_output:
        return res['x'], res['fun'], res['nit'], res['status'], res['message']
    else:
        return res['x']

def _minimize_slsqp(func, x0, args=(), jac=None, bounds=None,
                    constraints=(),
                    maxiter=100, ftol=1.0E-6, iprint=1, disp=False,
                    eps=_epsilon,
                    **unknown_options):
    """
    Minimize a scalar function of one or more variables using Sequential
    Least SQuares Programming (SLSQP).

    Options for the SLSQP algorithm are:
        ftol : float
            Precision goal for the value of f in the stopping criterion.
        eps : float
            Step size used for numerical approximation of the jacobian.
        disp : bool
            Set to True to print convergence messages. If False,
            `verbosity` is ignored and set to 0.
        maxiter : int
            Maximum number of iterations.

    This function is called by the `minimize` function with
    `method=SLSQP`. It is not supposed to be called directly.
    """
    _check_unknown_options(unknown_options)
    fprime = jac
    iter = maxiter
    acc = ftol
    epsilon = eps

    if not disp:
        iprint = 0

    # Constraints are triaged per type into a dictionnary of tuples
    if isinstance(constraints, dict):
        constraints = (constraints, )

    cons = {'eq': (), 'ineq': ()}
    for ic, con in enumerate(constraints):
        # check type
        try:
            ctype = con['type'].lower()
        except KeyError:
            raise KeyError('Constraint %d has no type defined.' % ic)
        except TypeError:
            raise TypeError('Constraints must be defined using a '
                            'dictionary.')
        except AttributeError:
            raise TypeError("Constraint's type must be a string.")
        else:
            if ctype not in ['eq', 'ineq']:
                raise ValueError("Unknown constraint type '%s'." % con['type'])

        # check function
        if 'fun' not in con:
            raise ValueError('Constraint %d has no function defined.' % ic)

        # check jacobian
        cjac = con.get('jac')
        if cjac is None:
            # approximate jacobian function
            def cjac(x, *args):
                return approx_jacobian(x, con['fun'], epsilon, *args)

        # update constraints' dictionary
        cons[ctype] += ({'fun' : con['fun'],
                         'jac' : cjac,
                         'args': con.get('args', ())}, )


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


    # Wrap func
    feval, func = wrap_function(func, args)

    # Wrap fprime, if provided, or approx_jacobian if not
    if fprime:
        geval, fprime = wrap_function(fprime, args)
    else:
        geval, fprime = wrap_function(approx_jacobian, (func, epsilon))

    # Transform x0 into an array.
    x = asfarray(x0).flatten()


    # Set the parameters that SLSQP will need
    # meq, mieq: number of equality and inequality constraints
    meq = sum(map(len, [atleast_1d(c['fun'](x, *c['args'])) for c in cons['eq']]))
    mieq = sum(map(len, [atleast_1d(c['fun'](x, *c['args'])) for c in cons['ineq']]))
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
    if bounds is None or len(bounds) == 0:
        xl, xu = array([-1.0E12]*n), array([1.0E12]*n)
    else:
        bnds = array(bounds, float)
        if bnds.shape[0] != n:
            raise IndexError('SLSQP Error: the length of bounds is not '
                             'compatible with that of x0.')

        bnderr = where(bnds[:, 0] > bnds[:, 1])[0]
        if bnderr.any():
            raise ValueError('SLSQP Error: lb > ub in bounds %s.' %
                             ', '.join(str(b) for b in bnderr))
        xl, xu = bnds[:, 0], bnds[:, 1]

        # filter -inf and inf values
        infbnd = isinf(bnds)
        xl[infbnd[:, 0]] = -1.0E12
        xu[infbnd[:, 1]] = 1.0E12

    # Initialize the iteration counter and the mode value
    mode = array(0,int)
    acc = array(acc,float)
    majiter = array(iter,int)
    majiter_prev = 0

    # Print the header if iprint >= 2
    if iprint >= 2:
        print("%5s %5s %16s %16s" % ("NIT","FC","OBJFUN","GNORM"))

    while 1:

        if mode == 0 or mode == 1: # objective and constraint evaluation requird

            # Compute objective function
            fx = func(x)
            # Compute the constraints
            if cons['eq']:
                c_eq  = concatenate([atleast_1d(con['fun'](x, *con['args']))
                                     for con in cons['eq']])
            else:
                c_eq = zeros(0)
            if cons['ineq']:
                c_ieq = concatenate([atleast_1d(con['fun'](x, *con['args']))
                                     for con in cons['ineq']])
            else:
                c_ieq = zeros(0)

            # Now combine c_eq and c_ieq into a single matrix
            c = concatenate((c_eq, c_ieq))

        if mode == 0 or mode == -1: # gradient evaluation required

            # Compute the derivatives of the objective function
            # For some reason SLSQP wants g dimensioned to n+1
            g = append(fprime(x),0.0)

            # Compute the normals of the constraints
            if cons['eq']:
                a_eq = vstack([con['jac'](x, *con['args'])
                               for con in cons['eq']])
            else: # no equality constraint
                a_eq = zeros((meq, n))

            if cons['ineq']:
                a_ieq = vstack([con['jac'](x, *con['args'])
                                for con in cons['ineq']])
            else: # no inequality constraint
                a_ieq = zeros((mieq, n))

            # Now combine a_eq and a_ieq into a single a matrix
            if m == 0: # no constraints
                a = zeros((la, n))
            else:
                a = vstack((a_eq, a_ieq))
            a = concatenate((a,zeros([la,1])),1)

        # Call SLSQP
        slsqp(m, meq, x, xl, xu, fx, c, g, a, acc, majiter, mode, w, jw)

        # Print the status of the current iterate if iprint > 2 and the
        # major iteration has incremented
        if iprint >= 2 and majiter > majiter_prev:
            print("%5i %5i % 16.6E % 16.6E" % (majiter,feval[0],
                                               fx,linalg.norm(g)))

        # If exit mode is not -1 or 1, slsqp has completed
        if abs(mode) != 1:
            break

        majiter_prev = int(majiter)

    # Optimization loop complete.  Print status if requested
    if iprint >= 1:
        print(exit_modes[int(mode)] + "    (Exit mode " + str(mode) + ')')
        print("            Current function value:", fx)
        print("            Iterations:", majiter)
        print("            Function evaluations:", feval[0])
        print("            Gradient evaluations:", geval[0])

    return Result(x=x, fun=fx, jac=g, nit=int(majiter), nfev=feval[0],
                  njev=geval[0], status=int(mode),
                  message=exit_modes[int(mode)], success=(mode == 0))

if __name__ == '__main__':

    # objective function
    def fun(x, r=[4, 2, 4, 2, 1]):
        """ Objective function """
        return exp(x[0]) * (r[0] * x[0]**2 + r[1] * x[1]**2 +
                            r[2] * x[0] * x[1] + r[3] * x[1] +
                            r[4])

    # bounds
    bnds = array([[-inf]*2, [inf]*2]).T
    bnds[:, 0] = [0.1, 0.2]

    # constraints
    def feqcon(x, b=1):
        """ Equality constraint """
        return array([x[0]**2 + x[1] - b])

    def jeqcon(x, b=1):
        """ Jacobian of equality constraint """
        return array([[2*x[0], 1]])

    def fieqcon(x, c=10):
        """ Inequality constraint """
        return array([x[0] * x[1] + c])

    def jieqcon(x, c=10):
        """ Jacobian of Inequality constraint """
        return array([[1, 1]])

    # constraints dictionaries
    cons=({'type': 'eq', 'fun' : feqcon, 'jac' : jeqcon, 'args': (1, )},
          {'type': 'ineq', 'fun' : fieqcon, 'jac' : jieqcon, 'args': (10,)})

    # Bounds constraint problem
    print(' Bounds constraints '.center(72, '-'))
    print(' * fmin_slsqp')
    x, f = fmin_slsqp(fun, array([-1, 1]), bounds=bnds, disp=1,
                      full_output=True)[:2]
    print(' * _minimize_slsqp')
    res = _minimize_slsqp(fun, array([-1, 1]), bounds=bnds,
                          **{'disp': True})

    # Equality and inequality constraints problem
    print(' Equality and inequality constraints '.center(72, '-'))
    print(' * fmin_slsqp')
    x, f = fmin_slsqp(fun, array([-1, 1]),
                      f_eqcons=feqcon, fprime_eqcons=jeqcon,
                      f_ieqcons=fieqcon, fprime_ieqcons=jieqcon,
                      disp=1, full_output=True)[:2]
    print(' * _minimize_slsqp')
    res = _minimize_slsqp(fun, array([-1, 1]), constraints=cons,
                          **{'disp': True})
