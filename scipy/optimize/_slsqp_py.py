"""
This module implements the algorithm taken from
"Sequential Least Squares Programming" optimization algorithm (SLSQP),
originally developed by Dieter Kraft in Fortran77.

See http://www.netlib.org/toms/733

Functions
---------
.. autosummary::
   :toctree: generated/

    approx_jacobian
    fmin_slsqp

"""

__all__ = ['approx_jacobian', 'fmin_slsqp']

import numpy as np
from scipy.optimize._slsqp import slsqp
from numpy import (zeros, array, linalg, append, concatenate, finfo,
                   sqrt, vstack, isfinite, atleast_1d)
from ._optimize import (OptimizeResult, _check_unknown_options,
                        _prepare_scalar_function)
from ._numdiff import approx_derivative
from ._constraints import old_bound_to_new, _arr_to_scalar, Bounds
from scipy._lib._array_api import atleast_nd, array_namespace

# deprecated imports to be removed in SciPy 1.13.0
from numpy import exp, inf  # noqa: F401


__docformat__ = "restructuredtext en"

_epsilon = sqrt(finfo(float).eps)


def _extract_bounds(B: Bounds | None, n: int):
    """
    A helper function that converts the information encoded by a Bounds instance to
    separate lower and upper bounds with the finite entry indices to be used in SLSQP.

    Parameters
    ----------
    bounds : Bounds, None
        Bounds object instance to be decoded.
    n : int
        Length of the independent variable vector x for the SLSQP problem.

    Returns
    -------
    LB: ndarray
        If any, the finite lower bound entries. Returns None if none found.
    LB_ind: list
        If any, the indices of finite entries. Returns [] if none found.
    UB: ndarray
        If any, the finite upper bound entries. Returns None if none found.
    UB_ind: list
        If any, the indices of finite entries. Returns [] if none found.

    """
    # There are too many cases that can go wrong hence check the common convenience
    # cases and bail out otherwise.

    # No bounds given
    if not B:
        return None, [], None, []

    # if it is a scalar inside array, take it out so that it enters the scalar branch
    if isinstance(B.lb, np.ndarray) and B.lb.size == 1:
        B.lb = B.lb.item()
    if isinstance(B.ub, np.ndarray) and B.ub.size == 1:
        B.ub = B.ub.item()

    # Is it a scalar setting like B.lb = 4.5 to be broadcasted
    if not isinstance(B.lb, np.ndarray):
        # Assuming scalar given; check if it is finite
        if np.isfinite(B.lb):
            # OK finite entry hence all variables are bounded.
            LB = np.ones(n, dtype=np.float64)*B.lb
            LB_ind = [i for i in range(n)]
        else:
            # infinity given -> no bounds
            LB = None
            LB_ind = []
    elif isinstance(B.lb, np.ndarray):
        # OK array given check if it is n-long, scalar case was handled above
        if len(B.lb) != n:
            raise ValueError("The lower bound array size is not compatible with the "
                             "variable vector 'x'")
        else:
            LB_ind = np.isfinite(B.lb).nonzero()[0].tolist()
            LB = B.lb[LB_ind] if LB_ind else None
    else:
        raise TypeError("The lower bound array size should either be a NumPy array"
                        f"or a real scalar. However {type(B.lb)} is given.")

    # Do the same for the upper bound
    if not isinstance(B.ub, np.ndarray):
        if np.isfinite(B.ub):
            UB = np.ones(n, dtype=np.float64)*B.ub
            UB_ind = [i for i in range(n)]
        else:
            UB = None
            UB_ind = []
    elif isinstance(B.ub, np.ndarray):
        if len(B.ub) != n:
            raise ValueError("The upper bound array size is not compatible with the "
                             "variable vector 'x'")
        else:
            UB_ind = np.isfinite(B.ub).nonzero()[0].tolist()
            UB = B.ub[UB_ind] if UB_ind else None

    else:
        raise TypeError("The upper bound array size should either be a NumPy array"
                        f"or a real scalar. However {type(B.ub)} is given.")

    return LB, LB_ind, UB, UB_ind


def approx_jacobian(x, func, epsilon, *args):
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
    # approx_derivative returns (m, n) == (lenf, lenx)
    jac = approx_derivative(func, x, method='2-point', abs_step=epsilon,
                            args=args)
    # if func returns a scalar jac.shape will be (lenx,). Make sure
    # it's at least a 2D array.
    return np.atleast_2d(jac)


def fmin_slsqp(func, x0, eqcons=(), f_eqcons=None, ieqcons=(), f_ieqcons=None,
               bounds=(), fprime=None, fprime_eqcons=None,
               fprime_ieqcons=None, args=(), iter=100, acc=1.0E-6,
               iprint=1, disp=None, full_output=0, epsilon=_epsilon,
               callback=None):
    """
    Minimize a function using Sequential Least Squares Programming

    Python interface function for the SLSQP Optimization subroutine
    originally implemented by Dieter Kraft.

    Parameters
    ----------
    func : callable f(x,*args)
        Objective function.  Must return a scalar.
    x0 : 1-D ndarray of float
        Initial guess for the independent variable(s).
    eqcons : list, optional
        A list of functions of length n such that
        eqcons[j](x,*args) == 0.0 in a successfully optimized
        problem.
    f_eqcons : callable f(x,*args), optional
        Returns a 1-D array in which each element must equal 0.0 in a
        successfully optimized problem. If f_eqcons is specified,
        eqcons is ignored.
    ieqcons : list, optional
        A list of functions of length n such that
        ieqcons[j](x,*args) >= 0.0 in a successfully optimized
        problem.
    f_ieqcons : callable f(x,*args), optional
        Returns a 1-D ndarray in which each element must be greater or
        equal to 0.0 in a successfully optimized problem. If
        f_ieqcons is specified, ieqcons is ignored.
    bounds : list, optional
        A list of tuples specifying the lower and upper bound
        for each independent variable [(xl0, xu0),(xl1, xu1),...]
        Infinite values will be interpreted as unbounded and discarded.
    fprime : callable `f(x,*args)`, optional
        A function that evaluates the partial derivatives of func.
    fprime_eqcons : callable `f(x,*args)`, optional
        A function of the form `f(x, *args)` that returns the m by n
        array of equality constraint normals. If not provided,
        the normals will be approximated. The array returned by
        fprime_eqcons should be sized as ( len(eqcons), len(x0) ).
    fprime_ieqcons : callable `f(x,*args)`, optional
        A function of the form `f(x, *args)` that returns the m by n
        array of inequality constraint normals. If not provided,
        the normals will be approximated. The array returned by
        fprime_ieqcons should be sized as ( len(ieqcons), len(x0) ).
    args : sequence, optional
        Additional arguments passed to func and fprime.
    iter : int, optional
        The maximum number of iterations.
    acc : float, optional
        Requested accuracy.
    iprint : int, optional
        The verbosity of fmin_slsqp :

        * iprint <= 0 : Silent operation
        * iprint == 1 : Print summary upon completion (default)
        * iprint >= 2 : Print status of each iterate and summary
    disp : int, optional
        Overrides the iprint interface (preferred).
    full_output : bool, optional
        If False, return only the minimizer of func (default).
        Otherwise, output final objective function and summary
        information.
    epsilon : float, optional
        The step size for finite-difference derivative estimates.
    callback : callable, optional
        Called after each iteration, as ``callback(x)``, where ``x`` is the
        current parameter vector.

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
         0 : Optimization terminated successfully
         1 : Function evaluation required (f & c)
         2 : More equality constraints than independent variables
         3 : More than 3*n iterations in LSQ subproblem
         4 : Inequality constraints incompatible
         5 : Singular matrix E in LSQ subproblem
         6 : Singular matrix C in LSQ subproblem
         7 : Rank-deficient equality constraint subproblem HFTI
         8 : Positive directional derivative for linesearch
         9 : Iteration limit reached

    Examples
    --------
    Examples are given :ref:`in the tutorial <tutorial-sqlsp>`.

    """
    if disp is not None:
        iprint = disp

    opts = {'maxiter': iter,
            'ftol': acc,
            'iprint': iprint,
            'disp': iprint != 0,
            'eps': epsilon,
            'callback': callback}

    # Build the constraints as a tuple of dictionaries
    cons = ()
    # 1. constraints of the 1st kind (eqcons, ieqcons); no Jacobian; take
    #    the same extra arguments as the objective function.
    cons += tuple({'type': 'eq', 'fun': c, 'args': args} for c in eqcons)
    cons += tuple({'type': 'ineq', 'fun': c, 'args': args} for c in ieqcons)
    # 2. constraints of the 2nd kind (f_eqcons, f_ieqcons) and their Jacobian
    #    (fprime_eqcons, fprime_ieqcons); also take the same extra arguments
    #    as the objective function.
    if f_eqcons:
        cons += ({'type': 'eq', 'fun': f_eqcons, 'jac': fprime_eqcons,
                  'args': args}, )
    if f_ieqcons:
        cons += ({'type': 'ineq', 'fun': f_ieqcons, 'jac': fprime_ieqcons,
                  'args': args}, )

    #
    # FIXME: Convert bounds to Bounds object
    #

    res = _minimize_slsqp(func, x0, args, jac=fprime, bounds=bounds,
                          constraints=cons, **opts)
    if full_output:
        return res['x'], res['fun'], res['nit'], res['status'], res['message']
    else:
        return res['x']


def _minimize_slsqp(func, x0, args=(), jac=None, bounds=None,
                    constraints=(),
                    maxiter=100, ftol=1.0E-6, iprint=1, disp=False,
                    eps=_epsilon, callback=None, finite_diff_rel_step=None,
                    **unknown_options):
    """
    Minimize a scalar function of one or more variables using Sequential
    Least Squares Programming (SLSQP).

    Options
    -------
    ftol : float
        Precision goal for the value of f in the stopping criterion.
    eps : float
        Step size used for numerical approximation of the Jacobian.
    disp : bool
        Set to True to print convergence messages. If False,
        `verbosity` is ignored and set to 0.
    maxiter : int
        Maximum number of iterations.
    finite_diff_rel_step : None or array_like, optional
        If `jac in ['2-point', '3-point', 'cs']` the relative step size to
        use for numerical approximation of `jac`. The absolute step
        size is computed as ``h = rel_step * sign(x) * max(1, abs(x))``,
        possibly adjusted to fit into the bounds. For ``method='3-point'``
        the sign of `h` is ignored. If None (default) then step is selected
        automatically.
    """
    _check_unknown_options(unknown_options)

    iter = maxiter - 1
    acc = ftol
    epsilon = eps

    exit_modes = {2: "More equality constraints than independent variables",
                  3: "More than 3*n iterations in LSQ subproblem",
                  4: "Inequality constraints incompatible",
                  5: "Singular matrix E in LSQ subproblem",
                  6: "Singular matrix C in LSQ subproblem",
                  7: "Rank-deficient equality constraint subproblem HFTI",
                  8: "Positive directional derivative for linesearch",
                  9: "Iteration limit reached"}

    if not disp:
        iprint = 0

    # Transform x0 into an array.
    xp = array_namespace(x0)
    x0 = atleast_nd(x0, ndim=1, xp=xp)
    dtype = xp.float64
    if xp.isdtype(x0.dtype, "real floating"):
        dtype = x0.dtype
    x = xp.reshape(xp.astype(x0, dtype), -1)

    LB, LB_ind, UB, UB_ind = _extract_bounds(bounds, len(x))

    # mbounds = The number of finite upper or lower bounds
    B_ind = LB_ind + UB_ind
    mlbounds = len(LB_ind)
    mbounds = len(B_ind)

    # Constraints are triaged per type into a dictionary of tuples
    if isinstance(constraints, dict):
        constraints = (constraints, )

    cons = {'eq': (), 'ineq': ()}
    for ic, con in enumerate(constraints):
        # check type
        try:
            ctype = con['type'].lower()
        except KeyError as e:
            raise KeyError(f'Constraint {ic} has no type defined.') from e
        except TypeError as e:
            raise TypeError('Constraints must be defined using a '
                            'dictionary.') from e
        except AttributeError as e:
            raise TypeError("Constraint's type must be a string.") from e
        else:
            if ctype not in ['eq', 'ineq']:
                raise ValueError(f"Unknown constraint type '{con['type']}'."
                                 " Known types are 'eq' and 'ineq'.")

        # check function
        if 'fun' not in con:
            raise ValueError('Constraint %d has no function defined.' % ic)

        # check Jacobian
        cjac = con.get('jac')
        if cjac is None:
            # approximate Jacobian function. The factory function is needed
            # to keep a reference to `fun`, see gh-4240.
            def cjac_factory(fun):
                def cjac(x, *args):
                    x = _check_clip_x(x, bounds)

                    if jac in ['2-point', '3-point', 'cs']:
                        return approx_derivative(fun, x, method=jac, args=args,
                                                 rel_step=finite_diff_rel_step,
                                                 bounds=bounds)
                    else:
                        return approx_derivative(fun, x, method='2-point',
                                                 abs_step=epsilon, args=args,
                                                 bounds=bounds)

                return cjac
            cjac = cjac_factory(con['fun'])

        # update constraints' dictionary
        cons[ctype] += ({'fun': con['fun'],
                         'jac': cjac,
                         'args': con.get('args', ())}, )

    # meq, mieq: number of equality and inequality constraints
    meq = sum(map(len, [atleast_1d(c['fun'](x, *c['args'])) for c in cons['eq']]))
    mineq = sum(map(len, [atleast_1d(c['fun'](x, *c['args'])) for c in cons['ineq']]))

    # m = The total number of constraints (excluding bounds)
    m = meq + mineq
    # n = The number of independent variables
    n = len(x)

    # ScalarFunction provides function and gradient evaluation
    sf = _prepare_scalar_function(func, x, jac=jac, args=args, epsilon=eps,
                                  finite_diff_rel_step=finite_diff_rel_step,
                                  bounds=bounds)
    # gh11403 SLSQP sometimes exceeds bounds by 1 or 2 ULP, make sure this
    # doesn't get sent to the func/grad evaluator.
    wrapped_fun = _clip_x_for_func(sf.fun, bounds)
    wrapped_grad = _clip_x_for_func(sf.grad, bounds)

    # Print the header if iprint >= 2
    if iprint >= 2:
        print("%5s %5s %16s %16s" % ("NIT", "FC", "OBJFUN", "GNORM"))

    # Call objective, constraints and gradients, there should be no func
    # evaluations here because it's cached from ScalarFunction
    fx = wrapped_fun(x)
    g = append(wrapped_grad(x), 0.0)

    # ===========================================
    # Prepare the least squares problem arguments
    # ===========================================

    # If any, evaluate all meq-many eq. constraints into (E)quality array
    # of Ex = f
    if meq:
        E = np.zeros([meq, n])
        f = np.zeros(n)
        for ind, eqc in enumerate(cons['eq']):
            E[ind, :] = eqc['jac'](x, *eqc['args'])
            f[ind] = eqc['fun'](x, *eqc['args'])
    else:
        E = np.empty([0, 0])
        f = np.array([])

    # Same with inequality constraints to Gx >= h, including bounds on x

    # Are there any inequalities?
    if mineq + mbounds > 0:
        # +1 is for the additional variable to avoid inconsistent linearization
        G = np.zeros([mineq + mbounds + 1, n])
        h = np.zeros([mineq + mbounds + 1])

        if mineq:
            for ind, ineqc in enumerate(cons['ineq']):
                G[ind, :] = con['jac'](x, *con['args'])
                h[ind] = ineqc['fun'](x, *ineqc['args'])

        # We are only interested in finite entries in the bounds which will
        # enter the inequality constraints as new rows of Gx >= h as the
        # following
        #
        #
        #
        #     [...              ]        [       ]
        #     [inequality const.]        [       ]
        #     [...              ]        [       ]
        #     [0 1 0 ... 0 0 0 0]        [ xl[1] ] <-----------------
        #     [0 0 0 ... 0 0 1 0] * X >= [ xl[k] ]     |- mlbounds  |
        #     [      ...        ]        [  ..   ] <---             |- mbounds
        #     [0-1 0 ... 0 0 0 0]        [-xu[1] ]                  |
        #     [0 0 0 ...-1 0 0 0]        [-xu[j] ] <----------------
        #
        #
        #           [G]             X >=     h
        #
        #
        # Infinite entries are discarded. Thus, we need the independent
        #  variable indices for the LHS and the value for the RHS
        # G will be constant however the right hand side will be changing
        # across iterations.

        if mbounds:
            G[[i for i in range(mineq, mineq+mlbounds)], LB_ind] = 1.
            G[[i for i in range(mineq+mlbounds, mineq+mbounds)], UB_ind] = -1.
            h[mineq:mineq+mlbounds] = LB - x[LB_ind]
            h[mineq+mlbounds:-1] = UB - x[UB_ind]

    # ====================
    # Start main iteration
    # ====================

    for iter in maxiter:
        # Call SLSQP





        # objective and constraint evaluation required
        fx = wrapped_fun(x)
        # Compute constraints
        if cons['eq']:
            c_eq = concatenate([atleast_1d(con['fun'](x, *con['args']))
                                for con in cons['eq']])

        if cons['ineq']:
            c_ieq = concatenate([atleast_1d(con['fun'](x, *con['args']))
                                 for con in cons['ineq']])

        # gradient evaluation required
        g = append(wrapped_grad(x), 0.0)
        # Compute the normals of the constraints
        if cons['eq']:
            a_eq = vstack([con['jac'](x, *con['args'])
                          for con in cons['eq']])

        if cons['ineq']:
            a_ieq = vstack([con['jac'](x, *con['args'])
                           for con in cons['ineq']])



            # Print the status of the current iterate if iprint > 2
            if iprint >= 2:
                print("%5i %5i % 16.6E % 16.6E" % (majiter, sf.nfev,
                                                   fx, linalg.norm(g)))

        # End of iteration do callbacks
        if callback is not None:
            callback(np.copy(x))

    else:  # all iterations are exhausted
        # TODO: Maxiter reached
        mode = 9

    # Optimization loop complete. Print status if requested
    if iprint >= 1:
        print(exit_modes[int(mode)] + "    (Exit mode " + str(mode) + ')')
        print("            Current function value:", fx)
        print("            Iterations:", iter)
        print("            Function evaluations:", sf.nfev)
        print("            Gradient evaluations:", sf.ngev)

    return OptimizeResult(x=x, fun=fx, jac=g[:-1], nit=int(iter),
                          nfev=sf.nfev, njev=sf.ngev, status=int(mode),
                          message=exit_modes[int(mode)], success=(mode == 0))
