"""
Interface to Constrained Optimization By Linear Approximation

Functions
---------
.. autosummary::
   :toctree: generated/

    fmin_cobyla

"""

import functools
from threading import RLock
from inspect import signature

import numpy as np
from ._optimize import (OptimizeResult, _check_unknown_options,
    _prepare_scalar_function)
from ._constraints import (NonlinearConstraint, LinearConstraint,
    old_constraint_to_new)

from ._prima._prima import minimize as _prima_minimize

__all__ = ['fmin_cobyla']

# Workaround as _prima_minimize is not threadsafe since it
# uses static variables to hold the objective, constraint,
# and callback functions, among other things.
_module_lock = RLock()
def synchronized(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with _module_lock:
            return func(*args, **kwargs)
    return wrapper

@synchronized
def fmin_cobyla(func, x0, cons, args=(), consargs=None, rhobeg=1.0,
                rhoend=1e-4, maxfun=1000, disp=None, catol=2e-4,
                *, callback=None):
    """
    Minimize a function using the Constrained Optimization By Linear
    Approximation (COBYLA) method. This method wraps a FORTRAN
    implementation of the algorithm.

    Parameters
    ----------
    func : callable
        Function to minimize. In the form func(x, \\*args).
    x0 : ndarray
        Initial guess.
    cons : sequence
        Constraint functions; must all be ``>=0`` (a single function
        if only 1 constraint). Each function takes the parameters `x`
        as its first argument, and it can return either a single number or
        an array or list of numbers.
    args : tuple, optional
        Extra arguments to pass to function.
    consargs : tuple, optional
        Extra arguments to pass to constraint functions (default of None means
        use same extra arguments as those passed to func).
        Use ``()`` for no extra arguments.
    rhobeg : float, optional
        Reasonable initial changes to the variables.
    rhoend : float, optional
        Final accuracy in the optimization (not precisely guaranteed). This
        is a lower bound on the size of the trust region.
    disp : {0, 1, 2, 3}, optional
        Controls the frequency of output; 0 implies no output.
    maxfun : int, optional
        Maximum number of function evaluations.
    catol : float, optional
        Absolute tolerance for constraint violations.
    callback : callable, optional
        Called after each iteration, as ``callback(x)``, where ``x`` is the
        current parameter vector.

    Returns
    -------
    x : ndarray
        The argument that minimises `f`.

    See also
    --------
    minimize: Interface to minimization algorithms for multivariate
        functions. See the 'COBYLA' `method` in particular.

    Notes
    -----
    This algorithm is based on linear approximations to the objective
    function and each constraint. We briefly describe the algorithm.

    Suppose the function is being minimized over k variables. At the
    jth iteration the algorithm has k+1 points v_1, ..., v_(k+1),
    an approximate solution x_j, and a radius RHO_j.
    (i.e., linear plus a constant) approximations to the objective
    function and constraint functions such that their function values
    agree with the linear approximation on the k+1 points v_1,.., v_(k+1).
    This gives a linear program to solve (where the linear approximations
    of the constraint functions are constrained to be non-negative).

    However, the linear approximations are likely only good
    approximations near the current simplex, so the linear program is
    given the further requirement that the solution, which
    will become x_(j+1), must be within RHO_j from x_j. RHO_j only
    decreases, never increases. The initial RHO_j is rhobeg and the
    final RHO_j is rhoend. In this way COBYLA's iterations behave
    like a trust region algorithm.

    Additionally, the linear program may be inconsistent, or the
    approximation may give poor improvement. For details about
    how these issues are resolved, as well as how the points v_i are
    updated, refer to the source code or the references below.


    References
    ----------
    Powell M.J.D. (1994), "A direct search optimization method that models
    the objective and constraint functions by linear interpolation.", in
    Advances in Optimization and Numerical Analysis, eds. S. Gomez and
    J-P Hennart, Kluwer Academic (Dordrecht), pp. 51-67

    Powell M.J.D. (1998), "Direct search algorithms for optimization
    calculations", Acta Numerica 7, 287-336

    Powell M.J.D. (2007), "A view of algorithms for optimization without
    derivatives", Cambridge University Technical Report DAMTP 2007/NA03


    Examples
    --------
    Minimize the objective function f(x,y) = x*y subject
    to the constraints x**2 + y**2 < 1 and y > 0::

        >>> def objective(x):
        ...     return x[0]*x[1]
        ...
        >>> def constr1(x):
        ...     return 1 - (x[0]**2 + x[1]**2)
        ...
        >>> def constr2(x):
        ...     return x[1]
        ...
        >>> from scipy.optimize import fmin_cobyla
        >>> fmin_cobyla(objective, [0.0, 0.1], [constr1, constr2], rhoend=1e-7)
        array([-0.70710685,  0.70710671])

    The exact solution is (-sqrt(2)/2, sqrt(2)/2).



    """
    err = "cons must be a sequence of callable functions or a single"\
          " callable function."
    try:
        len(cons)
    except TypeError as e:
        if callable(cons):
            cons = [cons]
        else:
            raise TypeError(err) from e
    else:
        for thisfunc in cons:
            if not callable(thisfunc):
                raise TypeError(err)

    if consargs is None:
        consargs = args

    # build constraints
    con = tuple(
        old_constraint_to_new(i, {'type': 'ineq', 'fun': c, 'args': consargs})
        for i, c in enumerate(cons)
    )

    # options
    opts = {'rhobeg': rhobeg,
            'tol': rhoend,
            'disp': disp,
            'maxiter': maxfun,
            'catol': catol,
            'callback': callback}

    sol = _minimize_cobyla(func, x0, args, constraints=con,
                           **opts)
    if disp and not sol['success']:
        print(f"COBYLA failed to find a solution: {sol.message}")
    return sol['x']


@synchronized
def _minimize_cobyla(fun, x0, args=(), constraints=(),
                     rhobeg=1.0, tol=1e-6, maxiter=1000,
                     disp=0, catol=np.sqrt(np.finfo(float).eps),
                     ftarget=-np.inf, callback=None, bounds=None,
                     **unknown_options):
    """
    Minimize a scalar function of one or more variables using the
    Constrained Optimization BY Linear Approximation (COBYLA) algorithm.

    Options
    -------
    rhobeg : float
        Reasonable initial changes to the variables.
    tol : float
        Final accuracy in the optimization (not precisely guaranteed).
        This is a lower bound on the size of the trust region.
    disp : int
        Controls the frequency of output:

            0. (default) There will be no printing
            1. A message will be printed to the screen at the end of iteration, showing
               the best vector of variables found and its objective function value
            2. in addition to 1, each new value of RHO is printed to the screen,
               with the best vector of variables so far and its objective function
               value.
            3. in addition to 2, each function evaluation with its variables will
               be printed to the screen.

    maxiter : int
        Maximum number of function evaluations.
    catol : float
        Tolerance (absolute) for constraint violations.
    ftarget : float
        Stop if the objective function is less than `ftarget`.

    """
    # import here to avoid circular imports
    from .._lib.prima._linear_constraints import (
        combine_multiple_linear_constraints, separate_LC_into_eq_and_ineq)
    from .._lib.prima._nonlinear_constraints import process_nl_constraints

    _check_unknown_options(unknown_options)
    maxfun = maxiter
    rhoend = tol
    iprint = disp if disp is not None else 0
    if iprint != 0 and iprint != 1 and iprint != 2 and iprint != 3:
        raise ValueError(f'disp argument to minimize must be 0, 1, 2, or 3,\
                          received {iprint}')

    # create the ScalarFunction, cobyla doesn't require derivative function
    def _jac(x, *args):
        return None

    sf = _prepare_scalar_function(fun, x0, args=args, jac=_jac)

    if callback is not None:
        sig = signature(callback)
        if set(sig.parameters) == {"intermediate_result"}:
            def wrapped_callback_intermediate(x, f, nf, tr, cstrv, nlconstrlist):
                intermediate_result = OptimizeResult(x=np.copy(x), fun=f, nfev=nf,
                                                     nit=tr, maxcv=cstrv)
                callback(intermediate_result=intermediate_result)
        else:
            def wrapped_callback_intermediate(x, f, nf, tr, cstrv, nlconstrlist):
                callback(np.copy(x))
        def wrapped_callback(x, f, nf, tr, cstrv, nlconstrlist):
            try:
                wrapped_callback_intermediate(x, f, nf, tr, cstrv, nlconstrlist)
                return False
            except StopIteration:
                return True
    else:
        wrapped_callback = None

    linear_constraints = []
    nonlinear_constraints = []
    for constraint in constraints:
        if isinstance(constraint, LinearConstraint):
            linear_constraints.append(constraint)
        elif isinstance(constraint, NonlinearConstraint):
            nonlinear_constraints.append(constraint)
        else:
            raise ValueError(f"Uknown constraint type '{type(constraint)}'.")

    # Combine the linear constraints
    if len(linear_constraints) > 1:
        linear_constraint = \
            combine_multiple_linear_constraints(linear_constraints)
    elif len(linear_constraints) == 1:
        linear_constraint = linear_constraints[0]
    else:
        linear_constraint = None
    
    # Then process the linear constraints into separate matrices for equality
    # and inequality.
    if linear_constraint is not None:
        A_eq, b_eq, A_ineq, b_ineq = \
            separate_LC_into_eq_and_ineq(linear_constraint)
    else:
        A_eq, b_eq, A_ineq, b_ineq = None, None, None, None

    options = {
        'ctol': catol,
        'ftarget': ftarget,
        'iprint': iprint,
        'maxfev': maxfun,
        'rhobeg': rhobeg,
        'rhoend': rhoend,
    }

    if len(nonlinear_constraints) > 0:
        # PRIMA's process_nl_constraints expects the constraints to be either a
        # number or to have attribute __len__, but did not expect that it is
        # possible to be both. For this case we want the constraint to be a number
        # so that it can be broadcast to the length of the output of the constraint
        # function.
        for nlc in nonlinear_constraints:
            if hasattr(nlc.lb, '__len__') and isinstance(nlc.lb, np.ndarray) \
                    and nlc.lb.shape == ():
                nlc.lb = nlc.lb.tolist()
            if hasattr(nlc.ub, '__len__') and isinstance(nlc.ub, np.ndarray) \
                    and nlc.ub.shape == ():
                nlc.ub = nlc.ub.tolist()
        
        nonlinear_constraint_function = \
            process_nl_constraints(nonlinear_constraints)

        # COBYLA requires knowledge of the number of nonlinear constraints. In
        # order to get this number we need to evaluate the constraint function
        # at x0. The constraint value at x0 (nlconstr0) is not discarded but
        # passed down to the Fortran backend, as its evaluation is assumed to
        # be expensive. We also evaluate the objective function at x0 and pass
        # the result # (f0) down to the Fortran backend, which expects
        # nlconstr0 and f0 to be provided in sync.

        f0 = sf.fun(x0)
        nlconstr0 = nonlinear_constraint_function(x0)
        options['f0'] = f0
        options['nlconstr0'] = nlconstr0
        options['m_nlcon'] = len(nlconstr0)
    else:
        nonlinear_constraint_function = None
        options['m_nlcon'] = 0

    result = _prima_minimize(
        sf.fun,
        x0,
        args=(),  # _prepare_scalar_function handles args for us
        method='cobyla',
        lb=bounds.lb if bounds else None,
        ub=bounds.ub if bounds else None,
        A_eq=A_eq,
        b_eq=b_eq,
        A_ineq=A_ineq,
        b_ineq=b_ineq,
        nonlinear_constraint_function=nonlinear_constraint_function,
        callback=wrapped_callback,
        options=options
    )

    return OptimizeResult(x=result.x,
                          status=result.status,
                          success=result.success,
                          message=result.message,
                          nfev=result.nfev,
                          fun=result.fun,
                          maxcv=result.maxcv)
