"""
Unified interfaces to minimization algorithms.

Functions
---------
- minimize : minimization of a function of several variables.
"""


__all__ = ['minimize', 'show_minimize_options']


from warnings import warn

# unconstrained minimization
from optimize import _minimize_neldermead, _minimize_powell, \
        _minimize_cg, _minimize_bfgs, _minimize_newtoncg, \
        MemoizeJac
from anneal import _minimize_anneal

# contrained minimization
from lbfgsb import _minimize_lbfgsb
from tnc import _minimize_tnc
from cobyla import _minimize_cobyla
from slsqp import _minimize_slsqp

def minimize(fun, x0, args=(), method='BFGS', jac=None, hess=None,
             hessp=None, bounds=None, constraints=(),
             options=dict(), full_output=False, callback=None,
             retall=False):
    """
    Minimization of scalar function of one or more variables.

    Parameters
    ----------
    fun : callable
        Objective function.
    x0 : ndarray
        Initial guess.
    args : tuple, optional
        Extra arguments passed to the objective function and its
        derivatives (Jacobian, Hessian).
    method : str, optional
        Type of solver.  Should be one of:
            {'Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG', 'Anneal',
             'L-BFGS-B', 'TNC', 'COBYLA', 'SLSQP'}.
    jac : bool or callable, optional
        Jacobian of objective function. Only for CG, BFGS, Newton-CG.
        If `jac` is a Boolean and is True, `fun` is assumed to return the
        value of Jacobian along with the objective function. If False, the
        Jacobian will be estimated numerically.
        `jac` can also be a callable returning the Jacobian of the
        objective. In this case, it must accept the same arguments as
        `fun`.
    hess, hessp : callable, optional
        Hessian of objective function or Hessian of objective function
        times an arbitrary vector p.  Only for Newton-CG.
        Only one of `hessp` or `hess` needs to be given.  If `hess` is
        provided, then `hessp` will be ignored.  If neither `hess` nor
        `hessp` is provided, then the hessian product will be approximated
        using finite differences on `jac`. `hessp` must compute the Hessian
        times an arbitrary vector.
    bounds : sequence, optional
        Bounds for variables (only for L-BFGS-B, TNC, COBYLA and SLSQP).
        ``(min, max)`` pairs for each element in ``x``, defining
        the bounds on that parameter. Use None for one of ``min`` or
        ``max`` when there is no bound in that direction.
    constraints : dict or sequence of dict, optional
        Constraints definition (only for COBYLA and SLSQP).
        Each constraint is defined in a dictionary with fields:
            type: str
                Constraint type: 'eq' for equality, 'ineq' for inequality.
            fun: callable
                The function defining the constraint.
            jac: callable, optional
                The Jacobian of `fun` (only for SLSQP).
            args: sequence, optional
                Extra arguments to be passed to the function and Jacobian.
        Equality constraint means that the constraint function result is to
        be zero whereas inequality means that it is to be non-negative.
        Note that COBYLA only supports inequality constraints.
    options : dict, optional
        A dictionary of solver options. All methods accept the following
        generic options:
            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.
        For method-specific options, see `show_minimize_options`.
    full_output : bool, optional
        If True, return optional outputs.  Default is False.
    callback : callable, optional
        Called after each iteration, as ``callback(xk)``, where ``xk`` is the
        current parameter vector.
    retall : bool, optional
        If True, return a list of the solution at each iteration.  This is only
        done if `full_output` is True.

    Returns
    -------
    xopt : ndarray
        The solution.
    info : dict
        A dictionary of optional outputs (depending on the chosen method)
        with the keys:
            solution : ndarray
                The solution (same as `xopt`).
            success : bool
                Boolean flag indicating if a solution was found.
            status : int
                An integer flag indicating the type of termination.  Its
                value depends on the underlying solver.  Refer to `message`
                for more information.
            message : str
                A string message giving information about the cause of the
                termination.
            fun, jac, hess : ndarray
                Values of objective function, Jacobian and Hessian (if
                available).
            nfev, njev, nhev: int
                Number of evaluations of the objective functions and of its
                jacobian and hessian.
            nit: int
                Number of iterations.
            direc: ndarray
                Current set of direction vectors for the Powell method.
            T : float
                Final temperature for simulated annealing.
            accept : int
                Number of tests accepted.
            allvecs : list
                Solution at each iteration (if ``retall == True``).

    Notes
    -----
    This section describes the available solvers that can be selected by the
    'method' parameter. The default method is *BFGS*.

    Unconstrained minimization
    ~~~~~~~~~~~~~~~~~~~~~~~~~~
    Method *Nelder-Mead* uses the Simplex algorithm [1]_, [2]_. This
    algorithm has been successful in many applications but other algorithms
    using the first and/or second derivatives information might be preferred
    for their better performances and robustness in general.

    Method *Powell* is a modification of Powell's method [3]_, [4]_ which
    is a conjugate direction method. It performs sequential one-dimensional
    minimizations along each vector of the directions set (`direc` field in
    `options` and `info`), which is updated at each iteration of the main
    minimization loop. The function need not be differentiable, and no
    derivatives are taken.

    Method *CG* uses a nonlinear conjugate gradient algorithm by Polak and
    Ribiere, a variant of the Fletcher-Reeves method described in [5]_ pp.
    120-122. Only the first derivatives are used.

    Method *BFGS* uses the quasi-Newton method of Broyden, Fletcher,
    Goldfarb, and Shanno (BFGS) [5]_ pp. 136. It uses the first derivatives
    only. BFGS has proven good performance even for non-smooth
    optimizations

    Method *Newton-CG* uses a Newton-CG algorithm [5]_ pp. 168 (also known
    as the truncated Newton method). It uses a CG method to the compute the
    search direction. See also *TNC* method for a box-constrained
    minimization with a similar algorithm.

    Method *Anneal* uses simulated annealing, which is a probabilistic
    metaheuristic algorithm for global optimization. It uses no derivative
    information from the function being optimized.

    Constrained minimization
    ~~~~~~~~~~~~~~~~~~~~~~~~
    Method *L-BFGS-B* uses the L-BFGS-B algorithm [6]_, [7]_ for bound
    constrained minimization.

    Method *TNC* uses a truncated Newton algorithm [5]_, [8]_ to minimize a
    function with variables subject to bounds. This algorithm is uses
    gradient information; it is also called Newton Conjugate-Gradient. It
    differs from the *Newton-CG* method described above as it wraps a C
    implementation and allows each variable to be given upper and lower
    bounds.

    Method *COBYLA* uses the Constrained Optimization BY Linear
    Approximation (COBYLA) method [9]_, [10]_, [11]_. The algorithm is
    based on linear approximations to the objective function and each
    constraint. The method wraps a FORTRAN implementation of the algorithm.

    Method *SLSQP* uses Sequential Least SQuares Programming to minimize a
    function of several variables with any combination of bounds, equality
    and inequality constraints. The method wraps the SLSQP Optimization
    subroutine originally implemented by Dieter Kraft [12]_.

    References
    ----------
    .. [1] Nelder, J A, and R Mead. 1965. A Simplex Method for Function
        Minimization. The Computer Journal 7: 308-13.
    .. [2] Wright M H. 1996. Direct search methods: Once scorned, now
        respectable, in Numerical Analysis 1995: Proceedings of the 1995
        Dundee Biennial Conference in Numerical Analysis (Eds. D F
        Griffiths and G A Watson). Addison Wesley Longman, Harlow, UK.
        191-208.
    .. [3] Powell, M J D. 1964. An efficient method for finding the minimum of
       a function of several variables without calculating derivatives. The
       Computer Journal 7: 155-162.
    .. [4] Press W, S A Teukolsky, W T Vetterling and B P Flannery.
       Numerical Recipes (any edition), Cambridge University Press.
    .. [5] Nocedal, J, and S J Wright. 2006. Numerical Optimization.
       Springer New York.
    .. [6] Byrd, R H and P Lu and J. Nocedal. 1995. A Limited Memory
       Algorithm for Bound Constrained Optimization. SIAM Journal on
       Scientific and Statistical Computing 16 (5): 1190-1208.
    .. [7] Zhu, C and R H Byrd and J Nocedal. 1997. L-BFGS-B: Algorithm
       778: L-BFGS-B, FORTRAN routines for large scale bound constrained
       optimization. ACM Transactions on Mathematical Software 23 (4):
       550-560.
    .. [8] Nash, S G. Newton-Type Minimization Via the Lanczos Method.
       1984. SIAM Journal of Numerical Analysis 21: 770-778.
    .. [9] Powell, M J D. A direct search optimization method that models
       the objective and constraint functions by linear interpolation.
       1994. Advances in Optimization and Numerical Analysis, eds. S. Gomez
       and J-P Hennart, Kluwer Academic (Dordrecht), 51-67.
    .. [10] Powell M J D. Direct search algorithms for optimization
       calculations. 1998. Acta Numerica 7: 287-336.
    .. [11] Powell M J D. A view of algorithms for optimization without
       derivatives. 2007.Cambridge University Technical Report DAMTP
       2007/NA03
    .. [12] Kraft, D. A software package for sequential quadratic
       programming. 1988. Tech. Rep. DFVLR-FB 88-28, DLR German Aerospace
       Center -- Institute for Flight Mechanics, Koln, Germany.

    Examples
    --------
    Let us consider the problem of minimizing the Rosenbrock function. This
    function (and its respective derivatives) is implemented in `rosen`
    (resp. `rosen_der`, `rosen_hess`) in the `scipy.optimize`.

    >>> from scipy.optimize import minimize, rosen, rosen_der

    A simple application of the *Nelder-Mead* method is:

    >>> x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
    >>> xopt = minimize(rosen, x0, method='Nelder-Mead')
    Optimization terminated successfully.
         Current function value: 0.000066
         Iterations: 141
         Function evaluations: 243
    >>> print xopt
    [ 1.  1.  1.  1.  1.]

    Now using the *BFGS* algorithm, using the first derivative and a few
    options:

    >>> xopt, info = minimize(rosen, x0, method='BFGS', jac=rosen_der,
    ...                       options={'gtol': 1e-6, 'disp': False},
    ...                       full_output=True)

    >>> print info['message']
    Optimization terminated successfully.
    >>> print info['solution']
    [ 1.  1.  1.  1.  1.]
    >>> print info['hess']
    [[ 0.00749589  0.01255155  0.02396251  0.04750988  0.09495377]
     [ 0.01255155  0.02510441  0.04794055  0.09502834  0.18996269]
     [ 0.02396251  0.04794055  0.09631614  0.19092151  0.38165151]
     [ 0.04750988  0.09502834  0.19092151  0.38341252  0.7664427 ]
     [ 0.09495377  0.18996269  0.38165151  0.7664427   1.53713523]]


    Next, consider a minimization problem with several constraints (namely
    Example 16.4 from [5]_). The objective function is:

    >>> fun = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2

    There are three constraints defined as:

    >>> cons = ({'type': 'ineq', 'fun': lambda x:  x[0] - 2 * x[1] + 2},
    ...         {'type': 'ineq', 'fun': lambda x: -x[0] - 2 * x[1] + 6},
    ...         {'type': 'ineq', 'fun': lambda x: -x[0] + 2 * x[1] + 2})

    And variables must be positive, hence the following bounds:

    >>> bnds = ((0, None), (0, None))

    The optimization problem is solved using the SLSQP method as:

    >>> xopt, info = minimize(fun, (2, 0), method='SLSQP', bounds=bnds,
    ...                       constraints=cons, full_output=True)

    It should converge to the theoretical solution (1.4 ,1.7).
    """
    meth = method.lower()
    # check if optional parameters are supported by the selected method
    # - jac
    if meth in ['nelder-mead', 'powell', 'anneal', 'cobyla'] and bool(jac):
        warn('Method %s does not use gradient information (jac).' % method,
             RuntimeWarning)
    # - hess
    if meth != 'newton-cg' and hess is not None:
        warn('Method %s does not use Hessian information (hess).' % method,
             RuntimeWarning)
    # - constraints or bounds
    if meth in ['nelder-mead', 'powell', 'cg', 'bfgs', 'newton-cg'] and \
        (bounds is not None or any(constraints)):
        warn('Method %s cannot handle constraints nor bounds.' % method,
             RuntimeWarning)
    if meth in ['l-bfgs-b', 'tnc'] and any(constraints):
        warn('Method %s cannot handle constraints.' % method,
             RuntimeWarning)
    if meth is 'cobyla' and bounds is not None:
        warn('Method %s cannot handle bounds.' % method,
             RuntimeWarning)
    # - callback
    if meth in ['anneal', 'l-bfgs-b', 'tnc', 'cobyla', 'slsqp'] and \
       callback is not None:
        warn('Method %s does not support callback.' % method,
             RuntimeWarning)
    # - retall
    if meth in ['anneal', 'l-bfgs-b', 'tnc', 'cobyla', 'slsqp'] and \
       retall:
        warn('Method %s does not support retall.' % method,
             RuntimeWarning)

    # fun also returns the jacobian
    if not callable(jac):
        if bool(jac):
            fun = MemoizeJac(fun)
            jac = fun.derivative
        else:
            jac = None

    if meth == 'nelder-mead':
        return _minimize_neldermead(fun, x0, args, options, full_output,
                                    retall, callback)
    elif meth == 'powell':
        return _minimize_powell(fun, x0, args, options, full_output,
                                retall, callback)
    elif meth == 'cg':
        return _minimize_cg(fun, x0, args, jac, options, full_output,
                            retall, callback)
    elif meth == 'bfgs':
        return _minimize_bfgs(fun, x0, args, jac, options, full_output,
                              retall, callback)
    elif meth == 'newton-cg':
        return _minimize_newtoncg(fun, x0, args, jac, hess, hessp, options,
                                  full_output, retall, callback)
    elif meth == 'anneal':
        return _minimize_anneal(fun, x0, args, options, full_output)
    elif meth == 'l-bfgs-b':
        return _minimize_lbfgsb(fun, x0, args, jac, bounds, options,
                                full_output)
    elif meth == 'tnc':
        return _minimize_tnc(fun, x0, args, jac, bounds, options,
                             full_output)
    elif meth == 'cobyla':
        return _minimize_cobyla(fun, x0, args, constraints, options,
                                full_output)
    elif meth == 'slsqp':
        return _minimize_slsqp(fun, x0, args, jac, bounds,
                               constraints, options, full_output)
    else:
        raise ValueError('Unknown solver %s' % method)


def show_minimize_options(method=None):
    """Show documentation for additional options of minimize's methods.

    These are method-specific options that can be supplied to `minimize` in the
    ``options`` dict.

    Parameters
    ----------
    method : str, optional
        If not given, shows all methods.  Otherwise, show only the options for
        the specified method.  Valid values are: 'BFGS', 'Newton-CG',
        'Nelder-Mead', 'Powell', 'CG', 'Anneal', 'L-BFGS-B', 'TNC',
        'COBYLA', 'SLSQP'.

    Notes
    -----
    * BFGS options:
        gtol : float
            Gradient norm must be less than `gtol` before successful
            termination.
        norm : float
            Order of norm (Inf is max, -Inf is min).
        eps : float or ndarray
            If `jac` is approximated, use this value for the step size.

    * Nelder-Mead options:
        xtol : float
            Relative error in solution `xopt` acceptable for convergence.
        ftol : float
            Relative error in ``fun(xopt)`` acceptable for convergence.
        maxfev : int
            Maximum number of function evaluations to make.

    * Newton-CG options:
        xtol : float
            Average relative error in solution `xopt` acceptable for
            convergence.
        eps : float or ndarray
            If `jac` is approximated, use this value for the step size.

    * CG options:
        gtol : float
            Gradient norm must be less than `gtol` before successful
            termination.
        norm : float
            Order of norm (Inf is max, -Inf is min).
        eps : float or ndarray
            If `jac` is approximated, use this value for the step size.

    * Powell options:
        xtol : float
            Relative error in solution `xopt` acceptable for convergence.
        ftol : float
            Relative error in ``fun(xopt)`` acceptable for convergence.
        maxfev : int
            Maximum number of function evaluations to make.
        direc : ndarray
            Initial set of direction vectors for the Powell method.

    * Anneal options:
        schedule : str
            Annealing schedule to use. One of: 'fast', 'cauchy' or
            'boltzmann'.
        T0 : float
            Initial Temperature (estimated as 1.2 times the largest
            cost-function deviation over random points in the range).
        Tf : float
            Final goal temperature.
        maxfev : int
            Maximum number of function evaluations to make.
        maxaccept : int
            Maximum changes to accept.
        boltzmann : float
            Boltzmann constant in acceptance test (increase for less
            stringent test at each temperature).
        learn_rate : float
            Scale constant for adjusting guesses.
        ftol : float
            Relative error in ``fun(x)`` acceptable for convergence.
        quench, m, n : float
            Parameters to alter fast_sa schedule.
        lower, upper : float or ndarray
            Lower and upper bounds on `x`.
        dwell : int
            The number of times to search the space at each temperature.

    * L-BFGS-B options:
        maxcor : int
            The maximum number of variable metric corrections used to
            define the limited memory matrix. (The limited memory BFGS
            method does not store the full hessian but uses this many terms
            in an approximation to it.)
        factr : float
            The iteration stops when ``(f^k -
            f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr * eps``, where ``eps``
            is the machine precision, which is automatically generated by
            the code. Typical values for `factr` are: 1e12 for low
            accuracy; 1e7 for moderate accuracy; 10.0 for extremely high
            accuracy.
        pgtol : float
            The iteration will stop when ``max{|proj g_i | i = 1, ..., n}
            <= pgtol`` where ``pg_i`` is the i-th component of the
            projected gradient.
        maxfev : int
            Maximum number of function evaluations.

    * TNC options:
        scale : list of floats
            Scaling factors to apply to each variable.  If None, the
            factors are up-low for interval bounded variables and
            1+|x] fo the others.  Defaults to None
        offset : float
            Value to substract from each variable.  If None, the
            offsets are (up+low)/2 for interval bounded variables
            and x for the others.
        maxCGit : int
            Maximum number of hessian*vector evaluations per main
            iteration.  If maxCGit == 0, the direction chosen is
            -gradient if maxCGit < 0, maxCGit is set to
            max(1,min(50,n/2)).  Defaults to -1.
        maxfev : int
            Maximum number of function evaluation.  if None, `maxfev` is
            set to max(100, 10*len(x0)).  Defaults to None.
        eta : float
            Severity of the line search. if < 0 or > 1, set to 0.25.
            Defaults to -1.
        stepmx : float
            Maximum step for the line search.  May be increased during
            call.  If too small, it will be set to 10.0.  Defaults to 0.
        accuracy : float
            Relative precision for finite difference calculations.  If
            <= machine_precision, set to sqrt(machine_precision).
            Defaults to 0.
        minfev : float
            Minimum function value estimate.  Defaults to 0.
        ftol : float
            Precision goal for the value of f in the stoping criterion.
            If ftol < 0.0, ftol is set to 0.0 defaults to -1.
        xtol : float
            Precision goal for the value of x in the stopping
            criterion (after applying x scaling factors).  If xtol <
            0.0, xtol is set to sqrt(machine_precision).  Defaults to
            -1.
        pgtol : float
            Precision goal for the value of the projected gradient in
            the stopping criterion (after applying x scaling factors).
            If pgtol < 0.0, pgtol is set to 1e-2 * sqrt(accuracy).
            Setting it to 0.0 is not recommended.  Defaults to -1.
        rescale : float
            Scaling factor (in log10) used to trigger f value
            rescaling.  If 0, rescale at each iteration.  If a large
            value, never rescale.  If < 0, rescale is set to 1.3.

    * COBYLA options:
        rhobeg : float
            Reasonable initial changes to the variables.
        rhoend : float
            Final accuracy in the optimization (not precisely guaranteed).
            This is a lower bound on the size of the trust region.
        maxfev : int
            Maximum number of function evaluations.

    * SLSQP options:
        eps : float
            Step size used for numerical approximation of the jacobian.
        maxiter : int
            Maximum number of iterations.
    """
    if method is None:
        notes_header = "Notes\n    -----"
        sections = show_minimize_options.__doc__.split(notes_header)[1:]
    else:
        sections = show_minimize_options.__doc__.split('*')[1:]
        sections = [s.strip() for s in sections]
        sections = [s for s in sections if s.lower().startswith(method.lower())]

    print '\n'.join(sections)

    return
