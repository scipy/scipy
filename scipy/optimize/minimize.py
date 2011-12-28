"""
Interfaces to minimization algorithms.

Functions
---------
- minimize : unconstrained minimization of a function of several variables.
"""


__all__ = ['minimize', '_minimize_neldermead', '_minimize_powell',
           '_minimize_cg', '_minimize_bfgs', '_minimize_newtoncg',
           '_minimize_anneal']


from warnings import warn

# unconstrained minimization
from optimize import _minimize_neldermead, _minimize_powell, \
        _minimize_cg, _minimize_bfgs, _minimize_newtoncg
from anneal import _minimize_anneal


def minimize(fun, x0, args=(), method='Nelder-Mead', jac=None, hess=None,
             hessp=None, options=dict(), full_output=False, callback=None,
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
            {'Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG', 'Anneal'}.
    jac : callable, optional
        Jacobian of objective function (if None, Jacobian will be
        estimated numerically). Only for CG, BFGS, Newton-CG.
    hess, hessp : callable, optional
        Hessian of objective function or Hessian of objective function
        times an arbitrary vector p.  Only for Newton-CG.
        Only one of `hessp` or `hess` needs to be given.  If `hess` is
        provided, then `hessp` will be ignored.  If neither `hess` nor
        `hessp` is provided, then the hessian product will be approximated
        using finite differences on `jac`.  `hessp` must compute the hessian
        times an arbitrary vector.  If it is not given, finite-differences
        on `jac` are used to compute it.
    options : dict, optional
        A dictionary of solver options. All methods accept the following
        generic options:
            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.
        For method-specific options, refer to the documentation of
        respective underlying functions `_minimize_METHOD` (where METHOD is
        one of 'neldermead', 'powell', 'cg', 'bfgs', 'newtoncg' or
        'anneal').
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
    'method' parameter.  Respective solver options are also listed here.

    Method *Nelder-Mead* uses the Simplex algorithm [1]_, [2]_. This
    algorithm has been successful in many applications but other algorithms
    using the first and/or second derivatives information might be preferred
    for their better performances and robustness in general.

    Relevant `options` are: `xtol`, `ftol`, `maxiter`, `maxfev`, `disp`.


    Method *Powell* is a modification of Powell's method [3]_, [4]_ which
    is a conjugate direction method. It performs sequential one-dimensional
    minimizations along each vector of the directions set (`direc` field in
    `options` and `info`), which is updated at each iteration of the main
    minimization loop. The function need not be differentiable, and no
    derivatives are taken.

    Relevant `options` are: `xtol`, `ftol`, `maxiter`, `maxfev`, `disp`.


    Method *CG* uses a nonlinear conjugate gradient algorithm by Polak and
    Ribiere, a variant of the Fletcher-Reeves method described in [5]_ pp.
    120-122. Only the first derivatives are used.

    Relevant `options` are: `gtol`, `norm`, `maxiter`, `eps`, `disp`.


    Method *BFGS* uses the quasi-Newton method of Broyden, Fletcher,
    Goldfarb, and Shanno (BFGS) [5]_ pp. 136. It uses the first derivatives
    only. BFGS has proven good performance even for non-smooth
    optimizations

    Relevant `options` are: `gtol`, `norm`, `maxiter`, `eps`, `disp`.


    Method *Newton-CG* uses a Newton-CG algorithm [5]_ pp. 168 (also known
    as the truncated Newton method). It uses a CG method to the compute the
    search direction. See also `fmin_tnc` for a box-constrained
    minimization with a similar algorithm.

    Relevant `options` are: `xtol`, `maxiter`, `eps`, `disp`.


    Method *Anneal* uses simulated annealing, which is a probabilistic
    metaheuristic algorithm for global optimization. It uses no derivative
    information from the function being optimized.

    Relevant `options` are: `schedule`, `T0`, `Tf`, `maxfev`, `maxaccept`,
    `maxiter`, `boltzmann`, `learn_rate`, `ftol`, `quench`, `m`, `n`,
    `lower`, `upper`, `dwell`.

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
    """
    if method.lower() == 'nelder-mead':
        return _minimize_neldermead(fun, x0, args, options, full_output,
                                    retall, callback)
    elif method.lower() == 'powell':
        return _minimize_powell(fun, x0, args, options, full_output,
                                retall, callback)
    elif method.lower() == 'cg':
        return _minimize_cg(fun, x0, args, jac, options, full_output,
                            retall, callback)
    elif method.lower() == 'bfgs':
        return _minimize_bfgs(fun, x0, args, jac, options, full_output,
                              retall, callback)
    elif method.lower() == 'newton-cg':
        return _minimize_newtoncg(fun, x0, args, jac, hess, hessp, options,
                                  full_output, retall, callback)
    elif method.lower() == 'anneal':
        if callback:
            warn("Method 'Anneal' does not support callback.",
                 RuntimeWarning)
        if retall:
            warn("Method 'Anneal' does not support retall option.",
                 RuntimeWarning)
        return _minimize_anneal(fun, x0, args, options, full_output)
    else:
        raise ValueError('Unknown solver %s' % method)
