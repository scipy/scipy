"""
Interface to minimization algorithms

Functions
---------
- minimize : unconstrained minimization of a function of several variables.
- con_minimize : constrained minimization of a function of several variables.
- minimize1d: minimization of a scalar function.

"""


__all__ = ['minimize']


from warnings import warn

# unconstrained minimization
from optimize import _minimize_neldermead, _minimize_powell, \
        _minimize_cg, _minimize_bfgs, _minimize_ncg


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
            {'Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG'}.
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
        A dictionary of solver options with the keys:
            disp : bool
                Set to True to print convergence messages.
            xtol : float
                Relative error in solution `xopt` acceptable for convergence.
            ftol : float
                Relative error in ``fun(xopt)`` acceptable for convergence.
            maxit : int
                Maximum number of iterations to perform.
            maxfev : int
                Maximum number of function evaluations to make.
            gtol : float
                Gradient norm must be less than `gtol` before successful
                termination.
            norm : float
                Order of norm (Inf is max, -Inf is min).
            eps : int or ndarray
                If `jac` is approximated, use this value for the step size.
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
        A dictionary of optional outputs with the keys:
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
                Current direction set.
            allvecs : list
                Solution at each iteration (if ``retall == True``).

    Notes
    -----
    This section describes the available solvers that can be selected by the
    'method' parameter.  Respective solver options are also listed here.

    By default, this function uses the *Nelder-Mead* simplex algorithm
    [1]_, [2]_, [3]_ to find the minimum of a function of one or more
    variables.  This algorithm has a long history of successful use in
    applications.  However, it will usually be slower than an algorithm that
    uses first or second derivative information.  In practice it can have poor
    performance for high-dimensional problems and is not robust when
    minimizing complicated functions.  Additionally, there currently is no
    complete theory describing when the algorithm will successfully
    converge to the global minimum, or how fast it will if it does.

    Relevant `options` are: `xtol`, `ftol`, `maxit`, `maxfev`, `disp`.


    Method *Powell* is a modification of Powell's method [4]_, [5]_ to find
    the minimum of a function of N variables.  Powell's method is a
    conjugate direction method.

    The algorithm has two loops.  The outer loop merely iterates over the
    inner loop.  The inner loop minimizes over each current direction in the
    direction set.  At the end of the inner loop, if certain conditions are
    met, the direction that gave the largest decrease is dropped and
    replaced with the difference between the current estimate of the solution
    and the estimate from the beginning of the inner loop.

    The technical conditions for replacing the direction of greatest
    increase amount to checking that::

      1. No further gain can be made along the direction of greatest increase
         from that iteration.
      2. The direction of greatest increase accounted for a large sufficient
         fraction of the decrease in the function value from that iteration of
         the inner loop.

    Relevant `options` are: `xtol`, `ftol`, `maxit`, `maxfev`, `disp`.


    Method *CG* uses a nonlinear conjugate gradient algorithm by Polak and
    Ribiere as described in [6]_ pp. 120-122.

    Relevant `options` are: `gtol`, `norm`, `maxit`, `eps`, `disp`.


    Method *BFGS* uses the quasi-Newton method of Broyden, Fletcher,
    Goldfarb, and Shanno (BFGS) [6]_ pp. 198

    Relevant `options` are: `gtol`, `norm`, `maxit`, `eps`, `disp`.


    Method *Newton-CG* uses a Newton-CG algorithm [6]_ pp.140. Newton-CG
    methods are also called truncated Newton methods. See also `fmincon`
    with method='TNC' for a box-constrained minimization with a similar
    algorithm.

    Relevant `options` are: `xtol`, `maxit`, `eps`, `disp`.

    References
    ----------
    .. [1] Nelder, J.A. and Mead, R. (1965), "A simplex method for function
       minimization", The Computer Journal, 7, pp. 308-313
    .. [2] Wright, M.H. (1996), "Direct Search Methods: Once Scorned, Now
       Respectable", in Numerical Analysis 1995, Proceedings of the
       1995 Dundee Biennial Conference in Numerical Analysis, D.F.
    .. [3] Griffiths and G.A. Watson (Eds.), Addison Wesley Longman,
       Harlow, UK, pp. 191-208.
    .. [4] Powell M.J.D. (1964) An efficient method for finding the minimum of a
       function of several variables without calculating derivatives,
       Computer Journal, 7 (2):155-162.
    .. [5] Press W., Teukolsky S.A., Vetterling W.T., and Flannery B.P.:
       Numerical Recipes (any edition), Cambridge University Press
    .. [6] Wright & Nocedal, 'Numerical Optimization', 1999.

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
        return _minimize_ncg(fun, x0, args, jac, hess, hessp, options,
                             full_output, retall, callback)
    else:
        raise ValueError('Unknown solver %s' % method)
