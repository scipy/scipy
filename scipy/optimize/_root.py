"""
Unified interfaces to root finding algorithms.

Functions
---------
- root : find a root of a vector function.
- root_scalar : find the root of a scalar function.
"""
from __future__ import division, print_function, absolute_import

__all__ = ['root', 'root_scalar']

import numpy as np

from scipy._lib.six import callable

from warnings import warn

from .optimize import MemoizeJac, OptimizeResult, _check_unknown_options
from .minpack import _root_hybr, leastsq
from ._spectral import _root_df_sane
from . import nonlin
from ._zeros import _brentq, _brenth, _ridder, _bisect


def root(fun, x0, args=(), method='hybr', jac=None, tol=None, callback=None,
         options=None):
    """
    Find a root of a vector function.

    Parameters
    ----------
    fun : callable
        A vector function to find a root of.
    x0 : ndarray
        Initial guess.
    args : tuple, optional
        Extra arguments passed to the objective function and its Jacobian.
    method : str, optional
        Type of solver.  Should be one of

            - 'hybr'             :ref:`(see here) <optimize.root-hybr>`
            - 'lm'               :ref:`(see here) <optimize.root-lm>`
            - 'broyden1'         :ref:`(see here) <optimize.root-broyden1>`
            - 'broyden2'         :ref:`(see here) <optimize.root-broyden2>`
            - 'anderson'         :ref:`(see here) <optimize.root-anderson>`
            - 'linearmixing'     :ref:`(see here) <optimize.root-linearmixing>`
            - 'diagbroyden'      :ref:`(see here) <optimize.root-diagbroyden>`
            - 'excitingmixing'   :ref:`(see here) <optimize.root-excitingmixing>`
            - 'krylov'           :ref:`(see here) <optimize.root-krylov>`
            - 'df-sane'          :ref:`(see here) <optimize.root-dfsane>`

    jac : bool or callable, optional
        If `jac` is a Boolean and is True, `fun` is assumed to return the
        value of Jacobian along with the objective function. If False, the
        Jacobian will be estimated numerically.
        `jac` can also be a callable returning the Jacobian of `fun`. In
        this case, it must accept the same arguments as `fun`.
    tol : float, optional
        Tolerance for termination. For detailed control, use solver-specific
        options.
    callback : function, optional
        Optional callback function. It is called on every iteration as
        ``callback(x, f)`` where `x` is the current solution and `f`
        the corresponding residual. For all methods but 'hybr' and 'lm'.
    options : dict, optional
        A dictionary of solver options. E.g. `xtol` or `maxiter`, see
        :obj:`show_options()` for details.

    Returns
    -------
    sol : OptimizeResult
        The solution represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the algorithm exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes.

    See also
    --------
    show_options : Additional options accepted by the solvers

    Notes
    -----
    This section describes the available solvers that can be selected by the
    'method' parameter. The default method is *hybr*.

    Method *hybr* uses a modification of the Powell hybrid method as
    implemented in MINPACK [1]_.

    Method *lm* solves the system of nonlinear equations in a least squares
    sense using a modification of the Levenberg-Marquardt algorithm as
    implemented in MINPACK [1]_.

    Method *df-sane* is a derivative-free spectral method. [3]_

    Methods *broyden1*, *broyden2*, *anderson*, *linearmixing*,
    *diagbroyden*, *excitingmixing*, *krylov* are inexact Newton methods,
    with backtracking or full line searches [2]_. Each method corresponds
    to a particular Jacobian approximations. See `nonlin` for details.

    - Method *broyden1* uses Broyden's first Jacobian approximation, it is
      known as Broyden's good method.
    - Method *broyden2* uses Broyden's second Jacobian approximation, it
      is known as Broyden's bad method.
    - Method *anderson* uses (extended) Anderson mixing.
    - Method *Krylov* uses Krylov approximation for inverse Jacobian. It
      is suitable for large-scale problem.
    - Method *diagbroyden* uses diagonal Broyden Jacobian approximation.
    - Method *linearmixing* uses a scalar Jacobian approximation.
    - Method *excitingmixing* uses a tuned diagonal Jacobian
      approximation.

    .. warning::

        The algorithms implemented for methods *diagbroyden*,
        *linearmixing* and *excitingmixing* may be useful for specific
        problems, but whether they will work may depend strongly on the
        problem.

    .. versionadded:: 0.11.0

    References
    ----------
    .. [1] More, Jorge J., Burton S. Garbow, and Kenneth E. Hillstrom.
       1980. User Guide for MINPACK-1.
    .. [2] C. T. Kelley. 1995. Iterative Methods for Linear and Nonlinear
        Equations. Society for Industrial and Applied Mathematics.
        <http://www.siam.org/books/kelley/>
    .. [3] W. La Cruz, J.M. Martinez, M. Raydan. Math. Comp. 75, 1429 (2006).

    Examples
    --------
    The following functions define a system of nonlinear equations and its
    jacobian.

    >>> def fun(x):
    ...     return [x[0]  + 0.5 * (x[0] - x[1])**3 - 1.0,
    ...             0.5 * (x[1] - x[0])**3 + x[1]]

    >>> def jac(x):
    ...     return np.array([[1 + 1.5 * (x[0] - x[1])**2,
    ...                       -1.5 * (x[0] - x[1])**2],
    ...                      [-1.5 * (x[1] - x[0])**2,
    ...                       1 + 1.5 * (x[1] - x[0])**2]])

    A solution can be obtained as follows.

    >>> from scipy import optimize
    >>> sol = optimize.root(fun, [0, 0], jac=jac, method='hybr')
    >>> sol.x
    array([ 0.8411639,  0.1588361])
    """
    if not isinstance(args, tuple):
        args = (args,)

    meth = method.lower()
    if options is None:
        options = {}

    if callback is not None and meth in ('hybr', 'lm'):
        warn('Method %s does not accept callback.' % method,
             RuntimeWarning)

    # fun also returns the jacobian
    if not callable(jac) and meth in ('hybr', 'lm'):
        if bool(jac):
            fun = MemoizeJac(fun)
            jac = fun.derivative
        else:
            jac = None

    # set default tolerances
    if tol is not None:
        options = dict(options)
        if meth in ('hybr', 'lm'):
            options.setdefault('xtol', tol)
        elif meth in ('df-sane',):
            options.setdefault('ftol', tol)
        elif meth in ('broyden1', 'broyden2', 'anderson', 'linearmixing',
                      'diagbroyden', 'excitingmixing', 'krylov'):
            options.setdefault('xtol', tol)
            options.setdefault('xatol', np.inf)
            options.setdefault('ftol', np.inf)
            options.setdefault('fatol', np.inf)

    if meth == 'hybr':
        sol = _root_hybr(fun, x0, args=args, jac=jac, **options)
    elif meth == 'lm':
        sol = _root_leastsq(fun, x0, args=args, jac=jac, **options)
    elif meth == 'df-sane':
        _warn_jac_unused(jac, method)
        sol = _root_df_sane(fun, x0, args=args, callback=callback,
                            **options)
    elif meth in ('broyden1', 'broyden2', 'anderson', 'linearmixing',
                  'diagbroyden', 'excitingmixing', 'krylov'):
        _warn_jac_unused(jac, method)
        sol = _root_nonlin_solve(fun, x0, args=args, jac=jac,
                                 _method=meth, _callback=callback,
                                 **options)
    else:
        raise ValueError('Unknown solver %s' % method)

    return sol


def _warn_jac_unused(jac, method):
    if jac is not None:
        warn('Method %s does not use the jacobian (jac).' % (method,),
             RuntimeWarning)


def _root_leastsq(func, x0, args=(), jac=None,
                  col_deriv=0, xtol=1.49012e-08, ftol=1.49012e-08,
                  gtol=0.0, maxiter=0, eps=0.0, factor=100, diag=None,
                  **unknown_options):
    """
    Solve for least squares with Levenberg-Marquardt

    Options
    -------
    col_deriv : bool
        non-zero to specify that the Jacobian function computes derivatives
        down the columns (faster, because there is no transpose operation).
    ftol : float
        Relative error desired in the sum of squares.
    xtol : float
        Relative error desired in the approximate solution.
    gtol : float
        Orthogonality desired between the function vector and the columns
        of the Jacobian.
    maxiter : int
        The maximum number of calls to the function. If zero, then
        100*(N+1) is the maximum where N is the number of elements in x0.
    epsfcn : float
        A suitable step length for the forward-difference approximation of
        the Jacobian (for Dfun=None). If epsfcn is less than the machine
        precision, it is assumed that the relative errors in the functions
        are of the order of the machine precision.
    factor : float
        A parameter determining the initial step bound
        (``factor * || diag * x||``). Should be in interval ``(0.1, 100)``.
    diag : sequence
        N positive entries that serve as a scale factors for the variables.
    """

    _check_unknown_options(unknown_options)
    x, cov_x, info, msg, ier = leastsq(func, x0, args=args, Dfun=jac,
                                       full_output=True,
                                       col_deriv=col_deriv, xtol=xtol,
                                       ftol=ftol, gtol=gtol,
                                       maxfev=maxiter, epsfcn=eps,
                                       factor=factor, diag=diag)
    sol = OptimizeResult(x=x, message=msg, status=ier,
                         success=ier in (1, 2, 3, 4), cov_x=cov_x,
                         fun=info.pop('fvec'))
    sol.update(info)
    return sol


def _root_nonlin_solve(func, x0, args=(), jac=None,
                       _callback=None, _method=None,
                       nit=None, disp=False, maxiter=None,
                       ftol=None, fatol=None, xtol=None, xatol=None,
                       tol_norm=None, line_search='armijo', jac_options=None,
                       **unknown_options):
    _check_unknown_options(unknown_options)

    f_tol = fatol
    f_rtol = ftol
    x_tol = xatol
    x_rtol = xtol
    verbose = disp
    if jac_options is None:
        jac_options = dict()

    jacobian = {'broyden1': nonlin.BroydenFirst,
                'broyden2': nonlin.BroydenSecond,
                'anderson': nonlin.Anderson,
                'linearmixing': nonlin.LinearMixing,
                'diagbroyden': nonlin.DiagBroyden,
                'excitingmixing': nonlin.ExcitingMixing,
                'krylov': nonlin.KrylovJacobian
                }[_method]

    if args:
        if jac:
            def f(x):
                return func(x, *args)[0]
        else:
            def f(x):
                return func(x, *args)
    else:
        f = func

    x, info = nonlin.nonlin_solve(f, x0, jacobian=jacobian(**jac_options),
                                  iter=nit, verbose=verbose,
                                  maxiter=maxiter, f_tol=f_tol,
                                  f_rtol=f_rtol, x_tol=x_tol,
                                  x_rtol=x_rtol, tol_norm=tol_norm,
                                  line_search=line_search,
                                  callback=_callback, full_output=True,
                                  raise_exception=False)
    sol = OptimizeResult(x=x)
    sol.update(info)
    return sol

def _root_broyden1_doc():
    """
    Options
    -------
    nit : int, optional
        Number of iterations to make. If omitted (default), make as many
        as required to meet tolerances.
    disp : bool, optional
        Print status to stdout on every iteration.
    maxiter : int, optional
        Maximum number of iterations to make. If more are needed to
        meet convergence, `NoConvergence` is raised.
    ftol : float, optional
        Relative tolerance for the residual. If omitted, not used.
    fatol : float, optional
        Absolute tolerance (in max-norm) for the residual.
        If omitted, default is 6e-6.
    xtol : float, optional
        Relative minimum step size. If omitted, not used.
    xatol : float, optional
        Absolute minimum step size, as determined from the Jacobian
        approximation. If the step size is smaller than this, optimization
        is terminated as successful. If omitted, not used.
    tol_norm : function(vector) -> scalar, optional
        Norm to use in convergence check. Default is the maximum norm.
    line_search : {None, 'armijo' (default), 'wolfe'}, optional
        Which type of a line search to use to determine the step size in
        the direction given by the Jacobian approximation. Defaults to
        'armijo'.
    jac_options : dict, optional
        Options for the respective Jacobian approximation.
            alpha : float, optional
                Initial guess for the Jacobian is (-1/alpha).
            reduction_method : str or tuple, optional
                Method used in ensuring that the rank of the Broyden
                matrix stays low. Can either be a string giving the
                name of the method, or a tuple of the form ``(method,
                param1, param2, ...)`` that gives the name of the
                method and values for additional parameters.

                Methods available:
                    - ``restart``: drop all matrix columns. Has no
                        extra parameters.
                    - ``simple``: drop oldest matrix column. Has no
                        extra parameters.
                    - ``svd``: keep only the most significant SVD
                        components.
                      Extra parameters:
                          - ``to_retain``: number of SVD components to
                              retain when rank reduction is done.
                              Default is ``max_rank - 2``.
            max_rank : int, optional
                Maximum rank for the Broyden matrix.
                Default is infinity (ie., no rank reduction).
    """
    pass

def _root_broyden2_doc():
    """
    Options
    -------
    nit : int, optional
        Number of iterations to make. If omitted (default), make as many
        as required to meet tolerances.
    disp : bool, optional
        Print status to stdout on every iteration.
    maxiter : int, optional
        Maximum number of iterations to make. If more are needed to
        meet convergence, `NoConvergence` is raised.
    ftol : float, optional
        Relative tolerance for the residual. If omitted, not used.
    fatol : float, optional
        Absolute tolerance (in max-norm) for the residual.
        If omitted, default is 6e-6.
    xtol : float, optional
        Relative minimum step size. If omitted, not used.
    xatol : float, optional
        Absolute minimum step size, as determined from the Jacobian
        approximation. If the step size is smaller than this, optimization
        is terminated as successful. If omitted, not used.
    tol_norm : function(vector) -> scalar, optional
        Norm to use in convergence check. Default is the maximum norm.
    line_search : {None, 'armijo' (default), 'wolfe'}, optional
        Which type of a line search to use to determine the step size in
        the direction given by the Jacobian approximation. Defaults to
        'armijo'.
    jac_options : dict, optional
        Options for the respective Jacobian approximation.

        alpha : float, optional
            Initial guess for the Jacobian is (-1/alpha).
        reduction_method : str or tuple, optional
            Method used in ensuring that the rank of the Broyden
            matrix stays low. Can either be a string giving the
            name of the method, or a tuple of the form ``(method,
            param1, param2, ...)`` that gives the name of the
            method and values for additional parameters.

            Methods available:
                - ``restart``: drop all matrix columns. Has no
                    extra parameters.
                - ``simple``: drop oldest matrix column. Has no
                    extra parameters.
                - ``svd``: keep only the most significant SVD
                    components.
                  Extra parameters:
                      - ``to_retain``: number of SVD components to
                          retain when rank reduction is done.
                          Default is ``max_rank - 2``.
        max_rank : int, optional
            Maximum rank for the Broyden matrix.
            Default is infinity (ie., no rank reduction).
    """
    pass

def _root_anderson_doc():
    """
    Options
    -------
    nit : int, optional
        Number of iterations to make. If omitted (default), make as many
        as required to meet tolerances.
    disp : bool, optional
        Print status to stdout on every iteration.
    maxiter : int, optional
        Maximum number of iterations to make. If more are needed to
        meet convergence, `NoConvergence` is raised.
    ftol : float, optional
        Relative tolerance for the residual. If omitted, not used.
    fatol : float, optional
        Absolute tolerance (in max-norm) for the residual.
        If omitted, default is 6e-6.
    xtol : float, optional
        Relative minimum step size. If omitted, not used.
    xatol : float, optional
        Absolute minimum step size, as determined from the Jacobian
        approximation. If the step size is smaller than this, optimization
        is terminated as successful. If omitted, not used.
    tol_norm : function(vector) -> scalar, optional
        Norm to use in convergence check. Default is the maximum norm.
    line_search : {None, 'armijo' (default), 'wolfe'}, optional
        Which type of a line search to use to determine the step size in
        the direction given by the Jacobian approximation. Defaults to
        'armijo'.
    jac_options : dict, optional
        Options for the respective Jacobian approximation.

        alpha : float, optional
            Initial guess for the Jacobian is (-1/alpha).
        M : float, optional
            Number of previous vectors to retain. Defaults to 5.
        w0 : float, optional
            Regularization parameter for numerical stability.
            Compared to unity, good values of the order of 0.01.
    """
    pass

def _root_linearmixing_doc():
    """
    Options
    -------
    nit : int, optional
        Number of iterations to make. If omitted (default), make as many
        as required to meet tolerances.
    disp : bool, optional
        Print status to stdout on every iteration.
    maxiter : int, optional
        Maximum number of iterations to make. If more are needed to
        meet convergence, ``NoConvergence`` is raised.
    ftol : float, optional
        Relative tolerance for the residual. If omitted, not used.
    fatol : float, optional
        Absolute tolerance (in max-norm) for the residual.
        If omitted, default is 6e-6.
    xtol : float, optional
        Relative minimum step size. If omitted, not used.
    xatol : float, optional
        Absolute minimum step size, as determined from the Jacobian
        approximation. If the step size is smaller than this, optimization
        is terminated as successful. If omitted, not used.
    tol_norm : function(vector) -> scalar, optional
        Norm to use in convergence check. Default is the maximum norm.
    line_search : {None, 'armijo' (default), 'wolfe'}, optional
        Which type of a line search to use to determine the step size in
        the direction given by the Jacobian approximation. Defaults to
        'armijo'.
    jac_options : dict, optional
        Options for the respective Jacobian approximation.

        alpha : float, optional
            initial guess for the jacobian is (-1/alpha).
    """
    pass

def _root_diagbroyden_doc():
    """
    Options
    -------
    nit : int, optional
        Number of iterations to make. If omitted (default), make as many
        as required to meet tolerances.
    disp : bool, optional
        Print status to stdout on every iteration.
    maxiter : int, optional
        Maximum number of iterations to make. If more are needed to
        meet convergence, `NoConvergence` is raised.
    ftol : float, optional
        Relative tolerance for the residual. If omitted, not used.
    fatol : float, optional
        Absolute tolerance (in max-norm) for the residual.
        If omitted, default is 6e-6.
    xtol : float, optional
        Relative minimum step size. If omitted, not used.
    xatol : float, optional
        Absolute minimum step size, as determined from the Jacobian
        approximation. If the step size is smaller than this, optimization
        is terminated as successful. If omitted, not used.
    tol_norm : function(vector) -> scalar, optional
        Norm to use in convergence check. Default is the maximum norm.
    line_search : {None, 'armijo' (default), 'wolfe'}, optional
        Which type of a line search to use to determine the step size in
        the direction given by the Jacobian approximation. Defaults to
        'armijo'.
    jac_options : dict, optional
        Options for the respective Jacobian approximation.

        alpha : float, optional
            initial guess for the jacobian is (-1/alpha).
    """
    pass

def _root_excitingmixing_doc():
    """
    Options
    -------
    nit : int, optional
        Number of iterations to make. If omitted (default), make as many
        as required to meet tolerances.
    disp : bool, optional
        Print status to stdout on every iteration.
    maxiter : int, optional
        Maximum number of iterations to make. If more are needed to
        meet convergence, `NoConvergence` is raised.
    ftol : float, optional
        Relative tolerance for the residual. If omitted, not used.
    fatol : float, optional
        Absolute tolerance (in max-norm) for the residual.
        If omitted, default is 6e-6.
    xtol : float, optional
        Relative minimum step size. If omitted, not used.
    xatol : float, optional
        Absolute minimum step size, as determined from the Jacobian
        approximation. If the step size is smaller than this, optimization
        is terminated as successful. If omitted, not used.
    tol_norm : function(vector) -> scalar, optional
        Norm to use in convergence check. Default is the maximum norm.
    line_search : {None, 'armijo' (default), 'wolfe'}, optional
        Which type of a line search to use to determine the step size in
        the direction given by the Jacobian approximation. Defaults to
        'armijo'.
    jac_options : dict, optional
        Options for the respective Jacobian approximation.

        alpha : float, optional
            Initial Jacobian approximation is (-1/alpha).
        alphamax : float, optional
            The entries of the diagonal Jacobian are kept in the range
            ``[alpha, alphamax]``.
    """
    pass

def _root_krylov_doc():
    """
    Options
    -------
    nit : int, optional
        Number of iterations to make. If omitted (default), make as many
        as required to meet tolerances.
    disp : bool, optional
        Print status to stdout on every iteration.
    maxiter : int, optional
        Maximum number of iterations to make. If more are needed to
        meet convergence, `NoConvergence` is raised.
    ftol : float, optional
        Relative tolerance for the residual. If omitted, not used.
    fatol : float, optional
        Absolute tolerance (in max-norm) for the residual.
        If omitted, default is 6e-6.
    xtol : float, optional
        Relative minimum step size. If omitted, not used.
    xatol : float, optional
        Absolute minimum step size, as determined from the Jacobian
        approximation. If the step size is smaller than this, optimization
        is terminated as successful. If omitted, not used.
    tol_norm : function(vector) -> scalar, optional
        Norm to use in convergence check. Default is the maximum norm.
    line_search : {None, 'armijo' (default), 'wolfe'}, optional
        Which type of a line search to use to determine the step size in
        the direction given by the Jacobian approximation. Defaults to
        'armijo'.
    jac_options : dict, optional
        Options for the respective Jacobian approximation.

        rdiff : float, optional
            Relative step size to use in numerical differentiation.
        method : {'lgmres', 'gmres', 'bicgstab', 'cgs', 'minres'} or function
            Krylov method to use to approximate the Jacobian.
            Can be a string, or a function implementing the same
            interface as the iterative solvers in
            `scipy.sparse.linalg`.

            The default is `scipy.sparse.linalg.lgmres`.
        inner_M : LinearOperator or InverseJacobian
            Preconditioner for the inner Krylov iteration.
            Note that you can use also inverse Jacobians as (adaptive)
            preconditioners. For example,

            >>> jac = BroydenFirst()
            >>> kjac = KrylovJacobian(inner_M=jac.inverse).

            If the preconditioner has a method named 'update', it will
            be called as ``update(x, f)`` after each nonlinear step,
            with ``x`` giving the current point, and ``f`` the current
            function value.
        inner_tol, inner_maxiter, ...
            Parameters to pass on to the "inner" Krylov solver.
            See `scipy.sparse.linalg.gmres` for details.
        outer_k : int, optional
            Size of the subspace kept across LGMRES nonlinear
            iterations.

            See `scipy.sparse.linalg.lgmres` for details.
    """
    pass


_rtol = np.finfo(float).eps * 2

def root_scalar(f, a, b, args=(), method='brentq', xtol=1e-12, rtol=_rtol,
                maxiter=100):
    """Find a root of a continuous function on a sign-changing interval, that
    is, an interval :math:`[a, b]` so that :math:`f(a)` and
    :math:`f(b)` have opposite signs.

    Parameters
    ----------
    f : function
        Python function returning a number. The function must be
        continuous and :math:`f(a)` and math:`f(b)` must have opposite
        signs.
    a : number
        One end of the sign-changing interval :math:`[a, b]`.
    b : number
        The other end of the sign-changing interval :math:`[a, b]`.
    args : tuple, optional
        Extra arguments for the function :math:`f`. In the code
        :math:`f` is called by ``apply(f, (x)+args)``.
    method : str, optional
        Type of solver. Should be one of

            - ``'brentq'`` :ref:`(see here) <root-scalar-brentq>`
            - ``'brenth'`` :ref:`(see here) <root-scalar-brenth>`
            - ``'ridder'`` :ref:`(see here) <root-scalar-ridder>`
            - ``'bisect'`` :ref:`(see here) <root-scalar-bisect>`

    xtol : number, optional
        The routine converges when a root is known to lie within
        `xtol` of the value returned. Should be nonnegative.  The
        routine modifies the value to take into account the relative
        precision of doubles.
    rtol : number, optional
        The routine converges when a root is known to lie within
        `rtol` times the value returned of the value returned. Should
        be greater than or equal to the default value of
        ``np.finfo(float).eps * 2``.
    maxiter : number, optional
        If convergence is not achieved in `maxiter` iterations, an error is
        raised.  Must be >= 0.

    Returns
    -------
    sol : OptimizeResult
        The solution represented as a ``OptimizeResult`` object.
        Attributes are: ``x``, a zero of :math:`f` between :math:`a`
        and :math:`b`, ``success``, a Boolean flag indicating if the
        algorithm exited successfully, ``message``, which describes
        the cause of termination, ``nit``, the number of iterations
        performed, and ``function_calls``, the number of calls to
        :math:`f`.

    See Also
    --------
    newton : One-dimensional root finding with just an initial value.

    Notes
    -----
    This section describes the available solvers that can be selected
    by the ``method`` parameter. The default is ``brentq``.

    .. _root-scalar-bisect:

    The method ``'bisect'`` uses bisection. At each iteration
    bisection begins with a sign-changing interval :math:`[a, b]`. A
    new sign-changing interval is chosen by taking :math:`m = (a +
    b)/2` and picking the sign-changing interval of :math:`[a, m]` and
    :math:`[m, b]`. (Note that exactly one of these intervals will be
    sign-changing unless :math:`f(m) = 0`, in which case the algorithm
    returns :math:`m`.) Since the width of the search interval is cut
    in half each time, it takes bisection about :math:`\log_2((b -
    a)/\delta)` iterations to achieve a tolerance of :math:`2\delta`.

    .. _root-scalar-brentq:

    The method ``'brentq'`` uses Brent's method, also known as the Van
    Winjngaarden-Dekker-Brent method. Brent's method seeks to speed up
    bisection by adding inverse linear and quadratic interpolation
    steps. At each iteration the method begins with a sign-changing
    interval :math:`[b, c]`, where :math:`b` is the current estimate
    for the root. The algorithm computes a point :math:`b''` via
    either interpolation (inverse quadratic interpolation if
    :math:`b`, :math:`c`, and the previous value of :math:`b` are
    distinct; inverse linear interpolation otherwise) or bisection. If
    :math:`|b - b''| > \delta` (where as in bisection the desired
    tolerance is :math:`2\delta`), then :math:`b''` is taken as the
    next :math:`b`, otherwise the next :math:`b` is taken to be the
    current value plus a step of size :math:`\delta` in the direction
    of the root. The power of Brent's method comes from balancing the
    use interpolation and bisection for computing :math:`b''` so that
    interpolation is used as much as possible in an attempt to speed
    up convergence while ensuring that enough bisection steps are
    taken that the algorithm is never much slower than
    bisection. Brent's algorithm will usually require less iterations
    than bisection for sufficiently smooth functions, though it can
    require more for non-smooth functions or functions with repeated
    roots; in the worst case scenario the number of iterations can be
    about :math:`[\log_2((b - a)/\delta)]^2`. The presentation here is
    taken from [Atkinson1989]_; for the full technial details see
    Brent's original work [Brent1973]_.

    .. _root-scalar-brenth:

    The method ``'brenth'`` is a variation on the classic Brent routine
    that uses hyperbolic extrapolation instead of inverse quadratic
    extrapolation. Generally on a par with ``'brentq'``, but not as
    heavily tested. It is a safe version of the secant method. The
    version here is by Chuck Harris.

    .. _root-scalar-ridder:

    The method ``'ridder'`` uses Ridders' method. Ridders' method is
    faster than bisection, but not generally as fast as the Brent
    rountines. The classic description and source of the algorithm is
    in [Ridders1979]_. A description can also be found in any recent
    edition of Numerical Recipes. The routine used here diverges
    slightly from standard presentations in order to be a bit more
    careful of tolerance.

    References
    ----------
    .. [Atkinson1989]
       Atkinson, K. E.
       *An Introduction to Numerical Analysis*
       University of Iowa: Wiley, pp. 91-94, 1989.
       Section 2.8: Brent's Rootfinding Algorithm

    .. [Brent1973]
       Brent, R. P.,
       *Algorithms for Minimization Without Derivatives*.
       Englewood Cliffs, NJ: Prentice-Hall, 1973. Ch. 3-4.

    .. [PressEtal1992]
       Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T.
       *Numerical Recipes in FORTRAN: The Art of Scientific Computing*, 2nd ed.
       Cambridge, England: Cambridge University Press, pp. 352-355, 1992.
       Section 9.3:  "Van Wijngaarden-Dekker-Brent Method."

    .. [Ridders1979]
       Ridders, C. F. J. "A New Algorithm for Computing a
       Single Root of a Real Continuous Function."
       IEEE Trans. Circuits Systems 26, 979-980, 1979.

    Examples
    --------
    The following finds a root of a function on `[-1, 1]` using bisection.

    >>> def f(x):
    ...     return (x - 0.1)**3

    >>> from scipy import optimize
    >>> sol = optimize.root_scalar(f, -1, 1, method='bisect')
    >>> sol.x
    0.0999999999994543

    """
    CONVERGED = 'converged'
    SIGNERR = 'sign error'
    CONVERR = 'convergence error'
    flag_map = {0: CONVERGED, -1: SIGNERR, -2: CONVERR}
    full_output = True
    disp = False

    if not isinstance(args, tuple):
        args = (args,)
    if xtol <= 0:
        raise ValueError("xtol too small (%g <= 0)" % xtol)
    if rtol < _rtol:
        raise ValueError("rtol too small (%g < %g)" % (rtol, _rtol))
    if method == 'brentq':
        res = _brentq(f, a, b, xtol, rtol, maxiter, args, full_output, disp)
    elif method == 'brenth':
        res = _brenth(f, a, b, xtol, rtol, maxiter, args, full_output, disp)
    elif method == 'ridder':
        res = _ridder(f, a, b, xtol, rtol, maxiter, args, full_output, disp)
    elif method == 'bisect':
        res = _bisect(f, a, b, xtol, rtol, maxiter, args, full_output, disp)
    else:
        raise ValueError('Unknown solver %s' % method)
    x, function_calls, nit, flag = res
    if flag == 0:
        success = True
    else:
        success = False
    message = flag_map[flag]
    sol = OptimizeResult(x=x, nit=nit, success=success, message=message)
    sol['function_calls'] = function_calls
    return sol
