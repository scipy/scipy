"""Generic interface for least-square minimization."""

from warnings import warn

import numpy as np
from numpy.linalg import norm

from . import _minpack
from .optimize import OptimizeResult
from ._numdiff import approx_derivative, group_columns
from ..sparse import issparse, csr_matrix, csc_matrix
from ..sparse.linalg import LinearOperator

from ._lsq_common import in_bounds, prepare_bounds
from ._lsq_trf import trf
from ._lsq_dogbox import dogbox

__all__ = ['least_squares']


EPS = np.finfo(float).eps


def check_tolerance(ftol, xtol, gtol):
    message = "{} is too low, setting to machine epsilon {}."
    if ftol < EPS:
        warn(message.format("`ftol`", EPS))
        ftol = EPS
    if xtol < EPS:
        warn(message.format("`xtol`", EPS))
        xtol = EPS
    if gtol < EPS:
        warn(message.format("`gtol`", EPS))
        gtol = EPS

    return ftol, xtol, gtol


TERMINATION_MESSAGES = {
    -1: "Improper input parameters status returned from `leastsq`",
    0: "The maximum number of function evaluations is exceeded.",
    1: "`gtol` termination condition is satisfied.",
    2: "`ftol` termination condition is satisfied.",
    3: "`xtol` termination condition is satisfied.",
    4: "Both `ftol` and `xtol` termination conditions are satisfied."
}


FROM_MINPACK_TO_COMMON = {
    0: -1,  # 0 improper input parameters for MINPACK.
    1: 2,
    2: 3,
    3: 4,
    4: 1,
    5: 0
    # There are 6, 7, 8 for too small tolerance parameters,
    # but we guard against them by checking ftol, xtol, gtol beforehand.
}


def call_minpack(fun, x0, jac, ftol, xtol, gtol, max_nfev, scaling, diff_step):
    n = x0.size

    if diff_step is None:
        epsfcn = np.finfo(float).eps
    else:
        epsfcn = diff_step**2

    if scaling == 'jac':
        scaling = None

    full_output = True
    col_deriv = False
    factor = 1.0

    if jac is None:
        if max_nfev is None:
            # n squared to account for Jacobian evaluations.
            max_nfev = 100 * n * (n + 1)
        x, info, status = _minpack._lmdif(
            fun, x0, (), full_output, ftol, xtol, gtol,
            max_nfev, epsfcn, factor, scaling)
    else:
        if max_nfev is None:
            max_nfev = 100 * n
        x, info, status = _minpack._lmder(
            fun, jac, x0, (), full_output, col_deriv,
            ftol, xtol, gtol, max_nfev, factor, scaling)

    f = info['fvec']

    if callable(jac):
        J = jac(x)
    else:
        J = np.atleast_2d(approx_derivative(fun, x))

    cost = 0.5 * np.dot(f, f)
    g = J.T.dot(f)
    g_norm = norm(g, ord=np.inf)

    nfev = info['nfev']
    njev = info.get('njev', None)

    status = FROM_MINPACK_TO_COMMON[status]
    active_mask = np.zeros_like(x0, dtype=int)

    return OptimizeResult(
        x=x, fun=f, jac=J, cost=cost, optimality=g_norm,
        active_mask=active_mask, nfev=nfev, njev=njev, status=status)


def check_scaling(scaling, x0):
    if scaling == 'jac':
        return scaling

    try:
        scaling = np.asarray(scaling, dtype=float)
        valid = np.all(np.isfinite(scaling)) and np.all(scaling > 0)
    except (ValueError, TypeError):
        valid = False

    if not valid:
        raise ValueError("`scaling` must be 'jac' or array_like with "
                         "positive numbers.")

    if scaling.ndim == 0:
        scaling = np.resize(scaling, x0.shape)

    if scaling.shape != x0.shape:
        raise ValueError("Inconsistent shapes between `scaling` and `x0`.")

    return scaling


def check_jac_sparsity(jac_sparsity, m, n):
    if jac_sparsity is None:
        return None

    if not issparse(jac_sparsity):
        jac_sparsity = np.atleast_2d(jac_sparsity)

    if jac_sparsity.shape != (m, n):
        raise ValueError("`jac_sparsity` has wrong shape.")

    return jac_sparsity, group_columns(jac_sparsity)


def least_squares(
        fun, x0, jac='2-point', bounds=(-np.inf, np.inf), method='trf',
        ftol=EPS**0.5, xtol=EPS**0.5, gtol=EPS**0.5, scaling=1.0,
        diff_step=None, tr_solver=None, tr_options={}, jac_sparsity=None,
        max_nfev=None, verbose=0, args=(), kwargs={}):
    """Minimize the sum of squares of nonlinear functions, subject to bound
    constraints on independent variables.

    Let f(x) maps from R^n to R^m, `least_squares` finds a local minimum of::

        F(x) = 0.5 * ||f(x)||**2 = 0.5 * sum(f_i(x)**2, i = 1, ..., m)
        lb <= x <= ub

    We call f(x) as a vector of residuals or simply residuals, and F(x) as a
    cost function or simply cost.

    Partial derivatives of f with respect to x form m-by-n matrix called
    Jacobian, where an element (i, j) equals the partial derivative of f[i]
    with respect to x[j].

    Parameters
    ----------
    fun : callable
        Function which computes the vector of residuals with the signature
        ``fun(x, *args, **kwargs)``, i.e., the minimization proceeds with
        respect to it's first argument. The argument ``x`` passed to this
        function is ndarray of shape (n,) (never a scalar, even for n=1).
        It must return a 1-d array_like of shape (m,) or a scalar.
    x0 : array_like with shape (n,) or float
        Initial guess on independent variables. If float, it will be treated
        as a 1-d array with one element.
    jac : {'2-point', '3-point', 'cs', callable}, optional
        Method of computing the Jacobian matrix. The keywords select the
        finite difference scheme for numerical estimation. The scheme '3-point'
        is more accurate, but requires twice as much operations compared to
        '2-point' (default). The scheme 'cs' uses complex steps, and while
        potentially the most accurate it is applicable only when `fun`
        correctly handles complex inputs and can be analytically continued to
        the complex plane.
        If callable, it is used as ``jac(x, *args, **kwargs)`` and should
        return a good approximation (or the exact value) for the Jacobian as
        an array_like (np.atleast_2d is applied), a sparse matrix or a
        `scipy.sparse.linalg.LinearOperator`, all with shape  (m, n).
    bounds : 2-tuple of array_like, optional
        Lower and upper bounds on independent variables. Defaults to no bounds.
        Each array must match the size of `x0` or be a scalar, in the latter
        case a bound will be the same for all variables. Use ``np.inf`` with
        an appropriate sign to disable bounds on all or some variables.
    method : {'trf', 'dogbox', 'lm'}, optional
        Algorithm to perform minimization.
            * 'trf' - Trust Region Reflective algorithm, particularly suitable
              for large sparse problems with bounds. Generally robust method.
            * 'dogbox' - dogleg algorithm with rectangular trust regions,
              typical use case is small problems with bounds. Not recommended
              to use in problems with rank-deficient Jacobian.
            * 'lm' - Levenberg-Marquardt algorithm as implemented in MINPACK.
              Doesn't handle bounds and sparse Jacobians. It is usually the
              most efficient method for small unconstrained problems.
        Default is 'trf'. See Notes for more information.
    ftol : float, optional
        Tolerance for termination by the change of the cost function.
        Default is the square root of machine epsilon. The optimization process
        is stopped when ``dF < ftol * F``, and there was an adequate agreement
        between a local quadratic model and the true model in the last step.
    xtol : float, optional
        Tolerance for termination by the change of the independent variables.
        Default is the square root of machine epsilon. The exact condition
        checked depends on the `method` used:
            * For 'trf' and 'dogbox': ``norm(dx) < xtol * (xtol + norm(x))``
            * For 'lm': ``Delta < xtol * norm(scaled_x)``, where ``Delta`` is
              a trust-region radius and ``scaled_x`` is the value of ``x``
              scaled according to `scaling` parameter (see below).
    gtol : float, optional
        Tolerance for termination by the norm of the gradient. Default is
        the square root of machine epsilon. The exact condition depends
        on a `method` used:
            * For 'trf' : ``norm(g_scaled, ord=np.inf) < gtol``, where
              ``g_scaled`` is the value of the gradient scaled to account for
              the presence of the bounds [STIR]_.
            * For 'dogbox' : ``norm(g_free, ord=np.inf) < gtol``, where
              ``g_free`` is the gradient with respect to the variables which
              are not in the optimal state on the boundary.
            * For 'lm' : the maximum absolute value of the cosine of angles
              between columns of the Jacobian and the residual vector is less
              than `gtol`, or the residual vector is zero.
    max_nfev : None or int, optional
        Maximum number of function evaluations before the termination.
        If None (default), the value is chosen automatically:
            * For 'trf' and 'dogbox' : 100 * n.
            * For 'lm':  100 * n if `jac` is callable and 100 * n * (n + 1)
              otherwise (because 'lm' counts function calls in Jacobian
              estimation).
    scaling : array_like or 'jac', optional
        Applies scaling to the variables to potentially improve the algorithm's
        convergence. Default is 1.0, which means no scaling. Scaling should be
        used to equalize the influence of each variable on the cost function.
        Alternatively you can think of `scaling` as diagonal elements of
        a matrix which determines the shape of a trust region. Use smaller
        values for variables which have larger characteristic scale compared
        to others. A scalar value won't affect the algorithm (except maybe
        fixing/introducing numerical issues and changing termination criteria).
        If 'jac', then scaling is proportional to the norms of columns of the
        Jacobian matrix. If the algorithm converges poorly on your problem
        try using this parameter.
    diff_step : None or array_like, optional
        Determines the relative step size for the finite difference
        approximation of the Jacobian. The actual step is computed as
        ``x * diff_step``. If None (default), then `diff_step` is taken to be
        a conventional "optimal" power of machine epsilon for the finite
        difference scheme used [NR]_.
    tr_solver : {None, 'exact', 'lsmr'}, optional
        Method for solving trust-region subproblems, relevant only for 'trf'
        and 'dogbox' methods.
            * 'exact' is suitable for not very large problems with dense
              Jacobian matrices. The computational complexity per iteration is
              comparable to a singular value decomposition of the Jacobian
              matrix.
            * 'lsmr' is suitable for problems with sparse and large Jacobian
              matrices. It uses the iterative procedure
              `scipy.sparse.linalg.lsmr` for finding a solution of a linear
              least-squares problem and only requires matrix-vector product
              evaluations.
        If None (default) the solver is chosen based on type of Jacobian
        returned on the first iteration.
    tr_options : dict, optional
        Keyword options passed to trust-region solver.
            * ``tr_solver='exact'``: `tr_options` are ignored.
            * ``tr_solver='lsmr'``: options for `scipy.sparse.linalg.lsmr`.
              Additionally  ``method='trf'`` supports  'regularize' option
              (bool, default is True) which adds a regularization term to the
              normal equations, which improves convergence if Jacobian is
              rank-deficient [Byrd]_ (eq. 3.4).
    jac_sparsity : {None, array_like, sparse matrix}, optional
        Defines the sparsity structure of the Jacobian matrix for finite
        differences. If the Jacobian has only few non-zeros in *each* row,
        providing the sparsity structure will greatly speed up the computations
        [Curtis]_. Should have shape (m, n). A zero entry means that a
        corresponding element in the Jacobian is identically zero. If provided,
        forces the use of 'lsmr' trust-region solver. If None (default) then
        dense differencing will be used. Has no effect for 'lm' method.
    verbose : {0, 1, 2}, optional
        Level of algorithm's verbosity:
            * 0 (default) - work silently.
            * 1 - display a termination report.
            * 2 - display progress during iterations (not supported by 'lm'
              method).
    args, kwargs : tuple and dict, optional
        Additional arguments passed to `fun` and `jac`. Both empty by default.
        The calling signature is ``fun(x, *args, **kwargs)`` and the same for
        `jac`.

    Returns
    -------
    `OptimizeResult` with the following fields defined:
    x : ndarray, shape (n,)
        Solution found.
    cost : float
        Value of the cost function at the solution.
    fun : ndarray, shape (m,)
        Vector of residuals at the solution.
    jac : ndarray, sparse matrix or LinearOperator, shape (m, n)
        Jacobian matrix at the solution. The type is the same as was used by
        the algorithm.
    optimality : float
        First-order optimality measure. In unconstrained problems it is always
        the uniform norm of the gradient. In constrained problems it is the
        quantity which was compared with `gtol` during iterations.
    active_mask : ndarray of int, shape (n,)
        Each component shows whether a corresponding constraint is active
        (that is, whether a variable is at the bound):
            *  0 - a constraint is not active.
            * -1 - a lower bound is active.
            *  1 - an upper bound is active.
        Might be somewhat arbitrary for 'trf' method as it does strictly
        feasible iterates and `active_mask` is determined within a tolerance
        threshold.
    nfev : int
        Number of function evaluations done. Methods 'trf' and 'dogbox' do not
        count function calls for numerical Jacobian approximation, as opposed
        to 'lm' method.
    njev : int or None
        Number of Jacobian evaluations done. If numerical Jacobian
        approximation is used in 'lm' method, it is set to None.
    status : int
        The reason for algorithm termination:
            * -1 - improper input parameters status returned from MINPACK.
            *  0 - the maximum number of function evaluations is exceeded.
            *  1 - `gtol` termination condition is satisfied.
            *  2 - `ftol` termination condition is satisfied.
            *  3 - `xtol` convergence test is satisfied.
            *  4 - Both `ftol` and `xtol` termination conditions are satisfied.
    message : string
        Verbal description of the termination reason.
    success : int
        True if one of the convergence criteria is satisfied.

    See Also
    --------
    leastsq : A legacy wrapper for the MINPACK implementation of the
              Levenberg-Marquadt algorithm.
    curve_fit : Least-squares minimization applied to a curve fitting problem.

    Notes
    -----
    Method 'lm' (Levenberg-Marquardt) calls a wrapper over least-squares
    algorithms implemented in MINPACK (lmder, lmdif). It runs
    Levenberg-Marquadrd algorithm formulated as a trust-region type algorithm.
    The implementation is based on paper [JJMore]_, it is very robust and
    efficient with a lot of smart tricks. It should be your first choice
    for unconstrained problems. Note that it doesn't support bounds. Also
    it doesn't work when m < n.

    Method 'trf' (Trust Region Reflective) is motivated by the process of
    solving a system of equations, which constitute the first-order optimality
    condition for a bound-constrained minimization problem as formulated in
    [STIR]_. The algorithm iteratively solves trust-region subproblems
    augmented by a special diagonal quadratic term and with trust-region shape
    determined by the distance from the bounds and the direction of the
    gradient. This enhancements help to avoid making steps directly into bounds
    and efficiently explore the whole space of variables. To further improve
    convergence, the algorithm considers search directions reflected from the
    bounds. To obey theoretical requirements, the algorithm keeps iterates
    strictly feasible. With dense Jacobians trust-region subproblems are
    solved by an exact method very similar to the one described in [JJMore]_
    (and implemented in MINPACK). The difference from the MINPACK
    implementation is that a singular value decomposition of a Jacobian
    matrix is done once per iteration, instead of a QR decomposition and series
    of Givens rotation eliminations. For large sparse Jacobians a 2-d subspace
    approach of solving trust-region subproblems is used [STIR]_, [Byrd]_.
    The subspace is spanned by a scaled gradient and an approximate
    Gauss-Newton solution delivered by `scipy.sparse.linalg.lsmr`. When no
    constraints are imposed the algorithm is very similar to MINPACK and has
    generally comparable performance. The algorithm works quite robust in
    unbounded and bounded problems, thus it is chosen as a default algorithm.

    Method 'dogbox' operates in a trust-region framework, but considers
    rectangular trust regions as opposed to conventional ellipsoids [Voglis]_.
    The intersection of a current trust region and initial bounds is again
    rectangular, so on each iteration a quadratic minimization problem subject
    to bound constraints is solved approximately by Powell's dogleg method
    [NumOpt]_. The required Gauss-Newton step can be computed exactly for
    dense Jacobians or approximately by `scipy.sparse.linalg.lsmr` for large
    sparse Jacobians. The algorithm is likely to exhibit slow convergence when
    the rank of Jacobian is less than the number of variables. The algorithm
    often outperforms 'trf' in bounded problems with a small number of
    variables.

    .. versionadded: 0.17.0

    References
    ----------
    .. [STIR] M. A. Branch, T. F. Coleman, and Y. Li, "A Subspace, Interior,
          and Conjugate Gradient Method for Large-Scale Bound-Constrained
          Minimization Problems," SIAM Journal on Scientific Computing,
          Vol. 21, Number 1, pp 1-23, 1999.
    .. [NR] William H. Press et. al. "Numerical Recipes. The Art of Scientific
            Computing. 3rd edition", Sec. 5.7.
    .. [Byrd] R. H. Byrd, R. B. Schnabel and G. A. Shultz, "Approximate
              solution of the trust region problem by minimization over
              two-dimensional subspaces", Math. Programming, 40, pp. 247-263,
              1988.
    .. [Curtis] A. Curtis, M. J. D. Powell, and J. Reid, "On the estimation of
                sparse Jacobian matrices", Journal of the Institute of
                Mathematics and its Applications, 13, pp. 117-120, 1974.
    .. [JJMore] J. J. More, "The Levenberg-Marquardt Algorithm: Implementation
                and Theory," Numerical Analysis, ed. G. A. Watson, Lecture
                Notes in Mathematics 630, Springer Verlag, pp. 105-116, 1977.
    .. [Voglis] C. Voglis and I. E. Lagaris, "A Rectangular Trust Region
                Dogleg Approach for Unconstrained and Bound Constrained
                Nonlinear Optimization", WSEAS International Conference on
                Applied Mathematics, Corfu, Greece, 2004.
    .. [NumOpt] J. Nocedal and S. J. Wright, "Numerical optimization,
                2nd edition", Chapter 4.

    Examples
    --------
    In this example we find a minimum of the Rosenbrock function without bounds
    on independed variables.

    >>> import numpy as np
    >>> def fun_rosenbrock(x):
    ...     return np.array([10 * (x[1] - x[0]**2), (1 - x[0])])

    Notice that we only provide the vector of the residuals. The algorithm
    constructs the cost function as a sum of squares of the residuals, which
    gives the Rosenbrock function. The exact minimum is at ``x = [1.0, 1.0]``.

    >>> from scipy.optimize import least_squares
    >>> x0_rosenbrock = np.array([2, 2])
    >>> res_1 = least_squares(fun_rosenbrock, x0_rosenbrock)
    >>> res_1.x
    array([ 1.,  1.])
    >>> res_1.cost
    2.4651903288156619e-30
    >>> res_1.optimality
    4.4408921315878507e-14

    We now constrain the variables, in such a way that the previous solution
    becomes infeasible. Specifically, we require that ``x[0] >= 1.5``, and
    ``x[1]`` left unconstrained. To this end, we specify the `bounds` parameter
    to `least_squares` in the form ``bounds=([-np.inf, 1.5], np.inf)``.

    We also provide the analytic Jacobian:

    >>> def jac_rosenbrock(x):
    ...     return np.array([
    ...         [-20 * x[0], 10],
    ...         [-1, 0]])

    Putting this all together, we see that the new solution lies on the bound:

    >>> res_2 = least_squares(fun_rosenbrock, x0_rosenbrock, jac_rosenbrock,
    ...                       bounds=([-np.inf, 1.5], np.inf))
    >>> res_2.x
    array([ 1.22437075,  1.5       ])
    >>> res_2.cost
    0.025213093946805685
    >>> res_2.optimality
    1.5885401433157753e-07

    Now we solve a system of equations (i.e., the cost function should be zero
    at a minimum) for a Broyden tridiagonal vector-valued function of 100000
    variables:

    >>> def fun_broyden(x):
    ...     f = (3 - x) * x + 1
    ...     f[1:] -= x[:-1]
    ...     f[:-1] -= 2 * x[1:]
    ...     return f

    The corresponding Jacobian matrix is sparse. We tell the algorithm to
    estimate it by finite differences and provide the sparsity structure of
    Jacobian to significantly speed up this process.

    >>> from scipy.sparse import lil_matrix
    >>> def sparsity_broyden(n):
    ...     sparsity = lil_matrix((n, n), dtype=int)
    ...     i = np.arange(n)
    ...     sparsity[i, i] = 1
    ...     i = np.arange(1, n)
    ...     sparsity[i, i - 1] = 1
    ...     i = np.arange(n - 1)
    ...     sparsity[i, i + 1] = 1
    ...     return sparsity
    ...
    >>> n = 100000
    >>> x0_broyden = -np.ones(n)
    ...
    >>> res_3 = least_squares(fun_broyden, x0_broyden,
    ...                       jac_sparsity=sparsity_broyden(n))
    >>> res_3.cost
    4.5687161966109073e-23
    >>> res_3.optimality
    1.1650454296851518e-11
    """
    if method not in ['trf', 'dogbox', 'lm']:
        raise ValueError("`method` must be 'trf', 'dogbox' or 'lm'.")

    if tr_solver not in [None, 'exact', 'lsmr']:
        raise ValueError("'tr_solver' must be None, 'exact' or 'lsmr'.")

    if verbose not in [0, 1, 2]:
        raise ValueError("`verbose` must be in [0, 1, 2].")

    if len(bounds) != 2:
        raise ValueError("`bounds` must contain 2 elements.")

    x0 = np.atleast_1d(x0).astype(float)

    if x0.ndim > 1:
        raise ValueError("`x0` must have at most 1 dimension.")

    lb, ub = prepare_bounds(bounds, x0)

    if lb.shape != x0.shape or ub.shape != x0.shape:
        raise ValueError("Inconsistent shapes between bounds and `x0`.")

    if np.any(lb >= ub):
        raise ValueError("Each lower bound mush be strictly less than each "
                         "upper bound.")

    if method == 'lm' and not np.all((lb == -np.inf) & (ub == np.inf)):
        raise ValueError("Method 'lm' doesn't support bounds.")

    if jac not in ['2-point', '3-point', 'cs'] and not callable(jac):
        raise ValueError("`jac` must be '2-point', '3-point', 'cs' or "
                         "callable.")

    scaling = check_scaling(scaling, x0)

    ftol, xtol, gtol = check_tolerance(ftol, xtol, gtol)

    if not in_bounds(x0, lb, ub):
        raise ValueError("`x0` is infeasible.")

    def fun_wrapped(x):
        return np.atleast_1d(fun(x, *args, **kwargs))

    f0 = fun_wrapped(x0)
    if f0.ndim != 1:
        raise RuntimeError("`fun` must return at most 1-d array_like.")

    n = x0.size
    m = f0.size
    if method == 'lm' and m < n:
        raise ValueError("Method 'lm' doesn't work when the number of "
                         "residuals is less than the number of variables.")

    if jac in ['2-point', '3-point', 'cs']:
        if method == 'lm':
            if jac_sparsity is not None:
                raise ValueError("Usage of method='lm' with `jac_sparsity` "
                                 "provided is forbidden.")

            if jac == '3-point':
                warn("jac='3-point' works equivalently to '2-point' "
                     "for method='lm'.")

            jac_wrapped = None
        else:
            if jac_sparsity is not None and tr_solver == 'exact':
                raise ValueError("Usage of tr_solver='exact' with "
                                 "`jac_sparsity` provided is forbidden.")

            jac_sparsity = check_jac_sparsity(jac_sparsity, m, n)

            def jac_wrapped(x, f):
                J = approx_derivative(fun, x, rel_step=diff_step, method=jac,
                                      f0=f, bounds=bounds, args=args,
                                      kwargs=kwargs, sparsity=jac_sparsity)

                if J.ndim != 2:
                    J = np.atleast_2d(J)

                return J
    else:
        def jac_wrapped(x, _=None):
            J = jac(x, *args, **kwargs)

            if issparse(J):
                J = csr_matrix(J)
            elif not isinstance(J, LinearOperator):
                J = np.atleast_2d(J)

            return J

    if jac_wrapped is not None:
        J0 = jac_wrapped(x0, f0)

        if J0.shape != (m, n):
            raise ValueError(
                "The return value of `jac` has wrong shape: expected {0}, "
                "actual {1}.".format((m, n), J0.shape))

        if not isinstance(J0, np.ndarray):
            if method == 'lm':
                raise ValueError("method='lm' works only with dense Jacobian.")

            if tr_solver == 'exact':
                raise ValueError(
                    "tr_solver='exact' works only with dense Jacobian.")

        if isinstance(J0, LinearOperator) and scaling == 'jac':
            raise ValueError("scaling='jac' can't be used when `jac` "
                             "returns LinearOperator.")

        if tr_solver is None:
            if isinstance(J0, np.ndarray):
                tr_solver = 'exact'
            else:
                tr_solver = 'lsmr'

    if method == 'lm':
        result = call_minpack(fun_wrapped, x0, jac_wrapped, ftol, xtol, gtol,
                              max_nfev, scaling, diff_step)

    if method == 'trf':
        result = trf(fun_wrapped, jac_wrapped, x0, f0, J0, lb, ub, ftol, xtol,
                     gtol, max_nfev, scaling, tr_solver, tr_options.copy(),
                     verbose)

    elif method == 'dogbox':
        result = dogbox(fun_wrapped, jac_wrapped, x0, f0, J0, lb, ub, ftol,
                        xtol, gtol, max_nfev, scaling, tr_solver, tr_options,
                        verbose)

    result.message = TERMINATION_MESSAGES[result.status]
    result.success = result.status > 0

    if verbose >= 1:
        print(result.message)
        print("Function evaluations: {0}, initial cost: {1:.4e}, final cost "
              "{2:.4e}, first-order optimality {3:.2e}."
              .format(result.nfev, 0.5 * np.dot(f0, f0),
                      result.cost, result.optimality))

    return result
