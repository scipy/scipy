"""Generic interface for least-square minimization."""

from warnings import warn

import numpy as np
from numpy.linalg import norm

from . import _minpack
from .optimize import OptimizeResult
from ._numdiff import approx_derivative, group_columns
from ..sparse import issparse, csr_matrix, csc_matrix
from ..sparse.linalg import LinearOperator

from ._lsq_bounds import in_bounds, prepare_bounds
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

    obj_value = np.dot(f, f)
    g = J.T.dot(f)
    g_norm = norm(g, ord=np.inf)

    nfev = info['nfev']
    njev = info.get('njev', None)

    status = FROM_MINPACK_TO_COMMON[status]
    active_mask = np.zeros_like(x0, dtype=int)

    return OptimizeResult(
        x=x, fun=f, jac=J, obj_value=obj_value, optimality=g_norm,
        active_mask=active_mask, nfev=nfev, njev=njev, status=status)


def check_scaling(scaling, x0):
    if scaling == 'jac':
        return scaling
    try:
        scaling = np.asarray(scaling, dtype=float)
    except ValueError:
        raise ValueError("`scaling` must be 'jac' or array_like with numbers.")

    if np.any(scaling <= 0):
        raise ValueError("`scaling` must contain only positive values.")

    if scaling.ndim == 0:
        scaling = np.resize(scaling, x0.shape)

    if scaling.shape != x0.shape:
        raise ValueError("Inconsistent shapes between `scaling` and `x0`.")

    return scaling


def least_squares(
        fun, x0, jac='2-point', bounds=(-np.inf, np.inf), method='trf',
        ftol=EPS**0.5, xtol=EPS**0.5, gtol=EPS**0.5, scaling=1.0,
        diff_step=None, tr_solver=None, tr_options={}, jac_sparsity=None,
        max_nfev=None, args=(), kwargs={}):
    """Minimize the sum of squares of nonlinear functions, subject to bound
    constraints on independent variables.

    Let f(x) maps from R^n to R^m, this function finds a local minimum of::

        F(x) = ||f(x)||**2 = sum(f_i(x)**2, i = 1, ..., m), lb <= x <= ub

    f(x) is called vector of residuals or simply residuals.

    Partial derivatives of f with respect to x form m-by-n matrix called
    Jacobian, where an element (i, j) equals to a partial derivative of f[i]
    with respect to x[j].

    Parameters
    ----------
    fun : callable
        Function which computes a vector of residuals. The argument x passed
        to this function is ndarray of shape (n,) (never a scalar, even
        if n=1). It must return 1-d array_like of shape (m,) or a scalar.
    x0 : array_like with shape (n,) or float
        Initial guess on independent variables. If float, for internal usage
        it will be converted to 1-d array.
    jac : '2-point', '3-point' or callable, optional
        Method of computing the Jacobian matrix. If set to '2-point' or
        '3-point', the Jacobian matrix is estimated by the corresponding
        finite difference scheme. The scheme '3-point' is more accurate, but
        requires twice as much operations compared to '2-point' (default).
        If callable, then it should return a reasonable approximation of
        Jacobian as array_like (np.atleast_2d is applied), sparse matrix or
        LinearOperator, all with shape  (m, n).
    bounds : 2-tuple of array_like, optional
        Lower and upper bounds on independent variables. Defaults to no bounds.
        Each array must match the size of `x0` or be a scalar, in the latter
        case a bound will be the same for all variables. Use ``np.inf`` with
        an appropriate sign to disable bounds on all or some variables.
    method : {'trf', 'dogbox', 'lm'}, optional
        Algorithm to perform optimization.
            * 'trf' - Trust Region Reflective algorithm, particularly suitable
              for large sparse problems with bounds. Generally robust method.
            * 'dogbox' - dogleg algorithm with rectangular trust regions,
              typical use case is small problems with bounds. Not recommended
              to use in rank-deficient problems.
            * 'lm' - Levenberg-Marquardt algorithm as implemented in MINPACK.
              Doesn't handle bounds and sparse Jacobians, but usually most
              efficient for small unconstrained problems.
        Default is 'trf'. See Notes for more information.
    ftol : float, optional
        Tolerance for termination by the change of the objective value.
        Default is square root of machine epsilon. The optimization process is
        stopped when ``dF < ftol * F``, and there was adequate agreement
        between a local quadratic model and the true model in the last step.
    xtol : float, optional
        Tolerance for termination by the change of the independent variables.
        Default is square root of machine epsilon. The exact condition checked
        depends on a `method` used:
            * For 'trf' and 'dogbox' : ``norm(dx) < xtol * (xtol + norm(x))``
            * For 'lm': ``Delta < xtol * norm(scaled_x)``, where Delta is a
              trust-region radius and scaled_x is a scaled value of x
              according to `scaling` parameter (see below).
    gtol : float, optional
        Tolerance for termination by the norm of gradient. Default is
        square root of machine epsilon. The exact condition depends
        on a `method` used:

            * For 'trf' : ``norm(g_scaled, ord=np.inf) < gtol``, where
              g_scaled is properly scaled gradient to account for the
              presence of bounds [STIR]_.
            * Fot 'dogbox' : ``norm(g_free, ord=np.inf) < gtol``, where
              g_free is the gradient with respect to the variables which
              aren't in the optimal state on the boundary.
            * For 'lm' : the maximum cosine of angles between columns of
              Jacobian and the residual vector is less than `gtol`, or the
              residual vector is zero.
    max_nfev : None or int, optional
        Maximum number of function evaluations before the termination.
        If None (default), the value is chosen automatically:
            * For 'trf' and 'dogbox' : 100 * n.
            * For 'lm':  100 * n if `jac` is callable and 100 * n * (n + 1)
              otherwise (because 'lm' counts function calls in Jacobian
              estimation).
    scaling : array_like or 'jac', optional
        Applies variables scaling to potentially improve algorithm convergence.
        Default is 1.0, which means no scaling. Scaling should be used to
        equalize the influence of each variable on the objective function.
        Alternatively you can think of `scaling` as diagonal elements of
        a matrix which determines the shape of a trust region. Use smaller
        values for variables which have bigger characteristic scale compared
        to others. A scalar value won't affect the algorithm (except maybe
        fixing/introducing numerical issues and changing termination criteria).
        If 'jac', then scaling is proportional to the norms of Jacobian
        columns. The latter option is often helpful in unconstrained problems,
        but not so much in constrained ones. From experience usage of 'jac'
        scaling is not recommended for bounded problems with 'trf' method.
    diff_step : None or array_like, optional
        Determines the step size for finite difference Jacobian approximation.
        The actual step is computed as ``x * diff_step``. If None (default),
        `diff_step` is assigned to a conventional "optimal" power of machine
        epsilon depending on a finite difference approximation method [NR]_.
    tr_solver : {None, 'exact', 'lsmr'}, optional
        Method for solving trust-region subproblems, relevant only for 'trf'
        and 'dogbox' methods.
            * 'exact' is suitable for not very large problems, which have
              dense Jacobian matrices. It requires work comparable to a single
              SVD of Jacobian per iteration.
            * 'lsmr' is suitable for problems with sparse and large Jacobian
              matrices. It uses iterative ``scipy.sparse.linalg.lsmr``
              procedure for finding a solution of linear least squares and
              requires only matrix-vector product evaluations.
        If None (default) the solver is chosen based on type of Jacobian
        returned on the first iteration.
    tr_options : dict, optional
        Keyword options passed to trust-region solver.
            * ``tr_solver='exact'`` : `tr_options` are ignored.
            * ``tr_solver='lsmr'`` : options for ``scipy.sparse.linalg.lsmr``.
              Additionally  ``method='trf'`` supports 'regularize' option
              (bool, default True) - add regularization term to the normal
              equation, which improves convergence with rank-deficient
              Jacobian [Byrd]_ (eq. 3.4).
    jac_sparsity : {None, array_like, sparse matrix}, optional
        Defines Jacobian sparsity structure for finite differencing. Provide
        this parameter to greatly speed up finite difference Jacobian
        estimation, if it has only few non-zeros in *each* row [Curtis]_.
        Should be array_like or sparse matrix with shape (m, n). A zero element
        means that a corresponding element in Jacobian is identically zero.
        Forces `tr_solver` to 'lsmr' if it wasn't set. If None (default) then
        dense differencing will be used. Has no effect for ``method='lm'``.
    args, kwargs : tuple and dict, optional
        Additional arguments passed to `fun` and `jac`. Both empty by default.
        The calling signature is ``fun(x, *args, **kwargs)`` and the same for
        `jac`.

    Returns
    -------
    OptimizeResult with the following fields defined.
    x : ndarray, shape (n,)
        Found solution.
    obj_value : float
        Sum of squares at the solution.
    fun : ndarray, shape (m,)
        Vector of residuals at the solution.
    jac : ndarray, sparse matrix or LinearOperator, shape (m, n)
        Jacobian matrix at the solution. The type is the same as was used by
        the algorithm.
    optimality : float
        First-order optimality measure. In unconstrained problems it is always
        the uniform norm of the gradient. In constrained problems this is the
        quantity which was compared with `gtol` during iterations.
    active_mask : ndarray of int, shape (n,)
        Each component shows whether a corresponding constraint is active
        (a variable is on the bound):
            *  0 - a constraint is not active.
            * -1 - a lower bound is active.
            *  1 - an upper bound is active.
        Might be somewhat arbitrary for 'trf' method as it does strictly
        feasible iterates and `active_mask` is determined with tolerance
        threshold.
    nfev : int
        Number of function evaluations done. Methods 'trf' and 'dogbox' don't
        count function calls for numerical Jacobian approximation, opposed to
        'lm' method.
    njev : int or None
        Number of Jacobian evaluations done. If numerical Jacobian
        approximation is used in 'lm' method, it is set to None.
    status : int
        Reason for algorithm termination:
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
    solving a system of equations, which constitutes the first-order optimality
    condition for a bound-constrained minimization problem as formulated in
    [STIR]_. The algorithm iteratively solves trust-region subproblems
    augmented by special diagonal quadratic term and with trust-region shape
    determined by the distance from the bounds and the direction of the
    gradient. This enhancements help not to take steps directly into bounds
    and explore the whole variable space. To improve convergence speed the
    reflected from the first bound search direction is considered. To obey
    theoretical requirements the algorithm keeps iterates strictly feasible.
    Trust-region subproblems are solved by exact method very similar to one
    described in [JJMore]_ and implemented in MINPACK, but with the help of
    one per iteration singular value decomposition of the Jacobian matrix.
    The algorithm's performance is generally comparable to MINPACK in
    unbounded problems. The algorithm works quite robust in unbounded and
    bounded problems, thus it is chosen as a default algorithm.

    Method 'dogbox' operates in a trust-region framework, but considers
    rectangular trust regions as opposed to conventional elliptical [Voglis]_.
    The intersection of a current trust region and initial bounds is again
    rectangular, so on each iteration a quadratic minimization problem subject
    to bounds is solved approximately by Powell's dogleg method [NumOpt]_.
    The algorithm is likely to exhibit slow convergence when the rank of
    Jacobian is less than the number of variables. The algorithm often
    outperform 'trf' in bounded problems with small number of variables.

    References
    ----------
    .. [STIR] Branch, M.A., T.F. Coleman, and Y. Li, "A Subspace, Interior,
          and Conjugate Gradient Method for Large-Scale Bound-Constrained
          Minimization Problems," SIAM Journal on Scientific Computing,
          Vol. 21, Number 1, pp 1-23, 1999.
    .. [NR] William H. Press et. al. "Numerical Recipes. The Art of Scientific
            Computing. 3rd edition", Sec. 5.7.
    .. [Byrd] R. H. Byrd, R. B. Schnabel and G. A. Shultz, "Approximate
              solution of the trust region problem by minimization over
              two-dimensional subspaces", Math. Programming, 40 (1988),
              pp. 247-263.
    .. [Curtis] A. Curtis, M. J. D. Powell, and J. Reid, "On the estimation of
                sparse Jacobian matrices", Journal of the Institute of
                Mathematics and its Applications, 13 (1974), pp. 117-120.
    .. [JJMore] More, J. J., "The Levenberg-Marquardt Algorithm: Implementation
                and Theory," Numerical Analysis, ed. G. A. Watson, Lecture Notes
                in Mathematics 630, Springer Verlag, pp. 105-116, 1977.
    .. [Voglis] C. Voglis and I. E. Lagaris, "A Rectangular Trust Region
                Dogleg Approach for Unconstrained and Bound Constrained
                Nonlinear Optimization", WSEAS International Conference on
                Applied Mathematics, Corfu, Greece, 2004.
    .. [NumOpt] J. Nocedal and S. J. Wright, "Numerical optimization,
                2nd edition", Chapter 4.
    """
    if method not in ['trf', 'dogbox', 'lm']:
        raise ValueError("`method` must be 'trf', 'dogbox' or 'lm'.")

    if tr_solver not in [None, 'exact', 'lsmr']:
        raise ValueError("'tr_solver' must be None, 'exact' or 'lsmr'.")

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

    if jac not in ['2-point', '3-point'] and not callable(jac):
        raise ValueError("`jac` must be '2-point', '3-point' or callable.")

    scaling = check_scaling(scaling, x0)

    ftol, xtol, gtol = check_tolerance(ftol, xtol, gtol)

    if not in_bounds(x0, lb, ub):
        raise ValueError("`x0` is infeasible.")

    def fun_wrapped(x):
        return np.atleast_1d(fun(x, *args, **kwargs))

    if jac in ['2-point', '3-point']:
        if jac_sparsity is not None:
            if method == 'lm':
                warn("`jac_sparsity` is ignored for method='lm', dense "
                     "differencing will be used.")
            else:
                structure = csc_matrix(jac_sparsity)
                groups = group_columns(structure)
                sparsity = (structure, groups)
        else:
            sparsity = None

        def jac_wrapped(x, f):
            J = approx_derivative(
                fun, x, rel_step=diff_step, method=jac, f0=f, bounds=bounds,
                args=args, kwargs=kwargs, sparsity=sparsity)

            if issparse(J) and tr_solver == 'exact':
                warn("Jacobian is converted to dense for tr_solver='exact' "
                     "even though `jac_sparsity` is provided. Consider using "
                     "'lsmr' trust-region solver or set `jac_sparsity` to "
                     "None.")
                J = J.toarray()
            elif not issparse(J):
                J = np.atleast_2d(J)

            return J
    else:
        def jac_wrapped(x, _=None):
            J = jac(x, *args, **kwargs)

            if issparse(J):
                if method == 'lm':
                    warn("Jacobian is converted to dense for "
                         "method='lm', consider using a different `method` "
                         "or return dense Jacobian.")
                    J = J.toarray()
                elif tr_solver == 'exact':
                    warn("Jacobian is converted to dense for "
                         "tr_solver='exact', consider using 'lsmr' "
                         "trust-region solver or return dense Jacobian.")
                    J = J.toarray()
                else:
                    J = csr_matrix(J)
            elif not isinstance(J, LinearOperator):
                J = np.atleast_2d(J)

            return J

    f0 = fun_wrapped(x0)
    if f0.ndim != 1:
        raise RuntimeError("`fun` must return at most 1-d array_like.")

    if method == 'lm' and f0.size < x0.size:
        raise ValueError("Method 'lm' doesn't work when the number of "
                         "variables is less than the number of residuals.")

    J0 = jac_wrapped(x0, f0)
    if len(J0.shape) != 2:
        raise RuntimeError("`jac` must return at most 2-d array_like, "
                           "sparse matrix or LinearOperator.")

    if f0.shape[0] != J0.shape[0]:
        raise RuntimeError("Inconsistent dimensions between the "
                           "returns of `fun` and `jac`.")

    if isinstance(J0, LinearOperator):
        if method == 'lm':
            raise ValueError("method='lm' can't be used when `jac` "
                             "returns LinearOperator.")

        if tr_solver == 'exact':
            raise ValueError("tr_solver='exact' can't be used when `jac` "
                             "returns LinearOperator.")

        if scaling == 'jac':
            raise ValueError("scaling='jac' can't be used when `jac` "
                             "return LinearOperator.")

    if tr_solver is None:
        if isinstance(J0, np.ndarray):
            tr_solver = 'exact'
        else:
            tr_solver = 'lsmr'

    if method == 'lm':
        if jac == '2-point':
            jac_wrapped = None
        elif jac == '3-point':
            jac_wrapped = None
            warn("jac='3-point' works equivalently to '2-point' "
                 "for 'lm' method.")

        if scaling == 'jac':
            scaling = None

        result = call_minpack(fun_wrapped, x0, jac_wrapped, ftol, xtol, gtol,
                              max_nfev, scaling, diff_step)

    if method == 'trf':
        result = trf(fun_wrapped, jac_wrapped, x0, f0, J0, lb, ub, ftol, xtol,
                     gtol, max_nfev, scaling, tr_solver, tr_options.copy())

    elif method == 'dogbox':
        result = dogbox(fun_wrapped, jac_wrapped, x0, f0, J0, lb, ub, ftol,
                        xtol, gtol, max_nfev, scaling, tr_solver, tr_options)

    result.message = TERMINATION_MESSAGES[result.status]
    result.success = result.status > 0
    return result
