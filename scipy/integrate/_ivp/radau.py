import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse import csc_matrix, issparse, eye, diags
from scipy.sparse.linalg import splu
from scipy.optimize._numdiff import group_columns
from .common import (validate_max_step, validate_tol, select_initial_step,
                     norm, num_jac, EPS, warn_extraneous,
                     validate_first_step)
from .base import DaeSolver, DenseOutput

S6 = 6 ** 0.5

# Butcher tableau. A is not used directly, see below.
C = np.array([(4 - S6) / 10, (4 + S6) / 10, 1])
E = np.array([-13 - 7 * S6, -13 + 7 * S6, -1]) / 3

# Eigendecomposition of A is done: A = T L T**-1. There is 1 real eigenvalue
# and a complex conjugate pair. They are written below.
MU_REAL = 3 + 3 ** (2 / 3) - 3 ** (1 / 3) # inverse of the real eigenvalue of A
MU_COMPLEX = (3 + 0.5 * (3 ** (1 / 3) - 3 ** (2 / 3))
              - 0.5j * (3 ** (5 / 6) + 3 ** (7 / 6)))

# These are transformation matrices.
T = np.array([
    [0.09443876248897524, -0.14125529502095421, 0.03002919410514742],
    [0.25021312296533332, 0.20412935229379994, -0.38294211275726192],
    [1, 1, 0]])
TI = np.array([
    [4.17871859155190428, 0.32768282076106237, 0.52337644549944951],
    [-4.17871859155190428, -0.32768282076106237, 0.47662355450055044],
    [0.50287263494578682, -2.57192694985560522, 0.59603920482822492]])

# These linear combinations are used in the algorithm.
TI_REAL = TI[0]
TI_COMPLEX = TI[1] + 1j * TI[2]

# Interpolator coefficients.
P = np.array([
    [13/3 + 7*S6/3, -23/3 - 22*S6/3, 10/3 + 5 * S6],
    [13/3 - 7*S6/3, -23/3 + 22*S6/3, 10/3 - 5 * S6],
    [1/3, -8/3, 10/3]])



def predict_factor(h_abs, h_abs_old, error_norm, error_norm_old):
    """Predict by which factor to increase/decrease the step size.

    The algorithm is described in [1]_.

    Parameters
    ----------
    h_abs, h_abs_old : float
        Current and previous values of the step size, `h_abs_old` can be None
        (see Notes).
    error_norm, error_norm_old : float
        Current and previous values of the error norm, `error_norm_old` can
        be None (see Notes).

    Returns
    -------
    factor : float
        Predicted factor.

    Notes
    -----
    If `h_abs_old` and `error_norm_old` are both not None then a two-step
    algorithm is used, otherwise a one-step algorithm is used.

    References
    ----------
    .. [1] E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential
           Equations II: Stiff and Differential-Algebraic Problems", Sec. IV.8.
    """
    if error_norm_old is None or h_abs_old is None or error_norm == 0:
        multiplier = 1
    else:
        multiplier = h_abs / h_abs_old * (error_norm_old / error_norm) ** 0.25

    with np.errstate(divide='ignore'):
        factor = min(1, multiplier) * error_norm ** -0.25

    return factor


class Radau(DaeSolver):
    """Implicit Runge-Kutta method of Radau IIA family of order 5.

    The implementation follows [1]_. The error is controlled with a
    third-order accurate embedded formula. A cubic polynomial which satisfies
    the collocation conditions is used for the dense output.
    Specific treatment of systems of differential-algebraic equations (DAEs) follows [3]_.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system: the time derivative of the state ``y``
        at time ``t``. The calling signature is ``fun(t, y)``, where ``t`` is a
        scalar and ``y`` is an ndarray with ``len(y) = len(y0)``. ``fun`` must
        return an array of the same shape as ``y``. See `vectorized` for more
        information.
    t0 : float
        Initial time.
    y0 : array_like, shape (n,)
        Initial state.
    t_bound : float
        Boundary time - the integration won't continue beyond it. It also
        determines the direction of the integration.
    first_step : float or None, optional
        Initial step size. Default is ``None`` which means that the algorithm
        should choose.
    max_step : float, optional
        Maximum allowed step size. Default is np.inf, i.e., the step size is not
        bounded and determined solely by the solver.
    rtol, atol : float and array_like, optional
        Relative and absolute tolerances. The solver keeps the local error
        estimates less than ``atol + rtol * abs(y)``. Here `rtol` controls a
        relative accuracy (number of correct digits), while `atol` controls
        absolute accuracy (number of correct decimal places). To achieve the
        desired `rtol`, set `atol` to be smaller than the smallest value that
        can be expected from ``rtol * abs(y)`` so that `rtol` dominates the
        allowable error. If `atol` is larger than ``rtol * abs(y)`` the
        number of correct digits is not guaranteed. Conversely, to achieve the
        desired `atol` set `rtol` such that ``rtol * abs(y)`` is always smaller
        than `atol`. If components of y have different scales, it might be
        beneficial to set different `atol` values for different components by
        passing array_like with shape (n,) for `atol`. Default values are
        1e-3 for `rtol` and 1e-6 for `atol`.
    jac : {None, array_like, sparse_matrix, callable}, optional
        Jacobian matrix of the right-hand side of the system with respect to
        y, required by this method. The Jacobian matrix has shape (n, n) and
        its element (i, j) is equal to ``d f_i / d y_j``.
        There are three ways to define the Jacobian:

            * If array_like or sparse_matrix, the Jacobian is assumed to
              be constant.
            * If callable, the Jacobian is assumed to depend on both
              t and y; it will be called as ``jac(t, y)`` as necessary.
              For the 'Radau' and 'BDF' methods, the return value might be a
              sparse matrix.
            * If None (default), the Jacobian will be approximated by
              finite differences.

        It is generally recommended to provide the Jacobian rather than
        relying on a finite-difference approximation.
    jac_sparsity : {None, array_like, sparse matrix}, optional
        Defines a sparsity structure of the Jacobian matrix for a
        finite-difference approximation. Its shape must be (n, n). This argument
        is ignored if `jac` is not `None`. If the Jacobian has only few non-zero
        elements in *each* row, providing the sparsity structure will greatly
        speed up the computations [2]_. A zero entry means that a corresponding
        element in the Jacobian is always zero. If None (default), the Jacobian
        is assumed to be dense.
    vectorized : bool, optional
        Whether `fun` can be called in a vectorized fashion. Default is False.

        If ``vectorized`` is False, `fun` will always be called with ``y`` of
        shape ``(n,)``, where ``n = len(y0)``.

        If ``vectorized`` is True, `fun` may be called with ``y`` of shape
        ``(n, k)``, where ``k`` is an integer. In this case, `fun` must behave
        such that ``fun(t, y)[:, i] == fun(t, y[:, i])`` (i.e. each column of
        the returned array is the time derivative of the state corresponding
        with a column of ``y``).

        Setting ``vectorized=True`` allows for faster finite difference
        approximation of the Jacobian by this method, but may result in slower
        execution overall in some circumstances (e.g. small ``len(y0)``).
    mass_matrix : {None, array_like, sparse_matrix}, shape (n,n), optional
        Defines the constant mass matrix M of the system, with shape (n,n).
        The problem considered is then of the form:
          ``M y' = fun(t,y)``
        This matrix may be singular, thus defining a problem of the differential-
        algebraic type (DAE), see [1]. The default value is None, correspoding to the identity matrix, i.e. all components are of the differential nature.
    var_index : {None, array_like}, shape (n,), optional
        In the case of a differential-algebraic system (DAE), i.e. a singular matrix is provided with the `mass_matrix`argument, this vector of integers of shape (n,) defines the algebraic index of each component. This index is 0 for differential components, and is larger than 0 for algebraic components. Note  that the Radau solver should be able to handle DAEs with components up to index 3. The values of the algebraic indices are used to apply various scalings (for the Newton iterations and the error estimate) to improve the robustness of the algorithm for DAEs.
    max_newton_ite : int, optional
        Defines the maximum number of Newton iterations that can be performed to solve each time step.

    bUsePredictiveNewtonStoppingCriterion : boolean, optional
        During the Newton iterations, their convergence is monitored. If this parameter is true (default value) and it is predicted that the Newton error level reached after performing all allowed  iteration (``max_newton_ite``) is too large, the iterations are stopped. This may greatly reduce the computational cost of the integration, by avoiding uneffective Newton iterations. Still, in particular when dealing with higher-index DAEs, this predictive monitoring may be too conservative and wrongly predict a poor convergence of the iterations. On certain problems, this may lead to a infinite loop of time step reduction. Setting this parameter to False may correct this issue.

    max_bad_ite: int, optional
      It is possible that one Newton iteration is such that residual or increment norm increases compared to the previous iteration. This is referred to as a bad iteration. Originally, the Newton iterations are stopped when such an event occurs, triggering a time step reduction. For certain problems, this may lead to unnecessary time reductions, potentially even leading to a failure of the integration. For some differential-algebraic equations (DAEs) of index higher than 1. It turns out that allowing one or more such bad iterations to occur without stopping the Newton iterations may lead to a successful convergence. Usually, setting `max_bad_ite` to 1 is a good all-around choice. It is set to 0 by default, or to 1 if a DAE problem of index higher than 1 (specified via ``var_index``). See [3, page TODO] for more details.

    scale_residuals : boolean, optional
      If a system of differential-algebraic equations (DAEs) is considered
      (see ``mass_matrix``), and if this parameter is true, the residuals
      of the implicit system (solved by a Newton method at each time step)
      are scaled component by component by h^i, i being the algebraic index
      of the component, and h the time step size.
      Defaults to true.
      This usually improves the robustness of the solver for DAEs of index higher than 1.
      See [3, page TODO] for more details.

    scale_newton_norm : boolean, optional
      If a system of differential-algebraic equations (DAEs) is considered
      (see ``mass_matrix``), and if this parameter is true, the convergence
      of the Newton iterations are evaluated based on a modified norm, where
      the components corresponding to algebraic variables are scaled by h^i,
      i being the algebraic index of the component, and h the time step size.
      Defaults to true.
      This usually improves the robustness of the solver for DAEs of index higher than 1.
      See [3, page TODO] for more details.

    scale_error : boolean, optional
      If a system of differential-algebraic equations (DAEs) is considered
      (see ``mass_matrix``), and if this parameter is true, the norm of the
      error estimate used for adapting the time step is computed with a
      specific scaling for the algberaic components. Each component of the
      error estimate is scaled by h^i, i being the algebraic index of the
      component, and h the time step size. Defaults to true. This lowers
      the estimated error on the algebraic variables, which is usually
      overestimated. This usually leads to larger step sizes for DAEs,
      while maintaining the desired accuracy on the differential components.

    zero_algebraic_error : boolean, optional
      If a system of differential-algebraic equations (DAEs) is considered
      (see ``mass_matrix``), and if this parameter is true, the norm of the
      error estimate is computed by considering only the differential
      components and discarding the error estimated on the algebraic
      components. Defaults to false. This may help with DAEs for which the
      time step seems to be unnecesseraily low, even when setting
      ``scale_error=True``. In any case, the error on the differential
      components remains below the prescribed accuracy.



    Attributes
    ----------
    n : int
        Number of equations.
    status : string
        Current status of the solver: 'running', 'finished' or 'failed'.
    t_bound : float
        Boundary time.
    direction : float
        Integration direction: +1 or -1.
    t : float
        Current time.
    y : ndarray
        Current state.
    t_old : float
        Previous time. None if no steps were made yet.
    step_size : float
        Size of the last successful step. None if no steps were made yet.
    nfev : int
        Number of evaluations of the right-hand side.
    njev : int
        Number of evaluations of the Jacobian.
    nlu : int
        Number of LU decompositions.
    nlusolve: int
        Number of linear systems solved with the LU decompositions
        (2 per Newton iteration, 1 or 2 per error estimation)
      

    References
    ----------
    .. [1] E. Hairer, G. Wanner, "Solving Ordinary Differential Equations II:
           Stiff and Differential-Algebraic Problems", Sec. IV.8.
    .. [2] A. Curtis, M. J. D. Powell, and J. Reid, "On the estimation of
           sparse Jacobian matrices", Journal of the Institute of Mathematics
           and its Applications, 13, pp. 117-120, 1974.
    .. [3] E. Hairer, C. Lubich, M. Roche, "The Numerical Solution of
           Differential-Algebraic Systems by Runge-Kutta Methods"
    """
    def __init__(self, fun, t0, y0, t_bound, max_step=np.inf,
                 rtol=1e-3, atol=1e-6, jac=None, jac_sparsity=None,
                 vectorized=False, first_step=None,
                 max_newton_ite=6, max_bad_ite=None,
                 mass_matrix=None, var_index=None,
                 bUsePredictiveNewtonStoppingCriterion=True,
                 scale_residuals=True, scale_newton_norm=True,
                 scale_error=True, zero_algebraic_error=False,
                 **extraneous):

        warn_extraneous(extraneous)
        super().__init__(fun, t0, y0, t_bound, vectorized)
        self.y_old = None
        self.max_step = validate_max_step(max_step)
        self.rtol, self.atol = validate_tol(rtol, atol, self.n)
        self.f = self.fun(self.t, self.y)
        # Select initial step assuming the same order which is used to control
        # the error.
        if first_step is None:
            if mass_matrix is None:
              self.h_abs = select_initial_step(
                  self.fun, self.t, self.y, t_bound, max_step, self.f, self.direction,
                  3, self.rtol, self.atol)
            else: # as in [1], default to 1e-6
                self.h_abs = self.direction * min(1e-6, abs(self.t_bound-self.t0))
        else:
            self.h_abs = validate_first_step(first_step, t0, t_bound)
        self.h_abs_old = None
        self.error_norm_old = None

        # Convergence tolerance for the Newton iterations
        # (relative to the integration error tolerance)
        self.newton_tol = max(10 * EPS / rtol, min(0.03, rtol ** 0.5))
        # see [1] and the original Radau5 Fortran code
        
        self.sol = None
        self.jac_factor = None
        self.jac, self.J = self._validate_jac(jac, jac_sparsity)
        self.nlusolve = 0

        if issparse(self.J):
            def lu(A):
                self.nlu += 1
                return splu(A)

            def solve_lu(LU, b):
                self.nlusolve += 1
                return LU.solve(b)

            I = eye(self.n, format='csc')
        else:
            def lu(A):
                self.nlu += 1
                return lu_factor(A, overwrite_a=True)

            def solve_lu(LU, b):
                self.nlusolve += 1
                return lu_solve(LU, b, overwrite_b=True)

            I = np.identity(self.n)

        self.lu = lu
        self.solve_lu = solve_lu
        self.I = I

        # DAE-specific treatment
        self.scale_residuals = scale_residuals
        self.scale_error = scale_error
        self.scale_newton_norm = scale_newton_norm
        self.zero_algebraic_error = zero_algebraic_error

        self.hscale = self.I
        self.mass_matrix = self._validate_mass_matrix(mass_matrix)
        if var_index is None: # vector of the algebraic index of each variable
            self.var_index = np.zeros((y0.size,)) # assume all differential
        else:
            assert isinstance(var_index, np.ndarray), '`var_index` must be an array'
            assert var_index.ndim == 1
            assert var_index.size == y0.size
            self.var_index = var_index

        if mass_matrix is None:
            assert (var_index is None or np.all(var_index==0)),
            """``var_index`` should not be specified, or only contains zeros,
            when ``mass_matrix`` is not specified"""
        
        self.var_exp = self.var_index - 1 # for DAE-specific scalings
        self.var_exp[ self.var_exp < 0 ] = 0 # for differential components

        if not ( max_bad_ite is None ):
            self.NMAX_BAD = max_bad_ite # Maximum number of bad Newton iterations per step
        else: # by default, if the DAE has index>1, allow one bad iteration
            if np.any(self.var_index>1):
                self.NMAX_BAD = 1
            else:
                self.NMAX_BAD = 0

        self.index_algebraic_vars = np.where( self.var_index != 0 )[0] #self.var_index != 0
        self.nvars_algebraic = self.index_algebraic_vars.size


        self.current_jac = True
        self.LU_real = None
        self.LU_complex = None
        self.Z = None

    def _validate_mass_matrix(self, mass_matrix):
        if mass_matrix is None:
            M = self.I
        elif callable(mass_matrix):
            raise ValueError("`mass_matrix` should be a constant matrix, but is"
                             " callable")
        else:
            if issparse(mass_matrix):
                M = csc_matrix(mass_matrix)
            else:
                M = np.asarray(mass_matrix, dtype=float)
            if M.shape != (self.n, self.n):
                raise ValueError("`mass_matrix` is expected to have shape {}, "
                                 "but actually has {}."
                                 .format((self.n, self.n), M.shape))
        return M

    def _validate_jac(self, jac, sparsity):
        t0 = self.t
        y0 = self.y

        if jac is None:
            if sparsity is not None:
                if issparse(sparsity):
                    sparsity = csc_matrix(sparsity)
                groups = group_columns(sparsity)
                sparsity = (sparsity, groups)

            def jac_wrapped(t, y, f):
                self.njev += 1
                J, self.jac_factor = num_jac(self.fun_vectorized, t, y, f,
                                             self.atol, self.jac_factor,
                                             sparsity)
                return J
            J = jac_wrapped(t0, y0, self.f)
        elif callable(jac):
            J = jac(t0, y0)
            self.njev = 1
            if issparse(J):
                J = csc_matrix(J)

                def jac_wrapped(t, y, _=None):
                    self.njev += 1
                    return csc_matrix(jac(t, y), dtype=float)

            else:
                J = np.asarray(J, dtype=float)

                def jac_wrapped(t, y, _=None):
                    self.njev += 1
                    return np.asarray(jac(t, y), dtype=float)

            if J.shape != (self.n, self.n):
                raise ValueError(f"`jac` is expected to have shape {(self.n, self.n)},"
                                 f" but actually has {J.shape}.")
        else:
            if issparse(jac):
                J = csc_matrix(jac)
            else:
                J = np.asarray(jac, dtype=float)

            if J.shape != (self.n, self.n):
                raise ValueError(f"`jac` is expected to have shape {(self.n, self.n)},"
                                 f" but actually has {J.shape}.")
            jac_wrapped = None

        return jac_wrapped, J

    def _step_impl(self):
        t = self.t
        y = self.y
        f = self.f
        n = y.size

        max_step = self.max_step
        atol = self.atol
        rtol = self.rtol

        min_step = max(1e-20, 10 * np.abs(np.nextafter(t, self.direction * np.inf) - t) )
        if self.h_abs > max_step:
            h_abs = max_step
            h_abs_old = None
            error_norm_old = None
        elif self.h_abs < min_step:
            h_abs = min_step
            h_abs_old = None
            error_norm_old = None
        else:
            h_abs = self.h_abs
            h_abs_old = self.h_abs_old
            error_norm_old = self.error_norm_old

        if abs(self.t_bound - (t + self.direction * h_abs)) < 1e-2*h_abs:
            # the next time step would be too small, hence we slightly
            # increase the step size to reach the final time directly
            h_abs = abs(self.t_bound - t)
            # require refactorization of the iteration matrices
            self.LU_real = None
            self.LU_complex = None

        J = self.J
        LU_real = self.LU_real
        LU_complex = self.LU_complex

        current_jac = self.current_jac
        jac = self.jac

        rejected = False
        step_accepted = False
        message = None
        while not step_accepted:
            if h_abs < min_step:
                return False, self.TOO_SMALL_STEP

            h = h_abs * self.direction
            t_new = t + h

            if self.direction * (t_new - self.t_bound) > 0:
                t_new = self.t_bound

            h = t_new - t
            h_abs = np.abs(h)

            # initial solution for the Newton solve
            if (self.sol is None):
                Z0 = np.zeros((3, y.shape[0]))
            else:
                Z0 = self.sol(t + h * C).T - y  # extrapolate using previous dense output

            newton_scale = atol + np.abs(y) * rtol
            if self.scale_newton_norm:
                # scale Newton increment for convergence measure (see [3], p. 95)
                newton_scale = newton_scale / (h**self.var_exp)

            converged = False
            while not converged:
                if LU_real is None or LU_complex is None:
                  if self.scale_residuals:
                    # residuals associated with high-index components are scaled
                    # this may help with the stability of the matrix decomposition
                    # (see [3], p97)
                    # self.hscale = np.diag(h**(-self.var_index))
                    if issparse(self.I):
                        self.hscale = diags(h**(-np.minimum(1,self.var_index)), offsets=0, format='csc')
                    else:
                        self.hscale = np.diag(h**(-np.minimum(1,self.var_index))) # only by h or 1
                  try:
                      LU_real    = self.lu( self.hscale @ (MU_REAL    / h * self.mass_matrix - J) )
                      LU_complex = self.lu( self.hscale @ (MU_COMPLEX / h * self.mass_matrix - J) )
                      self.nlu -= 1 # to match original Fortran code
                  except ValueError as e:
                    return False, 'LU decomposition failed ({})'.format(e)

                converged, n_iter, n_bad, Z, f_subs, rate = self.solve_collocation_system(
                    t, y, h, Z0, newton_scale, self.newton_tol,
                    LU_real, LU_complex, residual_scale=self.hscale)

                safety = self.safety_factor * (2 * self.NEWTON_MAXITER + 1) / (2 * self.NEWTON_MAXITER + n_iter)
                  
                if not converged:
                    if current_jac: # we only allow one Jacobian computation per time step
                        break

                    J = self.jac(t, y, f)
                    current_jac = True
                    LU_real = None
                    LU_complex = None

            if not converged:
                h_abs *= 0.5 # retry with lower time step
                LU_real = None # triggers refactorization
                LU_complex = None
                continue

            y_new = y + Z[-1]
            
            ZE = Z.T.dot(E) / h
            err_scale = atol + np.maximum(np.abs(y), np.abs(y_new)) * rtol
            if self.scale_error:
                # correct for the overestimation of the error on
                # algebraic variables, ideally multiply their errors by
                # (h ** index), see [3]
                err_scale = err_scale / (h**self.var_exp) # scale for algebraic variables

            # see [1], chapter IV.8, page 127
            # (akin to a single linearized backward Euler step, damping fast solution modes)
            error = self.solve_lu(LU_real, f + self.mass_matrix.dot(ZE))
            if self.zero_algebraic_error:
                # we exclude the algebraic components, otherwise
                # they artificially lower the error norm
                error[ self.index_algebraic_vars ] = 0.
                error_norm = np.linalg.norm(error / err_scale) / (n - self.nvars_algebraic) ** 0.5
            else:
                error_norm = norm(error / err_scale)
                
            if (rejected and error_norm > 1): # try with stabilised error estimate
                error = self.solve_lu(LU_real, self.fun(t, y + error) + self.mass_matrix.dot(ZE))
                # error is not corrected for algebraic variables
                if self.zero_algebraic_error:
                    # again, we exclude the algebraic components
                    error[ self.index_algebraic_vars ] = 0.
                    error_norm = np.linalg.norm(error / err_scale) / (n - self.nvars_algebraic) ** 0.5
                else:
                    error_norm = norm(error / err_scale)

            if error_norm > 1:
                factor = predict_factor(h_abs, h_abs_old,
                                        error_norm, error_norm_old)
                
                h_abs *= max(self.MIN_FACTOR, safety * factor)
                LU_real = None
                LU_complex = None
            else:
                step_accepted = True

        # Step is converged and accepted
        recompute_jac = jac is not None and n_iter > 2 and rate > 1e-3

        factor = predict_factor(h_abs, h_abs_old, error_norm, error_norm_old)
        factor = min(self.MAX_FACTOR, safety * factor)

        if not recompute_jac and factor < 1.2:
            # maintain previous step size to avoid refactorization cost
            factor = 1
        else:
            LU_real = None
            LU_complex = None

        f_new = self.fun(t_new, y_new)
        if recompute_jac:
            J = jac(t_new, y_new, f_new)
            current_jac = True
        elif jac is not None:
            current_jac = False

        self.h_abs_old = self.h_abs
        self.error_norm_old = error_norm

        self.h_abs = h_abs * factor

        self.y_old = y

        self.t = t_new
        self.y = y_new
        self.f = f_new

        self.Z = Z

        self.LU_real = LU_real
        self.LU_complex = LU_complex
        self.current_jac = current_jac
        self.J = J

        self.t_old = t
        self.sol = self._compute_dense_output()

        return step_accepted, message

    def solve_collocation_system(self, t, y, h, Z0, norm_scale, tol,
                                 LU_real, LU_complex, residual_scale):
        """Solve the collocation system.

        Parameters
        ----------
        t : float
            Current time.
        y : ndarray, shape (n,)
            Current state.
        h : float
            Step to try.
        Z0 : ndarray, shape (3, n)
            Initial guess for the solution. It determines new values of `y` at
            ``t + h * C`` as ``y + Z0``, where ``C`` is the Radau method constants.
        norm_scale : ndarray, shape (n)
            Problem tolerance scale, i.e. ``rtol * abs(y) + atol``.
        tol : float
            Tolerance to which solve the system. This value is compared with
            the normalized by `scale` error.
        LU_real, LU_complex
            LU decompositions of the system Jacobians.

        Returns
        -------
        converged : bool
            Whether iterations converged.
        n_iter : int
            Number of completed iterations.
        Z : ndarray, shape (3, n)
            Found solution.
        rate : float
            The rate of convergence.
        """
        n = y.shape[0]
        M_real = MU_REAL / h
        M_complex = MU_COMPLEX / h

        W = TI.dot(Z0) # state vector at each quadrature point (with complex transformation)
        Z = Z0 # state vector at each quadrature point

        F = np.empty((3, n))  # RHS evaluated at the quadrature points
        ch = h * C # quadrature time points

        dW_norm_old = None
        res_norm = None
        dW = np.empty_like(W)
        converged = False
        rate = None
        nbad_iter = 0

        for k in range(self.NEWTON_MAXITER):
            for i in range(3):
                F[i] = self.fun(t + ch[i], y + Z[i])

            if not np.all(np.isfinite(F)):
                break

            # compute residuals
            f_real    = F.T.dot(TI_REAL)    - M_real    * self.mass_matrix.dot( W[0] )
            f_complex = F.T.dot(TI_COMPLEX) - M_complex * self.mass_matrix.dot( W[1] + 1j * W[2] )

            # scale residuals
            f_real = residual_scale @ f_real
            f_complex = residual_scale @ f_complex

            # compute Newton increment
            dW_real    = self.solve_lu(LU_real,    f_real)
            dW_complex = self.solve_lu(LU_complex, f_complex)

            dW[0] = dW_real
            dW[1] = dW_complex.real
            dW[2] = dW_complex.imag

            dW_norm = norm(dW / norm_scale)

            W += dW
            Z = T.dot(W)

            if dW_norm_old is not None:
                rate = dW_norm / dW_norm_old
                dW_true = rate / (1 - rate) * dW_norm # estimated true error
                dW_true_max_ite = rate ** (self.NEWTON_MAXITER - k) / (1 - rate) * dW_norm # estimated true error after max number of iterations

            if dW_norm < tol:
                converged = True
                break

            if rate is not None:
                if rate >= 1: # Newton loop diverges
                    if rate<100: # divergence is not too extreme yet
                        if nbad_iter<self.NMAX_BAD:
                            # we accept a few number of bad iterations,
                            # which may be necessary for certain higher-index DAEs
                            nbad_iter+=1
                            continue
                if dW_true < tol:
                    converged = True
                    break
                if (dW_true_max_ite > tol) and self.bUsePredictiveNewtonStoppingCriterion :
                    # Newton will most likely not converge in the allowed number of iterations
                    if nbad_iter<self.NMAX_BAD:
                        nbad_iter+=1
                        continue

            dW_norm_old  = dW_norm
            res_norm_old = res_norm

        return converged, k + 1, nbad_iter, Z, F, rate

    def _compute_dense_output(self):
        Q = np.dot(self.Z.T, P)
        return RadauDenseOutput(self.t_old, self.t, self.y_old, Q)

    def _dense_output_impl(self):
        return self.sol


class RadauDenseOutput(DenseOutput):
    def __init__(self, t_old, t, y_old, Q):
        super().__init__(t_old, t)
        self.h = t - t_old
        self.Q = Q
        self.order = Q.shape[1] - 1
        self.y_old = y_old

    def _call_impl(self, t):
        x = (t - self.t_old) / self.h
        if t.ndim == 0:
            p = np.tile(x, self.order + 1)
            p = np.cumprod(p)
        else:
            p = np.tile(x, (self.order + 1, 1))
            p = np.cumprod(p, axis=0)
        # Here we don't multiply by h, not a mistake,
        # because we rely on a dimensionless time.
        y = np.dot(self.Q, p)
        if y.ndim == 2:
            y += self.y_old[:, None]
        else:
            y += self.y_old

        return y
