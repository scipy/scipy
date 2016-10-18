from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from .common import (validate_max_step, validate_tol, select_initial_step,
                     norm, EPS, num_jac, warn_extraneous)
from .base import OdeSolver, DenseOutput


MAX_ORDER = 5
NEWTON_MAXITER = 4
MIN_FACTOR = 0.2
MAX_FACTOR = 10


def compute_R(order, factor):
    I = np.arange(1, order + 1)[:, None]
    J = np.arange(1, order + 1)
    M = np.zeros((order + 1, order + 1))
    M[1:, 1:] = (I - 1 - factor * J) / I
    M[0] = 1
    return np.cumprod(M, axis=0)


def change_D(D, order, factor):
    R = compute_R(order, factor)
    U = compute_R(order, 1)
    RU = R.dot(U)
    D[:order + 1] = np.dot(RU.T, D[:order + 1])


def solve_corrector_system(fun, t_new, y_predict, c, psi, LU, scale, tol):
    d = 0
    y = y_predict.copy()
    dy_norm_old = None
    converged = False
    for k in range(NEWTON_MAXITER):
        f = fun(t_new, y)
        if not np.all(np.isfinite(f)):
            break

        dy = lu_solve(LU, c * f - psi - d, overwrite_b=True)
        dy_norm = norm(dy / scale)

        if dy_norm_old is None:
            rate = None
        else:
            rate = dy_norm / dy_norm_old

        if (rate is not None and (rate >= 1 or
                rate ** (NEWTON_MAXITER - k) / (1 - rate) * dy_norm > tol)):
            break

        y += dy
        d += dy

        if (dy_norm == 0 or
                rate is not None and rate / (1 - rate) * dy_norm < tol):
            converged = True
            break

        dy_norm_old = dy_norm

    return converged, k + 1, y, d


class BDF(OdeSolver):
    """Implicit method based on Backward Differentiation Formulas.

    This is a variable order method with the order varying automatically from
    1 to 5. The general framework of the BDF algorithm is described in [1]_.
    This class implements a quasi-constant step size approach as explained
    in [2]_. The error estimation strategy for the constant step BDF is derived
    in [3]_. An accuracy enhancement using modified formulas (NDF) [2]_ is also
    implemented.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system. The calling signature is ``fun(t, y)``.
        Here ``t`` is a scalar, and ``y`` is ndarray with shape (n,). It
        must return an array_like with shape (n,).
    t0 : float
        Initial time.
    y0 : array_like, shape (n,)
        Initial state.
    t_crit : float
        Boundary time --- the integration won't continue beyond it. It also
        determines the direction of the integration.
    max_step : float, optional
        Maximum allowed step size. Default is np.inf, i.e. the step is not
        bounded and determined solely by the solver.
    rtol, atol : float and array_like, optional
        Relative and absolute tolerances. The solver keeps the local error
        estimates less than ``atol + rtol * abs(y)``. Here `rtol` controls a
        relative accuracy (number of correct digits). But if a component of `y`
        is approximately below `atol` then the error only needs to fall within
        the same `atol` threshold, and the number of correct digits is not
        guaranteed. If components of y have different scales, it might be
        beneficial to set different `atol` values for different components by
        passing array_like with shape (n,) for `atol`. Default values are
        1e-3 for `rtol` and 1e-6 for `atol`.
    jac : array_like, callable or None, optional
        Jacobian matrix of the right-hand side of the system with respect to
        y. The Jacobian matrix has shape (n, n) and its element (i, j) is
        equal to d f_i / d y_j. There are 3 ways to define the Jacobian:

            * If array_like, then the Jacobian is assumed to be constant.
            * If callable, then the Jacobian is assumed to depend on both
              t and y, and will be called as ``jac(t, y)`` as necessary.
            * If None (default), then the Jacobian will be approximated by
              finite differences.

        It is generally recommended to provided the Jacobian rather than
        relying on finite difference approximation.

    Attributes
    ----------
    n : int
        Number of equations.
    status : string
        Current status of the solver: 'running', 'finished' or 'failed'.
    t_crit : float
        Boundary time.
    direction : float
        Integration direction: +1 or -1.
    t : float
        Current time.
    y : ndarray
        Current state.
    t_old : float
        End time of the last successful step. None if no steps were made yet.
    step_size : float
        Size of the last successful step. None if no steps were made yet.
    nfev : int
        Number of the system rhs evaluations.
    njev : int
        Number of the Jacobian evaluations.
    nlu : int
        Number of LU decompositions of the Jacobian.

    References
    ----------
    .. [1] G. D. Byrne, A. C. Hindmarsh, "A Polyalgorithm for the Numerical
           Solution of Ordinary Differential Equations", ACM Transactions on
           Mathematical Software, Vol. 1, No. 1, pp. 71-96, March 1975.
    .. [2] L. F. Shampine, M. W. Reichelt, "THE MATLAB ODE SUITE", SIAM J. SCI.
           COMPUTE., Vol. 18, No. 1, pp. 1-22, January 1997.
    .. [3] E. Hairer, G. Wanner, "Solving Ordinary Differential Equations I:
           Nonstiff Problems", Sec. III.2.
    """
    def __init__(self, fun, t0, y0, t_crit, max_step=np.inf,
                 rtol=1e-3, atol=1e-6, jac=None, **extraneous):
        warn_extraneous(extraneous)
        super(BDF, self).__init__(fun, t0, y0, t_crit)
        self.max_step = validate_max_step(max_step)
        self.rtol, self.atol = validate_tol(rtol, atol, self.n)
        f = self.fun(self.t, self.y)
        self.h_abs = select_initial_step(
            self.fun, self.t, self.y, f, self.direction,
            1, self.rtol, self.atol)
        self.h_abs_old = None
        self.error_norm_old = None

        self.newton_tol = max(10 * EPS / rtol, min(0.03, rtol ** 0.5))

        self.jac_factor = None
        self.jac, self.J = self._validate_jac(jac)

        kappa = np.array([0, -0.1850, -1/9, -0.0823, -0.0415, 0])
        self.gamma = np.hstack((0, np.cumsum(1 / np.arange(1, MAX_ORDER + 1))))
        self.alpha = (1 - kappa) * self.gamma
        self.error_const = kappa * self.gamma + 1 / np.arange(1, MAX_ORDER + 2)

        D = np.empty((MAX_ORDER + 3, self.n))
        D[0] = self.y
        D[1] = f * self.h_abs * self.direction
        self.D = D

        self.order = 1
        self.n_equal_steps = 0
        self.LU = None

    def _lu(self, *args, **kwargs):
        self.nlu += 1
        return lu_factor(*args, **kwargs)

    def _validate_jac(self, jac):
        fun = self._fun
        t0 = self.t
        y0 = self.y

        if jac is None:
            def jac_wrapped(t, y):
                self.njev += 1
                J, self.jac_factor = num_jac(fun, t, y, fun(t, y),
                                             self.atol, self.jac_factor)
                return J
            J = jac_wrapped(t0, y0)
        elif callable(jac):
            def jac_wrapped(t, y):
                self.njev += 1
                return np.asarray(jac(t, y), dtype=float)
            J = jac_wrapped(t0, y0)
            if J.shape != (self.n, self.n):
                raise ValueError(
                    "`jac` return is expected to have shape {}, but actually "
                    "has {}.".format((self.n, self.n), J.shape))
        else:
            J = np.asarray(jac, dtype=float)
            if J.shape != (self.n, self.n):
                raise ValueError("`jac` is expected to have shape {}, but "
                                 "actually has {}."
                                 .format((self.n, self.n), J.shape))
            jac_wrapped = None

        return jac_wrapped, J

    def _step_impl(self):
        t = self.t
        y = self.y
        D = self.D

        max_step = self.max_step
        min_step = 10 * np.abs(np.nextafter(t, self.direction * np.inf) - t)
        if self.h_abs > max_step:
            h_abs = max_step
            factor = max_step / self.h_abs
            change_D(D, self.order, factor)
            self.n_equal_steps = 0
        elif self.h_abs < min_step:
            h_abs = min_step
            factor = min_step / self.h_abs
            change_D(D, self.order, factor)
            self.n_equal_steps = 0
        else:
            h_abs = self.h_abs

        atol = self.atol
        rtol = self.rtol
        order = self.order

        alpha = self.alpha
        gamma = self.gamma
        error_const = self.error_const

        I = np.identity(self.n)
        J = self.J
        LU = self.LU
        current_jac = self.jac is None

        step_accepted = False
        while not step_accepted:
            if h_abs < min_step:
                return False, self.TOO_SMALL_STEP

            h = h_abs * self.direction
            t_new = t + h

            if self.direction * (t_new - self.t_crit) > 0:
                t_new = self.t_crit
                factor = np.abs(t_new - t) / h_abs
                change_D(D, order, factor)
                self.n_equal_steps = 0
                LU = None

            h = t_new - t
            h_abs = np.abs(h)

            t_new = t + h
            y_predict = np.sum(D[:order + 1], axis=0)

            scale = atol + rtol * np.abs(y_predict)
            psi = np.dot(D[1: order + 1].T, gamma[1: order + 1]) / alpha[order]

            converged = False
            c = h / alpha[order]
            while not converged:
                if LU is None:
                    LU = self._lu(I - c * J, overwrite_a=True)

                converged, n_iter, y, d = solve_corrector_system(
                    self.fun, t_new, y_predict, c, psi, LU, scale,
                    self.newton_tol)

                if not converged:
                    if not current_jac:
                        J = self.jac(t_new, y_predict)
                        LU = None
                        current_jac = True
                    else:
                        break

            if not converged:
                factor = 0.5
                h_abs *= factor
                change_D(D, order, factor)
                self.n_equal_steps = 0
                LU = None
                continue

            safety = 0.9 * (2 * NEWTON_MAXITER + 1) / (2 * NEWTON_MAXITER +
                                                       n_iter)

            scale = atol + rtol * np.abs(y)
            error = error_const[order] * d
            error_norm = norm(error / scale)

            if error_norm > 1:
                factor = max(MIN_FACTOR,
                             safety * error_norm ** (-1 / (order + 1)))
                h_abs *= factor
                change_D(D, order, factor)
                self.n_equal_steps = 0
                # As we didn't have problems with convergence, we don't
                # reset LU here.
            else:
                step_accepted = True

        self.n_equal_steps += 1

        self.t = t_new
        self.y = y

        self.h_abs = h_abs
        self.J = J
        self.LU = LU

        # Update divided differences. The principal relation here is
        # D^{j + 1} y_n = D^{j} y_n - D^{j} y_{n - 1}. Keep in mind that D
        # contained divided difference for previous interpolating
        # polynomial and d = D^{k + 1} y_n. Thus this elegant code follows.
        D[order + 2] = d - D[order + 1]
        D[order + 1] = d
        for i in reversed(range(order + 1)):
            D[i] += D[i + 1]

        if self.n_equal_steps < order + 1:
            return True, None

        if order > 1:
            error_m = error_const[order - 1] * D[order]
            error_m_norm = norm(error_m / scale)
        else:
            error_m_norm = np.inf

        if order < MAX_ORDER:
            error_p = error_const[order + 1] * D[order + 2]
            error_p_norm = norm(error_p / scale)
        else:
            error_p_norm = np.inf

        error_norms = np.array([error_m_norm, error_norm, error_p_norm])
        factors = error_norms ** (-1 / np.arange(order, order + 3))

        delta_order = np.argmax(factors) - 1
        order += delta_order
        self.order = order

        factor = min(MAX_FACTOR, safety * np.max(factors))
        self.h_abs *= factor
        change_D(D, order, factor)
        self.n_equal_steps = 0
        self.LU = None

        return True, None

    def _dense_output_impl(self):
        return BdfDenseOutput(self.t_old, self.t, self.h_abs * self.direction,
                              self.order, self.D[:self.order + 1].copy())


class BdfDenseOutput(DenseOutput):
    def __init__(self, t_old, t, h, order, D):
        super(BdfDenseOutput, self).__init__(t_old, t)
        self.order = order
        self.t_shift = self.t - h * np.arange(self.order)
        self.denom = h * (1 + np.arange(self.order))
        self.D = D

    def _call_impl(self, t):
        if t.ndim == 0:
            x = (t - self.t_shift) / self.denom
            p = np.cumprod(x)
        else:
            x = (t - self.t_shift[:, None]) / self.denom[:, None]
            p = np.cumprod(x, axis=0)

        y = np.dot(self.D[1:].T, p)
        if y.ndim == 1:
            y += self.D[0]
        else:
            y += self.D[0, :, None]

        return y
