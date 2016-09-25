from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from .common import select_initial_step, norm, EPS
from .base import OdeSolver, DenseOutput
from scipy.optimize._numdiff import approx_derivative

S6 = 6 ** 0.5

ORDER = 4

# Butcher tableau. A is not used directly, see below.
C = np.array([(4 - S6) / 10, (4 + S6) / 10, 1])
E = np.array([-13 - 7 * S6, -13 + 7 * S6, -1]) / 3

# Eigendecomposition of A is done: A = T L T**-1. There is 1 real eigenvalue
# and a complex conjugate pair. They are written below.
MU_REAL = 3 + 3 ** (2 / 3) - 3 ** (1 / 3)
MU_COMPLEX = 3 + 0.5 * (3 ** (1 / 3) - 3 ** (2 / 3)) - 0.5j * (
    3 ** (5 / 6) + 3 ** (7 / 6))

# These are transformation matrices.
T = np.array([
    [0.09443876248897524, -0.14125529502095421, 0.03002919410514742],
    [0.25021312296533332, 0.20412935229379994, -0.38294211275726192],
    [1, 1, 0]])
TI = np.array([
    [4.17871859155190428, 0.32768282076106237, 0.52337644549944951],
    [-4.17871859155190428, -0.32768282076106237, 0.47662355450055044],
    [0.50287263494578682, -2.57192694985560522, 0.59603920482822492]
])
# This linear combinations are used in the algorithm.
TI_REAL = TI[0]
TI_COMPLEX = TI[1] + 1j * TI[2]

# Interpolator coefficients.
P = np.array([
    [13/3 + 7*S6/3, -23/3 - 22*S6/3, 10/3 + 5 * S6],
    [13/3 - 7*S6/3, -23/3 + 22*S6/3, 10/3 - 5 * S6],
    [1/3, -8/3, 10/3]
])


NEWTON_MAXITER = 7  # Maximum number of Newton iterations.
MIN_FACTOR = 0.2  # Minimum allowed decrease in a step size.
MAX_FACTOR = 10  # Maximum allowed increase in a step size.


def solve_collocation_system(fun, t, y, h, Z0, scale, tol, LU_real,
                             LU_complex):
    """Solve the collocation system.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system.
    t : float
        Current value of the independent variable.
    y : ndarray, shape (n,)
        Current value of the dependent variable.
    h : float
        Step to try.
    J : ndarray, shape (n, n)
        Jacobian of `fun` with respect to `y`.
    Z0 : ndarray, shape (3, n)
        Initial guess for the solution.
    scale : float
        Problem tolerance scale, i.e. ``rtol * abs(y) + atol``.
    tol : float
        Tolerance to which solve the system.
    LU_real, LU_complex
        LU decompositions of the system Jacobian.

    Returns
    -------
    converged : bool
        Whether iterations converged.
    n_iter : int
        Number of completed iterations.
    Z : ndarray, shape (3, n)
        Found solution.
    f_new : ndarray, shape (3, n)
        Value of `fun(x + h, y(x + h))`.
    rate : float
        The rate of convergence.
    LU_real, LU_complex : tuple
        Computed LU decompositions.
    """
    n = y.shape[0]
    M_real = MU_REAL / h
    M_complex = MU_COMPLEX / h

    W = TI.dot(Z0)
    Z = Z0

    F = np.empty((3, n))
    ch = h * C

    dW_norm_old = None
    dW = np.empty_like(W)
    for k in range(NEWTON_MAXITER):
        for i in range(3):
            F[i] = fun(t + ch[i], y + Z[i])

        f_real = F.T.dot(TI_REAL) - M_real * W[0]
        f_complex = F.T.dot(TI_COMPLEX) - M_complex * (W[1] + 1j * W[2])

        dW_real = lu_solve(LU_real, f_real, overwrite_b=True)
        dW_complex = lu_solve(LU_complex, f_complex, overwrite_b=True)

        dW[0] = dW_real
        dW[1] = dW_complex.real
        dW[2] = dW_complex.imag

        dW_norm = norm(dW / scale)
        if dW_norm_old is not None:
            rate = dW_norm / dW_norm_old
        else:
            rate = None

        if (rate is not None and (rate > 1 or
                rate ** (NEWTON_MAXITER - k) / (1 - rate) * dW_norm > tol)):
            converged = False
            break

        W += dW
        Z = T.dot(W)

        if rate is not None and rate / (1 - rate) * dW_norm < tol:
            converged = True
            break

        dW_norm_old = dW_norm
    else:
        converged = False

    return converged, k + 1, Z, F[-1], rate


def predict_factor(h_abs, h_abs_old, error_norm, error_norm_old):
    """Predict by which factor to increase the step size.

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
    with np.errstate(divide='ignore'):
        if error_norm_old is None or h_abs_old is None or error_norm == 0:
            multiplier = 1
        else:
            multiplier = h_abs / h_abs_old * (error_norm_old /
                                              error_norm) ** (1/ORDER)

        factor = min(1, multiplier) * error_norm ** (-1/ORDER)

    return factor


class Radau(OdeSolver):
    """Implicit Runge-Kutta method of Radau IIA family of order 5.

    Implementation follows [1]_. The error is controlled for a 3rd order
    accurate embedded formula. A cubic polynomial which satisfies the
    collocation conditions is used for the dense output.

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
        `y`. The Jacobian matrix has shape (n, n) and its element (i, j) is
        equal to ``d f_i / d y_j``. There are 3 ways to define the Jacobian:

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
        Current status of the solver.
    t_crit : float
        Boundary time.
    direction : -1 or +1
        Integration direction.
    t : float
        Current time.
    y : ndarray, shape (n,)
        Current state.
    step_size : float or None
        Size of the last taken step. None if not steps were made yet.

    References
    ----------
    .. [1] E. Hairer, G. Wanner, "Solving Ordinary Differential Equations II:
           Stiff and Differential-Algebraic Problems", Sec. IV.8.
    """

    def __init__(self, fun, t0, y0, t_crit, rtol=1e-3, atol=1e-6, jac=None):
        super(Radau, self).__init__(fun, t0, y0, t_crit)
        self.t_old = None
        self.y_old = None
        self.rtol = rtol
        self.atol = atol
        self.f = self.fun(self.t, self.y)
        self.h_abs = select_initial_step(
            self.fun, self.t, self.y, self.f, self.direction,
            ORDER, self.rtol, self.atol)
        self.h_abs_old = None
        self.error_norm_old = None

        self.newton_tol = max(10 * EPS / rtol, min(0.03, rtol ** 0.5))
        self.sol = None

        self.jac, self.J = self._validate_jac(jac)
        self.current_jac = True
        self.LU_real = None
        self.LU_complex = None
        self.Z = None

    def _validate_jac(self, jac):
        fun = self.fun
        t0 = self.t
        y0 = self.y

        if jac is None:
            def jac_wrapped(t, y):
                return approx_derivative(lambda z: fun(t, z),
                                         y, method='2-point')
            J = jac_wrapped(t0, y0)
        elif callable(jac):
            def jac_wrapped(t, y):
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

    def _step_impl(self, max_step):
        t = self.t
        y = self.y
        f = self.f

        atol = self.atol
        rtol = self.rtol
        if self.h_abs > max_step:
            h_abs = max_step
            h_abs_old = None
            error_norm_old = None
        else:
            h_abs = self.h_abs
            h_abs_old = self.h_abs_old
            error_norm_old = self.error_norm_old

        J = self.J
        LU_real = self.LU_real
        LU_complex = self.LU_complex

        current_jac = self.current_jac
        jac = self.jac

        I = np.identity(self.n)

        rejected = False
        step_accepted = False
        message = None
        while not step_accepted:
            h = h_abs * self.direction
            t_new = t + h

            if self.direction * (t_new - self.t_crit) > 0:
                t_new = self.t_crit

            if t_new == t:  # h is less than spacing between numbers.
                return False, self.TOO_SMALL_STEP

            h = t_new - t
            h_abs = np.abs(h)

            if self.sol is None:
                Z0 = np.zeros((3, y.shape[0]))
            else:
                Z0 = self.sol(t + h * C).T - y

            scale = atol + np.abs(y) * rtol

            if LU_real is None or LU_complex is None:
                LU_real = lu_factor(MU_REAL / h * I - J, overwrite_a=True)
                LU_complex = lu_factor(MU_COMPLEX / h * I - J,
                                       overwrite_a=True)

            converged, n_iter, Z, f_new, rate = \
                solve_collocation_system(
                    self.fun, t, y, h, Z0, scale, self.newton_tol,
                    LU_real, LU_complex)

            if not converged:
                if not current_jac:
                    J = self.jac(t, y)
                    current_jac = True
                else:
                    h_abs *= 0.5
                LU_real = None
                LU_complex = None
                continue

            y_new = y + Z[-1]
            ZE = Z.T.dot(E) / h
            error = lu_solve(LU_real, f + ZE, overwrite_b=True)
            scale = atol + np.maximum(np.abs(y), np.abs(y_new)) * rtol
            error_norm = norm(error / scale)
            safety = 0.9 * (2 * NEWTON_MAXITER + 1) / (
                2 * NEWTON_MAXITER + n_iter)

            if rejected and error_norm > 1:
                error = lu_solve(LU_real, self.fun(t, y + error) + ZE,
                                 overwrite_b=True)
                error_norm = norm(error / scale)

            if error_norm > 1:
                factor = predict_factor(h_abs, h_abs_old,
                                        error_norm, error_norm_old)
                h_abs *= max(MIN_FACTOR, safety * factor)

                LU_real = None
                LU_complex = None
                rejected = True
            else:
                step_accepted = True

        recompute_jac = jac is not None and n_iter > 2 and rate > 1e-3

        factor = predict_factor(h_abs, h_abs_old, error_norm, error_norm_old)
        factor = min(MAX_FACTOR, safety * factor)

        if not recompute_jac and factor < 1.2:
            factor = 1
        else:
            LU_real = None
            LU_complex = None

        if recompute_jac:
            J = jac(t_new, y_new)
            current_jac = True
        elif jac is not None:
            current_jac = False

        self.h_abs_old = self.h_abs
        self.error_norm_old = error_norm

        h_abs *= factor
        self.h_abs = h_abs

        self.t_old = t
        self.y_old = y

        self.t = t_new
        self.y = y_new
        self.step_size = np.abs(t_new - t)

        self.f = f_new
        self.Z = Z

        self.LU_real = LU_real
        self.LU_complex = LU_complex
        self.current_jac = current_jac
        self.J = J

        self.sol = self._compute_dense_output()

        return step_accepted, message

    def _compute_dense_output(self):
        Q = np.dot(self.Z.T, P)
        return RadauDenseOutput(self.t_old, self.t, self.y_old, Q)

    def _dense_output_impl(self):
        return self.sol


class RadauDenseOutput(DenseOutput):
    def __init__(self, t_prev, t, y_old, Q):
        super(RadauDenseOutput, self).__init__(t_prev, t)
        self.h = t - t_prev
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
        # Here we don't multiply by h, not a mistake.
        y = np.dot(self.Q, p)
        if y.ndim == 2:
            y += self.y_old[:, None]
        else:
            y += self.y_old

        return y
