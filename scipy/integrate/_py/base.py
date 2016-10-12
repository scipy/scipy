from __future__ import division, print_function, absolute_import
import numpy as np


def check_arguments(fun, y0):
    y0 = np.asarray(y0, dtype=float)

    if y0.ndim != 1:
        raise ValueError("`y0` must be 1-dimensional.")

    def fun_wrapped(t, y):
        return np.asarray(fun(t, y))

    return fun_wrapped, y0


class OdeSolver(object):
    """Base class for ODE solvers.

    In order to implement a new solver you need to follow the guidelines:

        1. A constructor must accept parameters presented in the base class
           (listed below) followed by any other parameters specific to a
           solver.
        2. A solver must implement a private method
           `_step_impl(self, max_step=np.inf)` which propagates a solver one
           step further. It must return tuple ``(success, message)``, where
           ``success`` is a boolean indicating whether a step was successful,
           and ``message`` is a string containing description of a failure if
           a step failed or None otherwise.
        3. A solver must implement a private method `_dense_output_impl(self)`
           which returns `DenseOutput` object covering the last successful
           step.
        4. A solver must have attributes listed below in Attributes section.
        5. Use `fun(self, t, y)` method for the system rhs evaluation, this
           way the number of functions evaluations (`nfev`) will be tracked
           automatically.
        6. If a solver uses Jacobian and LU decompositions, it should track
           the number of Jacobian evaluations (`njev`) and the number of LU
           factorizations (`nlu`).
        7. By convention a function evaluations used to compute a finite
           difference approximation of the Jacobian should not be counted in
           `nfev`, thus use `_fun(self, t, y)` method when computing a finite
           difference approximation of the Jacobian.

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
    """
    TOO_SMALL_STEP = "Required step size is less than spacing between numbers."

    def __init__(self, fun, t0, y0, t_crit):
        self.t_old = None
        self.t = t0
        self._fun, self.y = check_arguments(fun, y0)
        self.t_crit = t_crit

        self.direction = np.sign(t_crit - t0) if t_crit != t0 else 1
        self.n = self.y.size
        self.status = 'running'
        self.step_size = None

        self.nfev = 0
        self.njev = 0
        self.nlu = 0

    def fun(self, t, y):
        self.nfev += 1
        return self._fun(t, y)

    def step(self, max_step=np.inf):
        """Perform one integration step.

        Parameters
        ----------
        max_step : float, optional
            Maximum allowed step size. Default is np.inf, i.e. the step is not
            bounded and determined solely by the solver.

        Returns
        -------
        message : string or None
            Report from the solver. Typically a reason for the failure if
            `self.status` is 'failed' after the step was taken or None
            otherwise.
        """
        if max_step <= 0:
            raise ValueError("`max_step` must be positive.")

        if self.status != 'running':
            raise RuntimeError("Attempt to step on a failed or finished "
                               "solver.")

        if self.n == 0 or self.t == self.t_crit:
            # Handle corner cases of empty solver and no integration
            t = self.t
            t_max = self.t + self.direction * max_step
            t_new = min(t_max, self.t_crit) if self.direction == 1 else max(t_max, self.t_crit)
            self.t_old = t
            self.t = t_new
            self.step_size = np.abs(t_new - t)

            message = None

            if not np.isfinite(self.t) or self.direction * (self.t - self.t_crit) >= 0:
                self.status = 'finished'
        else:
            success, message = self._step_impl(max_step)

            if not success:
                self.status = 'failed'
            elif self.direction * (self.t - self.t_crit) >= 0:
                self.status = 'finished'

        return message

    def dense_output(self):
        """Compute a local interpolant over the last successful step.

        Returns
        -------
        sol : `DenseOutput`
            Local interpolant over the last successful step.
        """
        if self.t_old is None:
            raise RuntimeError("Dense output is available after a successful "
                               "step was made.")

        if self.n == 0 or self.t == self.t_old:
            # Handle corner cases of empty solver and no integration
            return ConstantDenseOutput(self.y, self.t_old, self.t)
        else:
            return self._dense_output_impl()

    def _step_impl(self, max_step):
        raise NotImplementedError

    def _dense_output_impl(self):
        raise NotImplementedError


class DenseOutput(object):
    """Local interpolant over the last successful step made by an ODE solver.

    It interpolates between `t_min` and `t_max` (see Attributes below).
    Evaluation outside this interval is not forbidden, but the accuracy is not
    guaranteed.

    Attributes
    ----------
    t_min, t_max : float
        Time range of the interpolation.
    """
    def __init__(self, t_old, t):
        self.t_old = t_old
        self.t = t
        self.t_min = min(t, t_old)
        self.t_max = max(t, t_old)

    def __call__(self, t):
        """Evaluate the interpolant.

        Parameters
        ----------
        t : float or array_like with shape (n_points,)
            Points to evaluate the solution at.

        Returns
        -------
        y : ndarray, shape (n,) or (n, n_points)
            Computed values. Shape depends on whether `t` was a scalar or a
            1-d array.
        """
        t = np.asarray(t)
        if t.ndim > 1:
            raise ValueError("`t` must be float or 1-d array.")
        return self._call_impl(t)

    def _call_impl(self, t):
        raise NotImplementedError


class ConstantDenseOutput(DenseOutput):
    def __init__(self, value, t_old, t):
        super(ConstantDenseOutput, self).__init__(t_old, t)
        self.value = value

    def _call_impl(self, t):
        return self.value
