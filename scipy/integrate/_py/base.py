from __future__ import division, print_function, absolute_import
import numpy as np


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
        Size of the last successful step. None if no steps were made yet.
    """
    TOO_SMALL_STEP = "Required step size is less than spacing between numbers."

    def __init__(self, fun, t0, y0, t_crit):
        def fun_wrapped(t, y):
            return np.atleast_1d(fun(t, y))

        self.fun = fun_wrapped
        self.t = t0
        self.y = np.asarray(y0)
        self.t_crit = t_crit

        if self.t == self.t_crit:
            raise ValueError("Initial `t` coincides with `t_crit`.")

        self.direction = np.sign(t_crit - t0)
        self.n = self.y.shape[0]
        self.status = 'started'
        self.step_size = None

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
            Report from the solver. Typically a reason for failure if
            `solve.status` is 'failed' after the step was taken or None
            otherwise.
        """
        if max_step <= 0:
            raise ValueError("`max_step` must be positive.")

        if self.status in ['failed', 'finished']:
            raise RuntimeError("Attempt to step on a failed or finished "
                               "solver.")

        success, message = self._step_impl(max_step)

        if not success:
            self.status = 'failed'
        elif self.direction * (self.t - self.t_crit) >= 0:
            self.status = 'finished'
        else:
            self.status = 'running'

        return message

    def dense_output(self):
        """Compute a local interpolant over the last step.

        Returns
        -------
        sol : `DenseOutput`
            Local interpolant over the last step.
        """
        if self.status not in ['running', 'finished']:
            raise RuntimeError("Dense output is available after a successful "
                               "step was taken.")

        return self._dense_output_impl()

    def _step_impl(self, max_step):
        raise NotImplementedError

    def _dense_output_impl(self):
        raise NotImplementedError


class DenseOutput(object):
    """Local interpolator over the last step taken by an ODE solver.

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
        y : ndarray, shape (n_states,) or (n_states, n_points)
            Computed values. Shape depends on whether `t` was a scalar or a
            1-d array.
        """
        t = np.asarray(t)
        if t.ndim > 1:
            raise ValueError("`t` must be float or 1-d array.")
        return self._call_impl(t)

    def _call_impl(self, t):
        raise NotImplementedError
