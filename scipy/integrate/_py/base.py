from __future__ import division, print_function, absolute_import
import numpy as np


STATUS = ['started', 'running', 'failed', 'finished']


class OdeSolver(object):
    """Base class for ODE solvers."""
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
        status : {'running', 'failed', 'finished'}
            Solver status after the step was taken. Further steps are possible
            only if ``status='running'``.
        message : string or None
            Reason for failure if `status`  is 'failed', None otherwise.
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

        return self.status, message

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
