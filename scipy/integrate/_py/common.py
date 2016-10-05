from __future__ import division, print_function, absolute_import
import warnings as warn
import numpy as np


EPS = np.finfo(float).eps


def validate_tol(rtol, atol, n):
    """Validate tolerance values."""
    if rtol < 100 * EPS:
        warn("`rtol` is too low, setting to {}".format(100 * EPS))
        rtol = 100 * EPS

    atol = np.asarray(atol)
    if atol.ndim > 0 and atol.shape != (n,):
        raise ValueError("`atol` has wrong shape.")

    if np.any(atol < 0):
        raise ValueError("`atol` must be positive.")

    return rtol, atol


def norm(x):
    """Compute RMS norm."""
    return np.linalg.norm(x) / x.size ** 0.5


def select_initial_step(fun, t0, y0, f0, direction, order, rtol, atol):
    """Empirically select a good initial step.

    The algorithm is described in [1]_.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system.
    t0 : float
        Initial value of the independent variable.
    y0 : ndarray, shape (n,)
        Initial value of the dependent variable.
    f0 : ndarray, shape (n,)
        Initial value of the derivative, i. e. ``fun(t0, y0)``.
    direction : float
        Integration direction.
    order : float
        Method order.
    rtol : float
        Desired relative tolerance.
    atol : float
        Desired absolute tolerance.

    Returns
    -------
    h_abs : float
        Absolute value of the suggested initial step.

    References
    ----------
    .. [1] E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential
           Equations I: Nonstiff Problems", Sec. II.4.
    """
    if y0.size == 0:
        return np.inf

    scale = atol + np.abs(y0) * rtol
    d0 = norm(y0 / scale)
    d1 = norm(f0 / scale)
    if d0 < 1e-5 or d1 < 1e-5:
        h0 = 1e-6
    else:
        h0 = 0.01 * d0 / d1

    y1 = y0 + h0 * direction * f0
    f1 = fun(t0 + h0 * direction, y1)
    d2 = norm((f1 - f0) / scale) / h0

    if d1 <= 1e-15 and d2 <= 1e-15:
        h1 = max(1e-6, h0 * 1e-3)
    else:
        h1 = (0.01 / max(d1, d2)) ** (1 / (order + 1))

    return min(100 * h0, h1)


def _find_group_indices(a):
    c = None
    indices = []
    values = []
    for i, e in enumerate(a):
        if e != c:
            indices.append(i)
            values.append(e)
            c = e
    indices.append(a.shape[0])

    return indices, values


class OdeSolution(object):
    """Continuous ODE solution.

    It is organized as a collection of `DenseOutput` objects which represent
    local interpolants. It provides an algorithm to select a right interpolant
    for each given point.

    The interpolants cover the range between `t_min` and `t_max` (see
    Attributes below). Evaluation outside this interval is not forbidden, but
    the accuracy is not guaranteed.

    Attributes
    ----------
    t_min, t_max : float
        Time range of the interpolation.
    """
    def __init__(self, ts, interpolants):
        d = np.diff(ts)
        if not (np.all(d >= 0) or np.all(d <= 0)):
            raise ValueError("`ts` must be strictly increasing or decreasing.")

        ts = np.asarray(ts)
        self.n_segments = len(interpolants)
        if ts.shape != (self.n_segments + 1,):
            raise ValueError("Numbers of time stamps and interpolants "
                             "don't match.")

        self.ts = ts
        self.interpolants = interpolants
        if ts[-1] >= ts[0]:
            self.t_min = ts[0]
            self.t_max = ts[-1]
            self.ascending = True
            self.ts_sorted = ts
        else:
            self.t_min = ts[-1]
            self.t_max = ts[0]
            self.ascending = False
            self.ts_sorted = np.sort(ts)

    def _call_single(self, t):
        i = min(max(np.searchsorted(self.ts_sorted, t) - 1, 0),
                self.n_segments - 1)
        if not self.ascending:
            i = self.n_segments - 1 - i

        return self.interpolants[i](t)

    def __call__(self, t):
        """Evaluate the solution.

        Parameters
        ----------
        t : float or array_like with shape (n_points,)
            Points to evaluate at.

        Returns
        -------
        y : ndarray, shape (n_states,) or (n_states, n_points)
            Computed values. Shape depends on whether `t` was a scalar or a
            1-d array.
        """
        t = np.asarray(t)

        if t.ndim == 0:
            return self._call_single(t)

        order = np.argsort(t)
        reverse = np.empty_like(order)
        reverse[order] = np.arange(order.shape[0])

        t_sorted = t[order]

        index = np.searchsorted(self.ts_sorted, t_sorted)
        index -= 1
        index[index < 0] = 0
        index[index > self.n_segments - 1] = self.n_segments - 1

        if not self.ascending:
            index = self.n_segments - 1 - index

        t_index, sol_index = _find_group_indices(index)

        ys = []
        for s, i, j in zip(sol_index, t_index[:-1], t_index[1:]):
            y = self.interpolants[s](t_sorted[i: j])
            ys.append(y)

        ys = np.hstack(ys)
        ys = ys[:, reverse]

        return ys


def num_jac(fun, t, y, f, threshold, factor=None):
    """Finite differences Jacobian approximation tailored for ODE solvers.

    This function computes finite difference approximation to the Jacobian
    matrix of `fun` with respect to `y` using forward differences.
    The Jacobian matrix has shape (n, n) and its element (i, j) is equal to
    ``d f_i / d y_j``.

    A special feature of this function is the ability to correct the step
    size from iteration to iteration. The main idea is to keep the finite
    difference significantly separated from its round-off error which
    approximately equals ``EPS * np.abs(f)``. It reduces a possibility of a
    huge error and assures that the estimated derivative are reasonably close
    to the true values (i.e. the finite difference approximation is at least
    qualitatively reflects the structure of the true Jacobian).

    Parameters
    ----------
    fun : callable
        Right-hand side of the system.
    t : float
        Current time.
    y : ndarray, shape (n,)
        Current state.
    f : ndarray, shape (n,)
        Value of the right hand side at (t, y).
    threshold : float
        Threshold for `y` value used for computing the step size as
        ``factor * np.maximum(np.abs(y), threshold)``. Typically the value of
        absolute tolerance (atol) for a solver should be passed as `threshold`.
    factor : ndarray with shape (n,) or None
        Factor to use for computing the step size. Pass None for the very
        evaluation, then use the value returned from this function.

    Returns
    -------
    J : ndarray, shape (n, n)
        Jacobian matrix.
    step_size : ndarray, shape (n,)
        Suggested step size for the next evaluation.
    """
    y = np.asarray(y)
    n = y.shape[0]

    if factor is None:
        factor = np.ones(n) * EPS ** 0.5

    # Direct the step as ODE dictates, hoping that such step won't lead to
    # a problematic region.
    f_sign = 2 * (f >= 0).astype(float) - 1
    y_scale = f_sign * np.maximum(threshold, np.abs(y))
    h = factor * y_scale
    h = (y + h) - y

    # Make sure that the step is not 0 to start with. Not likely it will be
    # executed often.
    for i in np.nonzero(h == 0)[0]:
        while h[i] == 0:
            factor[i] *= 10
            h[i] = (y[i] + factor[i] * y_scale[i]) - y[i]

    h_vecs = np.diag(h)
    J = np.empty((n, n))
    for i in range(n):
        f_new = fun(t, y + h_vecs[i])
        diff = f_new - f
        max_ind = np.argmax(np.abs(diff))

        max_diff = np.abs(diff[max_ind])
        scale = max(np.abs(f_new[max_ind]), np.abs(f[max_ind]))

        # Difference is suspiciously small. Try a 10x step.
        if max_diff < EPS ** 0.875 * scale:
            new_factor = 10 * factor[i]
            h_new = (y[i] + new_factor * y_scale[i]) - y[i]
            h_vecs[i, i] = h_new
            f_new = fun(t, y + h_vecs[i])
            diff_new = f_new - f
            max_ind = np.argmax(np.abs(diff_new))
            max_diff_new = np.abs(diff_new[max_ind])
            scale_new = max(np.abs(f_new[max_ind]), np.abs(f[max_ind]))

            if max_diff * scale_new < scale * max_diff_new:
                factor[i] = new_factor
                h[i] = h_new
                diff = diff_new
                scale = scale_new
                max_diff = max_diff_new

        J[:, i] = diff / h[i]

        # Adjust factor based on the difference compared to the scale.
        if max_diff < EPS ** 0.75 * scale:
            factor[i] *= 10
        elif max_diff > EPS ** 0.25 * scale:
            factor[i] *= 0.1

    factor = np.maximum(factor, 1e3 * EPS)
    return J, factor
