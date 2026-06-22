"""Least-squares multivariate spline approximation."""

from __future__ import annotations

import warnings

import numpy as np

from ._lsq_multivariate_basis import (
    _basis_counts_from_knots,
    _construct_full_knots,
    _design_matrix,
)
from ._lsq_multivariate_grid import _grid_to_scattered
from ._lsq_multivariate_solver import _solve_lsq
from ._lsq_multivariate_validation import (
    _validate_derivative_orders,
    _validate_eval_domain,
    _validate_eval_input,
    _validate_input,
)


__all__ = ["LSQMultivariateSpline"]


class LSQMultivariateSpline:
    """Weighted least-squares multivariate spline approximation.

    .. versionadded:: 1.19.0

    Parameters
    ----------
    x : array_like, shape (n_samples, n_dimensions)
        Coordinates of the input data. Each row is one sample point; each
        column is one input dimension.
    y : array_like, shape (n_samples,) or (n_samples, n_outputs)
        Observed values at the sample coordinates.
    t : int or sequence of array_like or int
        Interior knot coordinates, or the number of polynomial pieces to use.
        A scalar integer applies the same automatic knot count to every input
        dimension.
    w : array_like, shape (n_samples,), optional
        Positive weights for weighted least-squares fitting.
    bbox : array_like, shape (n_dimensions, 2), optional
        Boundary of the approximation domain for each input dimension.
    k : int or sequence of int, optional
        Spline degree for each input dimension. A scalar applies the same
        degree to every dimension. Default is cubic, ``k=3``.
    eps : float, optional
        Rank threshold used by the least-squares solver.
    check_finite : bool, optional
        Whether to check input arrays for NaN or infinite values.
    sparse : bool, optional
        Whether to build and solve the design matrix using sparse matrices.
    smoothing : float, optional
        Nonnegative roughness penalty strength. The default, ``0.0``, gives
        an ordinary least-squares spline.
    penalty_order : int, optional
        Finite-difference order used for the coefficient roughness penalty.
        The default is ``2``.
    unsupported : {"raise", "warn", "ignore"}, optional
        How to handle tensor-product basis functions with no data support.
        The default, ``"raise"``, rejects the fit. ``"warn"`` emits a warning
        and continues with a least-squares solution. ``"ignore"`` continues
        silently.

    Notes
    -----
    The fitting model is a tensor-product B-spline basis. In two dimensions
    this is analogous to bivariate spline basis terms ``B_i(x0) * B_j(x1)``.
    In ``n`` dimensions the basis terms become products over all input axes.

    Examples
    --------
    Fit a one-dimensional cubic least-squares spline with automatic knots:

    >>> import numpy as np
    >>> from scipy.interpolate import LSQMultivariateSpline
    >>> x = np.linspace(0.0, 1.0, 20)
    >>> y = np.sin(6.0 * x)
    >>> spl = LSQMultivariateSpline(x, y, t=5, k=3)
    >>> spl([0.25, 0.5]).shape
    (2,)

    Fit gridded two-dimensional data:

    >>> x0 = np.linspace(-1.0, 1.0, 6)
    >>> x1 = np.linspace(-1.0, 1.0, 7)
    >>> x0_grid, x1_grid = np.meshgrid(x0, x1, indexing="ij")
    >>> z = x0_grid - 2.0 * x1_grid
    >>> spl = LSQMultivariateSpline.from_grid((x0, x1), z, t=[1, 1], k=1)
    >>> np.round(spl([[0.5, -0.5]]), 12)
    array([1.5])
    """

    def __init__(
        self,
        x,
        y,
        t,
        w=None,
        bbox=None,
        k=3,
        eps=None,
        check_finite=False,
        sparse=False,
        smoothing=0.0,
        penalty_order=2,
        unsupported="raise",
    ):
        self._input = _validate_input(
            x=x,
            y=y,
            t=t,
            w=w,
            bbox=bbox,
            k=k,
            eps=eps,
            check_finite=check_finite,
            sparse=sparse,
            smoothing=smoothing,
            penalty_order=penalty_order,
            unsupported=unsupported,
        )

        self._basis = None
        self._basis_counts = None
        self._coeffs = None
        self._coeffs_tensor = None
        self._residual = None

        self._fit()

    def __call__(self, x, nu=None):
        """Evaluate the spline or its derivatives at given positions."""
        return self.ev(x, nu=nu)

    @classmethod
    def from_grid(
        cls,
        axes,
        y,
        t,
        w=None,
        bbox=None,
        k=3,
        eps=None,
        check_finite=False,
        sparse=False,
        smoothing=0.0,
        penalty_order=2,
        unsupported="raise",
    ):
        """Create a spline from gridded data.

        Parameters
        ----------
        axes : sequence of array_like
            One 1D coordinate vector per input dimension.
        y : array_like
            Values on the tensor grid. Scalar-valued data should have shape
            ``(len(axis0), len(axis1), ...)``. Vector-valued data may either
            put the output dimension first or last.
        """
        x, y, w = _grid_to_scattered(axes, y, w=w)
        return cls(
            x=x,
            y=y,
            t=t,
            w=w,
            bbox=bbox,
            k=k,
            eps=eps,
            check_finite=check_finite,
            sparse=sparse,
            smoothing=smoothing,
            penalty_order=penalty_order,
            unsupported=unsupported,
        )

    def ev(self, x, nu=None):
        """Evaluate the spline at scattered points.

        Parameters
        ----------
        x : array_like, shape (n_eval, n_dimensions)
            Points where the spline should be evaluated.
        nu : int or sequence of int, optional
            Derivative order for each input dimension.
        """
        ndim = self._input.x.shape[1]
        x, was_scalar = _validate_eval_input(x, ndim=ndim)
        _validate_eval_domain(x, self._input.bbox)
        nu = _validate_derivative_orders(nu, ndim=ndim)

        s_matrix = _design_matrix(
            x,
            self._basis,
            self._input.k,
            nu=nu,
            sparse=self._input.sparse,
        )
        values = s_matrix @ self._coeffs

        if was_scalar:
            return values[0]

        return values

    def get_knots(self):
        """Return the knot vectors for each input dimension."""
        return self._basis

    def get_coeffs(self):
        """Return the fitted spline coefficients."""
        return self._coeffs.copy()

    def get_coeffs_tensor(self):
        """Return fitted coefficients in tensor-product basis shape."""
        return self._coeffs_tensor.copy()

    def get_residual(self):
        """Return the weighted sum of squared residuals."""
        return self._residual

    def _fit(self):
        """Fit the spline coefficients from the constructor inputs."""
        full_knots = _construct_full_knots(
            self._input.t,
            self._input.bbox,
            self._input.k,
        )
        basis_counts = _basis_counts_from_knots(full_knots, self._input.k)
        _check_fit_size(
            n_samples=self._input.x.shape[0],
            basis_counts=basis_counts,
            sparse=self._input.sparse,
            smoothing=self._input.smoothing,
        )
        s_matrix = _design_matrix(
            self._input.x,
            full_knots,
            self._input.k,
            sparse=self._input.sparse,
        )
        unsupported_count = _check_basis_support(
            s_matrix,
            unsupported=self._input.unsupported,
        )
        coeffs, residual = _solve_lsq(
            s_matrix,
            self._input.y,
            w=self._input.w,
            eps=self._input.eps,
            smoothing=self._input.smoothing,
            coefficient_shape=basis_counts,
            penalty_order=self._input.penalty_order,
            allow_rank_deficient=(
                unsupported_count > 0
                and self._input.unsupported != "raise"
            ),
            warn_rank_deficient=(
                unsupported_count > 0
                and self._input.unsupported == "warn"
            ),
        )

        self._basis = full_knots
        self._basis_counts = basis_counts
        self._coeffs = coeffs
        self._coeffs_tensor = _reshape_coefficients(coeffs, basis_counts)
        self._residual = residual


def _reshape_coefficients(coeffs, basis_counts):
    """Reshape flattened coefficients to tensor-product basis shape."""
    if coeffs.ndim == 1:
        return coeffs.reshape(basis_counts)

    return coeffs.reshape((*basis_counts, coeffs.shape[1]))


def _check_fit_size(n_samples, basis_counts, sparse, smoothing):
    """Validate that the requested fit is identifiable and reasonably sized."""
    n_coeffs = 1
    for basis_count in basis_counts:
        n_coeffs *= basis_count

    if smoothing == 0.0 and n_samples < n_coeffs:
        raise ValueError(
            "not enough data points for an unregularized least-squares "
            "spline; use fewer knots, lower degrees, or smoothing > 0."
        )

    max_dense_entries = 100_000_000
    if not sparse and n_samples * n_coeffs > max_dense_entries:
        raise ValueError(
            "dense design matrix would be too large; use sparse=True or "
            "fewer knots."
        )


def _check_basis_support(design_matrix, unsupported):
    """Validate that every tensor-product basis has sample support."""
    if hasattr(design_matrix, "getnnz"):
        supported = np.asarray(design_matrix.getnnz(axis=0)).ravel() > 0
    else:
        supported = np.any(design_matrix != 0, axis=0)

    unsupported_count = int(np.sum(~supported))
    if unsupported_count == 0:
        return unsupported_count

    message = (
        f"data points do not support {unsupported_count} tensor-product "
        "basis function(s); adjust knots, bbox, or sample locations."
    )
    if unsupported == "raise":
        raise ValueError(
            "data points do not support every tensor-product basis "
            "function; adjust knots, bbox, or sample locations."
        )

    if unsupported == "warn":
        warnings.warn(message, RuntimeWarning, stacklevel=3)

    return unsupported_count
