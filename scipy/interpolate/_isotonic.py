from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np

from scipy.optimize import isotonic_regression

if TYPE_CHECKING:
    import numpy.typing as npt


__all__ = ["IsotonicInterpolator"]


class IsotonicInterpolator():
    """Nonparametric isotonic interpolation via PAVA.

    This class provides nonparametric monotonically increasing interpolation
    via the pool adjacent violators algorithm (PAVA).

    This interpolator first sorts the values of ``y`` according to the order
    given by ``x``. It then fits the isotonic regression, see
    :func:`scipy.optimize.isotonic_regression` for details. Afterwards it
    constructs a linear interpolation of the solution of the isotonic
    regression.

    Parameters
    ----------
    x : (N,) array_like or None
        A 1-D array of real values according to which ``y`` will be sorted.
    y : (N,) array_like
        Response variable.
    weights : (N,) array_like or None
        Case weights.
    increasing : bool
        If true, do isotonic regression, i.e. monotonic increasing. If false, do
        antitonic regresission, i.e. monotonic decreasing.
        Default is True.
    check_finite : bool
        Whether to check that the input arrays contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
        Default is True.

    Attributes
    ----------
    x_ : array_like
        Increasing threshold values of ``x`` that define the blocks or pools
        or knot positions.
    y_ : array_like
        Estimated isotonic regression values for all blocks.

    Methods
    -------
    __call__

    See Also
    --------
    scipy.optimize.isotonic_regression : Solves the isotonic regression via PAVA.
    PchipInterpolator : PCHIP 1-D monotonic cubic interpolator.
    """
    def __init__(
            self,
            x: npt.ArrayLike,
            y: npt.ArrayLike,
            weights: npt.ArrayLike | None = None,
            *,
            increasing: bool = True,
            check_finite: bool = True,
    ):
        y = np.asarray(y)
        if y.ndim != 1:
            raise ValueError("The y array must have exactly one dimension.")
        if check_finite and not np.isfinite(y).all():
            raise ValueError("Array y must not contain infs or nans.")

        if x is not None:
            x = np.asarray(x)
            if x.ndim != 1:
                raise ValueError("The x array must have exactly one dimension.")
            if x.shape[0] != y.shape[0]:
                raise ValueError("The x and y arrays must have same length.")
            if check_finite and not np.isfinite(x).all():
                raise ValueError("Array x must not contain infs or nans.")
            # Note that sorting is often the performance bottleneck of
            # __init__.
            order = np.lexsort((y, x))  # The last key is the primary key.
            x = x[order]
            y = y[order]

        if weights is not None:
            weights = np.asarray(weights)
            if weights.ndim != 1:
                raise ValueError("The weights array must have exactly one dimension.")
            if weights.shape[0] != y.shape[0]:
                raise ValueError("The y and weights arrays must have same length.")
            if check_finite and not np.isfinite(weights).all():
                raise ValueError("Array weights must not contain infs or nans.")
            if x is not None:
                weights = weights[order]

        z, wz, indices = isotonic_regression(y, weights, increasing=increasing)

        # Construct knots at boundaries of blocks for interpolation.
        n_knots = np.sum(np.clip(np.diff(indices), None, 2))
        n_blocks = indices.shape[0] - 1
        self.x_ = np.empty(n_knots, dtype=x.dtype)
        self.y_ = np.empty(n_knots, dtype=z.dtype)
        j = 0
        for b in range(n_blocks):  # loop over blocks
            self.x_[j] = x[indices[b]]
            self.y_[j] = z[indices[b]]
            if indices[b + 1] > indices[b] + 1:
                # a block containing multiple elements
                self.x_[j + 1] = x[indices[b + 1] - 1]
                self.y_[j + 1] = z[indices[b + 1] - 1]
                j += 2
            else:
                # a single element block
                j += 1

    def __call__(self, x: npt.ArrayLike):
        return np.interp(x, self.x_, self.y_)
