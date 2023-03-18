from __future__ import annotations
from typing import TYPE_CHECKING

import numpy as np

from scipy.optimize import isotonic_regression

if TYPE_CHECKING:
    import numpy.typing as npt


__all__ = ["IsotonicInterpolator"]


# TODO: Should this inherit from _Interpolator1D?
class IsotonicInterpolator():
    """Nonparametric isotonic interpolation via PAVA.

    This class provides nonparametric monotonically increasing interpolation
    via the pool adjacent violators algorithm (PAVA).

    This interpolator first sorts the values of ``y`` according to the order
    given by ``x``. It then fits the isotonic regression, see
    ``scipy.optimize.isotonic_regression`` for details. Afterwards it
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

    Attributes
    ----------
    x_ : array_like
        Increasing threshold values of ``x`` that define the blocks or pools
        or knot positions.
    y_ : array_like
        Estimated isotonic regression values for all blocks.
    x_is_numeric : bool
        True if ``x_`` is numeric.

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
            increasing: bool = True,
    ):
        # TODO: Should we check np.isfinite(x, y, w)?
        y = np.asarray(y)
        if y.ndim != 1:
            raise ValueError("The y array must have exactly one dimension.")

        if x is not None:
            x = np.asarray(x)
            if x.ndim != 1:
                raise ValueError("The x array must have exactly one dimension.")
            if x.shape[0] != y.shape[0]:
                raise ValueError("The x and y arrays must have same length.")
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

        # TODO: Maybe this is overkill.
        self.x_is_numeric = np.issubdtype(self.x_.dtype, np.number)

    def __call__(self, x: npt.ArrayLike):
        if self.x_is_numeric:
            return np.interp(x, self.x_, self.y_)
        else:
            # TODO: Maybe this is overkill.
            # We do not interpolate linearly between blocks as for non numeric
            # data we have no notion of distance.
            x = np.asarray(x)
            idx = np.searchsorted(self.x_, x, side="right") - 1
            idx = np.clip(idx, 0, self.x_.shape[0] - 1)
            return self.y_[idx]
