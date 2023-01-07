from __future__ import annotations
from typing import Optional, TYPE_CHECKING

import numpy as np

from ._pava_pybind import pava

if TYPE_CHECKING:
    import numpy.typing as npt


__all__ = ["isotonic_regression", "IsotonicInterpolator"]


def isotonic_regression(
    y: npt.ArrayLike,
    weights: Optional[npt.ArrayLike] = None,
    increasing: bool = True,
):
    r"""Nonparametric isotonic regression.

    A monotonically increasing interpolant is calculated by the pool adjacent
    violators algorithm (PAVA), see [1]_. See the Notes section for more
    details.

    Parameters
    ----------
    y : (N,) array_like
        Response variable.
    weights : (N,) array_like or None
        Case weights.
    increasing : bool
        If True, fit monotonic increasing, i.e. isotonic, regression.
        If False, fit a monotonic decreasing, i.e. antitonic, regression.

    Returns
    -------
    x : (N,) array_like
        Isotonic regression solution, i.e. an increasing (or decresing) array
        of the same length than y.
    wx : (B,) array_like
        Sum of case weights for all blocks / pools B.
    r : (B+1,) array_like
        Array of indices with the start position of each block / pool B.
        For the j-th block, all values of ``x[r[j]:r[j+1]]`` are the same.

    Notes
    -----
    Given data :math:`y` and case weights :math:`w`, the isotonic regression
    solves the following optimazation problem:

    .. math::

        \operatorname{argmin}_{x_i} \sum_i w_i (y_i - x_i)^2 \quad
        \text{subject to } x_i \leq x_j \text{ whenever } i \leq j \,.
    
    For every input value :math:`y_i`, it generates an interpolated value
    :math:`x_i` which are increasing. This is accomplished by the PAVA.
    The solution consists of pools or blocks, i.e. neighboring elements of
    :math:`x`, e.g. :math:`x_i` and :math:`x_{i+1}`, that all have the same
    value.

    Most interestingly, the solution stays the same if the squared loss is
    replaced by the wide class of Bregman functions which are the unique
    class of strictly consistent scoring functions for the mean, see [2]_
    and references therein.

    References
    ----------
    .. [1] Busing, F. M. T. A. (2022).
           Monotone Regression: A Simple and Fast O(n) PAVA Implementation.
           Journal of Statistical Software, Code Snippets, 102(1), 1-25.
           :doi:`10.18637/jss.v102.c01`
    .. [2] Jordan, A.I., Mühlemann, A. & Ziegel, J.F.
           Characterizing the optimal solutions to the isotonic regression
           problem for identifiable functionals.
           Ann Inst Stat Math 74, 489-514 (2022).
           :doi:`10.1007/s10463-021-00808-0`
    """
    y = np.asarray(y)
    if weights is None:
        weights = np.ones_like(y)
    else:
        weights = np.asarray(weights)

        if not (y.ndim == weights.ndim and y.shape[0] == weights.shape[0]):
            raise ValueError(
                "Input arrays y and w must have one dimension of equal length."
            )
        if np.any(weights <= 0):
            raise ValueError("Weights w must be strictly positive.")

    order = np.s_[:] if increasing else np.s_[::-1]
    x = np.array(y[order], order="C", dtype=np.float64, copy=True)
    wx = np.array(weights[order], order="C", dtype=np.float64, copy=True)
    n = x.shape[0]
    r = np.full(shape=n + 1, fill_value=-1, dtype=np.intp)
    pava(x, wx, r)
    r = r[r >= 0]
    # Due to the pava implementation, after the last block index, there might
    # be smaller numbers appended to r, e.g. r = [0, 10, 8, 7] which should be
    # r = [0, 10].
    imax = np.argmax(r)
    if imax < r.shape[0]:
        r = r[:imax + 1]
    n = r.shape[0] - 1
    wx = wx[:n]
    if not increasing:
        x = x[::-1]
        wx = wx[::-1]
        r = r[-1] - r[::-1]
    return x, wx, r


# TODO: Should this inherit from _Interpolator1D?
class IsotonicInterpolator():
    """Nonparametric isotonic interpolation via PAVA.
    
    This class provides nonparametric monotonically increasing interpolation
    via the pool adjacent violators algorithm (PAVA).
    
    This interpolator first sorts the values of ``y`` according to the order
    given by ``x``. It then fits the isotonic regression, see
    ``isotonic_regression`` for details. Afterwards it constructs a linear
    interpolation of the solution of the isotonic regression.

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
    isotonic_regression : Solvevs the isotonic regression via PAVA.
    PchipInterpolator : PCHIP 1-D monotonic cubic interpolator.
    """
    def __init__(self, x, y, weights=None, increasing=True):
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

    def __call__(self, x):
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
