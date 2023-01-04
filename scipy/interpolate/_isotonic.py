from __future__ import annotations
from typing import Optional, TYPE_CHECKING

import numpy as np

from ._pava_pybind import pava

if TYPE_CHECKING:
    import numpy.typing as npt


__all__ = ["isotonic_regression"]


def isotonic_regression(
    y: npt.ArrayLike,
    w: Optional[npt.ArrayLike] = None,
    increasing: bool = True,
):
    """Nonparametric isotonic regression.
    
    Parameters
    ----------
    y : array_like
        Response variable.
    w : array_like or None
    increasing : bool
        If True, fit monotonic increasing, i.e. isotonic, regression. If False, fit
        a monotonic decreasing, i.e. antitonic, regression.

    Returns
    -------
    x : array_like
        Isotonic regression solution, i.e. an increasing (or decresing) array of the
        same length than y.

    Notes
    -----

    References
    ----------
    .. [1] Busing, F. M. T. A. (2022).
           Monotone Regression: A Simple and Fast O(n) PAVA Implementation.
           Journal of Statistical Software, Code Snippets, 102(1), 1-25.
           https://doi.org/10.18637/jss.v102.c01
    """
    y = np.asarray(y)
    if w is None:
        w = np.ones_like(y)
    else:
        w = np.asarray(w)
    
        if not (y.ndim == w.ndim and y.shape[0] == w.shape[0]):
            raise ValueError(
                "Input arrays y and w must have one dimension of equal length."
            )
        if np.any(w <= 0):
            raise ValueError("Weights w must be strictly positive.")

    order = np.s_[:] if increasing else np.s_[::-1]
    x = np.array(y[order], order="C", dtype=np.float64, copy=True)
    wx = np.array(w[order], order="C", dtype=np.float64, copy=True)
    n = x.shape[0]
    r = np.full(shape=n + 1, fill_value=-1, dtype=np.intp)
    pava(x, wx, r)
    r = r[r >= 0]
    n = r.shape[0] - 1
    wx = wx[:n]
    if not increasing:
       x = x[::-1]
       wx = wx[::-1]
       r = n - r[::-1]
    return x, wx, r
    