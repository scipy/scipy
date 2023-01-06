from __future__ import annotations

from typing import (
    Tuple,
    TYPE_CHECKING,
)

import numpy as np

if TYPE_CHECKING:
    import numpy.typing as npt


def _validate_bounds(
    l_bounds: npt.ArrayLike, u_bounds: npt.ArrayLike, x0: np.ndarray
) -> Tuple[np.ndarray, ...]:
    """Bounds input validation.

    Parameters
    ----------
    l_bounds, u_bounds : array_like (d,)
        Lower and upper bounds.
    x0 : np.ndarray
        Array to use for broadcasting.

    Returns
    -------
    l_bounds, u_bounds : array_like (d,)
        Lower and upper bounds.

    """
    try:
        lower = np.broadcast_to(l_bounds, x0.shape)
        upper = np.broadcast_to(u_bounds, x0.shape)
    except ValueError as exc:
        msg = (
            "The number of bounds is not compatible with the length of "
            f"``x0={x0}``."
        )
        raise ValueError(msg) from exc

    if np.any(lower > upper):
        msg = "An upper bound is less than the corresponding lower bound."
        raise ValueError(msg)

    return lower, upper
