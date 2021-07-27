from __future__ import annotations
from typing import Any, TYPE_CHECKING, Callable, Tuple
from .optimize import OptimizeResult
from ._constraints import Bounds
import numpy as np

if TYPE_CHECKING:
    import numpy.typing as npt

def _minimize_direct(
        func: Callable[[npt.ArrayLike, Tuple[Any]], float],
        bounds: Bounds,
        args: tuple,
        disp: bool,
        iatol: float,
        maxfun: int,
        maxiter: int,
        locally_biased: bool,
        f_min: float,
        f_min_per: float,
        vol_per: float,
        sigma_per: float
) -> OptimizeResult: ...
