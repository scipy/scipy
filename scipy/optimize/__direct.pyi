from __future__ import annotations
from _typeshed import NoneType
from typing import Any, Callable, Iterable, Tuple, TYPE_CHECKING, Union
from scipy.optimize import OptimizeResult
from ._constraints import Bounds

if TYPE_CHECKING:
    import numpy.typing as npt

def direct(
    func: Callable[[npt.ArrayLike, Tuple[Any]], float],
    bounds: Union[Iterable, Bounds],
    args: tuple,
    disp: bool,
    iatol: float,
    maxfun: int,
    maxiter: int,
    locally_biased: bool,
    f_min: float,
    f_min_per: float,
    vol_per: float,
    sigma_per: float,
    callback: Callable[[npt.ArrayLike], NoneType]
) -> OptimizeResult: ...
