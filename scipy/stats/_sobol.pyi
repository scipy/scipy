import numpy.typing as npt
from scipy._lib._util import IntNumber

def initialize_v(
    v : npt.ArrayLike, 
    dim : IntNumber
) -> None: ...

def _cscramble (
    dim : IntNumber,
    ltm : npt.ArrayLike,
    sv: npt.ArrayLike
) -> None: ...

def _fill_p_cumulative(
    p: npt.ArrayLike,
    p_cumulative: npt.ArrayLike
) -> None: ...

def _draw(
    n : IntNumber,
    num_gen: IntNumber,
    dim: IntNumber,
    sv: npt.ArrayLike,
    quasi: npt.ArrayLike,
    result: npt.ArrayLike
    ) -> None: ...

def _fast_forward(
    n: IntNumber,
    num_gen: IntNumber,
    dim: IntNumber,
    sv: npt.ArrayLike,
    quasi: npt.ArrayLike
    ) -> None: ...

def _categorize(
    draws: npt.ArrayLike,
    p_cumulative: npt.ArrayLike,
    result: npt.ArrayLike
    ) -> None: ...

def initialize_direction_numbers() -> None: ...

_MAXDIM: int = ...
_MAXBIT: int = ...

def _test_find_index(
    p_cumulative: npt.ArrayLike, 
    size: int, 
    value: float
    ) -> int: ...