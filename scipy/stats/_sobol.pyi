from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy.typing as npt

def initialize_v(
    v : npt.ArrayLike, 
    dim : int
) -> None: ...

def _cscramble (
    dim : int,
    ltm : npt.ArrayLike,
    sv: npt.ArrayLike
) -> None: ...

def _fill_p_cumulative(
    p: npt.ArrayLike,
    p_cumulative: npt.ArrayLike
) -> None: ...

def _categorize(
    draws: npt.ArrayLike,
    p_cumulative: npt.ArrayLike,
    result: npt.ArrayLike
) -> None: ...

def _draw(
    n : int,
    num_gen: int,
    dim: int,
    sv: npt.ArrayLike,
    quasi: npt.ArrayLike,
    result: npt.ArrayLike
    ) -> None: ...

def _fast_forward(
    n: int,
    num_gen: int,
    dim: int,
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
 



