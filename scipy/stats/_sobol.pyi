from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy.typing as npt

def _test_find_index(
    p_cumulative : npt.ArrayLike, 
    size : int, 
    value : float,
    ) -> int:...