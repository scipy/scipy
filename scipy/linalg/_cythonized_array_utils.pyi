from numpy.typing import NDArray
from typing import Any

__all__ = ['issymmetric', 'ishermitian']

def issymmetric(
    a: NDArray[Any],
    atol: None | float = ...,
    rtol: None | float = ...,
) -> bool: ...

def ishermitian(
    a: NDArray[Any],
    atol: None | float = ...,
    rtol: None | float = ...,
) -> bool: ...
