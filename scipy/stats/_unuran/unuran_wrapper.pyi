import numpy as np
from typing import Union, Any, Tuple, List, overload, Callable
from typing_extensions import Protocol
import numpy.typing as npt


SeedType = Union[np.random.RandomState, np.random.Generator, int,
                 npt.ArrayLike, np.random.SeedSequence]
ArrayLike0D = Union[bool, int, float, complex, str, bytes, np.generic]


__all__: List[str]


class UNURANError(RuntimeError):
    ...


class Method:
    @overload
    def rvs(self, size: None = ...) -> float | int: ...  # type: ignore[misc]
    @overload
    def rvs(self, size: int | Tuple[int, ...] = ...) -> np.ndarray: ...


class TDRDist(Protocol):
    @property
    def pdf(self) -> Callable[..., float]: ...
    @property
    def dpdf(self) -> Callable[..., float]: ...


class TransformedDensityRejection(Method):
    def __init__(self,
                 dist: TDRDist,
                 mode: None | float = ...,
                 center: None | float = ...,
                 params: Tuple[Any, ...] = ...,
                 domain: None | Tuple[float, float] = ...,
                 c: float = ...,
                 cpoints: int = ...,
                 variant: str = ...,
                 use_dars: bool = ...,
                 max_sqhratio: float = ...,
                 max_intervals: int = ...,
                 use_center: bool = ...,
                 use_mode: bool = ...,
                 guide_factor: float = ...,
                 seed: SeedType = ...) -> None: ...

    @property
    def sqhratio(self) -> float: ...
    @property
    def n_intervals(self) -> int: ...
    @property
    def squeeze_area(self) -> float: ...

    @overload
    def ppf_hat(self, u: ArrayLike0D) -> float: ...  # type: ignore[misc]
    @overload
    def ppf_hat(self, u: npt.ArrayLike) -> np.ndarray: ...


class DAUDist(Protocol):
    @property
    def pmf(self) -> Callable[..., float]: ...

class DiscreteAliasUrn(Method):
    def __init__(self,
                 pv: None | npt.ArrayLike = ...,
                 dist: None | DAUDist = ...,
                 params: Tuple[Any, ...] = ...,
                 domain: None | Tuple[float, float] = ...,
                 urn_factor: float = ...,
                 seed: SeedType = ...) -> None: ...
