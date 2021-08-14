import numpy as np
from typing import Union, Any, Tuple, List, overload, Callable
from typing_extensions import Protocol
import numpy.typing as npt
from scipy._lib._util import SeedType


UNURANSeedType = Union[SeedType, np.random.SeedSequence]
ArrayLike0D = Union[bool, int, float, complex, str, bytes, np.generic]


__all__: List[str]


class UNURANError(RuntimeError):
    ...


class Method:
    @overload
    def rvs(self, size: None = ...) -> float | int: ...  # type: ignore[misc]
    @overload
    def rvs(self, size: int | Tuple[int, ...] = ...) -> np.ndarray: ...
    def set_numpy_rng(self, numpy_rng: UNURANSeedType) -> None: ...


class TDRDist(Protocol):
    @property
    def pdf(self) -> Callable[..., float]: ...
    @property
    def dpdf(self) -> Callable[..., float]: ...
    @property
    def support(self) -> Tuple[float, float]: ...


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
                 max_squeeze_hat_ratio: float = ...,
                 max_intervals: int = ...,
                 guide_factor: float = ...,
                 numpy_rng: UNURANSeedType = ...) -> None: ...
    @property
    def squeeze_hat_ratio(self) -> float: ...
    @property
    def squeeze_area(self) -> float: ...
    @overload
    def ppf_hat(self, u: ArrayLike0D) -> float: ...  # type: ignore[misc]
    @overload
    def ppf_hat(self, u: npt.ArrayLike) -> np.ndarray: ...


class DAUDist(Protocol):
    @property
    def pmf(self) -> Callable[..., float]: ...
    @property
    def support(self) -> Tuple[float, float]: ...

class DiscreteAliasUrn(Method):
    def __init__(self,
                 dist: npt.ArrayLike | DAUDist,
                 params: Tuple[Any, ...] = ...,
                 domain: None | Tuple[float, float] = ...,
                 urn_factor: float = ...,
                 numpy_rng: SeedType = ...) -> None: ...
