import numpy as np
from typing import Union, Any, Protocol, Tuple, List
import numpy.typing as npt


SeedType = Union[np.random.RandomState, np.random.Generator, int,
                 npt.ArrayLike, np.random.SeedSequence]


__all__: List[str]


class Method:
    def rvs(self, size: None | int | tuple = ...) -> np.ndarray: ...


class TDRDist(Protocol):
    def pdf(self, x: float, *args) -> float: ...
    def dpdf(self, x: float, *args) -> float: ...


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

    def ppf_hat(self, u: npt.ArrayLike) -> float | np.ndarray: ...


class DAUDist(Protocol):
    def pmf(self, k: float, *args) -> float: ...

class DiscreteAliasUrn(Method):
    def __init__(self,
                 pv: None | npt.ArrayLike = ...,
                 dist: None | DAUDist = ...,
                 params: Tuple[Any, ...] = ...,
                 domain: None | Tuple[float, float] = ...,
                 urn_factor: float = ...,
                 seed: SeedType = ...) -> None: ...
