import numpy as np
from typing import Union, Any, Tuple, List, overload, Callable, NamedTuple
from typing_extensions import Protocol
import numpy.typing as npt
from scipy._lib._util import SeedType


ArrayLike0D = Union[bool, int, float, complex, str, bytes, np.generic]


__all__: List[str]


class UNURANError(RuntimeError):
    ...


class Method:
    @overload
    def rvs(self, size: None = ...) -> float | int: ...  # type: ignore[misc]
    @overload
    def rvs(self, size: int | Tuple[int, ...] = ...) -> np.ndarray: ...
    def set_random_state(self, random_state: SeedType) -> None: ...


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
                 *,
                 domain: None | Tuple[float, float] = ...,
                 c: float = ...,
                 construction_points: int | npt.ArrayLike = ...,
                 variant: str = ...,
                 use_dars: bool = ...,
                 max_squeeze_hat_ratio: float = ...,
                 max_intervals: int = ...,
                 guide_factor: float = ...,
                 random_state: SeedType = ...) -> None: ...
    @property
    def squeeze_hat_ratio(self) -> float: ...
    @property
    def squeeze_area(self) -> float: ...
    @overload
    def ppf_hat(self, u: ArrayLike0D) -> float: ...  # type: ignore[misc]
    @overload
    def ppf_hat(self, u: npt.ArrayLike) -> np.ndarray: ...


UError = NamedTuple('UError', [('max_error', float),
                               ('mean_absolute_error', float)])

class PINVDist(Protocol):
    @property
    def pdf(self) -> Callable[..., float]: ...
    @property
    def cdf(self) -> Callable[..., float]: ...


class NumericalInversePolynomial(Method):
    def __init__(self,
                 dist: PINVDist,
                 mode: None | float = ...,
                 center: None | float = ...,
                 *,
                 domain: None | Tuple[float, float] = ...,
                 order: int = ...,
                 u_resolution: float = ...,
                 max_intervals: int = ...,
                 keep_cdf: bool = ...,
                 random_state: SeedType = ...) -> None: ...
    @property
    def intervals(self) -> int: ...
    @overload
    def ppf(self, u: ArrayLike0D) -> float: ...  # type: ignore[misc]
    @overload
    def ppf(self, u: npt.ArrayLike) -> np.ndarray: ...
    @overload
    def cdf(self, x: ArrayLike0D) -> float: ...  # type: ignore[misc]
    @overload
    def cdf(self, x: npt.ArrayLike) -> np.ndarray: ...
    def u_error(self, sample_size: int = ...) -> UError: ...


class NROUDist(Protocol):
    @property
    def pdf(self) -> Callable[..., float]: ...
    @property
    def support(self) -> Tuple[float, float]: ...


class NaiveRatioUniforms(Method):
    def __init__(self,
                 dist: NROUDist,
                 *,
                 center: float = ...,
                 domain: None | Tuple[float, float] = ...,
                 r: float = ...,
                 u_min: None | float = ...,
                 u_max: None | float = ...,
                 v_max: None | float = ...,
                 random_state: SeedType = ...) -> None: ...


class DAUDist(Protocol):
    @property
    def pmf(self) -> Callable[..., float]: ...
    @property
    def support(self) -> Tuple[float, float]: ...

class DiscreteAliasUrn(Method):
    def __init__(self,
                 dist: npt.ArrayLike | DAUDist,
                 *,
                 domain: None | Tuple[float, float] = ...,
                 urn_factor: float = ...,
                 random_state: SeedType = ...) -> None: ...
