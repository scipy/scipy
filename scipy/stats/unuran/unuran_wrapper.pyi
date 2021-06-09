import numpy as np
from typing import TYPE_CHECKING, Union, Any, Callable, Tuple, Optional, List

if TYPE_CHECKING:
    import numpy.typing as npt


SeedType = Union[np.random.RandomState, np.random.Generator, int,
                 npt.ArrayLike]


__all__: List[str]


class Method:
    ...


class TransformedDensityRejection(Method):
    def __init__(self,
                 pdf: Callable, dpdf: Callable,
                 params: Tuple[Any] = ...,
                 domain: Optional[Tuple[float, float]] = ...,
                 c: float = ..., cpoints: int = ...,
                 variant: str = ...,
                 seed: SeedType = ...) -> None: ...

    def rvs(self, size: Optional[int] = ...) -> np.ndarray: ...


class DiscreteAliasUrn(Method):
    def __init__(self,
                 pv_or_pmf: Union[Callable, npt.ArrayLike],
                 params: Tuple[Any] = ...,
                 domain: Optional[Tuple[float, float]] = ...,
                 urnfactor: float = ...,
                 seed: SeedType = ...) -> None: ...

    def rvs(self, size: Optional[int] = ...) -> np.ndarray: ...
