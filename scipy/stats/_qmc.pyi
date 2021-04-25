# Prevents annotations from being evaluated during runtime,
# thus storing them as plain strings
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import (
    TYPE_CHECKING,
    Optional,
    Union,
    type_check_only,
)

import numpy as np

if TYPE_CHECKING:
    import numpy.typing as npt

_IntegerType = Union[int, np.integer]
_FloatingType = Union[float, np.floating]

__all__ = ['scale', 'discrepancy', 'QMCEngine', 'Sobol', 'Halton',
           'LatinHypercube', 'MultinomialQMC', 'MultivariateNormalQMC']


def check_random_state(
        seed: Optional[Union[_IntegerType, np.random.Generator,
                             np.random.RandomState]] = ...
) -> Union[np.random.Generator, np.random.RandomState]: ...


def scale(
        sample: npt.ArrayLike, l_bounds: npt.ArrayLike,
        u_bounds: npt.ArrayLike, reverse: bool = ...
) -> np.ndarray: ...


def discrepancy(
        sample: npt.ArrayLike, iterative: bool = ..., method: str = ...
) -> _FloatingType: ...


def update_discrepancy(
        x_new: npt.ArrayLike, sample: npt.ArrayLike,
        initial_disc: _FloatingType
) -> _FloatingType: ...


def primes_from_2_to(n: _IntegerType) -> np.ndarray: ...


def n_primes(n: _IntegerType) -> np.ndarray: ...


def van_der_corput(
        n: _IntegerType, base: _IntegerType = ...,
        start_index: _IntegerType = ..., scramble: bool = ...,
        seed: Optional[Union[_IntegerType, np.random.Generator]] = ...
) -> np.ndarray: ...


@type_check_only
class QMCEngine(ABC):
    @abstractmethod
    def __init__(
            self,
            d: _IntegerType,
            seed: Optional[Union[_IntegerType, np.random.Generator]] = ...
    ) -> None: ...

    @abstractmethod
    def random(self, n: _IntegerType = ...) -> np.ndarray: ...

    def reset(self) -> QMCEngine: ...

    def fast_forward(self, n: _IntegerType) -> QMCEngine: ...


class Halton(QMCEngine):
    def __init__(
            self, d: _IntegerType, scramble: bool = ...,
            seed: Optional[Union[_IntegerType, np.random.Generator]] = ...
    ) -> None: ...

    def random(self, n: _IntegerType = ...) -> np.ndarray: ...


class LatinHypercube(QMCEngine):
    def __init__(
            self, d: _IntegerType, centered: bool = ...,
            seed: Optional[Union[_IntegerType, np.random.Generator]] = ...
    ) -> None: ...

    def random(self, n: _IntegerType = ...) -> np.ndarray: ...

    def reset(self) -> LatinHypercube: ...


class Sobol(QMCEngine):
    def __init__(
            self, d: _IntegerType, scramble: bool = ...,
            seed: Optional[Union[_IntegerType, np.random.Generator]] = ...
    ) -> None: ...

    def _scramble(self) -> None: ...

    def random(self, n: _IntegerType = ...) -> np.ndarray: ...

    def random_base2(self, m: _IntegerType) -> np.ndarray: ...

    def reset(self) -> Sobol: ...

    def fast_forward(self, n: _IntegerType) -> Sobol: ...


class MultivariateNormalQMC(QMCEngine):
    def __init__(
            self, mean: npt.ArrayLike, cov: Optional[npt.ArrayLike] = ...,
            cov_root: Optional[npt.ArrayLike] = ...,
            inv_transform: bool = ...,
            engine: Optional[QMCEngine] = ...,
            seed: Optional[Union[_IntegerType, np.random.Generator]] = ...
    ) -> None: ...

    def random(self, n: _IntegerType = ...) -> np.ndarray: ...

    def _correlate(self, base_samples: np.ndarray) -> np.ndarray: ...

    def _standard_normal_samples(self,
                                 n: _IntegerType = ...) -> np.ndarray: ...


class MultinomialQMC(QMCEngine):
    def __init__(
            self, pvals: npt.ArrayLike, engine: Optional[QMCEngine] = ...,
            seed: Optional[Union[_IntegerType, np.random.Generator]] = ...
    ) -> None: ...

    def random(self, n: _IntegerType = ...) -> np.ndarray: ...
