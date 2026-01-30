"""Helper module with various typing-related utilities."""

import functools
from collections.abc import Callable
from typing import Any, Protocol, TypeVar, Iterable, type_check_only

import uarray

_T = TypeVar("_T")
_TT = TypeVar("_TT", bound=type[Any])
_T_co = TypeVar("_T_co", covariant=True)

@type_check_only
class _PySequence(Protocol[_T_co]):
    def __len__(self) -> int: ...
    def __getitem__(self, key: int, /) -> _T_co: ...

@type_check_only
class _SupportsUA(Protocol):
    @property
    def __ua_domain__(self) -> str | _PySequence[str]: ...
    @property
    def __ua_cache__(self) -> dict[Any, Any]: ...
    def __ua_convert__(
        self,
        dispatchables: tuple[uarray.Dispatchable[Any, Any], ...],
        coerce: bool,
        /,
    ) -> Iterable[Any]: ...

@type_check_only
class _PartialDispatchable(functools.partial[uarray.Dispatchable[Any, _TT]]):
    func: type[uarray.Dispatchable[Any, _TT]]
    args: tuple[_TT]
    def __call__(  # type: ignore[override]
        self,
        value: _T,
        coercible: bool = ...,
    ) -> uarray.Dispatchable[_T, _TT]: ...

_ReplacerFunc = Callable[
    [
        tuple[Any, ...],
        dict[str, Any],
        tuple[uarray.Dispatchable[Any, Any], ...],
    ],
    tuple[tuple[Any, ...], dict[str, Any]],
]

_PyGlobalDict = dict[
    str,
    tuple[
        tuple[_T | None, bool, bool],
        list[_T],
        bool,
    ],
]

_PyLocalDict = dict[
    str,
    tuple[
        list[_T],
        list[tuple[_T, bool, bool]],
    ],
]
