# cython: cpow=True

from __future__ import annotations
from typing import TYPE_CHECKING, Union, Sequence
import numpy as np

if TYPE_CHECKING:
    import numpy.typing as npt

_IntegerType = int | np.integer

cdef class Translation:
    """
    A class to represent a 3D translation vector.

    Attributes
    ----------
    vector : np.ndarray
        The translation vector.

    Methods
    -------
    from_vector(vector: npt.ArrayLike) -> Translation
        Creates a Translation object from a vector.
    as_vector() -> np.ndarray
        Returns the translation vector.
    apply(points: npt.ArrayLike) -> np.ndarray
        Applies the translation to a set of points.
    inv() -> Translation
        Returns the inverse of the translation.
    magnitude() -> float
        Returns the magnitude of the translation vector.
    approx_equal(other: Translation, atol: float = 1e-8) -> bool
        Checks if two translations are approximately equal.
    concatenate(translations: Sequence[Translation]) -> Translation
        Concatenates a sequence of translations.
    identity() -> Translation
        Returns the identity translation.
    random(random_state: Union[_IntegerType, np.random.Generator, np.random.RandomState, None] = None) -> Translation
        Returns a random translation.
    """

    cdef np.ndarray _vector

    def __init__(self, vector: npt.ArrayLike, bint copy=True) -> None:
        """
        Initializes the Translation object.

        Parameters
        ----------
        vector : npt.ArrayLike
            The translation vector.
        copy : bool, optional
            Whether to copy the input vector (default is True).
        """
        self._vector = np.array(vector, copy=copy)
        if self._vector.ndim != 1 or self._vector.size != 3:
            raise ValueError("Translation vector must be a 1D array of size 3")

    @property
    def vector(self) -> np.ndarray:
        """Returns the translation vector."""
        return self._vector

    def __len__(self) -> int:
        """Returns the length of the translation vector (always 3)."""
        return 3

    @classmethod
    def from_vector(cls, vector: npt.ArrayLike) -> Translation:
        """Creates a Translation object from a vector."""
        return cls(vector)

    def as_vector(self) -> np.ndarray:
        """Returns the translation vector."""
        return self._vector

    def apply(self, points: npt.ArrayLike) -> np.ndarray:
        """
        Applies the translation to a set of points.

        Parameters
        ----------
        points : npt.ArrayLike
            The points to translate.

        Returns
        -------
        np.ndarray
            The translated points.
        """
        points = np.asarray(points)
        if points.shape[-1] != 3:
            raise ValueError("Points must have shape (N, 3) or (3,)")
        return points + self._vector

    def __add__(self, other: Translation) -> Translation:
        """Adds two translations."""
        if not isinstance(other, Translation):
            raise TypeError("Can only add Translation to Translation")
        return Translation(self._vector + other.vector)

    def __sub__(self, other: Translation) -> Translation:
        """Subtracts one translation from another."""
        if not isinstance(other, Translation):
            raise TypeError("Can only subtract Translation from Translation")
        return Translation(self._vector - other.vector)

    def __mul__(self, scalar: float) -> Translation:
        """Multiplies the translation by a scalar."""
        return Translation(self._vector * scalar)

    def __truediv__(self, scalar: float) -> Translation:
        """Divides the translation by a scalar."""
        return Translation(self._vector / scalar)

    def inv(self) -> Translation:
        """Returns the inverse of the translation."""
        return Translation(-self._vector)

    def magnitude(self) -> float:
        """Returns the magnitude of the translation vector."""
        return np.linalg.norm(self._vector)

    def approx_equal(self, other: Translation, atol: float = 1e-8) -> bool:
        """Checks if two translations are approximately equal."""
        return np.allclose(self._vector, other.vector, atol=atol)

    @classmethod
    def concatenate(cls, translations: Sequence[Translation]) -> Translation:
        """Concatenates a sequence of translations."""
        total_vector = np.sum([t.vector for t in translations], axis=0)
        return cls(total_vector)

    @classmethod
    def identity(cls) -> Translation:
        """Returns the identity translation."""
        return cls(np.zeros(3))

    @classmethod
    def random(cls, random_state: Union[_IntegerType, np.random.Generator, np.random.RandomState, None] = None) -> Translation:
        """Returns a random translation."""
        rng = np.random.default_rng(random_state)
        return cls(rng.uniform(-1, 1, size=3))

    def __getitem__(self, indexer: Union[int, slice]) -> float:
        """Gets an item from the translation vector."""
        return self._vector[indexer]