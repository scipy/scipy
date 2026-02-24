from __future__ import annotations

import numpy as np

from scipy.sparse import issparse, sparray

from collections.abc import Callable
from typing import Protocol, TypeVar
from numbers import Number

#ToDo:
# - tests
# - examples
# - backwards compatibility with old LinearOperator (provide some utility to make it easier to manually upgrade)
# - fixing appearence of operators in the docs
# - fixing doc of AsLinearOperatorDunderProtocol as __aslinearoperator__ does not show up in the docs

__all__ = ['LinearOperator', 'MatrixLinearOperator', 'IdentityOperator', 'aslinearoperator']

def priority(func):
    """Priority decorator to request ``NotImplemented``.
    
    Protocol to provide interoperability between `LinearOperator` and unknown / custom classes.
    
    Consider the operation ``A @ B``, then Python will first try to call ``A.__matmul__(B)``. Only if this returns ``NotImplemented``, Python will then try to call ``B.__rmatmul__(A)``.
    
    If ``A`` is a `LinearOperator`, the object ``B`` might request a ``NotImplemented`` from ``A`` to take over the operation by adding the attribute ``__linop_priority__`` to ``B`` with a value higher than the ``__linop_priority__`` of ``A`` (which is 0.0 by default).
    """
    def priority_wrapper(self: LinearOperator, other):
        if self.__linop_priority__ < getattr(other,  '__linop_priority__', -1.0):
                return NotImplemented
        return func(self, other)
    return priority_wrapper


##############################################################
# Array Protocol
##############################################################

class MinimalArrayProtocol(Protocol):
    """Minimal array protocol to ensure compatibility with the `LinearOperator` interface."""
    ndim: int
    shape: tuple[int, ...]
    def transpose(self) -> MinimalArrayProtocol: ...
    
class ArrayProtocol(MinimalArrayProtocol):
    """Array protocol to ensure compatibility with the `LinearOperator` interface.
    
    Attributes
    ----------
    ndim : int
        Number of dimensions.
    shape : tuple[int, ...]
        Array dimensions.
    """
    
    def transpose(self) -> ArrayProtocol:
        """Transpose of the array.
        
        Returns
        -------
        T : `ArrayProtocol`
            Transposed array.
        """
        ...
    def conjugate(self) -> ArrayProtocol:
        """Optional: Complex conjugate of the array.
        
        This method is only needed in combination with complex-valued `LinearOperator` to provide default implementations for the
        
        1. Hermitian adjoint (if the transpose is available), 
        2. the transpose (if the Hermitian adjoint is available) and
        3. the complex conjugate.
        
        Thus, this method is not needed in combination with a real-valued `LinearOperator` or a complex-valued `LinearOperator` that does not use a default implementation of the adjoint, transpose or complex conjugate. 
        
        Returns
        -------
        T : `ArrayProtocol`
            Complex conjugate array.
        """
        ...

##############################################################
# Basic Linear Operator Interface
##############################################################

Array = TypeVar('Array', bound=MinimalArrayProtocol)

class LinearOperator[Array]:
    """Abstract interface for performing matrix-free matrix-vector products.
    
    Many iterative algorithms in `scipy.sparse.linalg`, e.g. `cg`, `gmres`, `lsmr`, `eigs`, `svds`, etc.
    do not need access to individual matrix entries, but only require the evaluation of matrix-vector products
    ``y = A @ x`` where ``x`` is a given dense vector (`numpy.ndarray`). This class serves as an abstract interface
    between these iterative algorithms and matrix-like objects that define a matrix-vector product.
    
    Although all algorithms within `scipy.sparse.linalg` take a `LinearOperator` as input and the user-provided implementations must be compatible with `numpy.ndarray`, the interface itself is compatible with any array-like object that satisfies the `ArrayProtocol`, e.g. `numpy.ndarray`, `scipy.sparse.sparray` or custom objects, as long as the user-provided implementations are compatible with that array-like object. For type annotations, the type of the compatible array-like objects of a `LinearOperator` can be denoted by ``LinearOperator[Array]`` where ``Array`` is the type of the array-like object, e.g. algorithms within `scipy.sparse.linalg` always require ``LinearOperator[np.ndarray]``.
    
    To construct a concrete `LinearOperator`, either pass appropriate callables to the constructor as described in the *Parameters* section below or subclass `LinearOperator`. A subclass representing a matrix ``A`` might implement the following attributes and methods:
    
    Mandatory:
        - ``shape``: tuple of two positive integers ``(m,n)`` denoting the shape of ``A``
        - ``dtype``: data type of ``A``
        - ``_matmul``: matrix-vector product ``A @ x``
    
    Optional:
        - Transpose: (implement one of the following)
            - ``_tmatmul``: transpose matrix-vector product ``A.T @ x``
            - ``_transpose``: transpose operator ``A.T``
        - Hermitian adjoint: (implement one of the following)
            - ``_ctmatmul``: Hermitian adjoint matrix-vector product ``A.H @ x``
            - ``_adjoint``: Hermitian adjoint operator ``A.H``
        - Complex conjugate: (implement one of the following)
            - ``_cmatmul``: complex conjugate matrix-vector product ``A.C @ x``
            - ``_conjugate``: complex conjugate operator ``A.C``
            
    For more information, see *Notes* below.
            
    Parameters
    ----------
    shape : tuple[int, int]
        matrix dimension ``(m,n)``
    matmul : Callable func(x)
        matrix-vector product, returning ``A @ x``
    tmatmul : Callable func(x) | None, optional
        transpose matrix-vector product, returning ``A.T @ x``
    ctmatmul : Callable func(x) | None, optional
        Hermitian adjoint matrix-vector product, returning ``A.H @ x``
    cmatmul : Callable func(x) | None, optional
        complex conjugate matrix-vector product, returning ``A.C @ x``
    dtype : np.dtype | None, optional
        data type of the matrix, if ``None``, the data type is determined as described below
        
    Attributes
    ----------
    ndim : int
        Number of dimensions (this is always 2).
    shape : tuple[int, int]
        Matrix dimension ``(m,n)``.
    dtype : `numpy.dtype`
        Data type of the matrix.
    T : `LinearOperator`
        Transpose linear operator.
    H : `LinearOperator`
        Hermitian adjoint linear operator.
    C : `LinearOperator`
        Complex conjugate linear operator.
    

    Notes
    -----
    
    - All user-provided implementations must properly handle the cases where ``x`` has shape ``(n,)`` (1D vector) and ``(n,k)`` (2D matrix).

    - It is recommended to call ``super().__init__(shape, dtype=dtype)`` in the subclass constructor to have the mandatory attributes validated. If the `LinearOperator` is compatible with `numpy.ndarray`, ``dtype`` might be ``None``. In this case, the ``dtype`` is determined by applying ``_matmul`` to ``np.zeros(n, dtype=np.int8)`` and using the promoted ``dtype`` of the output (`numpy.int8` is the smallest `numpy.dtype`).
    
    - To enable all functionalities of the ``LinearOperator``, it is sufficient to either provide an implementation of the transpose (recommended) `or` the Hermitian adjoint, while the other is then available by a default implementation using ``A.H @ x = (A.T @ x.conjugate()).conjugate()`` and ``A.T @ x = (A.H @ x.conjugate()).conjugate()``, respectively. Similarly, the complex conjugate is always available by a default implementation using ``A.C @ x = (A @ x.conjugate()).conjugate()``. The user might replace the default implementation by providing implementations for the corresponding operation.
    
    - To avoid constructing the same object multiple times, the transpose, adjoint and conjugate are cached as attributes ``_T``, ``_H`` and ``_C``, respectively and are constructed only once. If the content of the `LinearOperator` changes dynamically, the user might need to clear the cached attributes by setting them to `None`.
    
    The `LinearOperator` interface supports the following operations, where ``A`` and ``B`` are `LinearOperator` objects, ``alpha`` is a scalar and ``x`` is an array-like object compatible with the `LinearOperator`:
    
    - matrix-vector product: ``A @ x`` and ``x @ A``, see :meth:`LinearOperator.__matmul__` and :meth:`LinearOperator.__rmatmul__`
    - matrix-matrix product: ``A @ B`` and ``B @ A``, see :meth:`LinearOperator.__matmul__` and :meth:`LinearOperator.__rmatmul__`
    - scalar multiplication: ``alpha * A`` and ``A * alpha``, see :meth:`LinearOperator.__mul__` and :meth:`LinearOperator.__rmul__`
    - scalar division: ``A / alpha``, see :meth:`LinearOperator.__truediv__`
    - addition: ``A + B`` and ``B + A``, see :meth:`LinearOperator.__add__` and :meth:`LinearOperator.__radd__`
    - subtraction: ``A - B`` and ``B - A``, see :meth:`LinearOperator.__sub__` and :meth:`LinearOperator.__rsub__`
    - negation: ``-A``, see :meth:`LinearOperator.__neg__`
    - power: ``A ** p``, see :meth:`LinearOperator.__pow__`
        
    See Also
    --------
    MatrixLinearOperator : implicit construction of a `LinearOperator` from a matrix-like object
    aslinearoperator : function to convert an object to a `LinearOperator`
    

    """
    
    ndim: int = 2
    """Number of dimensions (this is always 2)."""
    shape: tuple[int, int]
    """Matrix dimension ``(m,n)``."""
    dtype: np.dtype
    """Data type of the matrix."""
    
    _T: LinearOperator[Array] | None = None
    _H: LinearOperator[Array] | None = None
    _C: LinearOperator[Array] | None = None
    
    __linop_priority__: float = 0.0
    
    def __init__(self, shape: tuple[int, int],
                       matmul:   Callable[[Array], Array] | None = None, 
                       tmatmul:  Callable[[Array], Array] | None = None,
                       ctmatmul: Callable[[Array], Array] | None = None, 
                       cmatmul:  Callable[[Array], Array] | None = None,
                       dtype: np.dtype | None = None):
        
        if len(shape) != 2 or not all(isinstance(s, int) and s > 0 for s in shape):
            raise ValueError(f'shape must be a tuple of two positive integers, got {shape}')
        self.shape = shape
        
        if matmul is not None:   self._matmul:   Callable[[Array], Array] = matmul
        if tmatmul is not None:  self._tmatmul:  Callable[[Array], Array] = tmatmul
        if ctmatmul is not None: self._ctmatmul: Callable[[Array], Array] = ctmatmul
        if cmatmul is not None:  self._cmatmul:  Callable[[Array], Array] = cmatmul
        
        if '_matmul' not in self.__dict__ and type(self)._matmul == LinearOperator._matmul:
            raise NotImplementedError("matmul is not defined. Please implement _matmul (subclassing LinearOperator) or provide the matmul argument (instantiating LinearOperator).")
        
        if dtype is None:
            self.dtype = self._get_dtype()
        else:
            self.dtype = np.dtype(dtype)
            

    ##################################################
    # Public API
    ##################################################
        
    def matmul(self, other: Array | LinearOperator[Array] | AsLinearOperatorDunderProtocol[Array]) -> Array | LinearOperator[Array]:
        """Matrix-vector or matrix-matrix product.
        
        Calculates ``y = self @ other`` where ``other`` is either
        
        1. an array-like object satisfying the `ArrayProtocol` (e.g. `numpy.ndarray` or `scipy.sparse.sparray`) which is compatible with this `LinearOperator`,
        2. another `LinearOperator` or 
        3. any object satisfying the `AsLinearOperatorDunderProtocol`.
        
        In the last two cases, ``y`` is a new `LinearOperator`, while in the first case, ``y`` is the result of the user-provided ``matmul`` implementation applied to ``other``.
        
        If ``self`` has shape ``(m,n)``, then ``other`` must have shape ``(n,)`` or ``(n,k)``, resulting in the output ``y`` having shape ``(m,)`` or ``(m,k)``, respectively.
        
        Parameters
        ----------
        other : `ArrayProtocol`, `LinearOperator`, `AsLinearOperatorDunderProtocol`
            The array-like object or linear operator to be multiplied.
            
        Returns
        -------
        y : `ArrayProtocol`, `LinearOperator`
        """
        
        if hasattr(other, '__aslinearoperator__'):
            other = aslinearoperator(other)
        
        if self.shape[1] != other.shape[0]:
            raise ValueError(f'shapes {self.shape} and {other.shape} are not compatible')
        
        if isinstance(other, LinearOperator):
            return _ProductLinearOperator[Array](self, other)
        else:
            return self._matmul(other)
        
    def transpose(self) -> LinearOperator[Array]:
        """Transpose linear operator."""
        if self._T is None: 
            self._T = self._transpose()
        return self._T
    
    def adjoint(self) -> LinearOperator[Array]:
        """Hermitian adjoint linear operator."""
        if self._H is None: 
            if self.dtype.kind == 'c':
                self._H = self._adjoint()
            else:
                self._H = self.transpose()
        return self._H
    
    def conjugate(self) -> LinearOperator[Array]:
        """Complex conjugate linear operator."""
        if self._C is None: 
            if self.dtype.kind == 'c':
                self._C = self._conjugate()
            else:
                self._C = self
        return self._C
    
    T: LinearOperator[Array] = property(transpose)
    H: LinearOperator[Array] = property(adjoint)
    C: LinearOperator[Array] = property(conjugate)
    
    ##################################################
    # Operator overloads and dunder methods
    ##################################################
    
    # @ operator
    @priority
    def __matmul__(self, other: Array | LinearOperator[Array]) -> Array | LinearOperator[Array]:
        f"""{LinearOperator.matmul.__doc__}
        
        Notes
        -----
        Respects the priority protocol to request ``NotImplemented``, see :func:`priority`.
        """
        return self.matmul(other)
        
    def __rmatmul__(self, other: Array | LinearOperator[Array]) -> Array | LinearOperator[Array]:
        """Right Matrix-vector or matrix-matrix product.
        
        Calculates ``y = other @ self`` where ``other`` is either
        
        1. an array-like object satisfying the `ArrayProtocol` (e.g. `numpy.ndarray` or `scipy.sparse.sparray`) which is compatible with this `LinearOperator`,
        2. another `LinearOperator` or
        3. any object satisfying the `AsLinearOperatorDunderProtocol`.
        
        In the last two cases, ``y`` is a new `LinearOperator`, while in the first case, ``y`` is an array-like object.
        
        If ``self`` has shape ``(m,n)``, then ``other`` must have shape ``(m,)`` or ``(k,m)``, resulting in the output ``y`` having shape ``(n,)`` or ``(k,n)``, respectively.
        
        Parameters
        ----------
        other : `ArrayProtocol`, `LinearOperator`, `AsLinearOperatorDunderProtocol`
            The array-like object or linear operator to be multiplied.
            
        Returns
        -------
        y : `ArrayProtocol`, `LinearOperator`
        """
        return self.matmul(other.transpose()).transpose()
    
    # * operator
    def __mul__(self, other: Number) -> LinearOperator[Array]:
        """Scalar multiplication.
        
        Calculates ``y = self * other`` where ``other`` is a scalar.
        
        Parameters
        ----------
        other : scalar
            The scalar to be multiplied.
            
        Returns
        -------
        y : `LinearOperator`
            The scaled linear operator.
            
        Notes
        -----
        Returns ``NotImplemented`` if ``other`` is not a scalar.
        """
        if np.isscalar(other):
            return _ScaledLinearOperator[Array](self, other)
        else:
            return NotImplemented
        
    def __rmul__(self, other: Number) -> LinearOperator[Array]:
        """Right Scalar multiplication.
        
        Calculates ``y = other * self``.
        
        Same as :meth:`LinearOperator.__mul__` since scalar multiplication is commutative.
        """
        return _ScaledLinearOperator[Array](self, other)
    
    # / operator
    def __truediv__(self, other: Number) -> LinearOperator[Array]:
        """Scalar division.
        
        Calculates ``y = self / other``.
        
        Same as :meth:`LinearOperator.__mul__` with the reciprocal of the scalar, i.e. ``y = self / other`` is equivalent to ``y = self * (1/other)``.
        """
        if np.isscalar(other):
            return _ScaledLinearOperator[Array](self, np.int8(1)/other)
        else:
            return NotImplemented
    
    # + operator    
    @priority
    def __add__(self, other: AsLinearOperatorProtocol[Array]) -> LinearOperator[Array]:
        """Addition of two linear operators.
        
        Calculates ``y = self + other``.
        
        Parameters
        ----------
        other : any object compatible with `aslinearoperator`
            The linear operator to be added.
            
        Returns
        -------
        y : `LinearOperator`
        
        Notes
        -----
        Respects the priority protocol to request ``NotImplemented``, see :func:`priority`.
        """
        return _SumLinearOperator[Array](self, aslinearoperator(other))
    
    def __radd__(self, other: AsLinearOperatorProtocol[Array]) -> LinearOperator[Array]:
        """Right addition of two linear operators.
        
        Calculates ``y = other + self``.
        
        Same as :meth:`LinearOperator.__add__` since addition is commutative.
        """
        return _SumLinearOperator[Array](aslinearoperator(other), self)
    
    # - operator
    def __neg__(self) -> LinearOperator[Array]:
        """Negation of a linear operator.
        
        Calculates ``y = -self``.
        
        Returns:
        --------
        y : `LinearOperator`
            The negated linear operator.
        """
        return _ScaledLinearOperator[Array](self, np.int8(-1))
    
    @priority
    def __sub__(self, other: AsLinearOperatorProtocol[Array]) -> LinearOperator[Array]:
        """Subtraction of two linear operators.
        
        Calculates ``y = self - other``.
        
        Parameters
        ----------
        other : any object compatible with `aslinearoperator`
            The linear operator to be subtracted.
            
        Returns
        -------
        y : `LinearOperator`
        
        Notes
        -----
        Respects the priority protocol to request ``NotImplemented``, see :func:`priority`.
        """
        return _SumLinearOperator[Array](self, -aslinearoperator(other))
    
    def __rsub__(self, other: AsLinearOperatorProtocol[Array]) -> LinearOperator[Array]:
        """Right subtraction of two linear operators.
        
        Calculates ``y = other - self``.
        
        Parameters
        ----------
        other : any object compatible with `aslinearoperator`
            The linear operator from which to subtract.
            
        Returns
        -------
        y : `LinearOperator`
        """
        
        return _SumLinearOperator[Array](-self, aslinearoperator(other))

    # ** operator
    def __pow__(self, other: int) -> LinearOperator[Array]:
        """Power or matrix-exponential of a linear operator.
        
        Calculates ``y = self ** other`` where ``other`` is a positive integer.
        
        Parameters
        ----------
        other : int
            The exponent to which the linear operator is raised. Must be a positive integer.
            
        Returns
        -------
        y : `LinearOperator`
            The linear operator raised to the power of ``other``.
        
        Notes
        -----
        Returns ``NotImplemented`` if ``other`` is not a positive integer.
        """
        if np.isscalar(other):
            if other < 0 or not np.issubdtype(type(other), np.integer):
                raise ValueError(f'only positive integers supported, got {other}')
            elif other == 0:
                return IdentityOperator[Array](shape=self.shape, dtype=self.dtype)
            elif other == 1:
                return self
            else:
                return _ProductLinearOperator[Array](*[self for _ in range(other)])
        else:
            return NotImplemented

    # one-line representation
    def __repr__(self):
        M,N = self.shape
        if self.dtype is None:
            dt = 'unspecified dtype'
        else:
            dt = 'dtype=' + str(self.dtype)

        return f'<{M}x{N} {self.__class__.__name__} with {dt}>'
    
    
    ##################################################
    # Internal API, might be overwritten by subclasses
    ##################################################
    
    def _matmul(self, other: Array) -> Array:
        """Matrix-vector or matrix-matrix product.
        
        This method must implement the action of this `LinearOperator` on an array-like object ``other`` compatible with this `LinearOperator`, i.e. it must implement the matrix-vector or matrix-matrix product ``y = self @ other``.
        
        Parameters
        ----------
        other : `ArrayProtocol`
            The array-like object to be multiplied. Must have shape ``(n,)`` or ``(n,k)`` if ``self`` has shape ``(m,n)``.
            
        Returns
        -------
        y : `ArrayProtocol`
             The result of the matrix-vector or matrix-matrix product, having shape ``(m,)`` or ``(m,k)``, respectively.
        """
        
        # should not be reached if super().__init__ is called in subclass __init__
        raise NotImplementedError("_matmul is not implemented, but this method is mandatory. Please implement _matmul or provide the matmul argument.")
        
    def _cmatmul(self, other: Array) -> Array:
        """Complex conjugate matrix-vector or matrix-matrix product.
        
        This is the default implementation of the action of the complex conjugate of this `LinearOperator` on an array-like object ``other`` compatible with this `LinearOperator`, i.e. it must implement the matrix-vector or matrix-matrix product ``y = self.C @ other``.
        
        Parameters
        ----------
        other : `ArrayProtocol`
            The array-like object to be multiplied. Must have shape ``(n,)`` or ``(n,k)`` if ``self`` has shape ``(m,n)``.
            
        Returns
        -------
        y : `ArrayProtocol`
             The result of the complex conjugate matrix-vector or matrix-matrix product, having shape ``(m,)`` or ``(m,k)``, respectively.
        """
        return self._matmul(other.conjugate()).conjugate()
        
    def _tmatmul(self, other: Array) -> Array:
        """Transpose matrix-vector or matrix-matrix product.
        
        This is the default implementation of the action of the transpose of this `LinearOperator` on an array-like object ``other`` compatible with this `LinearOperator`, i.e. it implements the matrix-vector or matrix-matrix product ``y = self.T @ other``.
        
        Parameters
        ----------
        other : `ArrayProtocol`
            The array-like object to be multiplied. Must have shape ``(m,)`` or ``(m,k)`` if ``self`` has shape ``(m,n)``.
            
        Returns
        -------
        y : `ArrayProtocol`
             The result of the transposed matrix-vector or matrix-matrix product, having shape ``(n,)`` or ``(n,k)``, respectively.
        """
        # check if _ctmatmul was a) set as attribute in __init__ or b) overwritten by subclass to avoid infinite recursion
        if '_ctmatmul' in self.__dict__ or type(self)._ctmatmul != LinearOperator._ctmatmul:
            return self._ctmatmul(other.conjugate()).conjugate()
        else:
            AH = self.adjoint()
            # check if adjoint is not the default implementation to avoid infinite recursion
            if type(AH) is not (_AdjointLinearOperator, _TransposeLinearOperator):
                return AH._matmul(other.conjugate()).conjugate()
            else:
                raise NotImplementedError("Please provide tmatmul or ctmatmul or overwrite _ctmatmul, _tmatmul, _transpose or _adjoint.")
    
    def _ctmatmul(self, other: Array) -> Array:
        """Hermitian adjoint matrix-vector or matrix-matrix product.
        
        This is the default implementation of the action of the Hermitian adjoint (i.e. transpose and conjugate) of this `LinearOperator` on an array-like object ``other`` compatible with this `LinearOperator`, i.e. it implements the matrix-vector or matrix-matrix product ``y = self.H @ other``.
        
        Parameters
        ----------
        other : `ArrayProtocol`
            The array-like object to be multiplied. Must have shape ``(m,)`` or ``(m,k)`` if ``self`` has shape ``(m,n)``.
            
        Returns
        -------
        y : `ArrayProtocol`
             The result of the Hermitian adjoint matrix-vector or matrix-matrix product, having shape ``(n,)`` or ``(n,k)``, respectively.
        """
        # check if _tmatmul was a) set as attribute in __init__ or b) overwritten by subclass to avoid infinite recursion
        if '_tmatmul' in self.__dict__ or type(self)._tmatmul != LinearOperator._tmatmul:
            return self._tmatmul(other.conjugate()).conjugate()
        else:
            AT = self.transpose()
            # check if transposed is not the default implementation to avoid infinite recursion
            if type(AT) is not _TransposeLinearOperator:
                return AT._matmul(other.conjugate()).conjugate()
            else:
                raise NotImplementedError("Please provide tmatmul or ctmatmul or overwrite _ctmatmul, _tmatmul, _transpose or _adjoint.")
            
    def _transpose(self) -> LinearOperator[Array]:
        """Assemble the transpose linear operator.
        
        This is the default implementation of the transpose linear operator ``self.T`` based on `_tmatmul`.
        
        Returns
        -------
        T : `LinearOperator`
             The transpose linear operator.
        """
        return _TransposeLinearOperator[Array](self)
    
    def _adjoint(self) -> LinearOperator[Array]:
        """Assemble the Hermitian adjoint linear operator.
        
        This is the default implementation of the Hermitian adjoint linear operator ``self.H`` based on `_ctmatmul`.
        
        Returns
        -------
        H : `LinearOperator`
             The Hermitian adjoint linear operator.
        """
        return _AdjointLinearOperator[Array](self)
    
    def _conjugate(self) -> LinearOperator[Array]:
        """Assemble the complex conjugate linear operator.
        
        This is the default implementation of the complex conjugate linear operator ``self.conjugate()`` based on `_cmatmul`.
        
        Returns
        -------
        C : `LinearOperator`
             The complex conjugate linear operator.
        """
        return _ConjugateLinearOperator[Array](self)
    
    def _get_dtype(self) -> np.dtype:
        """Determine dtype if not provided in constructor."""
        x = np.zeros(self.shape[-1], dtype=np.int8)
        try:
            Ax = self._matmul(x)
        except OverflowError:
            # Python large `int` promoted to `np.int64`or `np.int32`
            self.dtype = np.dtype(int)
        else:
            self.dtype = np.dtype(Ax.dtype)


##############################################################
# Helper classes for default implementations
##############################################################

######################################
# Adjoint, Transpose and Conjugate

class _TransposeLinearOperator[Array](LinearOperator[Array]):
    def __init__(self, linop: LinearOperator[Array]):
        LinearOperator.__init__(self, shape=(linop.shape[1], linop.shape[0]), dtype=linop.dtype)
        self._linop: LinearOperator[Array] = linop
        
    def _matmul(self, other: Array) -> Array:
        return self._linop._tmatmul(other)
        
    def _transpose(self) -> LinearOperator[Array]:
        return self._linop
    def _adjoint(self) -> LinearOperator[Array]:
        return self._linop.conjugate()
    def _conjugate(self) -> LinearOperator[Array]:
        return self._linop.adjoint()

class _AdjointLinearOperator[Array](LinearOperator[Array]):
    def __init__(self, linop: LinearOperator[Array]):
        LinearOperator.__init__(self, shape=(linop.shape[1], linop.shape[0]), dtype=linop.dtype)
        self._linop: LinearOperator[Array] = linop
        
    def _matmul(self, other: Array) -> Array:
        return self._linop._ctmatmul(other)
        
    def _transpose(self) -> LinearOperator[Array]:
        return self._linop.conjugate()
    def _adjoint(self) -> LinearOperator[Array]:
        return self._linop
    def _conjugate(self) -> LinearOperator[Array]:
        return self._linop.transpose()
 
class _ConjugateLinearOperator[Array](LinearOperator[Array]):
    def __init__(self, linop: LinearOperator[Array]):
        LinearOperator.__init__(self, shape=linop.shape, dtype=linop.dtype)
        self._linop: LinearOperator[Array] = linop
        
    def _matmul(self, other: Array) -> Array:
        return self._linop._cmatmul(other)
        
    def _transpose(self) -> LinearOperator[Array]:
        return self._linop.adjoint()
    def _adjoint(self) -> LinearOperator[Array]:
        return self._linop.transpose()
    def _conjugate(self) -> LinearOperator[Array]:
        return self._linop    
    
######################################
# Product
    
class _ProductLinearOperator[Array](LinearOperator[Array]):
    
    def __init__(self, *A: AsLinearOperatorProtocol[Array]):
        
        if len(A) < 2:
            raise ValueError(f'At least two linear operators required for product, got {len(A)}')
        
        A_ = []
        for a in A:
            if isinstance(a, _ProductLinearOperator):
                A_.extend(a._A)
            else:
                A_.append(aslinearoperator(a))
        A = A_
        
        for i in range(len(A) - 1):
             if A[i].shape[1] != A[i+1].shape[0]:
                raise ValueError(f'shapes {A[i].shape} and {A[i+1].shape} are not compatible for multiplication')
            
        shape = (A[0].shape[0], A[-1].shape[1])
        dtype = np.result_type(*[a.dtype for a in A])
        super().__init__(shape=shape, dtype=dtype)
        self._A: list[LinearOperator[Array]] = A
        
    def _matmul(self, other: Array) -> Array:
        result = other
        for A in reversed(self._A): result = A._matmul(result)
        return result
    
    def _tmatmul(self, other: Array) -> Array:
        result = other
        for A in self._A: result = A.transpose()._matmul(result)
        return result
    
    def _ctmatmul(self, other: Array) -> Array:
        result = other
        for A in self._A: result = A.adjoint()._matmul(result)
        return result
    
    def _cmatmul(self, other: Array) -> Array:
        result = other
        for A in reversed(self._A): result = A.conjugate()._matmul(result)
        return result

######################################
# Scaled

class _ScaledLinearOperator[Array](LinearOperator[Array]):
    
    def __init__(self, A: AsLinearOperatorProtocol[Array], alpha: Number):
        if isinstance(A, _ScaledLinearOperator):
            A = A._A
            alpha = alpha * A._alpha
        A = aslinearoperator(A)
        if not np.isscalar(alpha):
            raise TypeError(f'Scaling factor must be a number, got {type(alpha)}')
        super().__init__(shape=A.shape, dtype=np.result_type(alpha, A.dtype))
        self._A: LinearOperator[Array] = A
        self._alpha: Number = alpha
        
    def _matmul(self, other: Array) -> Array:
        return self._alpha * self._A._matmul(other)
    
    def _tmatmul(self, other: Array) -> Array:
        return self._alpha * self._A.transpose()._matmul(other)
    
    def _ctmatmul(self, other: Array) -> Array:
        return self._alpha.conjugate() * self._A.adjoint()._matmul(other)
    
    def _cmatmul(self, other: Array) -> Array:
        return self._alpha.conjugate() * self._A.conjugate()._matmul(other)

######################################
# Sum

class _SumLinearOperator[Array](LinearOperator[Array]):
    
    def __init__(self, *A: AsLinearOperatorProtocol[Array]):
        if len(A) < 2:
            raise ValueError(f'At least two linear operators required for addition, got {len(A)}')
        
        A_ = []
        for a in A:
            if isinstance(a, _SumLinearOperator):
                A_.extend(a._A)
            else:
                A_.append(aslinearoperator(a))
        A = A_
        
        shape = A[0].shape
        for i in range(1, len(A)):
             if A[i].shape != shape:
                raise ValueError(f'shapes {A[i].shape} and {shape} are not compatible for addition')
            
        dtype = np.result_type(*[a.dtype for a in A])
        super().__init__(shape=shape, dtype=dtype)
        self._A: list[LinearOperator[Array]] = A
        
    def _matmul(self, other: Array) -> Array:
        result = self._A[0]._matmul(other)
        for A in self._A[1:]: result = result + A._matmul(other)
        return result
    
    def _tmatmul(self, other: Array) -> Array:
        result = self._A[0].transpose()._matmul(other)
        for A in self._A[1:]: result = result + A.transpose()._matmul(other)
        return result
    
    def _ctmatmul(self, other: Array) -> Array:
        result = self._A[0].adjoint()._matmul(other)
        for A in self._A[1:]: result = result + A.adjoint()._matmul(other)
        return result
    
    def _cmatmul(self, other: Array) -> Array:
        result = self._A[0].conjugate()._matmul(other)
        for A in self._A[1:]: result = result + A.conjugate()._matmul(other)
        return result

##############################################################
# Identity Operator
##############################################################   

class IdentityOperator[Array](LinearOperator[Array]):
    """Identity operator.
    
    This is a subclass of `LinearOperator`, see the documentation there.
    
    Parameters
    ----------
    shape : tuple[int, int]
        The shape of the identity operator, must be a tuple of two positive integers.
    dtype : np.dtype, optional
        The data type of the identity operator. If not provided, it is set to `numpy.int8` by default as the smallest possible data type.
        
    Attributes
    ----------
    ndim : int
        Number of dimensions (this is always 2).
    shape : tuple[int, int]
        Matrix dimension ``(m,n)``.
    dtype : `numpy.dtype`
        Data type of the matrix.
    T : `LinearOperator`
        Transpose linear operator.
    H : `LinearOperator`
        Hermitian adjoint linear operator.
    C : `LinearOperator`
        Complex conjugate linear operator.
        
    See Also
    --------
    LinearOperator : base class for linear operators
    MatrixLinearOperator : implicit construction of a `LinearOperator` from a matrix-like object
    aslinearoperator : function to convert an object to a `LinearOperator`
    """
    def __init__(self, shape, dtype=np.int8):
        super().__init__(dtype, shape)

    def _matmul(self, other: Array) -> Array:
        return other
    
    def _tmatmul(self, other: Array) -> Array:
        return other
    
    def _ctmatmul(self, other: Array) -> Array:
        return other
    
    def _cmatmul(self, other: Array) -> Array:
        return other

##############################################################
# Matrix Protocol and MatrixLinearOperator
##############################################################       

######################################
# Matrix protocol

class MinimalMatrixProtocol[Array](Protocol):
    """Minimal matrix protocol for implicit construction of a `LinearOperator` from a matrix-like object."""
    
    ndim: int = 2
    shape: tuple[int, int]
    dtype: np.dtype
    
    def __matmul__(self, other: Array) -> Array:
        """Matrix-vector or matrix-matrix product.
        
        Calculates ``y = self @ other`` where ``other`` is an array-like object compatible with this matrix-like object, i.e. it must have shape ``(n,)`` or ``(n,k)`` if this matrix-like object has shape ``(m,n)``, resulting in the output ``y`` having shape ``(m,)`` or ``(m,k)``, respectively.
        
        Parameters
        ----------
        other : `ArrayProtocol`
            The array-like object to be multiplied. Must have shape ``(n,)`` or ``(n,k)`` if this matrix-like object has shape ``(m,n)``.
            
        Returns
        -------
        y : `ArrayProtocol`
            The result of the matrix-vector or matrix-matrix product, having shape ``(m,)`` or ``(m,k)``, respectively.
        """
        ...
    
    def transpose(self) -> MinimalMatrixProtocol[Array]: ...
    
class MatrixProtocol[Array](MinimalMatrixProtocol[Array]):
    """Matrix protocol for implicit construction of a `LinearOperator` from a matrix-like object.
    
    Attributes
    ----------
    ndim : int
        Number of dimensions (this is always 2).
    shape : tuple[int, int]
        Matrix dimension ``(m,n)``.
    dtype : `numpy.dtype`
        Data type of the matrix.
    """
    
    def transpose(self) -> MatrixProtocol[Array]:
        """Transpose matrix-like object.
        
        Returns
        -------
        T : `MatrixProtocol`
            The transpose matrix-like object.
        """
        ...
    
    def conjugate(self) -> MatrixProtocol[Array]:
        """Optional: Complex conjugate matrix-like object.
        
        If this method is not present, the default implementations of the Hermitian adjoint and complex conjugate are used as described in the `LinearOperator` class.
        
        Returns
        -------
        C : `MatrixProtocol`
            The complex conjugate matrix-like object.
        """
        ...

######################################
# MatrixLinearOperator

class MatrixLinearOperator[Array](LinearOperator[Array]):
    """Linear operator implicitly defined by a matrix-like object.
    
    This is a subclass of `LinearOperator`, see the documentation there.
    
    Parameters
    ----------
    M : `MatrixProtocol`
        The matrix-like object defining this linear operator. Must satisfy the `MatrixProtocol`.
        
    Attributes
    ----------
    ndim : int
        Number of dimensions (this is always 2).
    shape : tuple[int, int]
        Matrix dimension ``(m,n)``.
    dtype : `numpy.dtype`
        Data type of the matrix.
    T : `LinearOperator`
        Transpose linear operator.
    H : `LinearOperator`
        Hermitian adjoint linear operator.
    C : `LinearOperator`
        Complex conjugate linear operator.
        
    See Also
    --------
    LinearOperator : base class for linear operators
    aslinearoperator : function to convert an object to a `LinearOperator`
    """
    
    _MT: MinimalMatrixProtocol[Array] | None = None
    _MH: MatrixProtocol[Array] | False | None = None
    _MC: MatrixProtocol[Array] | False | None = None
    
    def __init__(self, M: MinimalMatrixProtocol[Array]):
        
        if not all(hasattr(M, attr) for attr in ['ndim', 'shape', 'dtype', '__matmul__', 'transpose']):
            raise TypeError(f'Object of type {type(M)} does not satisfy the matrix protocol, missing one of ndim, shape, dtype, __matmul__ or transpose.')
        
        if M.ndim != 2:
            raise ValueError(f'matrix must be 2D but got ndim = {M.ndim}')
        
        super().__init__(shape=M.shape, dtype=M.dtype)
        self._M: MinimalMatrixProtocol[Array] = M
        
    def _M_transpose(self) -> MinimalMatrixProtocol[Array]:
        if self._MT is None:
            self._MT = self._M.transpose()
        return self._MT
    
    def _M_adjoint(self) -> MatrixProtocol[Array] | False:
        if self._MH is None:
            MC = self._M_conjugate()
            if MC is not False and hasattr(MC, 'transpose'):
                self._MH = MC.transpose()
            else:
                self._MH = False
        return self._MH
    
    def _M_conjugate(self) -> MatrixProtocol[Array] | False:
        if self._MC is None:
            if hasattr(self._M, 'conjugate'):
                self._MC = self._M.conjugate()
            else:
                self._MC = False
        return self._MC
    
    def _matmul(self, x: Array) -> Array:
        return self._M @ x
    
    def _tmatmul(self, x: Array) -> Array:
        return self._M_transpose() @ x
    
    def _ctmatmul(self, x: Array) -> Array:
        MH = self._M_adjoint()
        return super()._ctmatmul(x) if MH is False else MH @ x
        
    def _cmatmul(self, x: Array) -> Array:
        MC = self._M_conjugate()
        return super()._cmatmul(x) if MC is False else MC @ x


##############################################################
# Conversion functionality: aslinearoperator
##############################################################

######################################
# Conversion protocols

class AsLinearOperatorDunderProtocol[Array](Protocol):
    """Protocol for objects that provide an explicit instruction on how to convert them to a `LinearOperator` via the `aslinearoperator` function.
    """
    
    def __aslinearoperator__(self) -> LinearOperator[Array]:
        """Convert this object to a `LinearOperator`.
        
        Returns
        -------
        linop : `LinearOperator`
            The linear operator representation of this object.
        """
        ...
    

type AsLinearOperatorProtocol[Array] =  LinearOperator[Array] | MinimalMatrixProtocol[Array] | AsLinearOperatorDunderProtocol[Array]
"""All objects that can be converted to a `LinearOperator` via the `aslinearoperator` function, namely `LinearOperator`, `MatrixProtocol` and `AsLinearOperatorDunderProtocol`."""

######################################
# aslinearoperator

def aslinearoperator[Array](A: AsLinearOperatorProtocol[Array]) -> LinearOperator[Array]:
    """Convert an object to a `LinearOperator`.
    
    Parameters
    ----------
    A : `LinearOperator`, `AsLinearOperatorDunderProtocol`, `MatrixProtocol`
        The object to be converted to a `LinearOperator`. Can be a
        
        - `LinearOperator`: no conversion needed
        - explicit conversion: any object that satisfies the `AsLinearOperatorDunderProtocol`, i.e. that implements the ``__aslinearoperator__()`` method
        - implicit conversion: any matrix-like object that satisfies the `MatrixProtocol`, e.g. `numpy.ndarray` or `scipy.sparse.sparray`.
        
    Returns
    -------
    linop : `LinearOperator`
        The linear operator representation of the input object. If the input was a `numpy.ndarray` or `scipy.sparse.sparray`, the output is a `LinearOperator` compatible with arrays of type `numpy.ndarray` and `scipy.sparse.sparray`, i.e. ``LinearOperator[np.ndarray | sp.sparse.sparray]``.
        
    See Also
    --------
    LinearOperator : base class for linear operators
    MatrixLinearOperator : implicit construction of a `LinearOperator` from a matrix-like object
    """
    
    if isinstance(A, LinearOperator):
        return A
    
    if hasattr(A, '__aslinearoperator__'):
        linop = A.__aslinearoperator__()
        if not isinstance(linop, LinearOperator):
            raise TypeError(f'__aslinearoperator__() did not return a LinearOperator, got {type(linop)}')
        return linop
    
    if isinstance(A, np.ndarray | np.matrix) or issparse(A):
        return MatrixLinearOperator[np.ndarray | sparray](A)
    
    else:
        try:
            return MatrixLinearOperator(A)
        except TypeError as e:
            if str(e).startswith('Object of type') and 'does not satisfy the matrix protocol' in str(e):
                raise TypeError(f'type {type(A)} not understood, expected {LinearOperator|np.ndarray|sparray|AsLinearOperatorDunderProtocol|MinimalMatrixProtocol}') from e
            else:
                raise