"""Base class for sparse matrice with a .data attribute

    subclasses must provide a _with_data() method that
    creates a new matrix with the same sparsity pattern
    as self but with a different data array

"""

from __future__ import division, print_function, absolute_import

__all__ = []

import numpy as np

from .base import spmatrix
from .sputils import isscalarlike


# TODO implement all relevant operations
# use .data.__methods__() instead of /=, *=, etc.
class _data_matrix(spmatrix):
    def __init__(self):
        spmatrix.__init__(self)

    def _get_dtype(self):
        return self.data.dtype

    def _set_dtype(self,newtype):
        self.data.dtype = newtype
    dtype = property(fget=_get_dtype,fset=_set_dtype)

    def __abs__(self):
        return self._with_data(abs(self.data))

    def _real(self):
        return self._with_data(self.data.real)

    def _imag(self):
        return self._with_data(self.data.imag)

    def __neg__(self):
        return self._with_data(-self.data)

    def __imul__(self, other):  # self *= other
        if isscalarlike(other):
            self.data *= other
            return self
        else:
            raise NotImplementedError

    def __itruediv__(self, other):  # self /= other
        if isscalarlike(other):
            recip = 1.0 / other
            self.data *= recip
            return self
        else:
            raise NotImplementedError

    def astype(self, t):
        return self._with_data(self.data.astype(t))

    def conj(self):
        return self._with_data(self.data.conj())

    def copy(self):
        return self._with_data(self.data.copy(), copy=True)

    ###########################
    # Multiplication handlers #
    ###########################

    def _mul_scalar(self, other):
        return self._with_data(self.data * other)


# Add the numpy unary ufuncs for which func(0) = 0 to _data_matrix.
for npfunc in [np.sin, np.tan, np.arcsin, np.arctan, np.sinh, np.tanh,
               np.arcsinh, np.arctanh, np.rint, np.sign, np.expm1, np.log1p,
               np.deg2rad, np.rad2deg, np.floor, np.ceil, np.trunc, np.sqrt]:
    name = npfunc.__name__

    def _create_method(op):
        def method(self):
            result = op(self.data)
            x = self._with_data(result, copy=True)
            return x

        method.__doc__ = ("Element-wise %s.\n\n"
                          "See numpy.%s for more information." % (name, name))
        method.__name__ = name

        return method

    setattr(_data_matrix, name, _create_method(npfunc))


class _minmax_mixin(object):
    """Mixin for min and max methods.

    These are not implemented for dia_matrix, hence the separate class.
    """

    def max(self):
        """Maximum of the elements of this matrix.

        This takes all elements into account, not just the non-zero ones.

        Returns
        -------
        amax : self.dtype
            Maximum element.
        """
        zero = self.dtype.type(0)
        if self.nnz == 0:
            return zero
        mx = np.max(self.data)
        if self.nnz != np.product(self.shape):
            mx = max(zero, mx)
        return mx

    def min(self):
        """Minimum of the elements of this matrix.

        This takes all elements into account, not just the non-zero ones.

        Returns
        -------
        amin : self.dtype
            Minimum element.
        """
        zero = self.dtype.type(0)
        if self.nnz == 0:
            return zero
        mn = np.min(self.data)
        if self.nnz != np.product(self.shape):
            mn = min(zero, mn)
        return mn
