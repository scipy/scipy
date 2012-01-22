"""Base class for sparse matrice with a .data attribute

    subclasses must provide a _with_data() method that
    creates a new matrix with the same sparsity pattern
    as self but with a different data array

"""

__all__ = []

import types
import numpy as np

from base import spmatrix
from sputils import isscalarlike


#TODO implement all relevant operations
#use .data.__methods__() instead of /=, *=, etc.
class _data_matrix(spmatrix):
    def __init__(self):
        spmatrix.__init__(self)

        # Add the numpy unary ufuncs for which func(0) = 0
        # to this class.
        for npfunc in [np.sin, np.tan, np.arcsin, np.arctan,
                       np.sinh, np.tanh, np.arcsinh, np.arctanh,
                       np.rint, np.sign, np.expm1, np.log1p,
                       np.deg2rad, np.rad2deg, np.floor, np.ceil,
                       np.trunc]:
            name = npfunc.__name__

            def create_func(op):
                def func(self):
                    result = op(self.data)
                    x = self._with_data(result, copy=True)
                    return x
                func.__doc__ = ("Element-wise %s\n\nSee numpy.%s "
                                "for more information." % (name, name))
                return func
            setattr(self, name, types.MethodType(create_func(npfunc),
                                                 self, self.__class__))

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

    def __imul__(self, other): #self *= other
        if isscalarlike(other):
            self.data *= other
            return self
        else:
            raise NotImplementedError

    def __itruediv__(self, other): #self /= other
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
