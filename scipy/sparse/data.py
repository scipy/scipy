"""Base class for sparse matrice with a .data attribute
    
    subclasses must provide a _with_data() method that
    creates a new matrix with the same sparsity pattern
    as self but with a different data array

"""

__all__ = []

from base import spmatrix

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
        return self._with_data(numpy.real(self.data),copy=False)

    def _imag(self):
        return self._with_data(numpy.imag(self.data),copy=False)
    
    def __neg__(self):
        return self._with_data(-self.data)

    def astype(self, t):
        return self._with_data(self.data.astype(t))

    def conj(self, copy=False):
        return self._with_data(self.data.conj(),copy=copy)
    
    def copy(self):
        return self._with_data(self.data.copy(),copy=True)

