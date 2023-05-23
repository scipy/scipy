from ._sputils import isintlike, isscalarlike


class spmatrix:
    """This class provides a base class for all sparse matrix classes.

    It cannot be instantiated.  Most of the work is provided by subclasses.
    """
    _is_array = False

    @property
    def _bsr_container(self):
        from ._bsr import bsr_matrix
        return bsr_matrix

    @property
    def _coo_container(self):
        from ._coo import coo_matrix
        return coo_matrix

    @property
    def _csc_container(self):
        from ._csc import csc_matrix
        return csc_matrix

    @property
    def _csr_container(self):
        from ._csr import csr_matrix
        return csr_matrix

    @property
    def _dps_container(self):
        from ._dps import dps_matrix
        return dps_matrix

    @property
    def _dia_container(self):
        from ._dia import dia_matrix
        return dia_matrix

    @property
    def _dok_container(self):
        from ._dok import dok_matrix
        return dok_matrix

    @property
    def _lil_container(self):
        from ._lil import lil_matrix
        return lil_matrix

    # Restore matrix multiplication
    def __mul__(self, other):
        return self._mul_dispatch(other)

    def __rmul__(self, other):
        return self._rmul_dispatch(other)

    # Restore matrix power
    def __pow__(self, other):
        M, N = self.shape
        if M != N:
            raise TypeError('sparse matrix is not square')

        if isintlike(other):
            other = int(other)
            if other < 0:
                raise ValueError('exponent must be >= 0')

            if other == 0:
                from ._construct import eye
                return eye(M, dtype=self.dtype)

            if other == 1:
                return self.copy()

            tmp = self.__pow__(other // 2)
            if other % 2:
                return self @ tmp @ tmp
            else:
                return tmp @ tmp

        if isscalarlike(other):
            raise ValueError('exponent must be an integer')
        return NotImplemented


def _array_doc_to_matrix(docstr):
    # For opimized builds with stripped docstrings
    if docstr is None:
        return None
    return (
        docstr.replace('sparse arrays', 'sparse matrices')
              .replace('sparse array', 'sparse matrix')
    )
