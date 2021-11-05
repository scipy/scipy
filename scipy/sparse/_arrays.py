from ._bsr import bsr_matrix
from ._coo import coo_matrix
from ._csc import csc_matrix
from ._csr import csr_matrix
from ._dia import dia_matrix
from ._dok import dok_matrix
from ._lil import lil_matrix


class sparray:
    _is_array = True

    @property
    def _bsr_container(self):
        return bsr_array

    @property
    def _coo_container(self):
        return coo_array

    @property
    def _csc_container(self):
        return csc_array

    @property
    def _csr_container(self):
        return csr_array

    @property
    def _dia_container(self):
        return dia_array

    @property
    def _dok_container(self):
        return dok_array

    @property
    def _lil_container(self):
        return lil_array

    # Restore elementwise multiplication
    def __mul__(self, *args, **kwargs):
        return self.multiply(*args, **kwargs)


class bsr_array(sparray, bsr_matrix):
    pass


class coo_array(sparray, coo_matrix):
    pass


class csc_array(sparray, csc_matrix):
    pass


class csr_array(sparray, csr_matrix):
    pass


class dia_array(sparray, dia_matrix):
    pass


class dok_array(sparray, dok_matrix):
    pass


class lil_array(sparray, lil_matrix):
    pass
