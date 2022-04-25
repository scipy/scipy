"""
Sparse array creation functions e.g. `eye`, `diags`, etc.
"""

from . import _construct
from ._arrays import (
    bsr_array,
    coo_array,
    csc_array,
    csr_array,
    dia_array,
    dok_array,
    lil_array,
    _matrix_doc_to_array,
)

__all__ = [
    "block_diag",
    "bmat",
    "diags",
    "eye",
    "hstack",
    "identity",
    "kron",
    "kronsum",
    "random",
    "spdiags",
    "vstack",
]


def __dir__():
    return __all__


# Map sparse format strings to the corresponding sparse array container
_fmt_to_sparray = {
    "bsr": bsr_array, "coo": coo_array, "csc": csc_array, "csr": csr_array,
    "dia": dia_array, "dok": dok_array, "lil": lil_array,
}


def block_diag(mats, format=None, dtype=None):
    M = _construct.block_diag(mats, format, dtype)
    return _fmt_to_sparray[M.format](M)


def bmat(blocks, format=None, dtype=None):
    M = _construct.bmat(blocks, format, dtype)
    return _fmt_to_sparray[M.format](M)


def diags(diagonals, offsets=0, shape=None, format=None, dtype=None):
    M = _construct.diags(diagonals, offsets, shape, format, dtype)
    return _fmt_to_sparray[M.format](M)


def eye(m, n=None, k=0, dtype=float, format=None):
    M = _construct.eye(m, n, k, dtype, format)
    return _fmt_to_sparray[M.format](M)


def hstack(blocks, format=None, dtype=None):
    M = _construct.hstack(blocks, format, dtype)
    return _fmt_to_sparray[M.format](M)


def identity(n, dtype='d', format=None):
    M = _construct.identity(n, dtype, format)
    return _fmt_to_sparray[M.format](M)


def kron(A, B, format=None):
    M = _construct.kron(A, B, format)
    return _fmt_to_sparray[M.format](M)


def kronsum(A, B, format=None):
    M = _construct.kronsum(A, B, format)
    return _fmt_to_sparray[M.format](M)


def random(m, n, density=0.01, format='coo', dtype=None, random_state=None,
           data_rvs=None):
    M = _construct.random(m, n, density, format, dtype, random_state, data_rvs)
    return _fmt_to_sparray[M.format](M)


def spdiags(data, diags, m, n, format=None):
    M = _construct.spdiags(data, diags, m, n, format)
    return _fmt_to_sparray[M.format](M)


def vstack(blocks, format=None, dtype=None):
    M = _construct.vstack(blocks, format, dtype)
    return _fmt_to_sparray[M.format](M)


block_diag.__doc__ = _matrix_doc_to_array(_construct.block_diag.__doc__)
bmat.__doc__ = _matrix_doc_to_array(_construct.bmat.__doc__)
diags.__doc__ = _matrix_doc_to_array(_construct.diags.__doc__)
eye.__doc__ = _matrix_doc_to_array(_construct.eye.__doc__)
hstack.__doc__ = _matrix_doc_to_array(_construct.hstack.__doc__)
identity.__doc__ = _matrix_doc_to_array(_construct.identity.__doc__)
kron.__doc__ = _matrix_doc_to_array(_construct.kron.__doc__)
kronsum.__doc__ = _matrix_doc_to_array(_construct.kronsum.__doc__)
random.__doc__ = _matrix_doc_to_array(_construct.random.__doc__)
spdiags.__doc__ = _matrix_doc_to_array(_construct.spdiags.__doc__)
vstack.__doc__ = _matrix_doc_to_array(_construct.vstack.__doc__)
