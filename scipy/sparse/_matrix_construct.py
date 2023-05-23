"""Functions to construct sparse matrices
"""

__all__ = ['spdiags', 'eye', 'identity', 'kron', 'kronsum',
           'hstack', 'vstack', 'bmat', 'rand', 'random', 'diags', 'block_diag']

from . import _construct
from ._bsr import bsr_matrix
from ._coo import coo_matrix
from ._csc import csc_matrix
from ._csr import csr_matrix
from ._dia import dia_matrix
from ._dok import dok_matrix
from ._lil import lil_matrix
from ._matrix import _array_doc_to_matrix

_fmt_to_spmatrix = {
    "bsr": bsr_matrix, "coo": coo_matrix, "csc": csc_matrix, "csr": csr_matrix,
    "dia": dia_matrix, "dok": dok_matrix, "lil": lil_matrix,
}


def block_diag(mats, format=None, dtype=None):
    A = _construct.block_diag(mats, format, dtype)
    return _fmt_to_spmatrix[A.format](A)


def bmat(blocks, format=None, dtype=None):
    A = _construct.bmat(blocks, format, dtype)
    return _fmt_to_spmatrix[A.format](A)


def diags(diagonals, offsets=0, shape=None, format=None, dtype=None):
    A = _construct.diags(diagonals, offsets, shape, format, dtype)
    return _fmt_to_spmatrix[A.format](A)


def eye(m, n=None, k=0, dtype=float, format=None):
    A = _construct.eye(m, n, k, dtype, format)
    return _fmt_to_spmatrix[A.format](A)


def hstack(blocks, format=None, dtype=None):
    A = _construct.hstack(blocks, format, dtype)
    return _fmt_to_spmatrix[A.format](A)


def identity(n, dtype='d', format=None):
    A = _construct.identity(n, dtype, format)
    return _fmt_to_spmatrix[A.format](A)


def kron(A, B, format=None):
    A = _construct.kron(A, B, format)
    return _fmt_to_spmatrix[A.format](A)


def kronsum(A, B, format=None):
    A = _construct.kronsum(A, B, format)
    return _fmt_to_spmatrix[A.format](A)


def random(m, n, density=0.01, format='coo', dtype=None, random_state=None,
           data_rvs=None):
    A = _construct.random(m, n, density, format, dtype, random_state, data_rvs)
    return _fmt_to_spmatrix[A.format](A)


def rand(m, n, density=0.01, format="coo", dtype=None, random_state=None):
    A = _construct.rand(m, n, density, format, dtype, random_state)
    return _fmt_to_spmatrix[A.format](A)


def spdiags(data, diags, m, n, format=None):
    A = _construct.spdiags(data, diags, m, n, format)
    return _fmt_to_spmatrix[A.format](A)


def vstack(blocks, format=None, dtype=None):
    A = _construct.vstack(blocks, format, dtype)
    return _fmt_to_spmatrix[A.format](A)

block_diag.__doc__ = _array_doc_to_matrix(_construct.block_diag.__doc__)
bmat.__doc__ = _array_doc_to_matrix(_construct.bmat.__doc__)
diags.__doc__ = _array_doc_to_matrix(_construct.diags.__doc__)
eye.__doc__ = _array_doc_to_matrix(_construct.eye.__doc__)
hstack.__doc__ = _array_doc_to_matrix(_construct.hstack.__doc__)
identity.__doc__ = _array_doc_to_matrix(_construct.identity.__doc__)
kron.__doc__ = _array_doc_to_matrix(_construct.kron.__doc__)
kronsum.__doc__ = _array_doc_to_matrix(_construct.kronsum.__doc__)
random.__doc__ = _array_doc_to_matrix(_construct.random.__doc__)
spdiags.__doc__ = _array_doc_to_matrix(_construct.spdiags.__doc__)
vstack.__doc__ = _array_doc_to_matrix(_construct.vstack.__doc__)
