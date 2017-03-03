"""
sparsetools is not a public module in scipy.sparse, but this file is
for backward compatibility if someone happens to use it.
"""
from numpy import deprecate

# This file shouldn't be imported by scipy --- Scipy code should use
# internally scipy.sparse._sparsetools


@deprecate(old_name="scipy.sparse.sparsetools",
           message=("scipy.sparse.sparsetools is a private module for scipy.sparse, "
                    "and should not be used."))
def _deprecated():
    pass

del deprecate

try:
    _deprecated()
except DeprecationWarning as e:
    # don't fail import if DeprecationWarnings raise error -- works around
    # the situation with Numpy's test framework
    pass

from ._sparsetools import *
import ._sparsetools

# the signatures of the matvec functions have changed in 0.18  These python
# wrappers prevent code that depended on their earlier signature from breaking.
# Note that coo_matvec is not wrapped as it now also needs to know the number
# of rows.

def csr_matvec(N, M, intptr, indices, data, x, y):
    _sparsetools.csr_matvec(N, M, indptr, indices, data, 1, x, 1, y)

def csr_matvecs(N, M, n_vecs, intptr, indices, data, x, y):
    _sparsetools.csr_matvecs(N, M, n_vecs, indptr, indices, data, 1, x, 1, y)

def csc_matvec(N, M, intptr, indices, data, x, y):
    _sparsetools.csc_matvec(N, M, indptr, indices, data, 1, x, 1, y)

def csc_matvecs(N, M, n_vecs, intptr, indices, data, x, y):
    _sparsetools.csc_matvecs(N, M, n_vecs, indptr, indices, data, 1, x, 1, y)

def bsr_matvec(N, M, R, C, intptr, indices, data, x, y):
    _sparsetools.bsr_matvec(N, M, R, C, indptr, indices, data, 1, x, 1, y)

def bsr_matvecs(N, M, n_vecs, R, C, intptr, indices, data, x, y):
    _sparsetools.bsr_matvecs(N, M, n_vecs, R, C, indptr, indices, data, 1, x, 1, y)

def dia_matvec(N, M, len_off, L, offsets, data, x, y):
    _sparsetools.dia_matvec(N, M, indptr, indices, data, 1, x, 1, y)
