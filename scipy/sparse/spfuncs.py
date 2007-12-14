""" Functions that operate on sparse matrices
"""

from numpy import empty

from base import isspmatrix
from csr import isspmatrix_csr
from csc import isspmatrix_csc
from sputils import upcast

import sparsetools

def extract_diagonal(A):
    """
    extract_diagonal(A) returns the main diagonal of A.
    """
    #TODO extract k-th diagonal
    if isspmatrix_csr(A) or isspmatrix_csc(A):
        fn = getattr(sparsetools, "extract_" + A.format + "_diagonal")
        y = empty( min(A.shape), dtype=upcast(A.dtype) )
        fn(A.shape[0],A.shape[1],A.indptr,A.indices,A.data,y)
        return y
    elif isspmatrix(A):
        return extract_diagonal(A.tocsr())
    else:
        raise ValueError,'expected sparse matrix'



