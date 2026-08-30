import numpy as np
import pytest

from scipy.sparse import (
        csr_array,
        # for deprecation testing
        bsr_matrix, csc_matrix, dia_matrix, lil_matrix,
        coo_matrix, csr_matrix, dok_matrix,
        isspmatrix_bsr, isspmatrix_csc, isspmatrix_dia, isspmatrix_lil,
        isspmatrix_coo, isspmatrix_csr, isspmatrix_dok,
        spdiags, diags, identity, eye, bmat, rand, random,
    )

# Test deprecations of spmatrix items
def test_isspmatrix_functions_deprecated():
    A = csr_array([[1, 0], [0, 1]])
    with pytest.deprecated_call(match=".isspmatrix_bsr. is being replaced"):
        isspmatrix_bsr(A)
    with pytest.deprecated_call(match=".isspmatrix_coo. is being replaced"):
        isspmatrix_coo(A)
    with pytest.deprecated_call(match=".isspmatrix_csc. is being replaced"):
        isspmatrix_csc(A)
    with pytest.deprecated_call(match=".isspmatrix_csr. is being replaced"):
        isspmatrix_csr(A)
    with pytest.deprecated_call(match=".isspmatrix_dia. is being replaced"):
        isspmatrix_dia(A)
    with pytest.deprecated_call(match=".isspmatrix_dok. is being replaced"):
        isspmatrix_dok(A)
    with pytest.deprecated_call(match=".isspmatrix_lil. is being replaced"):
        isspmatrix_lil(A)


def test_class_constructors_deprecated():
    A = np.array([[1, 0], [0, 1]])
    with pytest.deprecated_call(match=".*_matrix is being replaced"):
        bsr_matrix(A)
    with pytest.deprecated_call(match=".*_matrix is being replaced"):
        coo_matrix(A)
    with pytest.deprecated_call(match=".*_matrix is being replaced"):
        csc_matrix(A)
    with pytest.deprecated_call(match=".*_matrix is being replaced"):
        csr_matrix(A)
    with pytest.deprecated_call(match=".*_matrix is being replaced"):
        dia_matrix(A)
    with pytest.deprecated_call(match=".*_matrix is being replaced"):
        dok_matrix(A)
    with pytest.deprecated_call(match=".*_matrix is being replaced"):
        lil_matrix(A)


def test_construct_deprecated():
    # removed functions
    with pytest.deprecated_call(match=".* is being replaced"):
        spdiags([[1, 0, 0, 2], [0, 2, 2, 1]], [-1, 1])
    with pytest.deprecated_call(match=".* is being replaced"):
        diags([[1, 0, 0, 2], [0, 2, 2, 1]], [-1, 1], dtype=None)
    with pytest.deprecated_call(match=".* is being replaced"):
        identity(5)
    with pytest.deprecated_call(match=".* is being replaced"):
        eye(5)
    # construction of A, B also has a deprecation warning
    with pytest.deprecated_call(match=".* is being replaced"):
        A = B = coo_matrix([[1, 0], [0, 1]])
    with pytest.deprecated_call(match=".* is being replaced"):
        bmat([[A, B]])
    with pytest.deprecated_call(match=".* is being replaced"):
        rand(3, 4, density=0.2)
    with pytest.deprecated_call(match=".* is being replaced"):
        random(3, 4, density=0.2)
