import scipy as sp
import pytest


def test_array_api_deprecations():
    X = sp.sparse.csr_array([
        [1,2,3],
        [4,0,6]
    ])
    msg = "1.13.0"

    with pytest.deprecated_call(match=msg):
        X.get_shape()

    with pytest.deprecated_call(match=msg):
        X.set_shape((2,3))

    with pytest.deprecated_call(match=msg):
        X.asfptype()

    with pytest.deprecated_call(match=msg):
        X.getmaxprint()

    with pytest.deprecated_call(match=msg):
        X.getnnz()

    with pytest.deprecated_call(match=msg):
        X.getH()

    with pytest.deprecated_call(match=msg):
        X.getcol(1).todense()

    with pytest.deprecated_call(match=msg):
        X.getrow(1).todense()


def test_isspmatrix_deprecations():
    msg = "1.13.0"

    X = sp.sparse.csr_matrix([[1, 0], [0, 1]])
    with pytest.deprecated_call(match=msg):
        sp.sparse.isspmatrix(X)

    X = sp.sparse.bsr_matrix([[1, 0], [0, 1]])
    with pytest.deprecated_call(match=msg):
        sp.sparse.isspmatrix_bsr(X)

    X = sp.sparse.coo_matrix([[1, 0], [0, 1]])
    with pytest.deprecated_call(match=msg):
        sp.sparse.isspmatrix_coo(X)

    X = sp.sparse.csc_matrix([[1, 0], [0, 1]])
    with pytest.deprecated_call(match=msg):
        sp.sparse.isspmatrix_csc(X)

    X = sp.sparse.csr_matrix([[1, 0], [0, 1]])
    with pytest.deprecated_call(match=msg):
        sp.sparse.isspmatrix_csr(X)

    X = sp.sparse.dia_matrix([[1, 0], [0, 1]])
    with pytest.deprecated_call(match=msg):
        sp.sparse.isspmatrix_dia(X)

    X = sp.sparse.dok_matrix([[1, 0], [0, 1]])
    with pytest.deprecated_call(match=msg):
        sp.sparse.isspmatrix_dok(X)

    X = sp.sparse.lil_matrix([[1, 0], [0, 1]])
    with pytest.deprecated_call(match=msg):
        sp.sparse.isspmatrix_lil(X)
