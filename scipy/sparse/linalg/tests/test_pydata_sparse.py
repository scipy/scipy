import pytest

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as splin

from numpy.testing import assert_allclose

try:
    import sparse
except ImportError:
    sparse = None

pytestmark = pytest.mark.skipif(sparse is None,
                                reason="pydata/sparse not installed")


def test_isolve_gmres():
    # Several of the iterative solvers use the same
    # isolve.utils.make_system wrapper code, so test just one of them.
    np.random.seed(1234)

    A = sparse.COO(np.random.rand(5, 5))
    b = np.random.rand(5)

    x, info = splin.gmres(A, b, atol=1e-15)
    assert info == 0
    assert isinstance(x, np.ndarray)
    assert_allclose(A @ x, b)


def test_lsmr():
    np.random.seed(1234)

    A_dense = np.random.rand(5, 5)
    A = sparse.COO(A_dense)
    b = np.random.rand(5)

    res0 = splin.lsmr(A_dense, b)
    res = splin.lsmr(A, b)
    assert_allclose(res[0], res0[0])


def test_lsqr():
    np.random.seed(1234)

    A_dense = np.random.rand(5, 5)
    A = sparse.COO(A_dense)
    b = np.random.rand(5)

    res0 = splin.lsqr(A_dense, b)
    res = splin.lsqr(A, b)
    assert_allclose(res[0], res0[0])


def test_eigs():
    np.random.seed(1234)

    A_dense = np.random.rand(10, 10)
    A_dense = A_dense @ A_dense.T
    A = sparse.COO(A_dense)

    M_dense = np.random.rand(10, 10)
    M_dense = M_dense @ M_dense.T
    M = sparse.COO(M_dense)

    v0 = np.random.rand(10)

    w_dense, v_dense = splin.eigs(A_dense, k=3, v0=v0)
    w, v = splin.eigs(A, k=3, v0=v0)

    assert_allclose(w, w_dense)
    assert_allclose(v, v_dense)

    w_dense, v_dense = splin.eigs(A_dense, M=M_dense, k=3, v0=v0)
    w, v = splin.eigs(A, M=M_dense, k=3, v0=v0)

    assert_allclose(w, w_dense)
    assert_allclose(v, v_dense)

    w_dense, v_dense = splin.eigsh(A_dense, M=M_dense, k=3, v0=v0)
    w, v = splin.eigsh(A, M=M_dense, k=3, v0=v0)

    assert_allclose(w, w_dense)
    assert_allclose(v, v_dense)


def test_svds():
    np.random.seed(1234)

    A_dense = np.random.rand(10, 10)
    A_dense = A_dense @ A_dense.T
    A = sparse.COO(A_dense)

    v0 = np.random.rand(10)

    u0, s0, vt0 = splin.svds(A_dense, k=3, v0=v0)
    u, s, vt = splin.svds(A, k=3, v0=v0)

    assert_allclose(s, s0)
    assert_allclose(u, u0)
    assert_allclose(vt, vt0)


def test_lobpcg():
    np.random.seed(1234)

    A_dense = np.random.rand(10, 10)
    A_dense = A_dense @ A_dense.T
    A = sparse.COO(A_dense)

    X = np.random.rand(10, 1)

    w_dense, v_dense = splin.lobpcg(A_dense, X)
    w, v = splin.lobpcg(A, X)

    assert_allclose(w, w_dense)
    assert_allclose(v, v_dense)


def test_spsolve():
    np.random.seed(1234)

    A_dense = np.random.rand(10, 10) + 10*np.eye(10)
    A = sparse.COO(A_dense)

    b = np.random.rand(10)
    b2 = np.random.rand(10, 3)

    x0 = splin.spsolve(sp.csc_matrix(A_dense), b)
    x = splin.spsolve(A, b)
    assert isinstance(x, np.ndarray)
    assert_allclose(x, x0)

    x0 = splin.spsolve(sp.csc_matrix(A_dense), b)
    x = splin.spsolve(A, b, use_umfpack=True)
    assert isinstance(x, np.ndarray)
    assert_allclose(x, x0)

    x0 = splin.spsolve(sp.csc_matrix(A_dense), b2)
    x = splin.spsolve(A, b2)
    assert isinstance(x, np.ndarray)
    assert_allclose(x, x0)

    x0 = splin.spsolve(sp.csc_matrix(A_dense),
                       sp.csc_matrix(A_dense))
    x = splin.spsolve(A, A)
    assert isinstance(x, type(A))
    assert_allclose(x.todense(), x0.todense())


def test_splu():
    # NB. spilu follows same code paths, so no need to test separately
    np.random.seed(1234)

    A_dense = np.random.rand(10, 10) + 10*np.eye(10)
    A = sparse.COO(A_dense)
    lu = splin.splu(A)

    assert isinstance(lu.L, sparse.COO)
    assert isinstance(lu.U, sparse.COO)

    Pr = sparse.COO(sp.csc_matrix((np.ones(10), (lu.perm_r, np.arange(10)))))
    Pc = sparse.COO(sp.csc_matrix((np.ones(10), (np.arange(10), lu.perm_c))))
    A2 = Pr.T @ lu.L @ lu.U @ Pc.T

    assert_allclose(A2.todense(), A.todense())

    z = lu.solve(A.todense())
    assert_allclose(z, np.eye(10), atol=1e-10)


def test_spsolve_triangular():
    np.random.seed(1234)

    A_dense = np.tril(np.random.rand(10, 10) + 10*np.eye(10))
    A = sparse.COO(A_dense)
    b = np.random.rand(10)

    x = splin.spsolve_triangular(A, b)
    assert_allclose(A @ x, b)


def test_onenormest():
    np.random.seed(1234)

    A_dense = np.random.rand(10, 10)
    A = sparse.COO(A_dense)

    est0 = splin.onenormest(A_dense)
    est = splin.onenormest(A)
    assert_allclose(est, est0)


def test_inv():
    np.random.seed(1234)

    A_dense = np.random.rand(10, 10)
    A = sparse.COO(A_dense)

    x0 = splin.inv(sp.csc_matrix(A_dense))
    x = splin.inv(A)
    assert_allclose(x.todense(), x0.todense())


def test_expm():
    np.random.seed(1234)

    A_dense = np.random.rand(10, 10)
    A = sparse.COO(A_dense)

    x0 = splin.expm(sp.csc_matrix(A_dense))
    x = splin.expm(A)
    assert_allclose(x.todense(), x0.todense())


def test_expm_multiply():
    np.random.seed(1234)

    A_dense = np.random.rand(10, 10)
    A = sparse.COO(A_dense)
    b = np.random.rand(10)

    x0 = splin.expm_multiply(A_dense, b)
    x = splin.expm_multiply(A, b)
    assert_allclose(x, x0)
