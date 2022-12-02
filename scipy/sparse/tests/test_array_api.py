import pytest
import numpy as np
import numpy.testing as npt
import scipy.sparse
import scipy.sparse.linalg as spla

sparray_types = ('bsr', 'coo', 'csc', 'csr', 'dia', 'dok', 'lil')

sparray_classes = [
    getattr(scipy.sparse, f'{T}_array') for T in sparray_types
]

A = np.array([
    [0, 1, 2, 0],
    [2, 0, 0, 3],
    [1, 4, 0, 0]
])

B = np.array([
    [0, 1],
    [2, 0]
])

X = np.array([
    [1, 0, 0, 1],
    [2, 1, 2, 0],
    [0, 2, 1, 0],
    [0, 0, 1, 2]
], dtype=float)


sparrays = [sparray(A) for sparray in sparray_classes]
square_sparrays = [sparray(B) for sparray in sparray_classes]
eig_sparrays = [sparray(X) for sparray in sparray_classes]

parametrize_sparrays = pytest.mark.parametrize(
    "A", sparrays, ids=sparray_types
)
parametrize_square_sparrays = pytest.mark.parametrize(
    "B", square_sparrays, ids=sparray_types
)
parametrize_eig_sparrays = pytest.mark.parametrize(
    "X", eig_sparrays, ids=sparray_types
)


@parametrize_sparrays
def test_sum(A):
    assert not isinstance(A.sum(axis=0), np.matrix), \
        "Expected array, got matrix"
    assert A.sum(axis=0).shape == (4,)
    assert A.sum(axis=1).shape == (3,)


@parametrize_sparrays
def test_mean(A):
    assert not isinstance(A.mean(axis=1), np.matrix), \
        "Expected array, got matrix"


@parametrize_sparrays
def test_todense(A):
    assert not isinstance(A.todense(), np.matrix), \
        "Expected array, got matrix"


@parametrize_sparrays
def test_indexing(A):
    if A.__class__.__name__[:3] in ('dia', 'coo', 'bsr'):
        return

    with pytest.raises(NotImplementedError):
        A[1, :]

    with pytest.raises(NotImplementedError):
        A[:, 1]

    with pytest.raises(NotImplementedError):
        A[1, [1, 2]]

    with pytest.raises(NotImplementedError):
        A[[1, 2], 1]

    assert A[[0]]._is_array, "Expected sparse array, got sparse matrix"
    assert A[1, [[1, 2]]]._is_array, "Expected ndarray, got sparse array"
    assert A[[[1, 2]], 1]._is_array, "Expected ndarray, got sparse array"
    assert A[:, [1, 2]]._is_array, "Expected sparse array, got something else"


@parametrize_sparrays
def test_dense_addition(A):
    X = np.random.random(A.shape)
    assert not isinstance(A + X, np.matrix), "Expected array, got matrix"


@parametrize_sparrays
def test_sparse_addition(A):
    assert (A + A)._is_array, "Expected array, got matrix"


@parametrize_sparrays
def test_elementwise_mul(A):
    assert np.all((A * A).todense() == A.power(2).todense())


@parametrize_sparrays
def test_elementwise_rmul(A):
    with pytest.raises(TypeError):
        None * A

    with pytest.raises(ValueError):
        np.eye(3) * scipy.sparse.csr_array(np.arange(6).reshape(2, 3))

    assert np.all((2 * A) == (A.todense() * 2))

    assert np.all((A.todense() * A) == (A.todense() ** 2))


@parametrize_sparrays
def test_matmul(A):
    assert np.all((A @ A.T).todense() == A.dot(A.T).todense())


@parametrize_square_sparrays
def test_pow(B):
    assert (B**0)._is_array, "Expected array, got matrix"
    assert (B**2)._is_array, "Expected array, got matrix"


@parametrize_sparrays
def test_sparse_divide(A):
    assert isinstance(A / A, np.ndarray)


@parametrize_sparrays
def test_dense_divide(A):
    assert (A / 2)._is_array, "Expected array, got matrix"


@parametrize_sparrays
def test_no_A_attr(A):
    with pytest.warns(np.VisibleDeprecationWarning):
        A.A


@parametrize_sparrays
def test_no_H_attr(A):
    with pytest.warns(np.VisibleDeprecationWarning):
        A.H


@parametrize_sparrays
def test_getrow_getcol(A):
    assert A.getcol(0)._is_array
    assert A.getrow(0)._is_array


@parametrize_sparrays
def test_docstr(A):
    if A.__doc__ is None:
        return

    docstr = A.__doc__.lower()
    for phrase in ('matrix', 'matrices'):
        assert phrase not in docstr


# -- linalg --

@parametrize_sparrays
def test_as_linearoperator(A):
    L = spla.aslinearoperator(A)
    npt.assert_allclose(L * [1, 2, 3, 4], A @ [1, 2, 3, 4])


@parametrize_square_sparrays
def test_inv(B):
    if B.__class__.__name__[:3] != 'csc':
        return

    C = spla.inv(B)

    assert C._is_array
    npt.assert_allclose(C.todense(), np.linalg.inv(B.todense()))


@parametrize_square_sparrays
def test_expm(B):
    if B.__class__.__name__[:3] != 'csc':
        return

    Bmat = scipy.sparse.csc_matrix(B)

    C = spla.expm(B)

    assert C._is_array
    npt.assert_allclose(
        C.todense(),
        spla.expm(Bmat).todense()
    )


@parametrize_square_sparrays
def test_expm_multiply(B):
    if B.__class__.__name__[:3] != 'csc':
        return

    npt.assert_allclose(
        spla.expm_multiply(B, np.array([1, 2])),
        spla.expm(B) @ [1, 2]
    )


@parametrize_sparrays
def test_norm(A):
    C = spla.norm(A)
    npt.assert_allclose(C, np.linalg.norm(A.todense()))


@parametrize_square_sparrays
def test_onenormest(B):
    C = spla.onenormest(B)
    npt.assert_allclose(C, np.linalg.norm(B.todense(), 1))


@parametrize_square_sparrays
def test_spsolve(B):
    if B.__class__.__name__[:3] not in ('csc', 'csr'):
        return

    npt.assert_allclose(
        spla.spsolve(B, [1, 2]),
        np.linalg.solve(B.todense(), [1, 2])
    )


def test_spsolve_triangular():
    X = scipy.sparse.csr_array([
        [1, 0, 0, 0],
        [2, 1, 0, 0],
        [3, 2, 1, 0],
        [4, 3, 2, 1],
    ])
    spla.spsolve_triangular(X, [1, 2, 3, 4])


@parametrize_square_sparrays
def test_factorized(B):
    if B.__class__.__name__[:3] != 'csc':
        return

    LU = spla.factorized(B)
    npt.assert_allclose(
        LU(np.array([1, 2])),
        np.linalg.solve(B.todense(), [1, 2])
    )


@parametrize_square_sparrays
@pytest.mark.parametrize(
    "solver",
    ["bicg", "bicgstab", "cg", "cgs", "gmres", "lgmres", "minres", "qmr",
     "gcrotmk", "tfqmr"]
)
def test_solvers(B, solver):
    if solver == "minres":
        kwargs = {}
    else:
        kwargs = {'atol': 1e-5}

    x, info = getattr(spla, solver)(B, np.array([1, 2]), **kwargs)
    assert info >= 0  # no errors, even if perhaps did not converge fully
    npt.assert_allclose(x, [1, 1], atol=1e-1)


@parametrize_sparrays
@pytest.mark.parametrize(
    "solver",
    ["lsqr", "lsmr"]
)
def test_lstsqr(A, solver):
    x, *_ = getattr(spla, solver)(A, [1, 2, 3])
    npt.assert_allclose(A @ x, [1, 2, 3])


@parametrize_eig_sparrays
def test_eigs(X):
    e, v = spla.eigs(X, k=1)
    npt.assert_allclose(
        X @ v,
        e[0] * v
    )


@parametrize_eig_sparrays
def test_eigsh(X):
    X = X + X.T
    e, v = spla.eigsh(X, k=1)
    npt.assert_allclose(
        X @ v,
        e[0] * v
    )


@parametrize_eig_sparrays
def test_svds(X):
    u, s, vh = spla.svds(X, k=3)
    u2, s2, vh2 = np.linalg.svd(X.todense())
    s = np.sort(s)
    s2 = np.sort(s2[:3])
    npt.assert_allclose(s, s2, atol=1e-3)


def test_splu():
    X = scipy.sparse.csc_array([
        [1, 0, 0, 0],
        [2, 1, 0, 0],
        [3, 2, 1, 0],
        [4, 3, 2, 1],
    ])
    LU = spla.splu(X)
    npt.assert_allclose(LU.solve(np.array([1, 2, 3, 4])), [1, 0, 0, 0])


def test_spilu():
    X = scipy.sparse.csc_array([
        [1, 0, 0, 0],
        [2, 1, 0, 0],
        [3, 2, 1, 0],
        [4, 3, 2, 1],
    ])
    LU = spla.spilu(X)
    npt.assert_allclose(LU.solve(np.array([1, 2, 3, 4])), [1, 0, 0, 0])


@parametrize_sparrays
def test_power_operator(A):
    # https://github.com/scipy/scipy/issues/15948
    npt.assert_equal((A**2).todense(), (A.todense())**2)
