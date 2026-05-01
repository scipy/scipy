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
    assert isinstance(A.sum(axis=0), np.ndarray), "Expected ndarray"
    assert A.sum(axis=0).shape == (4,)
    assert A.sum(axis=1).shape == (3,)


@parametrize_sparrays
def test_mean(A):
    assert isinstance(A.mean(axis=1), np.ndarray), "Expected ndarray"


@parametrize_sparrays
def test_min_max(A):
    # Some formats don't support min/max operations, so we skip those here.
    if hasattr(A, 'min'):
        assert isinstance(A.min(axis=1), scipy.sparse.sparray), "Expected sparray min"
    if hasattr(A, 'max'):
        assert isinstance(A.max(axis=1), scipy.sparse.sparray), "Expected sparray max"
    if hasattr(A, 'argmin'):
        assert isinstance(A.argmin(axis=1), np.ndarray), "Expected ndarray argmin"
    if hasattr(A, 'argmax'):
        assert isinstance(A.argmax(axis=1), np.ndarray), "Expected ndarray argmax"


@parametrize_sparrays
def test_todense(A):
    assert isinstance(A.todense(), np.ndarray), "Expected ndarray"
    assert isinstance(A.toarray(), np.ndarray), "Expected ndarray"


@parametrize_sparrays
def test_indexing(A):
    if A.__class__.__name__[:3] in ('dia', 'coo', 'bsr'):
        return

    all_res = (
        A[1, :],
        A[:, 1],
        A[1, [1, 2]],
        A[[1, 2], 1],
        A[[0]],
        A[:, [1, 2]],
        A[[1, 2], :],
        A[1, [[1, 2]]],
        A[[[1, 2]], 1],
    )

    for res in all_res:
        assert isinstance(res, scipy.sparse.sparray), \
            f"Expected sparse array, got {res._class__.__name__}"


@parametrize_sparrays
def test_dense_addition(A):
    X = np.random.random(A.shape)
    assert isinstance(A + X, np.ndarray), "Expected ndarray"


@parametrize_sparrays
def test_sparse_addition(A):
    assert isinstance((A + A), scipy.sparse.sparray), "Expected sparray"


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


@parametrize_sparrays
def test_power_operator(A):
    assert isinstance((A**2), scipy.sparse.sparray), "Expected array, got matrix"

    # https://github.com/scipy/scipy/issues/15948
    npt.assert_equal((A**2).todense(), (A.todense())**2)

    # power of zero is all ones (dense) so helpful msg exception
    with pytest.raises(NotImplementedError, match="zero power"):
        A**0


@parametrize_sparrays
def test_sparse_divide(A):
    assert isinstance(A / A, np.ndarray)

@parametrize_sparrays
def test_sparse_dense_divide(A):
    with pytest.warns(RuntimeWarning):
        assert isinstance((A / A.todense()), scipy.sparse.sparray)

@parametrize_sparrays
def test_dense_divide(A):
    assert isinstance((A / 2), scipy.sparse.sparray), "Expected array, got matrix"


@parametrize_sparrays
def test_no_A_attr(A):
    with pytest.raises(AttributeError):
        A.A


@parametrize_sparrays
def test_no_H_attr(A):
    with pytest.raises(AttributeError):
        A.H


@parametrize_sparrays
def test_getrow_getcol(A):
    assert isinstance(A._getcol(0), scipy.sparse.sparray)
    assert isinstance(A._getrow(0), scipy.sparse.sparray)


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

    assert isinstance(C, scipy.sparse.sparray)
    npt.assert_allclose(C.todense(), np.linalg.inv(B.todense()))


@parametrize_square_sparrays
def test_expm(B):
    if B.__class__.__name__[:3] != 'csc':
        return

    Bmat = scipy.sparse.csc_array(B)
    C = spla.expm(B)

    assert isinstance(C, scipy.sparse.sparray)
    npt.assert_allclose(C.todense(), spla.expm(Bmat).todense())


@parametrize_square_sparrays
def test_expm_multiply(B):
    if B.__class__.__name__[:3] != 'csc':
        return

    npt.assert_allclose(spla.expm_multiply(B, np.array([1, 2])), spla.expm(B) @ [1, 2])


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


@pytest.mark.parametrize("fmt",["csr","csc"])
def test_spsolve_triangular(fmt):
    arr = [
        [1, 0, 0, 0],
        [2, 1, 0, 0],
        [3, 2, 1, 0],
        [4, 3, 2, 1],
    ]
    if fmt == "csr":
      X = scipy.sparse.csr_array(arr)
    else:
      X = scipy.sparse.csc_array(arr)
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
    npt.assert_allclose(
        LU.solve(np.array([1, 2, 3, 4])),
        np.asarray([1, 0, 0, 0], dtype=np.float64),
        rtol=1e-14, atol=3e-16
    )


def test_spilu():
    X = scipy.sparse.csc_array([
        [1, 0, 0, 0],
        [2, 1, 0, 0],
        [3, 2, 1, 0],
        [4, 3, 2, 1],
    ])
    LU = spla.spilu(X)
    npt.assert_allclose(
        LU.solve(np.array([1, 2, 3, 4])),
        np.asarray([1, 0, 0, 0], dtype=np.float64),
        rtol=1e-14, atol=3e-16
    )


@pytest.mark.parametrize(
    "cls,indices_attrs",
    [
        (
            scipy.sparse.csr_array,
            ["indices", "indptr"],
        ),
        (
            scipy.sparse.csc_array,
            ["indices", "indptr"],
        ),
        (
            scipy.sparse.coo_array,
            ["row", "col"],
        ),
    ]
)
@pytest.mark.parametrize("expected_dtype", [np.int64, np.int32])
def test_index_dtype_compressed(cls, indices_attrs, expected_dtype):
    input_array = scipy.sparse.coo_array(np.arange(9).reshape(3, 3))
    coo_tuple = (
        input_array.data,
        (
            input_array.row.astype(expected_dtype),
            input_array.col.astype(expected_dtype),
        )
    )

    result = cls(coo_tuple)
    for attr in indices_attrs:
        assert getattr(result, attr).dtype == expected_dtype

    result = cls(coo_tuple, shape=(3, 3))
    for attr in indices_attrs:
        assert getattr(result, attr).dtype == expected_dtype

    if issubclass(cls, scipy.sparse._compressed._cs_matrix):
        input_array_csr = input_array.tocsr()
        csr_tuple = (
            input_array_csr.data,
            input_array_csr.indices.astype(expected_dtype),
            input_array_csr.indptr.astype(expected_dtype),
        )

        result = cls(csr_tuple)
        for attr in indices_attrs:
            assert getattr(result, attr).dtype == expected_dtype

        result = cls(csr_tuple, shape=(3, 3))
        for attr in indices_attrs:
            assert getattr(result, attr).dtype == expected_dtype


def test_default_is_matrix_diags():
    assert not hasattr(scipy.sparse, "diags")
    m = scipy.sparse.diags_array([0.0, 1.0, 2.0])
    assert isinstance(m, scipy.sparse.sparray)


def test_default_is_matrix_eye():
    assert not hasattr(scipy.sparse, "eye")
    m = scipy.sparse.eye_array(3)
    assert isinstance(m, scipy.sparse.sparray)


def test_default_is_matrix_identity():
    assert not hasattr(scipy.sparse, "identity")
    m = scipy.sparse.eye_array(3)
    assert isinstance(m, scipy.sparse.sparray)


@pytest.mark.filterwarnings("ignore:.*switching.*sparse array int:DeprecationWarning")
def test_default_is_matrix_kron_dense():
    m = scipy.sparse.kron(
        np.array([[1, 2], [3, 4]]), np.array([[4, 3], [2, 1]])
    )
    assert isinstance(m, scipy.sparse.sparray)


@pytest.mark.filterwarnings("ignore:.*switching.*sparse array int:DeprecationWarning")
def test_default_is_matrix_kron_sparse():
    m = scipy.sparse.kron(
        np.array([[1, 2], [3, 4]]), np.array([[1, 0], [0, 0]])
    )
    assert isinstance(m, scipy.sparse.sparray)


@pytest.mark.filterwarnings("ignore:.*switching.*sparse array int:DeprecationWarning")
def test_default_is_matrix_kronsum():
    m = scipy.sparse.kronsum(
        np.array([[1, 0], [0, 1]]), np.array([[0, 1], [1, 0]])
    )
    assert isinstance(m, scipy.sparse.sparray)


def test_default_is_matrix_random():
    assert not hasattr(scipy.sparse, "random")
    m = scipy.sparse.random_array((3, 3))
    assert isinstance(m, scipy.sparse.sparray)



@pytest.mark.parametrize("fn", (scipy.sparse.hstack, scipy.sparse.vstack))
def test_default_is_matrix_stacks(fn):
    """Same idea as `test_default_construction_fn_matrices`, but for the
    stacking creation functions."""
    A = scipy.sparse.coo_array(np.eye(2))
    B = scipy.sparse.coo_array([[0, 1], [1, 0]])
    m = fn([A, B])
    assert isinstance(m, scipy.sparse.sparray)


def test_blocks_default_construction_fn_matrices():
    """Same idea as `test_default_construction_fn_matrices`, but for the block
    creation function"""
    A = scipy.sparse.coo_array(np.eye(2))
    B = scipy.sparse.coo_array([[2], [0]])
    C = scipy.sparse.coo_array([[3]])

    # block diag
    m = scipy.sparse.block_diag((A, B, C))
    assert isinstance(m, scipy.sparse.sparray)

    # block_array
    m = scipy.sparse.block_array([[A, None], [None, C]])
    assert isinstance(m, scipy.sparse.sparray)


def test_format_property():
    for fmt in sparray_types:
        arr_cls = getattr(scipy.sparse, f"{fmt}_array")
        M = arr_cls([[1, 2]])
        assert M.format == fmt
        assert M._format == fmt
        with pytest.raises(AttributeError):
            M.format = "qqq"


def test_issparse():
    a = scipy.sparse.eye_array(3)
    assert isinstance(a, scipy.sparse.sparray)
    assert scipy.sparse.issparse(a)
    assert not scipy.sparse.issparse(a.todense())
