import pytest
import numpy as np
import scipy.sparse

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

sparrays = [sparray(A) for sparray in sparray_classes]
square_sparrays = [sparray(B) for sparray in sparray_classes]

parametrize_sparrays = pytest.mark.parametrize(
    "A", sparrays, ids=sparray_types
)
parametrize_square_sparrays = pytest.mark.parametrize(
    "B", square_sparrays, ids=sparray_types
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

    assert A[[0]]._is_array, "Expected sparse array, got sparse matrix"
    assert A[1, [1, 2]]._is_array, "Expected ndarray, got sparse array"
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
