"""Test the sparse.linalg._cond1est module."""

import numpy as np
import pytest

from numpy.testing import assert_allclose
from scipy import sparse
from scipy.sparse.linalg import splu, cond1est

# TODO update rtol based on dtype?

rng = np.random.default_rng(565656)

DTYPE_PARAMS = [np.float32, np.float64, np.complex64, np.complex128]
DTYPE_IDS = ['float32', 'float64', 'complex64', 'complex128']

ORD_PARAMS = [1, np.inf]
ORD_IDS = ['ord_1', 'ord_inf']

N = 5  # arbitrary size for the matrices in tests


@pytest.fixture(params=DTYPE_PARAMS, ids=DTYPE_IDS, scope='class')
def dtype(request):
    """Fixture to provide the data type for tests.

    Parameters
    ----------
    request : pytest.FixtureRequest
        The request object for the fixture.

    Returns
    -------
    dtype : data-type
        The data type to be used in the tests.
    """
    return request.param


@pytest.fixture(params=ORD_PARAMS, ids=ORD_IDS, scope='class')
def norm_ord(request):
    """Fixture to provide the norm order for tests.

    Parameters
    ----------
    request : pytest.FixtureRequest
        The request object for the fixture.

    Returns
    -------
    norm_ord : int or np.inf
        The norm order to be used in the tests.
    """
    return request.param


@pytest.fixture
def return_dtype(dtype):
    """Fixture to return the appropriate data type for output."""
    return (np.float32 if (dtype in {np.float32, np.complex64}) else np.float64)


@pytest.fixture
def empty_matrix(dtype):
    """Define an empty CSC matrix."""
    return sparse.csc_array((0, 0), dtype=dtype)


@pytest.fixture
def zero_matrix(dtype):
    """Define a square CSC matrix with non-zero shape, but all zero entries."""
    return sparse.csc_array((N, N), dtype=dtype)


@pytest.fixture
def singleton_matrix(dtype):
    """Define a CSC matrix with a single non-zero element."""
    return sparse.csc_array([[2]], dtype=dtype)


@pytest.fixture
def array_1D(dtype):
    """Define a 1D array as a sparse matrix."""
    return sparse.coo_array(np.arange(N), dtype=dtype)


@pytest.fixture
def array_ND(dtype):
    """Define a non-singular ND array as a sparse matrix."""
    A = sparse.random_array(
        (N + 3, N, N),
        dtype=dtype,
        density=0.5,
        rng=rng,
        data_sampler=rng.normal
    ).toarray()

    # Make each slice non-singular
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            A[i, j, j] = N + i + j

    return sparse.coo_array(A)  # CSC format is only for 2D


@pytest.fixture
def array_rect(dtype):
    """Define a rectangular array as a sparse matrix."""
    return sparse.random_array(
        (N + 2, N),
        density=0.5,
        dtype=dtype,
        format='csc',
        rng=rng,
        data_sampler=rng.normal
    )


@pytest.fixture
def identity_matrix(dtype):
    """Define an identity matrix as a sparse matrix."""
    return sparse.csc_array(sparse.eye_array(N, dtype=dtype))


def generate_matrix(N, dtype, singular=None, density=0.5):
    """Generate a random sparse matrix of size N x N.

    Parameters
    ----------
    N : int
        Size of the matrix (N x N).
    dtype : data-type, optional
        Data type of the matrix elements (default is None, which uses
        float64).
    singular : None or str in {'exactly', 'nearly'}, optional
        If 'exactly', the matrix will be exactly singular (infinite
        condition number). If 'nearly', the matrix will be invertible, but
        ill-conditioned (with a very large condition number). The default is
        None, which creates an invertible, well-conditioned matrix.
    density : float, optional
        Density of the random matrix, between 0 and 1 (default is 0.5).

    Returns
    -------
    result : (N, N) sparse.csc_matrix
        A random sparse matrix in Compressed Sparse Column (CSC) format.
    """
    A = sparse.random_array(
        (N, N),
        density=density,
        format="lil",
        dtype=dtype,
        rng=rng,
        data_sampler=rng.normal
    )

    # Make the matrix non-singular
    A.setdiag(N * (1 + np.arange(N)))

    if singular == 'exactly':
        A[0] = 0
    elif singular == 'nearly':
        eps = np.finfo(dtype).eps
        A[0] = eps * eps

    return A.tocsc()


N_TRIALS = 10  # Number of trials for generating random matrices


@pytest.fixture(params=range(N_TRIALS), ids=lambda x: f"trial_{x}")
def random_nonsingular_matrix(request, dtype, N_max=100, d_scale=1):
    """Generate a random non-singular matrix of maximum size (N_max, N_max)."""
    N = rng.integers(1, N_max, endpoint=True)
    d = d_scale * rng.random()  # densities
    return generate_matrix(N, dtype=dtype, density=d)


@pytest.mark.parametrize("dtype", DTYPE_PARAMS, indirect=True, ids=DTYPE_IDS)
@pytest.mark.parametrize("norm_ord", ORD_PARAMS, indirect=True, ids=ORD_IDS)
class TestNormEstInv:
    def test_error_unsupported_norm(self, dtype, norm_ord):
        if norm_ord != 1:
            pytest.skip(reason="No need to test for every ord value.")

        A = generate_matrix(N, dtype)
        with pytest.raises(ValueError, match="ord must be 1 or np.inf"):
            splu(A).normest_inv(ord=2)

    def test_empty_matrix(self, empty_matrix, norm_ord):
        """Test that an empty matrix returns 0."""
        assert_allclose(splu(empty_matrix).normest_inv(ord=norm_ord), 0)

    def test_zero_matrix(self, zero_matrix, norm_ord):
        with pytest.raises(RuntimeError, match="Factor is exactly singular"):
            splu(zero_matrix).normest_inv(ord=norm_ord)

    def test_singleton_matrix(self, singleton_matrix, dtype, norm_ord, return_dtype):
        assert_allclose(
            splu(singleton_matrix).normest_inv(ord=norm_ord),
            np.array(0.5, dtype=return_dtype),
            strict=True
        )

    def test_identity_matrix(self, identity_matrix, dtype, norm_ord, return_dtype):
        assert_allclose(
            splu(identity_matrix).normest_inv(ord=norm_ord),
            np.array(1.0, dtype=return_dtype),
            strict=True
        )

    def test_exactly_singular_matrix(self, dtype, norm_ord):
        A = generate_matrix(N, dtype, singular='exactly')
        with pytest.raises(RuntimeError, match="Factor is exactly singular"):
            splu(A).normest_inv(ord=norm_ord)

    def test_nearly_singular_matrix(self, dtype, norm_ord):
        A = generate_matrix(N, dtype, singular='nearly')
        norm_A_inv = np.linalg.norm(np.linalg.inv(A.toarray()), ord=norm_ord)
        normest_A_inv = splu(A).normest_inv(ord=norm_ord)
        assert_allclose(normest_A_inv, norm_A_inv, strict=True, rtol=1e-6)

    def test_1D_array(self, array_1D, norm_ord):
        with pytest.raises(
            ValueError,
            match="Cannot convert. CSC format must be 2D. Got 1D"
        ):
            splu(array_1D).normest_inv(ord=norm_ord)

    def test_ND_array(self, array_ND, norm_ord):
        with pytest.raises(
            ValueError,
            match="Cannot convert. CSC format must be 2D. Got 3D"
        ):
            splu(array_ND).normest_inv(ord=norm_ord)

    def test_rectangular_array(self, array_rect, norm_ord):
        with pytest.raises(ValueError, match="can only factor square matrices"):
            splu(array_rect).normest_inv(ord=norm_ord)

    def test_random_nonsingular_matrix(
        self,
        random_nonsingular_matrix,
        norm_ord
    ):
        A = random_nonsingular_matrix
        norm_Ainv = np.linalg.norm(np.linalg.inv(A.toarray()), ord=norm_ord)
        normest_Ainv = splu(A).normest_inv(ord=norm_ord)
        assert_allclose(normest_Ainv, norm_Ainv, strict=True, rtol=1e-6)


class TestCond1Est:
    def test_dense_matrix(self, dtype):
        A = generate_matrix(N, dtype).toarray()
        with pytest.raises(TypeError, match="Input is not a sparse matrix."):
            cond1est(A)

    def test_empty_matrix(self, empty_matrix):
        """Test that an empty matrix raises an error."""
        with pytest.raises(
            ValueError,
            match="Condition number of an empty matrix is undefined"
        ):
            cond1est(empty_matrix)

    def test_zero_matrix(self, zero_matrix):
        assert(cond1est(zero_matrix) == np.inf)

    def test_singleton_matrix(self, singleton_matrix, dtype, return_dtype):
        assert_allclose(
            cond1est(singleton_matrix),
            np.array(1.0, dtype=return_dtype),
            strict=True
        )

    def test_identity_matrix(self, identity_matrix, dtype, return_dtype):
        # Check that we output the correct data type
        assert_allclose(
            cond1est(identity_matrix),
            np.array(1.0, dtype=return_dtype),
            strict=True
        )

    def test_exactly_singular_matrix(self, dtype):
        A = generate_matrix(N, dtype, singular='exactly')
        assert(cond1est(A) == np.inf)

    def test_nearly_singular_matrix(self, dtype):
        A = generate_matrix(N, dtype, singular='nearly')
        cond_A = np.linalg.cond(np.linalg.inv(A.toarray()), p=1)
        # NOTE there is a bug in np.linalg.cond that returns np.complex when
        # the matrix input is complex, so take the real part for comparison.
        assert_allclose(cond1est(A), cond_A.real, strict=True, rtol=1e-6)

    def test_1D_array(self, array_1D):
        with pytest.raises(
            ValueError,
            match="Input must be a 2-dimensional matrix."
        ):
            cond1est(array_1D)

    def test_ND_array(self, array_ND):
        with pytest.raises(
            ValueError,
            match="Input must be a 2-dimensional matrix."
        ):
            cond1est(array_ND)

    def test_rectangular_array(self, array_rect):
        with pytest.raises(ValueError, match="Matrix must be square."):
            cond1est(array_rect)

    def test_random_nonsingular_matrix(self, random_nonsingular_matrix):
        A = random_nonsingular_matrix
        cond_A = np.linalg.cond(A.toarray(), p=1)
        # NOTE there is a bug in np.linalg.cond that returns np.complex when
        # the matrix input is complex, so take the real part for comparison.
        assert_allclose(cond1est(A), cond_A.real, strict=True, rtol=1e-6)
