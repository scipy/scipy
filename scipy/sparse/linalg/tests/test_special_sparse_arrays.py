import pytest
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose

from scipy.sparse import diags, csgraph
from scipy.linalg import eigh

from scipy.sparse.linalg import LaplacianNd


class TestLaplacianNd:
    """
    LaplacianNd tests
    """
    INT_DTYPES = [np.int32, np.int64]
    REAL_DTYPES = [np.float32, np.float64]
    COMPLEX_DTYPES = [np.complex64, np.complex128]
    ALLDTYPES = INT_DTYPES + REAL_DTYPES + COMPLEX_DTYPES

    @pytest.mark.parametrize("bc", ['neumann', 'dirichlet', 'periodic'])
    def test_1d_specific_shape(self, bc):
        bc = "neumann"
        lap = LaplacianNd(grid_shape=(6, ), boundary_conditions=bc)
        lapa = lap.toarray()
        if bc == 'neumann':
            a = np.array(
                [
                    [-1, 1, 0, 0, 0, 0],
                    [1, -2, 1, 0, 0, 0],
                    [0, 1, -2, 1, 0, 0],
                    [0, 0, 1, -2, 1, 0],
                    [0, 0, 0, 1, -2, 1],
                    [0, 0, 0, 0, 1, -1],
                ]
            )
        elif bc == 'dirichlet':
            a = np.array(
                [
                    [-2, 1, 0, 0, 0, 0],
                    [1, -2, 1, 0, 0, 0],
                    [0, 1, -2, 1, 0, 0],
                    [0, 0, 1, -2, 1, 0],
                    [0, 0, 0, 1, -2, 1],
                    [0, 0, 0, 0, 1, -2],
                ]
            )
        else:
            a = np.array(
                [
                    [-2, 1, 0, 0, 0, 1],
                    [1, -2, 1, 0, 0, 0],
                    [0, 1, -2, 1, 0, 0],
                    [0, 0, 1, -2, 1, 0],
                    [0, 0, 0, 1, -2, 1],
                    [1, 0, 0, 0, 1, -2],
                ]
            )
        assert_array_equal(a, lapa)

    def test_1d_with_graph_laplacian(self):
        n = 6
        G = diags(np.ones(n - 1), 1, format="dia")
        Lf = csgraph.laplacian(G, symmetrized=True, form="function")
        La = csgraph.laplacian(G, symmetrized=True, form="array")
        grid_shape = (n,)
        bc = "neumann"
        lap = LaplacianNd(grid_shape, boundary_conditions=bc)
        assert_array_equal(lap(np.eye(n)), -Lf(np.eye(n)))
        assert_array_equal(lap.toarray(), -La.toarray())
        # https://github.com/numpy/numpy/issues/24351
        assert_array_equal(lap.tosparse().toarray(), -La.toarray())

    @pytest.mark.parametrize("grid_shape", [(6, ), (2, 3), (2, 3, 4)])
    @pytest.mark.parametrize("bc", ['neumann', 'dirichlet', 'periodic'])
    def test_eigenvalues(self, grid_shape, bc):
        lap = LaplacianNd(grid_shape, boundary_conditions=bc, dtype=np.float64)
        eigenvalues = lap.eigs
        L = lap.toarray()
        eigvals = eigh(L, eigvals_only=True)
        dtype = eigenvalues.dtype
        n = np.prod(grid_shape)
        atol = n * n * max(np.finfo(dtype).eps, np.finfo(np.double).eps)
        assert_allclose(eigenvalues, eigvals, atol=atol)

    @pytest.mark.parametrize("grid_shape", [(6, ), (2, 3), (2, 3, 4)])
    @pytest.mark.parametrize("bc", ['neumann', 'dirichlet', 'periodic'])
    def test_toarray_tosparse_consistency(self, grid_shape, bc):
        lap = LaplacianNd(grid_shape, boundary_conditions=bc)
        n = np.prod(grid_shape)
        assert_array_equal(lap.toarray(), lap(np.eye(n)))
        assert_array_equal(lap.tosparse().toarray(), lap.toarray())

    @pytest.mark.parametrize("dtype", ALLDTYPES)
    @pytest.mark.parametrize("grid_shape", [(6, ), (2, 3), (2, 3, 4)])
    @pytest.mark.parametrize("bc", ['neumann', 'dirichlet', 'periodic'])
    def test_input_output_dtype_shape_consistency(self, grid_shape, bc, dtype):
        lap = LaplacianNd(grid_shape, boundary_conditions=bc, dtype=dtype)
        n = np.prod(grid_shape)
        assert lap.shape == (n, n)
        assert lap.dtype == dtype

    @pytest.mark.parametrize("dtype", ALLDTYPES)
    @pytest.mark.parametrize("grid_shape", [(6, ), (2, 3), (2, 3, 4)])
    @pytest.mark.parametrize("bc", ['neumann', 'dirichlet', 'periodic'])
    def test_matmat(self, grid_shape, bc, dtype):
        lap = LaplacianNd(grid_shape, boundary_conditions=bc)
        n = np.prod(grid_shape)
        x0 = np.arange(n)
        x1 = x0.reshape((-1, 1))
        x2 = np.arange(2 * n).reshape((n, 2))
        input_set = [x0, x1, x2]
        for x in input_set:
            y = lap.dot(x.astype(dtype))
            assert x.shape == y.shape
            assert y.dtype == dtype
