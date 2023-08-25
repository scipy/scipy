import pytest
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose

from scipy.sparse import diags, csgraph
from scipy.linalg import eigh

from scipy.sparse.linalg import LaplacianNd
from scipy.sparse.linalg._special_sparse_arrays import Sakurai

INT_DTYPES = [np.int8, np.int16, np.int32, np.int64]
REAL_DTYPES = [np.float32, np.float64]
COMPLEX_DTYPES = [np.complex64, np.complex128]

class TestLaplacianNd:
    """
    LaplacianNd tests
    """
    ALLDTYPES = INT_DTYPES + REAL_DTYPES + COMPLEX_DTYPES

    @pytest.mark.parametrize('bc', ['neumann', 'dirichlet', 'periodic'])
    def test_1d_specific_shape(self, bc):
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
        G = diags(np.ones(n - 1), 1, format='dia')
        Lf = csgraph.laplacian(G, symmetrized=True, form='function')
        La = csgraph.laplacian(G, symmetrized=True, form='array')
        grid_shape = (n,)
        bc = 'neumann'
        lap = LaplacianNd(grid_shape, boundary_conditions=bc)
        assert_array_equal(lap(np.eye(n)), -Lf(np.eye(n)))
        assert_array_equal(lap.toarray(), -La.toarray())
        # https://github.com/numpy/numpy/issues/24351
        assert_array_equal(lap.tosparse().toarray(), -La.toarray())

    @pytest.mark.parametrize('grid_shape', [(6, ), (2, 3), (2, 3, 4)])
    @pytest.mark.parametrize('bc', ['neumann', 'dirichlet', 'periodic'])
    def test_eigenvalues(self, grid_shape, bc):
        lap = LaplacianNd(grid_shape, boundary_conditions=bc, dtype=np.float64)
        L = lap.toarray()
        eigvals = eigh(L, eigvals_only=True)
        n = np.prod(grid_shape)
        eigenvalues = lap.eigenvalues()
        dtype = eigenvalues.dtype
        atol = n * n * np.finfo(dtype).eps
        # test the default ``m = None``
        assert_allclose(eigenvalues, eigvals, atol=atol)
        # test every ``m > 0``
        for m in np.arange(1, n + 1):
            assert_array_equal(lap.eigenvalues(m), eigenvalues[-m:])

    @pytest.mark.parametrize('grid_shape', [(6, ), (2, 3), (2, 3, 4)])
    @pytest.mark.parametrize('bc', ['neumann', 'dirichlet', 'periodic'])
    def test_eigenvectors(self, grid_shape, bc):
        lap = LaplacianNd(grid_shape, boundary_conditions=bc, dtype=np.float64)
        n = np.prod(grid_shape)
        eigenvalues = lap.eigenvalues()
        eigenvectors = lap.eigenvectors()
        dtype = eigenvectors.dtype
        atol = n * n * max(np.finfo(dtype).eps, np.finfo(np.double).eps)
        # test the default ``m = None`` every individual eigenvector
        for i in np.arange(n):
            r = lap.toarray() @ eigenvectors[:, i] - eigenvectors[:, i] * eigenvalues[i]
            assert_allclose(r, np.zeros_like(r), atol=atol)
        # test every ``m > 0``
        for m in np.arange(1, n + 1):
            e = lap.eigenvalues(m)
            ev = lap.eigenvectors(m)
            r = lap.toarray() @ ev - ev @ np.diag(e)
            assert_allclose(r, np.zeros_like(r), atol=atol)

    @pytest.mark.parametrize('grid_shape', [(6, ), (2, 3), (2, 3, 4)])
    @pytest.mark.parametrize('bc', ['neumann', 'dirichlet', 'periodic'])
    def test_toarray_tosparse_consistency(self, grid_shape, bc):
        lap = LaplacianNd(grid_shape, boundary_conditions=bc)
        n = np.prod(grid_shape)
        assert_array_equal(lap.toarray(), lap(np.eye(n)))
        assert_array_equal(lap.tosparse().toarray(), lap.toarray())

    @pytest.mark.parametrize('dtype', ALLDTYPES)
    @pytest.mark.parametrize('grid_shape', [(6, ), (2, 3), (2, 3, 4)])
    @pytest.mark.parametrize('bc', ['neumann', 'dirichlet', 'periodic'])
    def test_linearoperator_shape_dtype(self, grid_shape, bc, dtype):
        lap = LaplacianNd(grid_shape, boundary_conditions=bc, dtype=dtype)
        n = np.prod(grid_shape)
        assert lap.shape == (n, n)
        assert lap.dtype == dtype
        assert_array_equal(
            LaplacianNd(
                grid_shape, boundary_conditions=bc, dtype=dtype
            ).toarray(),
            LaplacianNd(grid_shape, boundary_conditions=bc)
            .toarray()
            .astype(dtype),
        )
        assert_array_equal(
            LaplacianNd(grid_shape, boundary_conditions=bc, dtype=dtype)
            .tosparse()
            .toarray(),
            LaplacianNd(grid_shape, boundary_conditions=bc)
            .tosparse()
            .toarray()
            .astype(dtype),
        )

    @pytest.mark.parametrize('dtype', ALLDTYPES)
    @pytest.mark.parametrize('grid_shape', [(6, ), (2, 3), (2, 3, 4)])
    @pytest.mark.parametrize('bc', ['neumann', 'dirichlet', 'periodic'])
    def test_dot(self, grid_shape, bc, dtype):
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

    def test_boundary_conditions_value_error(self):
        with pytest.raises(ValueError, match="Unknown value 'robin'"):
            LaplacianNd(grid_shape=(6, ), boundary_conditions='robin')

            
class TestSakurai:
    """
    Sakurai tests
    """
    ALLDTYPES = [np.int64] + REAL_DTYPES + COMPLEX_DTYPES
    def test_specific_shape(self):
        sak = Sakurai(6)
        assert_array_equal(sak.toarray(), sak(np.eye(6)))
        a = np.array(
            [
                [ 5, -4,  1,  0,  0,  0],
                [-4,  6, -4,  1,  0,  0],
                [ 1, -4,  6, -4,  1,  0],
                [ 0,  1, -4,  6, -4,  1],
                [ 0,  0,  1, -4,  6, -4],
                [ 0,  0,  0,  1, -4,  5]
            ]
        )

        np.array_equal(a, sak.toarray())
        np.array_equal(sak.tosparse().toarray(), sak.toarray())
        ab = np.array(
            [
                [ 1,  1,  1,  1,  1,  1],
                [-4, -4, -4, -4, -4, -4],
                [ 5,  6,  6,  6,  6,  5]
            ]
        )
        np.array_equal(ab, sak.tobanded())
        e = np.array(
                [0.03922866, 0.56703972, 2.41789479, 5.97822974,
                 10.54287655, 14.45473055]
            )
        np.array_equal(e, sak.eigenvalues)

    @pytest.mark.parametrize('dtype', ALLDTYPES)
    def test_dot(self, dtype):
        n = 5
        sak = Sakurai(n)
        x0 = np.arange(n)
        x1 = x0.reshape((-1, 1))
        x2 = np.arange(2 * n).reshape((n, 2))
        input_set = [x0, x1, x2]
        for x in input_set:
            y = sak.dot(x.astype(dtype))
            assert x.shape == y.shape
            assert y.dtype == dtype
            x = x.reshape((-1, 1)) if x.ndim == 1 else pass
            yy = sak.toarray() @ x
            assert yy.dtype == dtype
            np.array_equal(y, yy)
