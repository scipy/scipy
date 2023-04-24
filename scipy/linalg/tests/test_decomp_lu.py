import pytest
from pytest import raises as assert_raises

import numpy as np
from scipy.linalg import lu, lu_factor, lu_solve, get_lapack_funcs, solve
from numpy.testing import assert_allclose, assert_array_equal


class TestLU:
    def setup_method(self):
        self.rng = np.random.default_rng(1682281250228846)

    @pytest.mark.parametrize('shape', [[2, 2], [2, 4], [4, 2],
                                       [20, 20], [20, 4], [4, 20],
                                       [3, 2, 23, 23], [3, 2, 17, 5],
                                       [3, 2, 11, 59]])
    def test_simple_det_shapes_real_complex(self, shape):
        a = self.rng.uniform(-10., 10., size=shape)
        p, l, u = lu(a)
        assert_allclose(a, p @ l @ u)
        pl, u = lu(a, permute_l=True)
        assert_allclose(a, pl @ u)

        b = self.rng.uniform(-10., 10., size=shape)*1j
        b += self.rng.uniform(-10, 10, size=shape)
        pl, u = lu(b, permute_l=True)
        assert_allclose(b, pl @ u)

    @pytest.mark.parametrize('shape', [[2, 2], [2, 4], [4, 2], [20, 20],
                                       [20, 4], [4, 20]])
    def test_simple_det_shapes_real_complex_2d_indices(self, shape):
        a = self.rng.uniform(-10., 10., size=shape)
        p, l, u = lu(a, p_indices=True)
        assert_allclose(a, l[p, :] @ u)

    def test_1by1_input_output(self):
        a = self.rng.random([4, 5, 1, 1], dtype=np.float32)
        p, l, u = lu(a, p_indices=True)
        assert_allclose(p, np.zeros(shape=(4, 5, 1), dtype=int))
        assert_allclose(l, np.ones(shape=(4, 5, 1, 1), dtype=np.float32))
        assert_allclose(u, a)

        a = self.rng.random([4, 5, 1, 1], dtype=np.float32)
        p, l, u = lu(a)
        assert_allclose(p, np.ones(shape=(4, 5, 1, 1), dtype=np.float32))
        assert_allclose(l, np.ones(shape=(4, 5, 1, 1), dtype=np.float32))
        assert_allclose(u, a)

        pl, u = lu(a, permute_l=True)
        assert_allclose(pl, np.ones(shape=(4, 5, 1, 1), dtype=np.float32))
        assert_allclose(u, a)

        a = self.rng.random([4, 5, 1, 1], dtype=np.float32)*np.complex64(1.j)
        p, l, u = lu(a)
        assert_allclose(p, np.ones(shape=(4, 5, 1, 1), dtype=np.complex64))
        assert_allclose(l, np.ones(shape=(4, 5, 1, 1), dtype=np.complex64))
        assert_allclose(u, a)

    def test_empty_edge_cases(self):
        a = np.empty([0, 0])
        p, l, u = lu(a)
        assert_allclose(p, np.empty(shape=(0, 0), dtype=np.float64))
        assert_allclose(l, np.empty(shape=(0, 0), dtype=np.float64))
        assert_allclose(u, np.empty(shape=(0, 0), dtype=np.float64))

        a = np.empty([0, 3], dtype=np.float16)
        p, l, u = lu(a)
        assert_allclose(p, np.empty(shape=(0, 0), dtype=np.float32))
        assert_allclose(l, np.empty(shape=(0, 0), dtype=np.float32))
        assert_allclose(u, np.empty(shape=(0, 3), dtype=np.float32))

        a = np.empty([3, 0], dtype=np.complex64)
        p, l, u = lu(a)
        assert_allclose(p, np.empty(shape=(0, 0), dtype=np.float32))
        assert_allclose(l, np.empty(shape=(3, 0), dtype=np.complex64))
        assert_allclose(u, np.empty(shape=(0, 0), dtype=np.complex64))
        p, l, u = lu(a, p_indices=True)
        assert_allclose(p, np.empty(shape=(0,), dtype=int))
        assert_allclose(l, np.empty(shape=(3, 0), dtype=np.complex64))
        assert_allclose(u, np.empty(shape=(0, 0), dtype=np.complex64))
        pl, u = lu(a, permute_l=True)
        assert_allclose(pl, np.empty(shape=(3, 0), dtype=np.complex64))
        assert_allclose(u, np.empty(shape=(0, 0), dtype=np.complex64))

        a = np.empty([3, 0, 0], dtype=np.complex64)
        p, l, u = lu(a)
        assert_allclose(p, np.empty(shape=(3, 0, 0), dtype=np.float32))
        assert_allclose(l, np.empty(shape=(3, 0, 0), dtype=np.complex64))
        assert_allclose(u, np.empty(shape=(3, 0, 0), dtype=np.complex64))

        a = np.empty([0, 0, 3])
        p, l, u = lu(a)
        assert_allclose(p, np.empty(shape=(0, 0, 0)))
        assert_allclose(l, np.empty(shape=(0, 0, 0)))
        assert_allclose(u, np.empty(shape=(0, 0, 3)))

        with assert_raises(ValueError, match='at least two-dimensional'):
            lu(np.array([]))

        a = np.array([[]])
        p, l, u = lu(a)
        assert_allclose(p, np.empty(shape=(0, 0)))
        assert_allclose(l, np.empty(shape=(1, 0)))
        assert_allclose(u, np.empty(shape=(0, 0)))

        a = np.array([[[]]])
        p, l, u = lu(a)
        assert_allclose(p, np.empty(shape=(1, 0, 0)))
        assert_allclose(l, np.empty(shape=(1, 1, 0)))
        assert_allclose(u, np.empty(shape=(1, 0, 0)))


class TestLUFactor:
    def setup_method(self):
        self.rng = np.random.default_rng(1682281250228846)

        self.a = np.array([[1, 2, 3], [1, 2, 3], [2, 5, 6]])
        self.ca = np.array([[1, 2, 3], [1, 2, 3], [2, 5j, 6]])
        # Those matrices are more robust to detect problems in permutation
        # matrices than the ones above
        self.b = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        self.cb = np.array([[1j, 2j, 3j], [4j, 5j, 6j], [7j, 8j, 9j]])

        # Reectangular matrices
        self.hrect = np.array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 12, 12]])
        self.chrect = np.array([[1, 2, 3, 4], [5, 6, 7, 8],
                                [9, 10, 12, 12]]) * 1.j

        self.vrect = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 12, 12]])
        self.cvrect = 1.j * np.array([[1, 2, 3],
                                      [4, 5, 6],
                                      [7, 8, 9],
                                      [10, 12, 12]])

        # Medium sizes matrices
        self.med = self.rng.random((30, 40))
        self.cmed = self.rng.random((30, 40)) + 1.j*self.rng.random((30, 40))

    def _test_common_lu_factor(self, data):
        l_and_u1, piv1 = lu_factor(data)
        (getrf,) = get_lapack_funcs(("getrf",), (data,))
        l_and_u2, piv2, _ = getrf(data, overwrite_a=False)
        assert_allclose(l_and_u1, l_and_u2)
        assert_allclose(piv1, piv2)

    # Simple tests.
    # For lu_factor gives a LinAlgWarning because these matrices are singular
    def test_hrectangular(self):
        self._test_common_lu_factor(self.hrect)

    def test_vrectangular(self):
        self._test_common_lu_factor(self.vrect)

    def test_hrectangular_complex(self):
        self._test_common_lu_factor(self.chrect)

    def test_vrectangular_complex(self):
        self._test_common_lu_factor(self.cvrect)

    # Bigger matrices
    def test_medium1(self):
        """Check lu decomposition on medium size, rectangular matrix."""
        self._test_common_lu_factor(self.med)

    def test_medium1_complex(self):
        """Check lu decomposition on medium size, rectangular matrix."""
        self._test_common_lu_factor(self.cmed)

    def test_check_finite(self):
        p, l, u = lu(self.a, check_finite=False)
        assert_allclose(p @ l @ u, self.a)

    def test_simple_known(self):
        # Ticket #1458
        for order in ['C', 'F']:
            A = np.array([[2, 1], [0, 1.]], order=order)
            LU, P = lu_factor(A)
            assert_allclose(LU, np.array([[2, 1], [0, 1]]))
            assert_array_equal(P, np.array([0, 1]))


class TestLUSolve:
    def setup_method(self):
        self.rng = np.random.default_rng(1682281250228846)

    def test_lu(self):
        a0 = self.rng.random((10, 10))
        b = self.rng.random((10,))

        for order in ['C', 'F']:
            a = np.array(a0, order=order)
            x1 = solve(a, b)
            lu_a = lu_factor(a)
            x2 = lu_solve(lu_a, b)
            assert_allclose(x1, x2)

    def test_check_finite(self):
        a = self.rng.random((10, 10))
        b = self.rng.random((10,))
        x1 = solve(a, b)
        lu_a = lu_factor(a, check_finite=False)
        x2 = lu_solve(lu_a, b, check_finite=False)
        assert_allclose(x1, x2)
