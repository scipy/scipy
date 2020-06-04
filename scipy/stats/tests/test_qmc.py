import copy
import numpy as np
from numpy.testing import assert_allclose, assert_equal, assert_almost_equal
from scipy.stats import qmc

import pytest


class TestUtils(object):
    def test_discrepancy(self):
        space_0 = [[0.1, 0.5], [0.2, 0.4], [0.3, 0.3], [0.4, 0.2], [0.5, 0.1]]
        space_1 = [[1, 3], [2, 6], [3, 2], [4, 5], [5, 1], [6, 4]]
        space_2 = [[1, 5], [2, 4], [3, 3], [4, 2], [5, 1], [6, 6]]

        corners = np.array([[0.5, 0.5], [6.5, 6.5]])

        assert_allclose(qmc.discrepancy(space_0), 0.1353, atol=1e-4)

        # From Fang et al. Design and modeling for computer experiments, 2006
        assert_allclose(qmc.discrepancy(space_1, corners), 0.0081, atol=1e-4)
        assert_allclose(qmc.discrepancy(space_2, corners), 0.0105, atol=1e-4)

    def test_update_discrepancy(self):
        space_0 = [[0.1, 0.5], [0.2, 0.4], [0.3, 0.3], [0.4, 0.2], [0.5, 0.1]]
        space_1 = [[1, 3], [2, 6], [3, 2], [4, 5], [5, 1], [6, 4]]

        # without bounds
        disc_init = qmc.discrepancy(space_0[:-1], iterative=True)
        disc_iter = qmc._update_discrepancy(space_0[-1], space_0[:-1],
                                            disc_init)

        assert_allclose(disc_iter, 0.1353, atol=1e-4)

        # with bounds
        corners = np.array([[0.5, 0.5], [6.5, 6.5]])

        disc_init = qmc.discrepancy(space_1[:-1], corners, iterative=True)
        disc_iter = qmc._update_discrepancy(space_1[-1], space_1[:-1],
                                            disc_init, bounds=corners)

        assert_allclose(disc_iter, 0.0081, atol=1e-4)

    def test_perm_discrepancy(self):
        doe_init = np.array([[1, 3], [2, 6], [3, 2], [4, 5], [5, 1], [6, 4]])
        corners = np.array([[0.5, 0.5], [6.5, 6.5]])
        disc_init = qmc.discrepancy(doe_init, corners)

        row_1, row_2, col = 5, 2, 1

        doe = copy.deepcopy(doe_init)
        doe[row_1, col], doe[row_2, col] = doe[row_2, col], doe[row_1, col]

        disc_valid = qmc.discrepancy(doe, corners)
        disc_perm = qmc._perturb_discrepancy(doe_init, row_1, row_2, col,
                                             disc_init, corners)

        assert_allclose(disc_valid, disc_perm)

    def test_n_primes(self):
        n_primes = qmc.n_primes(10)
        assert n_primes[-1] == 29

        n_primes = qmc.n_primes(168)
        assert n_primes[-1] == 997

        n_primes = qmc.n_primes(350)
        assert n_primes[-1] == 2357

    def test_primes(self):
        primes = qmc.primes_from_2_to(50)
        out = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
        assert_allclose(primes, out)


class TestQMC(object):
    def test_van_der_corput(self):
        sample = qmc.van_der_corput(10)
        out = [0., 0.5, 0.25, 0.75, 0.125, 0.625, 0.375, 0.875, 0.0625, 0.5625]
        assert_almost_equal(sample, out)

        sample = qmc.van_der_corput(5, start_index=3)
        out = [0.75, 0.125, 0.625, 0.375, 0.875]
        assert_almost_equal(sample, out)

    def test_halton(self):
        # without bounds
        sample = qmc.halton(dim=2, n_samples=5)

        out = np.array([[1/2, 1/3], [1/4, 2/3], [3/4, 1/9], [1/8, 4/9], [5/8, 7/9]])
        assert_almost_equal(sample, out, decimal=1)

        # with bounds
        corners = np.array([[0, 2], [10, 5]])
        sample = qmc.halton(dim=2, n_samples=5, bounds=corners)

        out = np.array([[5., 3.], [2.5, 4.], [7.5, 2.3], [1.25, 3.3], [6.25, 4.3]])
        assert_almost_equal(sample, out, decimal=1)

        # continuing
        sample = qmc.halton(dim=2, n_samples=3, bounds=corners, start_index=2)
        out = np.array([[7.5, 2.3], [1.25, 3.3], [6.25, 4.3]])
        assert_almost_equal(sample, out, decimal=1)


class TestLHS(object):
    def test_lhs(self):
        np.random.seed(123456)

        corners = np.array([[0, 2], [10, 5]])

        sample = qmc.latin_hypercube(dim=2, n_samples=5, bounds=corners)
        out = np.array([[5.7, 3.2], [5.5, 3.9], [5.2, 3.6],
                        [5.1, 3.3], [5.8, 4.1]])
        assert_almost_equal(sample, out, decimal=1)

        sample = qmc.latin_hypercube(dim=2, n_samples=5, centered=True)
        out = np.array([[0.1, 0.5], [0.3, 0.1], [0.7, 0.1],
                        [0.1, 0.1], [0.3, 0.7]])
        assert_almost_equal(sample, out, decimal=1)

    def test_orthogonal_lhs(self):
        np.random.seed(123456)

        corners = np.array([[0, 2], [10, 5]])

        sample = qmc.orthogonal_latin_hypercube(2, 5, bounds=corners)
        out = np.array([[3.933, 2.670], [7.794, 4.031], [4.520, 2.129],
                        [0.253, 4.976], [8.753, 3.249]])
        assert_almost_equal(sample, out, decimal=1)

        # Checking independency of the random numbers generated
        n_samples = 500
        sample = qmc.orthogonal_latin_hypercube(dim=2, n_samples=n_samples)
        min_b = 50  # number of bins
        bins = np.linspace(0, 1, min(min_b, n_samples) + 1)
        hist = np.histogram(sample[:, 0], bins=bins)
        out = np.array([n_samples / min_b] * min_b)
        assert_equal(hist[0], out)

        hist = np.histogram(sample[:, 1], bins=bins)
        assert_equal(hist[0], out)

    @pytest.mark.xfail(raises=AssertionError, reason='Global optimization')
    def test_optimal_design(self):
        np.random.seed(123456)

        start_design = qmc.orthogonal_latin_hypercube(2, 5)
        sample = qmc.optimal_design(dim=2, n_samples=5,
                                    start_design=start_design)
        out = np.array([[0.025, 0.223], [0.779, 0.677], [0.452, 0.043],
                        [0.393, 0.992], [0.875, 0.416]])
        assert_almost_equal(sample, out, decimal=1)

        corners = np.array([[0, 2], [10, 5]])
        sample = qmc.optimal_design(2, 5, bounds=corners)
        out = np.array([[5.189, 4.604], [3.553, 2.344], [6.275, 3.947],
                        [0.457, 3.554], [9.705, 2.636]])
        assert_almost_equal(sample, out, decimal=1)

        sample = qmc.optimal_design(2, 5, niter=2)
        out = np.array([[0.681, 0.231], [0.007, 0.719], [0.372, 0.101],
                        [0.550, 0.456], [0.868, 0.845]])
        assert_almost_equal(sample, out, decimal=1)

        sample = qmc.optimal_design(2, 5, bounds=corners, force=True)
        out = np.array([[8.610, 4.303], [5.318, 3.498], [7.323, 2.288],
                        [1.135, 2.657], [3.561, 4.938]])
        assert_almost_equal(sample, out, decimal=1)

        sample = qmc.optimal_design(2, 5, bounds=corners, optimization=False)
        out = np.array([[1.052, 4.218], [2.477, 2.987], [7.616, 4.527],
                        [9.134, 3.393], [4.064, 2.430]])
        assert_almost_equal(sample, out, decimal=1)

        sample = qmc.optimal_design(2, 5, bounds=corners, optimization=False, niter=2)
        print(sample)
        out = np.array([[7.902, 2.166], [4.915, 2.741], [3.797, 3.365],
                        [0.602, 4.896], [9.880, 4.215]])
        assert_almost_equal(sample, out, decimal=1)
