import copy
import numpy as np
from numpy.testing import assert_allclose, assert_equal, assert_almost_equal
from scipy.stats import qmc

import pytest


class TestUtils(object):
    def test_scale(self):
        space = [[0, 0], [1, 1], [0.5, 0.5]]
        corners = np.array([[-2, 0], [6, 5]])
        out = [[-2, 0], [6, 5], [2, 2.5]]

        scaled_space = qmc.scale(space, bounds=corners)

        assert_allclose(scaled_space, out)

        scaled_back_space = qmc.scale(scaled_space, bounds=corners, reverse=True)
        assert_allclose(scaled_back_space, space)

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
        engine = qmc.Halton(k_vars=2)
        sample = engine.random(n_samples=5)

        out = np.array([[1 / 2, 1 / 3], [1 / 4, 2 / 3], [3 / 4, 1 / 9], [1 / 8, 4 / 9], [5 / 8, 7 / 9]])
        assert_almost_equal(sample, out, decimal=1)

        assert engine.num_generated == 5

        # with bounds
        corners = np.array([[0, 2], [10, 5]])
        sample = qmc.scale(sample, bounds=corners)

        out = np.array([[5., 3.], [2.5, 4.], [7.5, 2.3], [1.25, 3.3], [6.25, 4.3]])
        assert_almost_equal(sample, out, decimal=1)

        # continuing
        engine = qmc.Halton(k_vars=2)
        engine.fast_forward(2)
        sample = engine.random(n_samples=3)
        sample = qmc.scale(sample, bounds=corners)

        out = np.array([[7.5, 2.3], [1.25, 3.3], [6.25, 4.3]])
        assert_almost_equal(sample, out, decimal=1)


class TestLHS(object):
    def test_lhs(self):
        np.random.seed(123456)

        corners = np.array([[0, 2], [10, 5]])

        lhs = qmc.LatinHypercube(k_vars=2)
        sample = lhs.random(n_samples=5)
        sample = qmc.scale(sample, bounds=corners)
        out = np.array([[5.7, 3.2], [5.5, 3.9], [5.2, 3.6],
                        [5.1, 3.3], [5.8, 4.1]])
        assert_almost_equal(sample, out, decimal=1)

        lhs = qmc.LatinHypercube(k_vars=2, centered=True)
        sample = lhs.random(n_samples=5)
        out = np.array([[0.1, 0.5], [0.3, 0.1], [0.7, 0.1],
                        [0.1, 0.1], [0.3, 0.7]])
        assert_almost_equal(sample, out, decimal=1)

    def test_orthogonal_lhs(self):
        np.random.seed(123456)

        corners = np.array([[0, 2], [10, 5]])

        olhs = qmc.OrthogonalLatinHypercube(k_vars=2)
        sample = olhs.random(n_samples=5)
        sample = qmc.scale(sample, bounds=corners)
        out = np.array([[3.933, 2.670], [7.794, 4.031], [4.520, 2.129],
                        [0.253, 4.976], [8.753, 3.249]])
        assert_almost_equal(sample, out, decimal=1)

        # Checking independency of the random numbers generated
        n_samples = 500
        sample = olhs.random(n_samples=n_samples)
        min_b = 50  # number of bins
        bins = np.linspace(0, 1, min(min_b, n_samples) + 1)
        hist = np.histogram(sample[:, 0], bins=bins)
        out = np.array([n_samples / min_b] * min_b)
        assert_equal(hist[0], out)

        hist = np.histogram(sample[:, 1], bins=bins)
        assert_equal(hist[0], out)

    def test_optimal_design(self):
        seed = 123456

        olhs = qmc.OrthogonalLatinHypercube(k_vars=2, seed=seed)
        sample_ref = olhs.random(n_samples=20)
        disc_ref = qmc.discrepancy(sample_ref)

        optimal_1 = qmc.OptimalDesign(k_vars=2, start_design=sample_ref, seed=seed)
        sample_1 = optimal_1.random(n_samples=20)
        disc_1 = qmc.discrepancy(sample_1)

        assert disc_1 < disc_ref

        optimal_2 = qmc.OptimalDesign(k_vars=2, start_design=sample_ref, niter=2, seed=seed)
        sample_2 = optimal_2.random(n_samples=20)
        disc_2 = qmc.discrepancy(sample_2)
        assert disc_2 < disc_1

        # resample is like doing another iteration
        sample_1 = optimal_1.random(n_samples=20)
        assert_allclose(sample_1, sample_2)

        optimal_3 = qmc.OptimalDesign(k_vars=2, start_design=sample_ref, force=True)
        sample_3 = optimal_3.random(n_samples=20)
        disc_3 = qmc.discrepancy(sample_3)
        assert disc_3 < disc_ref

        # no optimization
        optimal_4 = qmc.OptimalDesign(k_vars=2, start_design=sample_ref,
                                      optimization=False, niter=100, seed=seed)
        sample_4 = optimal_4.random(n_samples=20)
        disc_4 = qmc.discrepancy(sample_4)

        assert disc_4 < disc_ref
