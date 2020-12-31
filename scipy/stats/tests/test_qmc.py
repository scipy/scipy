import copy

import pytest
import numpy as np
from numpy.testing import (assert_allclose, assert_almost_equal, assert_,
                           assert_equal, assert_array_almost_equal,
                           assert_array_equal)
from scipy.stats import shapiro
from scipy.optimize import basinhopping

from scipy.stats._sobol import _test_find_index
from scipy.stats import qmc
from scipy.stats._qmc import (van_der_corput, n_primes, primes_from_2_to,
                              _perturb_discrepancy, _update_discrepancy)


class TestUtils(object):
    def test_scale(self):
        space = [[0, 0], [1, 1], [0.5, 0.5]]
        corners = np.array([[-2, 0], [6, 5]])
        out = [[-2, 0], [6, 5], [2, 2.5]]

        scaled_space = qmc.scale(space, bounds=corners)

        assert_allclose(scaled_space, out)

        scaled_back_space = qmc.scale(scaled_space, bounds=corners,
                                      reverse=True)
        assert_allclose(scaled_back_space, space)

    def test_discrepancy(self):
        space_1 = np.array([[1, 3], [2, 6], [3, 2], [4, 5], [5, 1], [6, 4]])
        space_1 = (2.0 * space_1 - 1.0) / (2.0 * 6.0)
        space_2 = np.array([[1, 5], [2, 4], [3, 3], [4, 2], [5, 1], [6, 6]])
        space_2 = (2.0 * space_2 - 1.0) / (2.0 * 6.0)

        # From Fang et al. Design and modeling for computer experiments, 2006
        assert_allclose(qmc.discrepancy(space_1), 0.0081, atol=1e-4)
        assert_allclose(qmc.discrepancy(space_2), 0.0105, atol=1e-4)

        # From Zhou Y.-D. et al. Mixture discrepancy for quasi-random point
        # sets. Journal of Complexity, 29 (3-4), pp. 283-301, 2013.
        sample = np.array([[2, 1, 1, 2, 2, 2],
                           [1, 2, 2, 2, 2, 2],
                           [2, 1, 1, 1, 1, 1],
                           [1, 1, 1, 1, 2, 2],
                           [1, 2, 2, 2, 1, 1],
                           [2, 2, 2, 2, 1, 1],
                           [2, 2, 2, 1, 2, 2]])
        sample = (2.0 * sample - 1.0) / (2.0 * 2.0)

        assert_allclose(qmc.discrepancy(sample, method='MD'), 2.5000,
                        atol=1e-4)
        assert_allclose(qmc.discrepancy(sample, method='WD'), 1.3680,
                        atol=1e-4)
        assert_allclose(qmc.discrepancy(sample, method='CD'), 0.3172,
                        atol=1e-4)
        assert_allclose(qmc.discrepancy(sample, method='star'), 0.037451,
                        atol=1e-4)

        with pytest.raises(ValueError, match=r"toto is not a valid method."):
            qmc.discrepancy(sample, method='toto')

    def test_update_discrepancy(self):
        space_1 = np.array([[1, 3], [2, 6], [3, 2], [4, 5], [5, 1], [6, 4]])
        space_1 = (2.0 * space_1 - 1.0) / (2.0 * 6.0)

        disc_init = qmc.discrepancy(space_1[:-1], iterative=True)
        disc_iter = _update_discrepancy(space_1[-1], space_1[:-1],
                                        disc_init)

        assert_allclose(disc_iter, 0.0081, atol=1e-4)

    def test_perm_discrepancy(self):
        doe_init = np.array([[1, 3], [2, 6], [3, 2], [4, 5], [5, 1], [6, 4]])
        doe_init = (2.0 * doe_init - 1.0) / (2.0 * 6.0)

        disc_init = qmc.discrepancy(doe_init)

        row_1, row_2, col = 5, 2, 1

        doe = copy.deepcopy(doe_init)
        doe[row_1, col], doe[row_2, col] = doe[row_2, col], doe[row_1, col]

        disc_valid = qmc.discrepancy(doe)
        disc_perm = _perturb_discrepancy(doe_init, row_1, row_2, col,
                                         disc_init)

        assert_allclose(disc_valid, disc_perm)

    def test_n_primes(self):
        primes = n_primes(10)
        assert primes[-1] == 29

        primes = n_primes(168)
        assert primes[-1] == 997

        primes = n_primes(350)
        assert primes[-1] == 2357

    def test_primes(self):
        primes = primes_from_2_to(50)
        out = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
        assert_allclose(primes, out)


class TestQMC(object):
    def test_van_der_corput(self):
        seed = np.random.RandomState(12345)
        sample = van_der_corput(10, seed=seed)
        out = [0., 0.5, 0.25, 0.75, 0.125, 0.625, 0.375, 0.875, 0.0625, 0.5625]
        assert_almost_equal(sample, out)

        sample = van_der_corput(7, start_index=3, seed=seed)
        assert_almost_equal(sample, out[3:])

    def test_van_der_corput_scramble(self):
        seed = np.random.RandomState(123456)
        out = van_der_corput(10, scramble=True, seed=seed)

        seed = np.random.RandomState(123456)
        sample = van_der_corput(7, start_index=3, scramble=True,
                                seed=seed)
        assert_almost_equal(sample, out[3:])

    def test_halton(self):
        seed = np.random.RandomState(123456)
        engine = qmc.Halton(d=2, scramble=False, seed=seed)
        sample = engine.random(n=3)

        out = np.array([[0, 0], [1 / 2, 1 / 3], [1 / 4, 2 / 3], [3 / 4, 1 / 9],
                        [1 / 8, 4 / 9], [5 / 8, 7 / 9]])
        assert_almost_equal(sample, out[:3], decimal=1)

        assert engine.num_generated == 3

        # continuing
        sample = engine.random(n=3)
        assert_almost_equal(sample, out[3:], decimal=1)

        # reset
        engine.reset()
        sample = engine.random(n=6)
        assert_almost_equal(sample, out, decimal=1)

        # continuing with fast_forward
        engine = qmc.Halton(d=2, scramble=False, seed=seed)
        engine.fast_forward(2)
        sample = engine.random(n=4)

        assert_almost_equal(sample, out[2:], decimal=1)

    def test_halton_scramble(self):
        seed = np.random.RandomState(123456)
        engine = qmc.Halton(d=2, scramble=True, seed=seed)
        out = engine.random(n=10)

        seed = np.random.RandomState(123456)
        engine = qmc.Halton(d=2, scramble=True, seed=seed)
        sample = engine.random(n=3)
        assert_almost_equal(sample, out[:3], decimal=1)

        # continuing
        sample = engine.random(n=7)
        assert_almost_equal(sample, out[3:], decimal=1)

        # reset
        engine.reset()
        sample = engine.random(n=10)
        assert_almost_equal(sample, out, decimal=1)

        # continuing with fast_forward
        seed = np.random.RandomState(123456)
        engine = qmc.Halton(d=2, scramble=True, seed=seed)
        engine.fast_forward(2)
        sample = engine.random(n=8)

        assert_almost_equal(sample, out[2:], decimal=1)


class TestLHS:
    def test_lhs(self):
        seed = np.random.RandomState(123456)
        corners = np.array([[0, 2], [10, 5]])

        lhs = qmc.LatinHypercube(d=2, seed=seed)
        sample = lhs.random(n=5)
        sample = qmc.scale(sample, bounds=corners)
        out = np.array([[5.7, 3.2], [5.5, 3.9], [5.2, 3.6],
                        [5.1, 3.3], [5.8, 4.1]])
        assert_almost_equal(sample, out, decimal=1)

        lhs = qmc.LatinHypercube(d=2, centered=True, seed=seed)
        sample = lhs.random(n=5)
        out = np.array([[0.1, 0.5], [0.3, 0.1], [0.7, 0.1],
                        [0.1, 0.1], [0.3, 0.7]])
        assert_almost_equal(sample, out, decimal=1)

    def test_orthogonal_lhs(self):
        seed = np.random.RandomState(123456)
        corners = np.array([[0, 2], [10, 5]])

        olhs = qmc.OrthogonalLatinHypercube(d=2, seed=seed)
        sample = olhs.random(n=5)
        sample = qmc.scale(sample, bounds=corners)
        out = np.array([[3.933, 2.670], [7.794, 4.031], [4.520, 2.129],
                        [0.253, 4.976], [8.753, 3.249]])
        assert_almost_equal(sample, out, decimal=1)

        # Checking independency of the random numbers generated
        n_samples = 500
        sample = olhs.random(n=n_samples)
        min_b = 50  # number of bins
        bins = np.linspace(0, 1, min(min_b, n_samples) + 1)
        hist = np.histogram(sample[:, 0], bins=bins)
        out = np.array([n_samples / min_b] * min_b)
        assert_equal(hist[0], out)

        hist = np.histogram(sample[:, 1], bins=bins)
        assert_equal(hist[0], out)

    def test_optimal_design(self):
        # base discrepancy as a reference for testing OptimalDesign is better
        seed = np.random.RandomState(123456)
        olhs = qmc.OrthogonalLatinHypercube(d=2, seed=seed)
        sample_ref = olhs.random(n=20)
        disc_ref = qmc.discrepancy(sample_ref)

        # all defaults
        seed = np.random.RandomState(123456)
        optimal_ = qmc.OptimalDesign(d=2, seed=seed)
        sample_ = optimal_.random(n=20)
        disc_ = qmc.discrepancy(sample_)

        assert disc_ < disc_ref

        # using an initial sample
        seed = np.random.RandomState(123456)
        optimal_1 = qmc.OptimalDesign(d=2, start_design=sample_ref, seed=seed)
        sample_1 = optimal_1.random(n=20)
        disc_1 = qmc.discrepancy(sample_1)

        assert disc_1 < disc_ref

        # 5 iterations is better than 1
        seed = np.random.RandomState(123456)
        optimal_2 = qmc.OptimalDesign(d=2, start_design=sample_ref, niter=5,
                                      seed=seed)
        sample_2 = optimal_2.random(n=20)
        disc_2 = qmc.discrepancy(sample_2)
        assert disc_2 < disc_1

        # another optimization method
        def method(func, x0, bounds):
            seed = np.random.RandomState(123456)
            minimizer_kwargs = {"method": "L-BFGS-B", "bounds": bounds}
            _ = basinhopping(func, x0, niter=100,
                             minimizer_kwargs=minimizer_kwargs,
                             seed=seed)
        seed = np.random.RandomState(123456)
        optimal_3 = qmc.OptimalDesign(d=2, start_design=sample_ref,
                                      method=method, seed=seed)
        sample_3 = optimal_3.random(n=20)
        disc_3 = qmc.discrepancy(sample_3)
        assert disc_3 < disc_ref


class TestMultinomialQMC:
    def test_MultinomialNegativePs(self):
        p = np.array([0.12, 0.26, -0.05, 0.35, 0.22])
        with pytest.raises(ValueError, match=r"Elements of pvals must "
                                             r"be non-negative."):
            qmc.multinomial_qmc(10, p)

    def test_MultinomialSumOfPTooLarge(self):
        p = np.array([0.12, 0.26, 0.1, 0.35, 0.22])
        with pytest.raises(ValueError, match=r"Elements of pvals must sum "
                                             r"to 1."):
            qmc.multinomial_qmc(10, p)

    @pytest.mark.filterwarnings('ignore::UserWarning')
    def test_MultinomialBasicDraw(self):
        seed = np.random.RandomState(12345)
        p = np.array([0.12, 0.26, 0.05, 0.35, 0.22])
        expected = np.array([12, 25, 6, 35, 22])
        assert_array_equal(qmc.multinomial_qmc(100, p, seed=seed), expected)

    def test_MultinomialDistribution(self):
        seed = np.random.RandomState(12345)
        p = np.array([0.12, 0.26, 0.05, 0.35, 0.22])
        draws = qmc.multinomial_qmc(8192, p, seed=seed)
        assert_array_almost_equal(draws / np.sum(draws), p, decimal=4)

    def test_FindIndex(self):
        p_cumulative = np.array([0.1, 0.4, 0.45, 0.6, 0.75, 0.9, 0.99, 1.0])
        size = len(p_cumulative)
        assert_equal(_test_find_index(p_cumulative, size, 0.0), 0)
        assert_equal(_test_find_index(p_cumulative, size, 0.4), 2)
        assert_equal(_test_find_index(p_cumulative, size, 0.44999), 2)
        assert_equal(_test_find_index(p_cumulative, size, 0.45001), 3)
        assert_equal(_test_find_index(p_cumulative, size, 1.0), size - 1)

    @pytest.mark.filterwarnings('ignore::UserWarning')
    def test_other_engine(self):
        # same as test_MultinomialBasicDraw with different engine
        seed = np.random.RandomState(12345)
        p = np.array([0.12, 0.26, 0.05, 0.35, 0.22])
        expected = np.array([12, 25, 6, 35, 22])
        engine = qmc.Sobol(1, scramble=True, seed=seed)
        assert_array_equal(qmc.multinomial_qmc(100, p, engine=engine,
                                               seed=seed),
                           expected)


class TestNormalQMC:
    def test_NormalQMC(self):
        # d = 1
        seed = np.random.RandomState(123456)
        engine = qmc.NormalQMC(d=1, seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 1))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 1))
        # d = 2
        engine = qmc.NormalQMC(d=2, seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 2))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 2))

    def test_NormalQMCInvTransform(self):
        # d = 1
        seed = np.random.RandomState(123456)
        engine = qmc.NormalQMC(d=1, inv_transform=True, seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 1))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 1))
        # d = 2
        engine = qmc.NormalQMC(d=2, inv_transform=True, seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 2))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 2))

    def test_other_engine(self):
        seed = np.random.RandomState(123456)
        engine = qmc.NormalQMC(d=2, engine=qmc.Sobol(d=2, scramble=False,
                                                     seed=seed),
                               inv_transform=True, seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 2))

    def test_NormalQMCSeeded(self):
        # test even dimension
        seed = np.random.RandomState(12345)
        engine = qmc.NormalQMC(d=2, inv_transform=False, seed=seed)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [[-0.943472, 0.405116], [-0.63099602, -1.32950772]]
        )
        assert_array_almost_equal(samples, samples_expected)

        # test odd dimension
        seed = np.random.RandomState(12345)
        engine = qmc.NormalQMC(d=3, inv_transform=False, seed=seed)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [
                [-0.943472, 0.405116, 0.268828],
                [1.83169884, -1.40473647, 0.24334828],
            ]
        )
        assert_array_almost_equal(samples, samples_expected)

    def test_NormalQMCSeededInvTransform(self):
        # test even dimension
        seed = np.random.RandomState(12345)
        engine = qmc.NormalQMC(d=2, seed=seed, inv_transform=True)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [[0.228309, -0.162516], [-0.41622922, 0.46622792]]
        )
        assert_array_almost_equal(samples, samples_expected)

        # test odd dimension
        seed = np.random.RandomState(12345)
        engine = qmc.NormalQMC(d=3, seed=seed, inv_transform=True)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [
                [0.228309, -0.162516, 0.167352],
                [-1.40525266, 1.37652443, -0.8519666],
            ]
        )
        assert_array_almost_equal(samples, samples_expected)

    def test_NormalQMCShapiro(self):
        seed = np.random.RandomState(12345)
        engine = qmc.NormalQMC(d=2, seed=seed)
        samples = engine.random(n=256)
        assert_(all(np.abs(samples.mean(axis=0)) < 1e-2))
        assert_(all(np.abs(samples.std(axis=0) - 1) < 1e-2))
        # perform Shapiro-Wilk test for normality
        for i in (0, 1):
            _, pval = shapiro(samples[:, i])
            assert_(pval > 0.9)
        # make sure samples are uncorrelated
        cov = np.cov(samples.transpose())
        assert_(np.abs(cov[0, 1]) < 1e-2)

    def test_NormalQMCShapiroInvTransform(self):
        seed = np.random.RandomState(12345)
        engine = qmc.NormalQMC(d=2, seed=seed, inv_transform=True)
        samples = engine.random(n=256)
        assert_(all(np.abs(samples.mean(axis=0)) < 1e-2))
        assert_(all(np.abs(samples.std(axis=0) - 1) < 1e-2))
        # perform Shapiro-Wilk test for normality
        for i in (0, 1):
            _, pval = shapiro(samples[:, i])
            assert_(pval > 0.9)
        # make sure samples are uncorrelated
        cov = np.cov(samples.transpose())
        assert_(np.abs(cov[0, 1]) < 1e-2)


class TestMultivariateNormalQMC:
    def test_MultivariateNormalQMCNonPSD(self):
        # try with non-psd, non-pd cov and expect an assertion error
        with pytest.raises(ValueError, match=r"Covariance matrix not PSD."):
            seed = np.random.RandomState(123456)
            qmc.MultivariateNormalQMC([0, 0], [[1, 2], [2, 1]], seed=seed)

    def test_MultivariateNormalQMCNonPD(self):
        # try with non-pd but psd cov; should work
        seed = np.random.RandomState(123456)
        engine = qmc.MultivariateNormalQMC(
            [0, 0, 0], [[1, 0, 1], [0, 1, 1], [1, 1, 2]],
            seed=seed
        )
        assert_(engine._corr_matrix is not None)

    def test_MultivariateNormalQMCSymmetric(self):
        # try with non-symmetric cov and expect an error
        with pytest.raises(ValueError, match=r"Covariance matrix is not "
                                             r"symmetric."):
            seed = np.random.RandomState(123456)
            qmc.MultivariateNormalQMC([0, 0], [[1, 0], [2, 1]], seed=seed)

    def test_MultivariateNormalQMCDim(self):
        # incompatible dimension of mean/cov
        with pytest.raises(ValueError, match=r"Dimension mismatch between "
                                             r"mean and covariance."):
            seed = np.random.RandomState(123456)
            qmc.MultivariateNormalQMC([0], [[1, 0], [0, 1]], seed=seed)

    def test_MultivariateNormalQMC(self):
        # d = 1 scalar
        seed = np.random.RandomState(123456)
        engine = qmc.MultivariateNormalQMC(mean=0, cov=5, seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 1))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 1))

        # d = 2 list
        engine = qmc.MultivariateNormalQMC(mean=[0, 1], cov=[[1, 0], [0, 1]],
                                           seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 2))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 2))

        # d = 3 np.array
        mean = np.array([0, 1, 2])
        cov = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        engine = qmc.MultivariateNormalQMC(mean, cov, seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 3))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 3))

    def test_MultivariateNormalQMCInvTransform(self):
        # d = 1 scalar
        seed = np.random.RandomState(123456)
        engine = qmc.MultivariateNormalQMC(mean=0, cov=5, inv_transform=True,
                                           seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 1))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 1))

        # d = 2 list
        engine = qmc.MultivariateNormalQMC(
            mean=[0, 1], cov=[[1, 0], [0, 1]], inv_transform=True,
            seed=seed
        )
        samples = engine.random()
        assert_equal(samples.shape, (1, 2))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 2))

        # d = 3 np.array
        mean = np.array([0, 1, 2])
        cov = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        engine = qmc.MultivariateNormalQMC(mean, cov, inv_transform=True,
                                           seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 3))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 3))

    def test_MultivariateNormalQMCSeeded(self):
        # test even dimension
        seed = np.random.RandomState(12345)
        np.random.seed(54321)
        a = np.random.randn(2, 2)
        A = a @ a.transpose() + np.diag(np.random.rand(2))
        engine = qmc.MultivariateNormalQMC(np.array([0, 0]), A,
                                           inv_transform=False, seed=seed)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [[-1.010703, -0.324223], [-0.67595995, -2.27437872]]
        )
        assert_array_almost_equal(samples, samples_expected)

        # test odd dimension
        seed = np.random.RandomState(12345)
        np.random.seed(54321)
        a = np.random.randn(3, 3)
        A = a @ a.transpose() + np.diag(np.random.rand(3))
        engine = qmc.MultivariateNormalQMC(np.array([0, 0, 0]), A,
                                           inv_transform=False, seed=seed)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [
                [-1.056834, 2.493251, 0.114556],
                [2.05178452, -6.35744194, 0.67944512],
            ]
        )
        assert_array_almost_equal(samples, samples_expected)

    def test_MultivariateNormalQMCSeededInvTransform(self):
        # test even dimension
        seed = np.random.RandomState(12345)
        np.random.seed(54321)
        a = np.random.randn(2, 2)
        A = a @ a.transpose() + np.diag(np.random.rand(2))
        engine = qmc.MultivariateNormalQMC(
            np.array([0, 0]), A, seed=seed, inv_transform=True
        )
        samples = engine.random(n=2)
        samples_expected = np.array(
            [[0.244578, -0.004441], [-0.44588916, 0.22657776]]
        )
        assert_array_almost_equal(samples, samples_expected)

        # test odd dimension
        seed = np.random.RandomState(12345)
        np.random.seed(54321)
        a = np.random.randn(3, 3)
        A = a @ a.transpose() + np.diag(np.random.rand(3))
        engine = qmc.MultivariateNormalQMC(
            np.array([0, 0, 0]), A, seed=seed, inv_transform=True
        )
        samples = engine.random(n=2)
        samples_expected = np.array(
            [
                [0.255741, -0.761559, 0.234236],
                [-1.5740992, 5.61057598, -1.28218525],
            ]
        )
        assert_array_almost_equal(samples, samples_expected)

    def test_MultivariateNormalQMCShapiro(self):
        # test the standard case
        seed = np.random.RandomState(12345)
        engine = qmc.MultivariateNormalQMC(
            mean=[0, 0], cov=[[1, 0], [0, 1]], seed=seed
        )
        samples = engine.random(n=256)
        assert_(all(np.abs(samples.mean(axis=0)) < 1e-2))
        assert_(all(np.abs(samples.std(axis=0) - 1) < 1e-2))
        # perform Shapiro-Wilk test for normality
        for i in (0, 1):
            _, pval = shapiro(samples[:, i])
            assert_(pval > 0.9)
        # make sure samples are uncorrelated
        cov = np.cov(samples.transpose())
        assert_(np.abs(cov[0, 1]) < 1e-2)

        # test the correlated, non-zero mean case
        seed = np.random.RandomState(12345)
        engine = qmc.MultivariateNormalQMC(
            mean=[1.0, 2.0], cov=[[1.5, 0.5], [0.5, 1.5]], seed=seed
        )
        samples = engine.random(n=256)
        assert_(all(np.abs(samples.mean(axis=0) - [1, 2]) < 1e-2))
        assert_(all(np.abs(samples.std(axis=0) - np.sqrt(1.5)) < 1e-2))
        # perform Shapiro-Wilk test for normality
        for i in (0, 1):
            _, pval = shapiro(samples[:, i])
            assert_(pval > 0.9)
        # check covariance
        cov = np.cov(samples.transpose())
        assert_(np.abs(cov[0, 1] - 0.5) < 1e-2)

    def test_MultivariateNormalQMCShapiroInvTransform(self):
        # test the standard case
        seed = np.random.RandomState(12345)
        engine = qmc.MultivariateNormalQMC(
            mean=[0, 0], cov=[[1, 0], [0, 1]], seed=seed, inv_transform=True
        )
        samples = engine.random(n=256)
        assert_(all(np.abs(samples.mean(axis=0)) < 1e-2))
        assert_(all(np.abs(samples.std(axis=0) - 1) < 1e-2))
        # perform Shapiro-Wilk test for normality
        for i in (0, 1):
            _, pval = shapiro(samples[:, i])
            assert_(pval > 0.9)
        # make sure samples are uncorrelated
        cov = np.cov(samples.transpose())
        assert_(np.abs(cov[0, 1]) < 1e-2)

        # test the correlated, non-zero mean case
        seed = np.random.RandomState(12345)
        engine = qmc.MultivariateNormalQMC(
            mean=[1.0, 2.0],
            cov=[[1.5, 0.5], [0.5, 1.5]],
            seed=seed,
            inv_transform=True,
        )
        samples = engine.random(n=256)
        assert_(all(np.abs(samples.mean(axis=0) - [1, 2]) < 1e-2))
        assert_(all(np.abs(samples.std(axis=0) - np.sqrt(1.5)) < 1e-2))
        # perform Shapiro-Wilk test for normality
        for i in (0, 1):
            _, pval = shapiro(samples[:, i])
            assert_(pval > 0.9)
        # check covariance
        cov = np.cov(samples.transpose())
        assert_(np.abs(cov[0, 1] - 0.5) < 1e-2)

    def test_MultivariateNormalQMCDegenerate(self):
        # X, Y iid standard Normal and Z = X + Y, random vector (X, Y, Z)
        seed = np.random.RandomState(12345)
        engine = qmc.MultivariateNormalQMC(
            mean=[0.0, 0.0, 0.0],
            cov=[[1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 2.0]],
            seed=seed,
        )
        samples = engine.random(n=512)
        assert_(all(np.abs(samples.mean(axis=0)) < 1e-2))
        assert_(np.abs(np.std(samples[:, 0]) - 1) < 1e-2)
        assert_(np.abs(np.std(samples[:, 1]) - 1) < 1e-2)
        assert_(np.abs(np.std(samples[:, 2]) - np.sqrt(2)) < 1e-2)
        for i in (0, 1, 2):
            _, pval = shapiro(samples[:, i])
            assert_(pval > 0.8)
        cov = np.cov(samples.transpose())
        assert_(np.abs(cov[0, 1]) < 1e-2)
        assert_(np.abs(cov[0, 2] - 1) < 1e-2)
        # check to see if X + Y = Z almost exactly
        assert_(
            all(np.abs(samples[:, 0] + samples[:, 1] - samples[:, 2]) < 1e-5)
        )
