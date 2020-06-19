from collections import Counter

import numpy as np
from scipy.stats._sobol import _test_find_index
from scipy.stats.qmc import (Sobol, MultivariateNormalQMC,
                             NormalQMC, multinomial_qmc)
from scipy.stats import shapiro

from numpy.testing import (assert_, assert_equal, assert_array_almost_equal,
                           assert_array_equal)
from pytest import raises as assert_raises


class TestSobol:
    # set maxDiff to None to show all differences when tests fail
    maxDiff = None

    def setUp(self):
        engine_unscrambled_1d = Sobol(1)
        self.draws_unscrambled_1d = engine_unscrambled_1d.random(10)
        engine_unscrambled_3d = Sobol(3)
        self.draws_unscrambled_3d = engine_unscrambled_3d.random(10)
        engine_scrambled_1d = Sobol(1, scramble=True, seed=12345)
        self.draws_scrambled_1d = engine_scrambled_1d.random(10)
        engine_scrambled_3d = Sobol(3, scramble=True, seed=12345)
        self.draws_scrambled_3d = engine_scrambled_3d.random(10)

    def test_Unscrambled1DSobol(self):
        self.setUp()
        expected = [
            0.5, 0.75, 0.25, 0.375, 0.875, 0.625, 0.125, 0.1875, 0.6875, 0.9375,
        ]
        assert_equal(self.draws_unscrambled_1d.shape[0], 10)
        assert_equal(self.draws_unscrambled_1d.shape[1], 1)
        assert_array_equal(self.draws_unscrambled_1d.flatten(), np.array(expected))

    def test_Unscrambled3DSobol(self):
        self.setUp()
        expected_dim3 = [
            0.5,
            0.75,
            0.25,
            0.625,
            0.125,
            0.375,
            0.875,
            0.3125,
            0.8125,
            0.5625,
        ]
        assert_equal(self.draws_unscrambled_3d.shape[0], 10)
        assert_equal(self.draws_unscrambled_3d.shape[1], 3)
        assert_array_equal(self.draws_unscrambled_3d[:, 2], np.array(expected_dim3))
        assert_array_equal(
            self.draws_unscrambled_3d[:, 0], self.draws_unscrambled_1d.flatten()
        )

    def test_Unscrambled3DAsyncSobol(self):
        self.setUp()
        engine_unscrambled_3d = Sobol(3)
        draws = np.vstack([engine_unscrambled_3d.random() for i in range(10)])
        assert_array_equal(self.draws_unscrambled_3d, draws)

    def test_UnscrambledFastForwardAndResetSobol(self):
        self.setUp()
        engine_unscrambled_3d = Sobol(3).fast_forward(5)
        draws = engine_unscrambled_3d.random(5)
        assert_array_equal(self.draws_unscrambled_3d[5:10, :], draws)

        engine_unscrambled_3d.reset()
        even_draws = []
        for i in range(10):
            if i % 2 == 0:
                even_draws.append(engine_unscrambled_3d.random())
            else:
                engine_unscrambled_3d.fast_forward(1)
        assert_array_equal(
            self.draws_unscrambled_3d[[i for i in range(10) if i % 2 == 0]],
            np.vstack(even_draws),
        )

    def test_UnscrambledHighDimSobol(self):
        engine = Sobol(1111)
        count1 = Counter(engine.random().flatten().tolist())
        count2 = Counter(engine.random().flatten().tolist())
        count3 = Counter(engine.random().flatten().tolist())
        assert_equal(count1, Counter({0.5: 1111}))
        assert_equal(count2, Counter({0.25: 580, 0.75: 531}))
        assert_equal(count3, Counter({0.25: 531, 0.75: 580}))

    def test_UnscrambledSobolBounds(self):
        engine = Sobol(1111)
        draws = engine.random(1000)
        assert_(np.all(draws >= 0))
        assert_(np.all(draws <= 1))

    def test_UnscrambledDistributionSobol(self):
        engine = Sobol(1111)
        draws = engine.random(1000)
        assert_array_almost_equal(
            np.mean(draws, axis=0), np.repeat(0.5, 1111), decimal=2
        )
        assert_array_almost_equal(
            np.percentile(draws, 25, axis=0), np.repeat(0.25, 1111), decimal=2
        )
        assert_array_almost_equal(
            np.percentile(draws, 75, axis=0), np.repeat(0.75, 1111), decimal=2
        )

    def test_Scrambled1DSobol(self):
        self.setUp()
        expected = [
            0.46784395,
            0.03562005,
            0.91319746,
            0.86014303,
            0.23796839,
            0.25856809,
            0.63636296,
            0.69455189,
            0.316758,
            0.18673652,
        ]
        print(self.draws_scrambled_1d.flatten())
        assert_equal(self.draws_scrambled_1d.shape[0], 10)
        assert_equal(self.draws_scrambled_1d.shape[1], 1)
        assert_array_almost_equal(
            self.draws_scrambled_1d.flatten(), np.array(expected)
        )

    def test_Scrambled3DSobol(self):
        self.setUp()
        expected_dim3 = [
            0.19711632,
            0.43653634,
            0.79965184,
            0.08670237,
            0.70811484,
            0.90994149,
            0.29499525,
            0.83833538,
            0.46057166,
            0.15769824,
        ]
        assert_equal(self.draws_scrambled_3d.shape[0], 10)
        assert_equal(self.draws_scrambled_3d.shape[1], 3)
        assert_array_almost_equal(
            self.draws_scrambled_3d[:, 2], np.array(expected_dim3), decimal=5
        )

    def test_Scrambled3DAsyncSobol(self):
        self.setUp()
        engine_unscrambled_3d = Sobol(3)
        draws = np.vstack([engine_unscrambled_3d.random() for i in range(10)])
        assert_array_equal(self.draws_unscrambled_3d, draws)

    def test_ScrambledSobolBounds(self):
        engine = Sobol(100, scramble=True)
        draws = engine.random(1000)
        assert_(np.all(draws >= 0))
        assert_(np.all(draws <= 1))

    def test_ScrambledFastForwardAndResetSobol(self):
        self.setUp()
        engine_scrambled_3d = Sobol(3, scramble=True, seed=12345).fast_forward(5)
        draws = engine_scrambled_3d.random(5)
        assert_array_equal(self.draws_scrambled_3d[5:10], draws)

        engine_scrambled_3d.reset()
        even_draws = []
        for i in range(10):
            if i % 2 == 0:
                even_draws.append(engine_scrambled_3d.random())
            else:
                engine_scrambled_3d.fast_forward(1)
        assert_array_equal(
            self.draws_scrambled_3d[[i for i in range(10) if i % 2 == 0]],
            np.vstack(even_draws),
        )

    def test_ScrambledDistributionSobol(self):
        engine = Sobol(10, scramble=True, seed=12345)
        draws = engine.random(1000)
        assert_array_almost_equal(
            np.mean(draws, axis=0), np.repeat(0.5, 10), decimal=2
        )
        assert_array_almost_equal(
            np.percentile(draws, 25, axis=0), np.repeat(0.25, 10), decimal=2
        )
        assert_array_almost_equal(
            np.percentile(draws, 75, axis=0), np.repeat(0.75, 10), decimal=2
        )

    def test_0Dim(self):
        engine = Sobol(0)
        draws = engine.random(5)
        assert_array_equal(np.empty((5, 0)), draws)


class TestMultinomialQMC:
    def test_MultinomialNegativePs(self):
        p = np.array([0.12, 0.26, -0.05, 0.35, 0.22])
        assert_raises(ValueError, multinomial_qmc, 10, p)

    def test_MultinomialSumOfPTooLarge(self):
        p = np.array([0.12, 0.26, 0.1, 0.35, 0.22])
        assert_raises(ValueError, multinomial_qmc, 10, p)

    def test_MultinomialBasicDraw(self):
        p = np.array([0.12, 0.26, 0.05, 0.35, 0.22])
        expected = np.array([12, 25, 6, 34, 23])
        assert_array_equal(multinomial_qmc(100, p, seed=12345), expected)

    def test_MultinomialDistribution(self):
        p = np.array([0.12, 0.26, 0.05, 0.35, 0.22])
        draws = multinomial_qmc(10000, p, seed=12345)
        assert_array_almost_equal(draws / np.sum(draws), p, decimal=4)

    def test_FindIndex(self):
        p_cumulative = np.array([0.1, 0.4, 0.45, 0.6, 0.75, 0.9, 0.99, 1.0])
        size = len(p_cumulative)
        assert_equal(_test_find_index(p_cumulative, size, 0.0), 0)
        assert_equal(_test_find_index(p_cumulative, size, 0.4), 2)
        assert_equal(_test_find_index(p_cumulative, size, 0.44999), 2)
        assert_equal(_test_find_index(p_cumulative, size, 0.45001), 3)
        assert_equal(_test_find_index(p_cumulative, size, 1.0), size - 1)


class TestNormalQMC:
    def test_NormalQMC(self):
        # d = 1
        engine = NormalQMC(k_vars=1)
        samples = engine.random()
        assert_equal(samples.shape, (1, 1))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 1))
        # d = 2
        engine = NormalQMC(k_vars=2)
        samples = engine.random()
        assert_equal(samples.shape, (1, 2))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 2))

    def test_NormalQMCInvTransform(self):
        # d = 1
        engine = NormalQMC(k_vars=1, inv_transform=True)
        samples = engine.random()
        assert_equal(samples.shape, (1, 1))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 1))
        # d = 2
        engine = NormalQMC(k_vars=2, inv_transform=True)
        samples = engine.random()
        assert_equal(samples.shape, (1, 2))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 2))

    def test_NormalQMCSeeded(self):
        # test even dimension
        engine = NormalQMC(k_vars=2, seed=12345)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [[-0.63099602, -1.32950772], [0.29625805, 1.86425618]]
        )
        assert_array_almost_equal(samples, samples_expected)
        # test odd dimension
        engine = NormalQMC(k_vars=3, seed=12345)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [
                [1.83169884, -1.40473647, 0.24334828],
                [0.36596099, 1.2987395, -1.47556275],
            ]
        )
        assert_array_almost_equal(samples, samples_expected)

    def test_NormalQMCSeededInvTransform(self):
        # test even dimension
        engine = NormalQMC(k_vars=2, seed=12345, inv_transform=True)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [[-0.41622922, 0.46622792], [-0.96063897, -0.75568963]]
        )
        assert_array_almost_equal(samples, samples_expected)
        # test odd dimension
        engine = NormalQMC(k_vars=3, seed=12345, inv_transform=True)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [
                [-1.40525266, 1.37652443, -0.8519666],
                [-0.166497, -2.3153681, -0.15975676],
            ]
        )
        assert_array_almost_equal(samples, samples_expected)

    def test_NormalQMCShapiro(self):
        engine = NormalQMC(k_vars=2, seed=12345)
        samples = engine.random(n=250)
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
        engine = NormalQMC(k_vars=2, seed=12345, inv_transform=True)
        samples = engine.random(n=250)
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
        assert_raises(
            ValueError, MultivariateNormalQMC, [0, 0], [[1, 2], [2, 1]]
        )

    def test_MultivariateNormalQMCNonPD(self):
        # try with non-pd but psd cov; should work
        engine = MultivariateNormalQMC(
            [0, 0, 0], [[1, 0, 1], [0, 1, 1], [1, 1, 2]]
        )
        assert_(engine._corr_matrix is not None)

    def test_MultivariateNormalQMCSymmetric(self):
        # try with non-symmetric cov and expect an error
        assert_raises(
            ValueError, MultivariateNormalQMC, [0, 0], [[1, 0], [2, 1]]
        )

    def test_MultivariateNormalQMC(self):
        # d = 1 scalar
        engine = MultivariateNormalQMC(mean=0, cov=5)
        samples = engine.random()
        assert_equal(samples.shape, (1, 1))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 1))

        # d = 2 list
        engine = MultivariateNormalQMC(mean=[0, 1], cov=[[1, 0], [0, 1]])
        samples = engine.random()
        assert_equal(samples.shape, (1, 2))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 2))

        # d = 3 np.array
        mean = np.array([0, 1, 2])
        cov = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        engine = MultivariateNormalQMC(mean, cov)
        samples = engine.random()
        assert_equal(samples.shape, (1, 3))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 3))

    def test_MultivariateNormalQMCInvTransform(self):
        # d = 1 scalar
        engine = MultivariateNormalQMC(mean=0, cov=5, inv_transform=True)
        samples = engine.random()
        assert_equal(samples.shape, (1, 1))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 1))

        # d = 2 list
        engine = MultivariateNormalQMC(
            mean=[0, 1], cov=[[1, 0], [0, 1]], inv_transform=True
        )
        samples = engine.random()
        assert_equal(samples.shape, (1, 2))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 2))

        # d = 3 np.array
        mean = np.array([0, 1, 2])
        cov = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        engine = MultivariateNormalQMC(mean, cov, inv_transform=True)
        samples = engine.random()
        assert_equal(samples.shape, (1, 3))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 3))

    def test_MultivariateNormalQMCSeeded(self):
        # test even dimension
        np.random.seed(54321)
        a = np.random.randn(2, 2)
        A = a @ a.transpose() + np.diag(np.random.rand(2))
        engine = MultivariateNormalQMC(np.array([0, 0]), A, seed=12345)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [[-0.67595995, -2.27437872], [0.317369, 2.66203577]]
        )
        assert_array_almost_equal(samples, samples_expected)

        # test odd dimension
        np.random.seed(54321)
        a = np.random.randn(3, 3)
        A = a @ a.transpose() + np.diag(np.random.rand(3))
        engine = MultivariateNormalQMC(np.array([0, 0, 0]), A, seed=12345)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [
                [2.05178452, -6.35744194, 0.67944512],
                [0.40993262, 2.60517697, -1.69415825],
            ]
        )
        assert_array_almost_equal(samples, samples_expected)

    def test_MultivariateNormalQMCSeededInvTransform(self):
        # test even dimension
        np.random.seed(54321)
        a = np.random.randn(2, 2)
        A = a @ a.transpose() + np.diag(np.random.rand(2))
        engine = MultivariateNormalQMC(
            np.array([0, 0]), A, seed=12345, inv_transform=True
        )
        samples = engine.random(n=2)
        samples_expected = np.array(
            [[-0.44588916, 0.22657776], [-1.02909281, -1.83193033]]
        )
        assert_array_almost_equal(samples, samples_expected)

        # test odd dimension
        np.random.seed(54321)
        a = np.random.randn(3, 3)
        A = a @ a.transpose() + np.diag(np.random.rand(3))
        engine = MultivariateNormalQMC(
            np.array([0, 0, 0]), A, seed=12345, inv_transform=True
        )
        samples = engine.random(n=2)
        samples_expected = np.array(
            [
                [-1.5740992, 5.61057598, -1.28218525],
                [-0.18650226, -5.41662685, 0.023199],
            ]
        )
        assert_array_almost_equal(samples, samples_expected)

    def test_MultivariateNormalQMCShapiro(self):
        # test the standard case
        engine = MultivariateNormalQMC(
            mean=[0, 0], cov=[[1, 0], [0, 1]], seed=12345
        )
        samples = engine.random(n=250)
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
        engine = MultivariateNormalQMC(
            mean=[1.0, 2.0], cov=[[1.5, 0.5], [0.5, 1.5]], seed=12345
        )
        samples = engine.random(n=250)
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
        engine = MultivariateNormalQMC(
            mean=[0, 0], cov=[[1, 0], [0, 1]], seed=12345, inv_transform=True
        )
        samples = engine.random(n=250)
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
        engine = MultivariateNormalQMC(
            mean=[1.0, 2.0],
            cov=[[1.5, 0.5], [0.5, 1.5]],
            seed=12345,
            inv_transform=True,
        )
        samples = engine.random(n=250)
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
        engine = MultivariateNormalQMC(
            mean=[0.0, 0.0, 0.0],
            cov=[[1.0, 0.0, 1.0], [0.0, 1.0, 1.0], [1.0, 1.0, 2.0]],
            seed=12345,
        )
        samples = engine.random(n=2000)
        assert_(all(np.abs(samples.mean(axis=0)) < 1e-2))
        assert_(np.abs(np.std(samples[:, 0]) - 1) < 1e-2)
        assert_(np.abs(np.std(samples[:, 1]) - 1) < 1e-2)
        assert_(np.abs(np.std(samples[:, 2]) - np.sqrt(2)) < 1e-2)
        for i in (0, 1, 2):
            _, pval = shapiro(samples[:, i])
            assert_(pval > 0.9)
        cov = np.cov(samples.transpose())
        assert_(np.abs(cov[0, 1]) < 1e-2)
        assert_(np.abs(cov[0, 2] - 1) < 1e-2)
        # check to see if X + Y = Z almost exactly
        assert_(
            all(np.abs(samples[:, 0] + samples[:, 1] - samples[:, 2]) < 1e-5)
        )
