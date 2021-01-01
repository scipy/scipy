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
                              _perturb_discrepancy, _update_discrepancy, QMCEngine)


class QMCEngineTests:
    qmce = ...
    can_scramble = ...
    unscramble_1d = ...
    scramble_1d = ...
    unscramble_nd = ...
    scramble_nd = ...

    scramble = [True, False]

    def engine(self, scramble: bool, **kwargs) -> QMCEngine:
        seed = np.random.RandomState(123456)
        if self.can_scramble:
            return self.qmce(scramble=scramble, seed=seed, **kwargs)
        else:
            if scramble:
                pytest.skip()
            else:
                return self.qmce(seed=seed, **kwargs)

    def reference(self, d: int, scramble: bool) -> np.ndarray:
        if scramble:
            return self.scramble_1d if d == 1 else self.scramble_nd
        else:
            return self.unscramble_1d if d == 1 else self.unscramble_nd

    @pytest.mark.parametrize("scramble", scramble)
    def test_0dim(self, scramble):
        engine = self.engine(d=0, scramble=scramble)
        sample = engine.random(4)
        assert_array_equal(np.empty((4, 0)), sample)

    @pytest.mark.parametrize("scramble", scramble)
    def test_bounds(self, scramble):
        engine = self.engine(d=100, scramble=scramble)
        sample = engine.random(1024)
        assert_(np.all(sample >= 0))
        assert_(np.all(sample <= 1))

    @pytest.mark.parametrize("scramble", scramble)
    def test_sample(self, scramble):
        ref_sample = self.reference(d=2, scramble=scramble)
        engine = self.engine(d=2, scramble=scramble)
        sample = engine.random(n=len(ref_sample))

        assert_almost_equal(sample, ref_sample, decimal=1)
        assert engine.num_generated == len(ref_sample)

    @pytest.mark.parametrize("scramble", scramble)
    def test_continuing(self, scramble):
        ref_sample = self.reference(d=2, scramble=scramble)
        engine = self.engine(d=2, scramble=scramble)

        _ = engine.random(n=len(ref_sample) // 2)
        sample = engine.random(n=len(ref_sample) // 2)
        assert_almost_equal(sample, ref_sample[len(ref_sample) // 2:],
                            decimal=1)

    @pytest.mark.parametrize("scramble", scramble)
    def test_reset(self, scramble):
        ref_sample = self.reference(d=2, scramble=scramble)
        engine = self.engine(d=2, scramble=scramble)

        _ = engine.random(n=len(ref_sample) // 2)

        engine.reset()
        assert engine.num_generated == 0

        sample = engine.random(n=len(ref_sample))
        assert_almost_equal(sample, ref_sample, decimal=1)

    @pytest.mark.parametrize("scramble", scramble)
    def test_fast_forward(self, scramble):
        ref_sample = self.reference(d=2, scramble=scramble)
        engine = self.engine(d=2, scramble=scramble)

        engine.fast_forward(4)
        sample = engine.random(n=4)

        assert_almost_equal(sample, ref_sample[4:], decimal=1)


class TestHalton(QMCEngineTests):
    qmce = qmc.Halton
    can_scramble = True
    # theoretical values known from Van der Corput
    unscramble_nd = np.array([[0, 0], [1 / 2, 1 / 3],
                              [1 / 4, 2 / 3], [3 / 4, 1 / 9],
                              [1 / 8, 4 / 9], [5 / 8, 7 / 9],
                              [3 / 8, 2 / 9], [7 / 8, 5 / 9]])
    # theoretical values unknown: convergence properties checked
    scramble_nd = np.array([[0.34229571, 0.89178423],
                            [0.84229571, 0.07696942],
                            [0.21729571, 0.41030275],
                            [0.71729571, 0.74363609],
                            [0.46729571, 0.18808053],
                            [0.96729571, 0.52141386],
                            [0.06104571, 0.8547472],
                            [0.56104571, 0.29919164]])


class TestLHS(QMCEngineTests):
    qmce = qmc.LatinHypercube
    can_scramble = False
    unscramble_nd = np.array([[0.73412877, 0.50416027],
                              [0.5924405, 0.51284543],
                              [0.57790629, 0.70797228],
                              [0.44357794, 0.64496811],
                              [0.23461223, 0.55712172],
                              [0.45337347, 0.4440004],
                              [0.73381992, 0.01751516],
                              [0.52245145, 0.33099331]])

    def test_continuing(self, *args):
        pytest.skip("Not applicable: not a sequence.")

    def test_fast_forward(self, *args):
        pytest.skip("Not applicable: not a sequence.")

    def test_sample_centered(self):
        engine = self.engine(d=2, scramble=False, centered=True)
        sample = engine.random(n=5)
        out = np.array([[0.3, 0.5],
                        [0.5, 0.3],
                        [0.1, 0.7],
                        [0.7, 0.7],
                        [0.7, 0.1]])
        assert_almost_equal(sample, out, decimal=1)


class TestOLHS(QMCEngineTests):
    qmce = qmc.OrthogonalLatinHypercube
    can_scramble = False
    unscramble_nd = np.array([[0.01587123, 0.01618008],
                              [0.24583973, 0.35254855],
                              [0.66702772, 0.82434795],
                              [0.80642206, 0.89219419],
                              [0.2825595, 0.41900669],
                              [0.98003189, 0.52861091],
                              [0.54709371, 0.23248484],
                              [0.48715457, 0.72209797]])

    def test_continuing(self, *args):
        pytest.skip("Not applicable: not a sequence.")

    def test_fast_forward(self, *args):
        pytest.skip("Not applicable: not a sequence.")

    def test_iid(self):
        # Checking independency of the random numbers generated
        engine = self.engine(d=2, scramble=False)
        n_samples = 500
        sample = engine.random(n=n_samples)
        min_b = 50  # number of bins
        bins = np.linspace(0, 1, min(min_b, n_samples) + 1)
        hist = np.histogram(sample[:, 0], bins=bins)
        out = np.array([n_samples / min_b] * min_b)
        assert_equal(hist[0], out)

        hist = np.histogram(sample[:, 1], bins=bins)
        assert_equal(hist[0], out)


class TestOptimalDesign(QMCEngineTests):
    qmce = qmc.OptimalDesign
    can_scramble = False
    unscramble_nd = np.array([[0.24583973, 0.01618008],
                              [0.01587123, 0.82434795],
                              [0.66702772, 0.35254855],
                              [0.80642206, 0.89219419],
                              [0.2825595, 0.41900669],
                              [0.98003189, 0.52861091],
                              [0.54709371, 0.23248484],
                              [0.48715457, 0.72209797]])

    def test_continuing(self, *args):
        pytest.skip("Not applicable: not a sequence.")

    def test_fast_forward(self, *args):
        pytest.skip("Not applicable: not a sequence.")

    def test_discrepancy_hierarchy(self):
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


class TestUtils:
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


class TestVDC:
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


class TestMultinomialQMC:
    def test_MultinomialNegativePs(self):
        p = np.array([0.12, 0.26, -0.05, 0.35, 0.22])
        with pytest.raises(ValueError, match=r"Elements of pvals must "
                                             r"be non-negative."):
            qmc.MultinomialQMC(p)

    def test_MultinomialSumOfPTooLarge(self):
        p = np.array([0.12, 0.26, 0.1, 0.35, 0.22])
        with pytest.raises(ValueError, match=r"Elements of pvals must sum "
                                             r"to 1."):
            qmc.MultinomialQMC(p)

    @pytest.mark.filterwarnings('ignore::UserWarning')
    def test_MultinomialBasicDraw(self):
        seed = np.random.RandomState(12345)
        p = np.array([0.12, 0.26, 0.05, 0.35, 0.22])
        expected = np.array([12, 25, 6, 35, 22])
        engine = qmc.MultinomialQMC(p, seed=seed)
        assert_array_equal(engine.random(100), expected)

    def test_MultinomialDistribution(self):
        seed = np.random.RandomState(12345)
        p = np.array([0.12, 0.26, 0.05, 0.35, 0.22])
        engine = qmc.MultinomialQMC(p, seed=seed)
        draws = engine.random(8192)
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
        base_engine = qmc.Sobol(1, scramble=True, seed=seed)
        engine = qmc.MultinomialQMC(p, engine=base_engine, seed=seed)
        assert_array_equal(engine.random(100), expected)

    def test_reset(self):
        p = np.array([0.12, 0.26, 0.05, 0.35, 0.22])
        engine = qmc.MultinomialQMC(p)
        samples = engine.random(2)
        engine.reset()
        samples_reset = engine.random(2)
        assert_array_equal(samples, samples_reset)


class TestNormalQMC:
    def test_NormalQMC(self):
        # d = 1
        seed = np.random.RandomState(123456)
        engine = qmc.MultivariateNormalQMC(mean=np.zeros(1), seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 1))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 1))
        # d = 2
        engine = qmc.MultivariateNormalQMC(mean=np.zeros(2), seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 2))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 2))

    def test_NormalQMCInvTransform(self):
        # d = 1
        seed = np.random.RandomState(123456)
        engine = qmc.MultivariateNormalQMC(
            mean=np.zeros(1), inv_transform=True, seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 1))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 1))
        # d = 2
        engine = qmc.MultivariateNormalQMC(
            mean=np.zeros(2), inv_transform=True, seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 2))
        samples = engine.random(n=5)
        assert_equal(samples.shape, (5, 2))

    def test_other_engine(self):
        seed = np.random.RandomState(123456)
        base_engine = qmc.Sobol(d=2, scramble=False, seed=seed)
        engine = qmc.MultivariateNormalQMC(mean=np.zeros(2), engine=base_engine,
                               inv_transform=True, seed=seed)
        samples = engine.random()
        assert_equal(samples.shape, (1, 2))

    def test_NormalQMCSeeded(self):
        # test even dimension
        seed = np.random.RandomState(12345)
        engine = qmc.MultivariateNormalQMC(
            mean=np.zeros(2), inv_transform=False, seed=seed)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [[-0.943472, 0.405116], [-0.63099602, -1.32950772]]
        )
        assert_array_almost_equal(samples, samples_expected)

        # test odd dimension
        seed = np.random.RandomState(12345)
        engine = qmc.MultivariateNormalQMC(
            mean=np.zeros(3), inv_transform=False, seed=seed)
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
        engine = qmc.MultivariateNormalQMC(
            mean=np.zeros(2), seed=seed, inv_transform=True)
        samples = engine.random(n=2)
        samples_expected = np.array(
            [[0.228309, -0.162516], [-0.41622922, 0.46622792]]
        )
        assert_array_almost_equal(samples, samples_expected)

        # test odd dimension
        seed = np.random.RandomState(12345)
        engine = qmc.MultivariateNormalQMC(
            mean=np.zeros(3), seed=seed, inv_transform=True)
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
        engine = qmc.MultivariateNormalQMC(mean=np.zeros(2), seed=seed)
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
        engine = qmc.MultivariateNormalQMC(
            mean=np.zeros(2), seed=seed, inv_transform=True)
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

    def test_reset(self):
        engine = qmc.MultivariateNormalQMC(mean=np.zeros(1))
        samples = engine.random(2)
        engine.reset()
        samples_reset = engine.random(2)
        assert_array_equal(samples, samples_reset)


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
