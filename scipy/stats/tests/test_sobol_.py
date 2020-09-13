from collections import Counter

import numpy as np
import pytest

from scipy.stats.qmc import Sobol

from numpy.testing import (assert_, assert_equal, assert_array_almost_equal,
                           assert_array_equal)
from pytest import raises as assert_raises
import scipy.stats.qmc as qmc


class TestSobol:
    # set maxDiff to None to show all differences when tests fail
    maxDiff = None

    def setUp(self):
        with pytest.warns(UserWarning):
            engine_unscrambled_1d = Sobol(1, scramble=False)
            self.draws_unscrambled_1d = engine_unscrambled_1d.random(10)
            engine_unscrambled_3d = Sobol(3, scramble=False)
            self.draws_unscrambled_3d = engine_unscrambled_3d.random(10)
            engine_scrambled_1d = Sobol(1, scramble=True, seed=12345)
            self.draws_scrambled_1d = engine_scrambled_1d.random(10)
            engine_scrambled_3d = Sobol(3, scramble=True, seed=12345)
            self.draws_scrambled_3d = engine_scrambled_3d.random(10)

    def test_warning(self):
        with pytest.warns(UserWarning):
            engine = Sobol(1)
            engine.random(10)

    def test_random_base2(self):
        self.setUp()
        engine = Sobol(1, scramble=False)
        sample = engine.random_base2(2)
        assert_array_equal(self.draws_unscrambled_1d[:4],
                           sample)

        # resampling still having N=2**n
        sample = engine.random_base2(2)
        assert_array_equal(self.draws_unscrambled_1d[4:8],
                           sample)

        # resampling again but leading to N!=2**n
        assert_raises(ValueError, engine.random_base2, 2)

    def test_raise(self):
        assert_raises(ValueError, Sobol, Sobol.MAXDIM + 1)

    def test_Unscrambled1DSobol(self):
        self.setUp()
        expected = [
            0, 0.5, 0.75, 0.25, 0.375, 0.875, 0.625, 0.125, 0.1875, 0.6875
        ]
        assert_equal(self.draws_unscrambled_1d.shape[0], 10)
        assert_equal(self.draws_unscrambled_1d.shape[1], 1)
        assert_array_equal(self.draws_unscrambled_1d.flatten(),
                           np.array(expected))

    def test_Unscrambled3DSobol(self):
        self.setUp()
        expected_dim3 = [0, 0.5, 0.25, 0.75, 0.625, 0.125, 0.875, 0.375,
                         0.9375, 0.4375]
        assert_equal(self.draws_unscrambled_3d.shape[0], 10)
        assert_equal(self.draws_unscrambled_3d.shape[1], 3)
        assert_array_equal(self.draws_unscrambled_3d[:, 2],
                           np.array(expected_dim3))
        assert_array_equal(
            self.draws_unscrambled_3d[:, 0],
            self.draws_unscrambled_1d.flatten()
        )

    def test_Unscrambled3DAsyncSobol(self):
        self.setUp()
        engine_unscrambled_3d = Sobol(3, scramble=False)
        draws = np.vstack([engine_unscrambled_3d.random() for i in range(10)])
        assert_array_equal(self.draws_unscrambled_3d, draws)

    def test_UnscrambledFastForwardAndResetSobol(self):
        self.setUp()
        engine_unscrambled_3d = Sobol(3, scramble=False).fast_forward(5)
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
        engine = Sobol(1111, scramble=False)
        count1 = Counter(engine.random().flatten().tolist())
        count2 = Counter(engine.random().flatten().tolist())
        count3 = Counter(engine.random().flatten().tolist())
        assert_equal(count1, Counter({0.0: 1111}))
        assert_equal(count2, Counter({0.5: 1111}))
        assert_equal(count3, Counter({0.25: 557, 0.75: 554}))

    def test_UnscrambledSobolBounds(self):
        engine = Sobol(1111, scramble=False)
        draws = engine.random(1024)
        assert_(np.all(draws >= 0))
        assert_(np.all(draws <= 1))

    def test_UnscrambledDistributionSobol(self):
        engine = Sobol(1111, scramble=False)
        draws = engine.random(1024)
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
            0.590297,
            0.46784395,
            0.03562005,
            0.91319746,
            0.86014303,
            0.23796839,
            0.25856809,
            0.63636296,
            0.69455189,
            0.316758
        ]
        print(self.draws_scrambled_1d.flatten())
        assert_equal(self.draws_scrambled_1d.shape[0], 10)
        assert_equal(self.draws_scrambled_1d.shape[1], 1)
        assert_array_almost_equal(
            self.draws_scrambled_1d.flatten(), np.array(expected)
        )

    def test_Scrambled3DSobol(self):
        self.setUp()
        expected_dim3 = [0.56645, 0.19712, 0.79965, 0.43654, 0.0867, 0.70811,
                         0.295, 0.90994, 0.31903, 0.94863]
        assert_equal(self.draws_scrambled_3d.shape[0], 10)
        assert_equal(self.draws_scrambled_3d.shape[1], 3)
        assert_array_almost_equal(
            self.draws_scrambled_3d[:, 2], np.array(expected_dim3), decimal=5
        )

    def test_Scrambled3DAsyncSobol(self):
        self.setUp()
        engine_unscrambled_3d = Sobol(3, scramble=False)
        draws = np.vstack([engine_unscrambled_3d.random() for i in range(10)])
        assert_array_equal(self.draws_unscrambled_3d, draws)

    def test_ScrambledSobolBounds(self):
        engine = Sobol(100, scramble=True)
        draws = engine.random(1024)
        assert_(np.all(draws >= 0))
        assert_(np.all(draws <= 1))

    def test_ScrambledFastForwardAndResetSobol(self):
        self.setUp()
        engine_scrambled_3d = Sobol(3, scramble=True,
                                    seed=12345).fast_forward(5)
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
        draws = engine.random(512)
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
        engine = Sobol(0, scramble=False)
        draws = engine.random(4)
        assert_array_equal(np.empty((4, 0)), draws)

    def test_discrepancy(self):
        engine_sobol = Sobol(10, scramble=False)
        sample_sobol = engine_sobol.random(128)

        engine_olhs = qmc.OrthogonalLatinHypercube(10)
        sample_olhs = engine_olhs.random(128)

        assert qmc.discrepancy(sample_sobol) < qmc.discrepancy(sample_olhs)

        engine_halton = qmc.Halton(10)
        sample_halton = engine_halton.random(128)

        assert qmc.discrepancy(sample_sobol) < qmc.discrepancy(sample_halton)
