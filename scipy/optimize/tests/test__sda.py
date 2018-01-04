# Simulated dual annealing unit tests implementation.
# Copyright (c) 2018 Sylvain Gubian <sylvain.gubian@pmi.com>,
# Yang Xiang <yang.xiang@pmi.com>
# Author: Sylvain Gubian, PMP S.A.
"""
Unit tests for the simulted dual annealing global optimizer
"""
from scipy.optimize import sda
from scipy.optimize._sda import VisitingDistribution
from scipy.optimize._sda import ObjectiveFunWrapper
from scipy.optimize._sda import EnergyState
import numpy as np
from numpy.testing import (assert_equal, TestCase, assert_allclose,
                           assert_almost_equal, assert_raises,
                           assert_array_less)
from scipy._lib._util import check_random_state


class TestSDA(TestCase):

    def setUp(self):
        # Using Rastrigin function for performing tests
        self.func = lambda x: np.sum(x * x - 10 * np.cos(
            2 * np.pi * x)) + 10 * np.size(x)
        # A function that returns always a big value for initialization tests
        self.weirdfunc = lambda x: 2e16
        self.ld_bounds = [(-5.12, 5.12)] * 2
        self.hd_bounds = self.ld_bounds * 4
        self.nbtestvalues = 5000
        self.high_temperature = 5230
        self.low_temperature = 0.1
        self.qv = 2.62
        self.seed = 1234
        self.rs = check_random_state(self.seed)

    def tearDown(self):
        pass

    def test_low_dim(self):
        ret = sda(self.func, None, self.ld_bounds, seed=self.seed)
        assert_allclose(ret.fun, 0., atol=1e-12)

    def test__visiting_stepping(self):
        lu = list(zip(*self.ld_bounds))
        lower = np.array(lu[0])
        upper = np.array(lu[1])
        dim = lower.size
        vd = VisitingDistribution(lower, upper, self.qv, self.rs)
        values = np.zeros(dim)
        x_step_low = vd.visiting(values, 0, self.high_temperature)
        # Make sure that only the first component is changed
        assert_equal(np.not_equal(x_step_low, 0), True)
        values = np.zeros(dim)
        x_step_high = vd.visiting(values, dim, self.high_temperature)
        # Make sure that component other than at dim has changed
        assert_equal(np.not_equal(x_step_high[0], 0), True)

    def test__visiting_dist_high_temperature(self):
        lu = list(zip(*self.ld_bounds))
        lower = np.array(lu[0])
        upper = np.array(lu[1])
        vd = VisitingDistribution(lower, upper, self.qv, self.rs)
        values = np.zeros(self.nbtestvalues)
        for i in np.arange(self.nbtestvalues):
            values[i] = vd.visit_fn(self.high_temperature)
        # Visiting distribution is a distorted version of Cauchy-Lorentz
        # distribution, and as no 1st and higher moments (no mean defined,
        # no variance defined).
        # Check that big tails values are generated
        assert_array_less(np.min(values), 1e-10)
        assert_array_less(1e+10, np.max(values))

    def test__reset(self):
        owf = ObjectiveFunWrapper(self.ld_bounds, self.weirdfunc)
        lu = list(zip(*self.ld_bounds))
        lower = np.array(lu[0])
        upper = np.array(lu[1])
        es = EnergyState(lower, upper)
        assert_raises(ValueError, es.reset, *(owf, check_random_state(None)))

    def test_high_dim(self):
        ret = sda(self.func, None, self.hd_bounds)
        assert_allclose(ret.fun, 0., atol=1e-12)

    def test__gaussian(self):
        lu = list(zip(*self.ld_bounds))
        lower = np.array(lu[0])
        upper = np.array(lu[1])
        vd = VisitingDistribution(lower, upper, self.qv, self.rs)
        values = np.zeros(self.nbtestvalues)
        for i in np.arange(self.nbtestvalues):
            values[i] = vd.gaussian_fn(1)
        assert_almost_equal(np.mean(values), 0, 1)
        assert_almost_equal(np.median(values), 0, 1)
        assert_almost_equal(np.std(values), 1., 1)

    def test_max_reinit(self):
        assert_raises(ValueError, sda, *(self.weirdfunc, None,
                                           self.ld_bounds))

    def test_reproduce(self):
        seed = 1234
        res1 = sda(self.func, None, self.ld_bounds, seed=seed)
        res2 = sda(self.func, None, self.ld_bounds, seed=seed)
        res3 = sda(self.func, None, self.ld_bounds, seed=seed)
        # If we have reproducible results, x components found has to
        # be exactly the same, which is not the case with no seeding
        assert_equal(res1.x, res2.x)
        assert_equal(res1.x, res3.x)

    def test_bounds_integrity(self):
        wrong_bounds = [(-5.12, 5.12), (1, 0), (5.12, 5.12)]
        assert_raises(ValueError, sda, *(self.func, None, wrong_bounds))
