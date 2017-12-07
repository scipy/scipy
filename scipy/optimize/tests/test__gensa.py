"""
Unit tests for the Generalized Simulated Annealing global optimizer
"""
from scipy.optimize import _gensa
from scipy.optimize._gensa import GenSARunner
import numpy as np
from numpy.testing import (assert_equal, TestCase, assert_allclose,
                           assert_almost_equal, assert_raises,
                           assert_array_less)


class TestGenSA(TestCase):

    def setUp(self):
        # Using Rastrigin function for performing tests
        self.func = lambda x: np.sum(x * x - 10 * np.cos(
            2 * np.pi * x)) + 10 * np.size(x)
        self.weirdfunc = lambda x: 1e15
        self.ld_bounds = [(-5.12, 5.12)] * 2
        self.hd_bounds = self.ld_bounds * 5
        self.nbtestvalues = 5000
        self.defautgr = (self.func, None, self.ld_bounds)

    def tearDown(self):
        pass

    def test_low_dim(self):
        ret = _gensa.gensa(self.func, None, self.ld_bounds)
        assert_allclose(ret.fun, 0., atol=1e-12)

    def test__visiting_dist(self):
        gr = GenSARunner(*(self.defautgr))
        values = np.zeros(self.nbtestvalues)
        for i in np.arange(self.nbtestvalues):
            values[i] = gr._visita()
        # Visiting distribution is a distorted version of Cauchy-Lorentz
        # distribution, and as no 1st and higher moments (no mean defined,
        # no variance defined).
        # Check that big tails values are generated
        assert_array_less(np.min(values), 1e-10)
        assert_array_less(1e+10, np.max(values))

    def test_high_dim(self):
        ret = _gensa.gensa(self.func, None, self.hd_bounds)
        assert_allclose(ret.fun, 0., atol=1e-12)

    def test__smooth_search(self):
        gr = GenSARunner(*(self.defautgr))
        gr._xbuffer = np.array([0.05, 0.05])
        gr._smooth_search()
        assert_allclose(gr._xbuffer, 0, atol=1e-8)
        assert_equal(gr._fvalue, 0)

    def test__yygas(self):
        gr = GenSARunner(*(self.defautgr))
        values = np.zeros(self.nbtestvalues)
        for i in np.arange(self.nbtestvalues):
            values[i] = gr._yygas()
        assert_almost_equal(np.mean(values), 0, 1)
        assert_almost_equal(np.median(values), 0, 1)
        assert_almost_equal(np.std(values), 1., 1)

    def test_max_reinit(self):
        assert_raises(ValueError, _gensa.gensa, *(self.weirdfunc,
            None, self.ld_bounds))

    def test_reproduce(self):
        seed = 1234
        res1 = _gensa.gensa(self.func, None, self.ld_bounds, seed=seed)
        res2 = _gensa.gensa(self.func, None, self.ld_bounds, seed=seed)
        res3 = _gensa.gensa(self.func, None, self.ld_bounds, seed=seed)
        # If we have reproducible results, x components found has to
        # be exactly the same, which is not the case with no seeding
        assert_equal(res1.x, res2.x)
        assert_equal(res1.x, res3.x)

    def test_bounds_integrity(self):
        wrong_bounds = [(-5.12, 5.12), (1, 0), (5.12, 5,12)]
        assert_raises(ValueError, _gensa.gensa, *(self.func,
            None, wrong_bounds))

