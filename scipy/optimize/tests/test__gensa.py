"""
Unit tests for the Generalized Simulated Annealing global optimizer
"""
from scipy.optimize import _gensa
from scipy.optimize._gensa import GenSARunner
from scipy.optimize._gensa import GenSARunnerException
import numpy as np
from numpy.testing import (assert_equal, TestCase, assert_allclose,
                           run_module_suite, assert_almost_equal,
                           assert_string_equal, assert_raises, assert_)


class TestGenSA(TestCase):

    def setUp(self):
        # Using Rastrigin function for performing tests
        self.func = lambda x: np.sum(x * x - 10 * np.cos(2 * np.pi * x))\
                + 10 * np.size(x)
        self.weirdfunc = lambda x: 1e15
        self.ld_lowerb = [-5.12] * 2
        self.ld_upperb = [5.12] * 2
        self.hd_lowerb = self.ld_lowerb * 5
        self.hd_upperb = self.ld_upperb * 5
        self.nbtestvalues = 5000
        self.defautgr = (self.func, None, self.ld_lowerb, self.ld_upperb)

    def tearDown(self):
        pass

    def test_low_dim(self):
        ret = _gensa.gensa(self.func, None, self.ld_lowerb, self.ld_upperb)
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
        assert(np.min(values) < 1e-10)
        assert(np.max(values) > 1e+10)

    def test_high_dim(self):
        ret = _gensa.gensa(self.func, None, self.hd_lowerb, self.hd_upperb)
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
        assert_raises(GenSARunnerException, _gensa.gensa, *(self.weirdfunc,
            None, self.ld_lowerb, self.ld_upperb))


    def test__check_stopping_cond(self):
        gr = GenSARunner(self.func, None, self.hd_lowerb, self.hd_upperb)
        gr._maxtime = 0.5
        gr._extensive = True
        gr.initialize()
        gr.start_search()
        assert(gr._message[0].startswith('Time'))

