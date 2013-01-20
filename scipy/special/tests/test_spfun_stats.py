from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_array_equal, TestCase, run_module_suite, \
    assert_array_almost_equal_nulp

from scipy.special import gammaln, multigammaln


class TestMultiGammaLn(TestCase):
    def test1(self):
        np.random.seed(1234)
        a = np.abs(np.random.randn())
        assert_array_equal(multigammaln(a, 1), gammaln(a))

    def test_ararg(self):
        d = 5
        np.random.seed(1234)
        a = np.abs(np.random.randn(3, 2)) + d

        tr = multigammaln(a, d)
        assert_array_equal(tr.shape, a.shape)
        for i in range(a.size):
            assert_array_almost_equal_nulp(tr.ravel()[i],
                                           multigammaln(a.ravel()[i], d))

        d = 5
        a = np.abs(np.random.randn(1, 2)) + d

        tr = multigammaln(a, d)
        assert_array_equal(tr.shape, a.shape)
        for i in range(a.size):
            assert_array_equal(tr.ravel()[i], multigammaln(a.ravel()[i], d))

    def test_bararg(self):
        try:
            multigammaln(0.5, 1.2)
            raise Exception("Expected this call to fail")
        except ValueError:
            pass


if __name__ == '__main__':
    run_module_suite()
