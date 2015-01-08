from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.optimize import (rosen, differential_evolution, basinhopping,
                            minimize, rosen_der, rosen_hess, rosen_hess_prod)
from scipy.optimize.optimize import _status_message
from numpy.testing import (TestCase, run_module_suite, assert_)
import warnings

halt_message = _status_message['halted']
methods = ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG', 'L-BFGS-B',
           'TNC', 'COBYLA', 'SLSQP', 'dogleg', 'trust-ncg']

class Callback(object):
    def __init__(self):
        self.iter = 0

    def callback(self, x0, *args, **kwds):
        self.iter += 1
        if self.iter == 2:
            return True
        else:
            return False


class TestCallbackHaltsMinimizer(TestCase):
    #check that returning True from the callback stops all the minimizers

    def setUp(self):
        self.x0 = [5, 5]
        self.func = rosen
        self.jac = rosen_der
        self.hess = rosen_hess
        self.bounds = [(0., 10.), (0., 10.)]
        self.Tcallback = Callback()
        self.callback = self.Tcallback.callback

    def test_diffev(self):
        res = differential_evolution(self.func, self.bounds,
                                     callback=self.callback)
        assert_(res.message == halt_message)
        assert_(self.Tcallback.iter == 2)

    def test_basin(self):
        res = basinhopping(self.func, self.x0,
                                     callback=self.callback)
        assert_(res.message == halt_message)
        assert_(self.Tcallback.iter == 2)

    def test_minimizer(self):
        for method in methods:
            #TODO TNC
            if method in ['TNC', 'COBYLA']:
                continue

            self.Tcallback.iter = 0
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = minimize(self.func, self.x0, method=method,
                               callback=self.callback, jac=self.jac,
                               hess=rosen_hess, hessp=rosen_hess_prod)
            assert_(res.message == halt_message)
            assert_(self.Tcallback.iter == 2)


if __name__ == '__main__':
    run_module_suite()