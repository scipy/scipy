"""
Unit tests for trust-region optimization routines.

To run it in its simplest form::
  nosetests test_optimize.py

"""
from __future__ import division, print_function, absolute_import

import numpy as np
import scipy.optimize
from numpy.testing import (TestCase, assert_, assert_equal, assert_allclose,
        run_module_suite)


class Accumulator:
    """
    This is for testing callbacks.
    """
    def __init__(self):
        self.count = 0
        self.accum = None
    def __call__(self, x):
        self.count += 1
        if self.accum is None:
            self.accum = np.array(x)
        else:
            self.accum += x


class TestTrustRegionSolvers(TestCase):

    def test_dogleg_accuracy(self):
        # test the accuracy and the retall option
        x0 = [-1.2, 1.0]
        x_opt = [1.0, 1.0]
        x_final, allvecs = scipy.optimize.fmin_dogleg(
                scipy.optimize.rosen,
                x0,
                scipy.optimize.rosen_der,
                scipy.optimize.rosen_hess,
                gtol=1e-8,
                retall=True,
                )
        assert_allclose(x0, allvecs[0])
        assert_allclose(x_final, allvecs[-1])
        assert_allclose(x_final, x_opt)

    def test_dogleg_callback(self):
        # test the callback mechanism and the maxiter and retall options
        accumulator = Accumulator()
        maxiter = 5
        xopt, allvecs = scipy.optimize.fmin_dogleg(
                scipy.optimize.rosen,
                [-1.2, 1.0],
                scipy.optimize.rosen_der,
                scipy.optimize.rosen_hess,
                callback=accumulator,
                retall=True,
                maxiter=maxiter,
                )
        assert_equal(accumulator.count, maxiter)
        assert_equal(len(allvecs), maxiter+1)
        assert_allclose(xopt, allvecs[-1])
        assert_allclose(sum(allvecs[1:]), accumulator.accum)

    def test_solver_concordance(self):
        # Assert that dogleg uses fewer iterations than ncg on the Rosenbrock
        # test function, although this does not necessarily mean
        # that dogleg is faster or better than ncg even for this function
        # and especially not for other test functions.
        f = scipy.optimize.rosen
        g = scipy.optimize.rosen_der
        h = scipy.optimize.rosen_hess
        x_opt = [1.0, 1.0]
        easy_guess = [2.0, 2.0]
        hard_guess = [-1.2, 1.0]
        for x0 in (easy_guess, hard_guess):
            x_dogleg, allvecs_dogleg = scipy.optimize.fmin_dogleg(
                    f, x0, fprime=g, fhess=h, gtol=1e-8, retall=True)
            x_trust_ncg, allvecs_trust_ncg = scipy.optimize.fmin_trust_ncg(
                    f, x0, fprime=g, fhess=h, gtol=1e-8, retall=True)
            x_ncg, allvecs_ncg = scipy.optimize.fmin_ncg(
                    f, x0, fprime=g, fhess=h, avextol=1e-8, retall=True)
            assert_allclose(x_opt, x_dogleg)
            assert_allclose(x_opt, x_trust_ncg)
            assert_allclose(x_opt, x_ncg)
            assert_(len(allvecs_dogleg) < len(allvecs_ncg))

    def test_trust_ncg_hessp(self):
        x_opt = [1.0, 1.0]
        easy_guess = [2.0, 2.0]
        hard_guess = [-1.2, 1.0]
        for x0 in (easy_guess, hard_guess):
            x_trust_ncg = scipy.optimize.fmin_trust_ncg(
                    scipy.optimize.rosen,
                    x0,
                    fprime=scipy.optimize.rosen_der,
                    fhessp=scipy.optimize.rosen_hess_prod,
                    gtol=1e-8)
            assert_allclose(x_opt, x_trust_ncg)

    def test_dogleg_return_options(self):
        f = scipy.optimize.rosen
        g = scipy.optimize.rosen_der
        h = scipy.optimize.rosen_hess
        x0 = [2.0, 2.0]
        x_opt = [1.0, 1.0]
        # by default the output should be the optimal x
        out = scipy.optimize.fmin_dogleg(f, x0, fprime=g, fhess=h, gtol=1e-8)
        assert_allclose(out, x_opt)
        # this is also the case when both full_output=False and retall=False
        out = scipy.optimize.fmin_dogleg(f, x0, fprime=g, fhess=h, gtol=1e-8,
                full_output=False, retall=False)
        assert_allclose(out, x_opt)
        # check full_output=False and retall=True
        out = scipy.optimize.fmin_dogleg(f, x0, fprime=g, fhess=h, gtol=1e-8,
                full_output=False, retall=True)
        assert_equal(len(out), 2)
        assert_allclose(out[0], x_opt)
        # check full_output=True and retall=False
        out = scipy.optimize.fmin_dogleg(f, x0, fprime=g, fhess=h, gtol=1e-8,
                full_output=True, retall=False)
        assert_equal(len(out), 6)
        # check full_output=True and retall=True
        out = scipy.optimize.fmin_dogleg(f, x0, fprime=g, fhess=h, gtol=1e-8,
                full_output=True, retall=True)
        assert_equal(len(out), 7)


if __name__ == '__main__':
    run_module_suite()

