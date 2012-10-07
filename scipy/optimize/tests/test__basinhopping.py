"""
Unit tests for the basin hopping global minimization algorithm.
"""
import copy

from numpy.testing import TestCase, run_module_suite, \
    assert_almost_equal, assert_, dec
import numpy as np
from numpy import cos, sin

from scipy.optimize import basinhopping, basinhopping_advanced, minimize
from scipy.optimize._basinhopping import _AdaptiveStepsize, _RandomDisplacement


def func1d(x):
    f = cos(14.5 * x - 0.3) + (x + 0.2) * x
    df = np.array(-14.5 * sin(14.5 * x - 0.3) + 2. * x + 0.2)
    return f, df


def func1d_nograd(x):
    f = cos(14.5 * x - 0.3) + (x + 0.2) * x
    df = np.array(-14.5 * sin(14.5 * x - 0.3) + 2. * x + 0.2)
    return f, df


def func2d_nograd(x):
    f = cos(14.5 * x[0] - 0.3) + (x[1] + 0.2) * x[1] + (x[0] + 0.2) * x[0]
    return f


def func2d(x):
    f = cos(14.5 * x[0] - 0.3) + (x[1] + 0.2) * x[1] + (x[0] + 0.2) * x[0]
    df = np.zeros(2)
    df[0] = -14.5 * sin(14.5 * x[0] - 0.3) + 2. * x[0] + 0.2
    df[1] = 2. * x[1] + 0.2
    return f, df


class Minimizer(object):
    def __init__(self, func, **kwargs):
        self.kwargs = kwargs
        self.func = func

    def __call__(self, x0, **newkwargs):
        #combine the two kwargs
        kwargs = dict(newkwargs.items() + self.kwargs.items())
        res = minimize(self.func, x0, **kwargs)
        return res


class TestBasinHopping(TestCase):
    """ Tests for basinhopping """
    def setUp(self):
        """ Tests setup.

        run tests based on the 1-D and 2-D functions described above.  These
        are the same functions as used in the anneal algorithm with some
        gradients added.
        """
        self.x0 = (1.0, [1.0, 1.0])
        self.sol = (-0.195, np.array([-0.195, -0.1]))
        self.upper = (3., [3., 3.])
        self.lower = (-3., [-3., -3.])
        self.tol = 3  # number of decimal places

        self.maxiter = 100
        self.disp = False

        # fix random seed
        np.random.seed(1234)

        self.kwargs = {"method": "L-BFGS-B", "jac": True}
        self.kwargs_nograd = {"method": "L-BFGS-B"}

    def test_1d_grad(self):
        """test 1d minimizations with gradient"""
        i = 0
        res = basinhopping(self.x0[i], func1d, minimizer_kwargs=self.kwargs,
                           maxiter=self.maxiter, disp=self.disp)
        assert_almost_equal(res.x, self.sol[i], self.tol)

    def test_2d(self):
        """test 2d minimizations with gradient"""
        i = 1
        res = basinhopping(self.x0[i], func2d, minimizer_kwargs=self.kwargs,
                           maxiter=self.maxiter, disp=self.disp)
        assert_almost_equal(res.x, self.sol[i], self.tol)

    def test_2d_nograd(self):
        """test 2d minimizations without gradient"""
        i = 1
        res = basinhopping(self.x0[i], func2d_nograd,
                           minimizer_kwargs=self.kwargs_nograd,
                           maxiter=self.maxiter, disp=self.disp)
        assert_almost_equal(res.x, self.sol[i], self.tol)

    def test_pass_minimizer(self):
        """test 2d minimizations with user defined minimizer"""
        i = 1
        minimizer = Minimizer(func2d, **self.kwargs)
        res = basinhopping(self.x0[i], minimizer=minimizer,
                           maxiter=self.maxiter, disp=self.disp)
        assert_almost_equal(res.x, self.sol[i], self.tol)

    def test_all_minimizers(self):
        """test 2d minimizations with gradient"""
        i = 1
        methods = ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG',
                   'L-BFGS-B', 'TNC', 'COBYLA', 'SLSQP']
        minimizer_kwargs = copy.copy(self.kwargs)
        for method in methods:
            minimizer_kwargs["method"] = method
            res = basinhopping(self.x0[i], func2d,
                               minimizer_kwargs=self.kwargs,
                               maxiter=self.maxiter, disp=self.disp)
            assert_almost_equal(res.x, self.sol[i], self.tol)

    #below here we are testing basinhopping_advanced

    def test_bh_advanced(self):
        """test 2d minimizations with gradient"""
        i = 1
        res = basinhopping_advanced(self.x0[i], func2d,
                                    minimizer_kwargs=self.kwargs,
                                    maxiter=self.maxiter, disp=self.disp)
        assert_almost_equal(res.x, self.sol[i], self.tol)

    def test_pass_takestep(self):
        class MyTakeStep(_RandomDisplacement):
            """use a copy of displace, but have it set a special parameter to
            make sure it's actually being used."""
            def __init__(self):
                self.been_called = False
                return super(MyTakeStep, self).__init__()

            def __call__(self, x):
                self.been_called = True
                return super(MyTakeStep, self).__call__(x)

        takestep = MyTakeStep()
        initial_step_size = takestep.stepsize
        i = 1
        res = basinhopping_advanced(self.x0[i], func2d,
                                    minimizer_kwargs=self.kwargs,
                                    maxiter=self.maxiter, disp=self.disp,
                                    take_step=takestep)
        assert_almost_equal(res.x, self.sol[i], self.tol)
        assert_(takestep.been_called)
        #make sure that the built in adaptive step size has been used
        assert_(initial_step_size != takestep.stepsize)

    def test_pass_takestep2(self):
        class MyTakeStep(_RandomDisplacement):
            """use a copy of displace, but have it set a special parameter to
            make sure it's actually being used.

            this time add a function report which overrides the default
            adaptive step size routine.
            """
            def __init__(self):
                self.been_called = False
                return super(MyTakeStep, self).__init__()

            def __call__(self, x):
                self.been_called = True
                return super(MyTakeStep, self).__call__(x)

            def report(self, accept, **kwargs):
                return

        takestep = MyTakeStep()
        initial_step_size = takestep.stepsize
        i = 1
        res = basinhopping_advanced(self.x0[i], func2d,
                                    minimizer_kwargs=self.kwargs,
                                    maxiter=self.maxiter, disp=self.disp,
                                    take_step=takestep)
        assert_almost_equal(res.x, self.sol[i], self.tol)
        assert_(takestep.been_called)
        #make sure that the built in adaptive step size has **not** been used
        assert_(initial_step_size == takestep.stepsize)

    def test_pass_accept_test(self):
        class AcceptTest(_RandomDisplacement):
            """pass a custom accept test

            This does nothing but make sure it's being used and ensure all the
            possible return values are accepted
            """
            def __init__(self):
                self.been_called = False
                self.ncalls = 0

            def __call__(self, **kwargs):
                self.been_called = True
                self.ncalls += 1
                if self.ncalls == 1:
                    return False
                elif self.ncalls == 2:
                    return 'force accept'
                else:
                    return True

        accept_test = AcceptTest()
        i = 1
        #there's no point in running it more than a few steps.
        res = basinhopping_advanced(self.x0[i], func2d,
                                    minimizer_kwargs=self.kwargs, maxiter=10,
                                    disp=self.disp, accept_test=accept_test)
        assert_(accept_test.been_called)

    def test_pass_callback(self):
        class CallBack(_RandomDisplacement):
            """pass a custom callback function

            This makes sure it's being used.  It also returns True after 10
            steps to ensure that it's stopping early.

            """
            def __init__(self):
                self.been_called = False
                self.ncalls = 0

            def __call__(self, x, f, accepted):
                self.been_called = True
                self.ncalls += 1
                if self.ncalls == 10:
                    return True

        callback = CallBack()
        i = 1
        #there's no point in running it more than a few steps.
        res = basinhopping_advanced(self.x0[i], func2d,
                                    minimizer_kwargs=self.kwargs, maxiter=30,
                                    disp=self.disp, callback=callback)
        assert_(callback.been_called)
        assert_("callback" in res.message[0])
        assert_(res.nit == 10)


if __name__ == "__main__":
    run_module_suite()
