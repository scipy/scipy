#!/usr/bin/env python
"""
Unit tests for the shuffled complex evolution global minimization algorithm.

"""
import unittest
import numpy as np
from scipy.optimize._constraints import Bounds
from scipy.optimize import rosen

from pytest import raises as assert_raises, warns
from numpy.testing import assert_equal, assert_allclose

from scipy.optimize._shuffled_complex_evolution import (
    ShuffledComplexEvolutionSolver, _strtobool)
from scipy.optimize import shuffled_complex_evolution


def rosenbrock_args(x, a, b=100):
    """
    Rosenbrock, global optimum: f(a, a^2) = 0.0.

    """
    f = (a - x[0])**2 + b * (x[1] - x[0]**2)**2

    return f


def rosenbrock_kwargs(x, a=1, b=100):
    """
    Rosenbrock, global optimum: f(a, a^2) = 0.0.

    """
    f = (a - x[0])**2 + b * (x[1] - x[0]**2)**2

    return f


def deb03(x):
    r"""
    Deb 3 objective function from scipy.benchmarks

    This is a multimodal minimization problem defined as follows:

    .. math::

        f_{\text{Deb02}}(x) = - \frac{1}{N} \sum_{i=1}^n \sin^6 \left[ 5 \pi
        \left ( x_i^{3/4} - 0.05 \right) \right ]


    Here, :math:`n` represents the number of dimensions and :math:`x_i \in
    [0, 1]` for :math:`i = 1, ..., n`.

    *Global optimum*: :math:`f(x) = 0.0`. The number of global minima is
    :math:`5^n` that are evenly spaced in the function landscape, where
    :math:`n` represents the dimension of the problem.

    .. [1] Jamil, M. & Yang, X.-S. A Literature Survey of Benchmark Functions
    For Global Optimization Problems Int. Journal of Mathematical Modelling
    and Numerical Optimisation, 2013, 4, 150-194.

    """
    return -(1.0 / len(x)) * np.sum(np.sin(5 * np.pi * (x**0.75 - 0.05))**6)


# Some lines are not covered:
#   120-121: second args is empty dict
#   919: if npt < 1
#   942: redo random number if exactly 0 in sampling='open'
#   957-960: if lb and ub <=0 with sampling='log'
class TestShuffledComplexEvolutionSolver(unittest.TestCase):

    def setUp(self):
        self.x0 = [0.5, 0.1]
        self.lower_bounds = [0., 0.]
        self.upper_bounds = np.array([2., 2.])
        self.rosenx = [1., 1.]

    def negative_rosen(self, x):
        return -rosen(x)

    def test_defaults(self):
        # test that defaults are set correctly
        solver = ShuffledComplexEvolutionSolver(
            rosen, self.x0, self.lower_bounds, self.upper_bounds)
        assert_equal(solver.sampling, 'half-open')
        assert_equal(solver.maxfev, 6496)
        assert_equal(solver.n_check, 10)
        assert_equal(solver.f_tol, 0.0001)
        assert_equal(solver.n_complex, 2)
        assert_equal(solver.n_point_complex, 2 * len(self.lower_bounds) + 1)
        assert_equal(solver.n_point_subcomplex, len(self.lower_bounds) + 1)
        assert_equal(solver.n_eval_complex_per_shuffle,
                     2 * len(self.lower_bounds) + 1)
        assert_equal(solver.min_n_complex, 2)
        assert_equal(solver.p_tol, 0.001)
        assert_equal(solver.alpha, 0.8)
        assert_equal(solver.beta, 0.45)
        assert_equal(solver.maximize, False)
        assert_equal(solver.disp, 0)
        assert_equal(solver.polish, True)

        solver = ShuffledComplexEvolutionSolver(
            rosen, self.x0, self.lower_bounds, self.upper_bounds,
            sampling='left-half-open',
            maxfev=100, n_check=1, f_tol=0.001,
            n_complex=20, n_point_complex=10, n_point_subcomplex=20,
            n_eval_complex_per_shuffle=10, min_n_complex=20,
            p_tol=0.01,
            alpha=0.9, beta=0.55, maximize=True, disp=1,
            polish=False)
        assert_equal(solver.sampling, 'left-half-open')
        assert_equal(solver.maxfev, 100)
        assert_equal(solver.n_check, 1)
        assert_equal(solver.f_tol, 0.001)
        assert_equal(solver.n_complex, 20)
        assert_equal(solver.n_point_complex, 10)
        assert_equal(solver.n_point_subcomplex, 20)
        assert_equal(solver.n_eval_complex_per_shuffle, 10)
        assert_equal(solver.min_n_complex, 20)
        assert_equal(solver.p_tol, 0.01)
        assert_equal(solver.alpha, 0.9)
        assert_equal(solver.beta, 0.55)
        assert_equal(solver.maximize, True)
        assert_equal(solver.disp, 1)
        assert_equal(solver.polish, False)

        solver = ShuffledComplexEvolutionSolver(
            rosen, self.x0, self.lower_bounds, self.upper_bounds,
            sampling='left-half-open',
            maxfev=100, n_check=1, f_tol=0.001,
            n_complex=20, n_point_complex=10, n_point_subcomplex=20,
            n_eval_complex_per_shuffle=10, min_n_complex=20,
            p_tol=0.01,
            alpha=0.9, beta=0.55, maximize=True, disp=1,
            polish=True)
        assert_equal(solver.sampling, 'left-half-open')
        assert_equal(solver.maxfev, 100)
        assert_equal(solver.n_check, 1)
        assert_equal(solver.f_tol, 0.001)
        assert_equal(solver.n_complex, 20)
        assert_equal(solver.n_point_complex, 10)
        assert_equal(solver.n_point_subcomplex, 20)
        assert_equal(solver.n_eval_complex_per_shuffle, 10)
        assert_equal(solver.min_n_complex, 20)
        assert_equal(solver.p_tol, 0.01)
        assert_equal(solver.alpha, 0.9)
        assert_equal(solver.beta, 0.55)
        assert_equal(solver.maximize, True)
        assert_equal(solver.disp, 1)
        assert_equal(solver.polish, True)

    def test_ShuffledComplexEvolutionSolver(self):
        solver = ShuffledComplexEvolutionSolver(
            rosen, self.x0, self.lower_bounds, self.upper_bounds, maxfev=100,
            polish=False)
        result = solver.solve()
        assert_equal(result.fun, rosen(result.x))

        solver = ShuffledComplexEvolutionSolver(
            rosen, self.x0, self.lower_bounds, self.upper_bounds, maxfev=100,
            polish=True)
        result = solver.solve()
        assert_equal(result.fun, rosen(result.x))

    def test_shuffled_complex_evolution(self):
        # standard
        result = shuffled_complex_evolution(rosen, self.x0, self.lower_bounds,
                                            self.upper_bounds)
        assert_equal(result.fun, rosen(result.x))

        # bounds
        result = shuffled_complex_evolution(rosen, self.x0, self.lower_bounds,
                                            self.upper_bounds)
        assert_allclose(result.x, self.rosenx, atol=1e-3)
        result = shuffled_complex_evolution(rosen, self.x0,
                                            list(zip(self.lower_bounds,
                                                     self.upper_bounds)),
                                            polish=False)
        assert_allclose(result.x, self.rosenx, atol=1e-3)
        result = shuffled_complex_evolution(rosen, self.x0,
                                            Bounds(self.lower_bounds,
                                                   self.upper_bounds))
        assert_allclose(result.x, self.rosenx, atol=1e-3)
        result = shuffled_complex_evolution(rosen, self.x0,
                                            self.lower_bounds[0],
                                            self.upper_bounds[0],
                                            polish=False)
        assert_allclose(result.x, self.rosenx, atol=1e-3)
        result = shuffled_complex_evolution(rosen, self.x0,
                                            self.lower_bounds[0:1],
                                            self.upper_bounds[0:1])
        assert_allclose(result.x, self.rosenx, atol=1e-3)
        # degenerated bounds
        x0 = [0.999, 0.5, 0.1]
        lower_bounds = [2., 0., 0.]
        upper_bounds = [2., 2., 2.]
        rosenx = [x0[0], 1., 1.]
        result = shuffled_complex_evolution(rosen, x0, lower_bounds,
                                            upper_bounds)
        assert_allclose(result.x, rosenx, atol=1e-2)
        x0 = [0.999, 0.5, 0.1]
        lower_bounds = [3., 0., 0.]
        upper_bounds = [2., 2., 2.]
        rosenx = [x0[0], 1., 1.]
        with warns(UserWarning):
            result = shuffled_complex_evolution(rosen, x0, lower_bounds,
                                                upper_bounds, polish=False)
        assert_allclose(result.x, rosenx, atol=1e-2)

        # seed
        result = shuffled_complex_evolution(rosen, self.x0, self.lower_bounds,
                                            self.upper_bounds,
                                            polish=False, seed=1)
        result2 = shuffled_complex_evolution(rosen, self.x0, self.lower_bounds,
                                             self.upper_bounds, polish=False,
                                             seed=1)
        assert_equal(result.x, result2.x)
        assert_equal(result.nfev, result2.nfev)

        # disp
        result = shuffled_complex_evolution(rosen, self.x0, self.lower_bounds,
                                            self.upper_bounds, disp=0,
                                            polish=False)
        assert_equal(result.fun, rosen(result.x))
        result = shuffled_complex_evolution(rosen, self.x0, self.lower_bounds,
                                            self.upper_bounds, disp=1,
                                            polish=False)
        assert_equal(result.fun, rosen(result.x))
        result = shuffled_complex_evolution(rosen, self.x0, self.lower_bounds,
                                            self.upper_bounds, disp=2,
                                            polish=False)
        assert_equal(result.fun, rosen(result.x))
        result = shuffled_complex_evolution(rosen, self.x0, self.lower_bounds,
                                            self.upper_bounds, disp=1,
                                            maxfev=10, polish=False)
        assert_equal(result.fun, rosen(result.x))
        result = shuffled_complex_evolution(rosen, self.x0, self.lower_bounds,
                                            self.upper_bounds, disp=2,
                                            maxfev=10, polish=False)
        assert_equal(result.fun, rosen(result.x))
        result = shuffled_complex_evolution(
            rosen, [0.5] * 5, [0.] * 5, [2.] * 5, disp=1,
            f_tol=10)
        # print(result.message)
        assert_equal(result.fun, rosen(result.x))

        # sampling
        result = shuffled_complex_evolution(
            rosen, self.x0, self.lower_bounds, self.upper_bounds,
            sampling=['half-open'] * len(self.lower_bounds))
        assert_allclose(result.x, self.rosenx, atol=1e-3)
        smpls = ['half-open', 'right-half-open', 'left-half-open', 'open',
                 'log']
        for sm in smpls:
            result = shuffled_complex_evolution(rosen, self.x0,
                                                self.lower_bounds,
                                                self.upper_bounds, sampling=sm)
            assert_allclose(result.x, self.rosenx, atol=1e-3)
        lower_bounds = [-1, -1.]
        result = shuffled_complex_evolution(rosen, self.x0, lower_bounds,
                                            self.upper_bounds, sampling='log')
        assert_allclose(result.x, self.rosenx, atol=1e-3)
        lower_bounds = [0., 0.]
        result = shuffled_complex_evolution(rosen, self.x0, lower_bounds,
                                            self.upper_bounds, sampling='log')
        assert_allclose(result.x, self.rosenx, atol=1e-3)

        # maximize
        result = shuffled_complex_evolution(self.negative_rosen, self.x0,
                                            self.lower_bounds,
                                            self.upper_bounds, maxfev=100,
                                            polish=True, maximize=True)
        assert_equal(result.fun, self.negative_rosen(result.x))

    def test_errors(self):
        # test that the bounds checking works
        func = rosen
        x0 = [1.]
        bounds = [(-3)]
        with assert_raises(TypeError, match='object is not iterable'):
            shuffled_complex_evolution(func, x0, bounds)
        bounds = [(-1, 1), (-1, 1)]
        with assert_raises(ValueError, match='unknown sampling option'):
            shuffled_complex_evolution(func, x0, bounds, sampling='unknown')
        # test correct bool string
        with assert_raises(ValueError, match='invalid truth value'):
            _strtobool('Ja')
        # test no initial population found
        func = deb03
        x0 = [-0.5, -0.5]
        bounds = [(-1, 0), (-1, 0.)]  # should be (0, 1) to work
        # deb03 with values < 0 gives NaN.
        # Hence SCE raises ValueError if standard numpy error settings
        # because SCE cannot produce initial population.
        # It raises RuntimeWarning if np.seterr('invalid'='raise')
        # because of the NaN.
        # with assert_raises(ValueError,
        #                    match='Did not succeed to produce initial'):
        print('Numpy errors', np.geterr())
        with assert_raises(RuntimeWarning,
                           match='invalid value encountered'):
            shuffled_complex_evolution(func, x0, bounds, disp=1)

    def test_shuffled_complex_evolution_args(self):
        # args
        a = 0.5
        result = shuffled_complex_evolution(rosenbrock_args, self.x0,
                                            self.lower_bounds,
                                            self.upper_bounds,
                                            args=(a,))
        assert_equal(result.fun, rosenbrock_args(result.x, a))
        assert_allclose(result.x, (a, a**2), atol=1e-3)
        # args and kwargs
        a = 0.5
        b = 200
        result = shuffled_complex_evolution(rosenbrock_args, self.x0,
                                            self.lower_bounds,
                                            self.upper_bounds,
                                            args=(a,), kwargs={'b': b},
                                            polish=False)
        assert_equal(result.fun, rosenbrock_args(result.x, a, b=b))
        assert_allclose(result.x, (a, a**2), atol=1e-3)
        # kwargs
        a = 0.5
        b = 200
        result = shuffled_complex_evolution(rosenbrock_kwargs, self.x0,
                                            self.lower_bounds,
                                            self.upper_bounds,
                                            kwargs={'a': a, 'b': b})
        assert_equal(result.fun, rosenbrock_kwargs(result.x, a=a, b=b))
        assert_allclose(result.x, (a, a**2), atol=1e-3)


if __name__ == "__main__":
    unittest.main()
