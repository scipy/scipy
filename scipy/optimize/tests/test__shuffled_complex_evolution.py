#!/usr/bin/env python
"""
Unit tests for the shuffled complex evolution global minimization algorithm.
"""
import unittest
import os
import numpy as np
from scipy.optimize._constraints import Bounds
from scipy.optimize import rosen

from pytest import raises as assert_raises, warns
from numpy.testing import assert_equal, assert_almost_equal

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
        self.lb = [0., 0.]
        self.ub = np.array([2., 2.])
        self.rosenx = [1., 1.]

    def negative_rosen(self, x):
        return -rosen(x)

    def test_defaults(self):
        # test that defaults are set correctly
        solver = ShuffledComplexEvolutionSolver(
            rosen, self.x0, self.lb, self.ub)
        assert_equal(solver.sampling, 'half-open')
        assert_equal(solver.maxn, 1000)
        assert_equal(solver.kstop, 10)
        assert_equal(solver.pcento, 0.0001)
        assert_equal(solver.ngs, 2)
        assert_equal(solver.npg, 2 * len(self.lb) + 1)
        assert_equal(solver.nps, len(self.lb) + 1)
        assert_equal(solver.nspl, 2 * len(self.lb) + 1)
        assert_equal(solver.mings, 2)
        assert_equal(solver.peps, 0.001)
        assert_equal(solver.alpha, 0.8)
        assert_equal(solver.beta, 0.45)
        assert_equal(solver.maxit, False)
        assert_equal(solver.printit, 2)
        assert_equal(solver.polish, True)
        assert_equal(solver.restartfile1, '')
        assert_equal(solver.restartfile2, '')

        solver = ShuffledComplexEvolutionSolver(
            rosen, self.x0, self.lb, self.ub,
            sampling='left-half-open',
            maxn=100, kstop=1, pcento=0.001,
            ngs=20, npg=10, nps=20, nspl=10, mings=20,
            peps=0.01,
            alpha=0.9, beta=0.55, maxit=True, printit=1,
            polish=False,
            restart=False, restartfile1='sce.restart.npz',
            restartfile2='sce.restart.txt')
        assert_equal(solver.sampling, 'left-half-open')
        assert_equal(solver.maxn, 100)
        assert_equal(solver.kstop, 1)
        assert_equal(solver.pcento, 0.001)
        assert_equal(solver.ngs, 20)
        assert_equal(solver.npg, 10)
        assert_equal(solver.nps, 20)
        assert_equal(solver.nspl, 10)
        assert_equal(solver.mings, 20)
        assert_equal(solver.peps, 0.01)
        assert_equal(solver.alpha, 0.9)
        assert_equal(solver.beta, 0.55)
        assert_equal(solver.maxit, True)
        assert_equal(solver.printit, 1)
        assert_equal(solver.polish, False)
        assert_equal(solver.restartfile1, 'sce.restart.npz')
        assert_equal(solver.restartfile2, 'sce.restart.txt')

        solver = ShuffledComplexEvolutionSolver(
            rosen, self.x0, self.lb, self.ub,
            sampling='left-half-open',
            maxn=100, kstop=1, pcento=0.001,
            ngs=20, npg=10, nps=20, nspl=10, mings=20,
            peps=0.01,
            alpha=0.9, beta=0.55, maxit=True, printit=1,
            polish=False,
            restart=True, restartfile1='sce.restart.npz',
            restartfile2='sce.restart.txt')
        assert_equal(solver.sampling, 'left-half-open')
        assert_equal(solver.maxn, 100)
        assert_equal(solver.kstop, 1)
        assert_equal(solver.pcento, 0.001)
        assert_equal(solver.ngs, 20)
        assert_equal(solver.npg, 10)
        assert_equal(solver.nps, 20)
        assert_equal(solver.nspl, 10)
        assert_equal(solver.mings, 20)
        assert_equal(solver.peps, 0.01)
        assert_equal(solver.alpha, 0.9)
        assert_equal(solver.beta, 0.55)
        assert_equal(solver.maxit, True)
        assert_equal(solver.printit, 1)
        assert_equal(solver.polish, False)
        assert_equal(solver.restartfile1, 'sce.restart.npz')
        assert_equal(solver.restartfile2, 'sce.restart.txt')

        toremove = [solver.restartfile1, solver.restartfile2]
        for ff in toremove:
            if os.path.exists(ff):
                os.remove(ff)

    def test_ShuffledComplexEvolutionSolver(self):
        solver = ShuffledComplexEvolutionSolver(
            rosen, self.x0, self.lb, self.ub, maxn=100,
            polish=False)
        result = solver.solve()
        assert_equal(result.fun, rosen(result.x))

        solver = ShuffledComplexEvolutionSolver(
            rosen, self.x0, self.lb, self.ub, maxn=100,
            polish=True)
        result = solver.solve()
        assert_equal(result.fun, rosen(result.x))

    def test_shuffled_complex_evolution(self):
        # standard
        result = shuffled_complex_evolution(rosen, self.x0, self.lb, self.ub)
        assert_equal(result.fun, rosen(result.x))

        # bounds
        result = shuffled_complex_evolution(rosen, self.x0, self.lb, self.ub)
        assert_almost_equal(result.x, self.rosenx, decimal=3)
        result = shuffled_complex_evolution(rosen, self.x0,
                                            list(zip(self.lb, self.ub)))
        assert_almost_equal(result.x, self.rosenx, decimal=3)
        result = shuffled_complex_evolution(rosen, self.x0,
                                            Bounds(self.lb, self.ub))
        assert_almost_equal(result.x, self.rosenx, decimal=3)
        result = shuffled_complex_evolution(rosen, self.x0, self.lb[0],
                                            self.ub[0])
        assert_almost_equal(result.x, self.rosenx, decimal=3)
        result = shuffled_complex_evolution(rosen, self.x0, self.lb[0:1],
                                            self.ub[0:1])
        assert_almost_equal(result.x, self.rosenx, decimal=3)
        # degenerated bounds
        x0 = [0.999, 0.5, 0.1]
        lb = [2., 0., 0.]
        ub = [2., 2., 2.]
        rosenx = [x0[0], 1., 1.]
        result = shuffled_complex_evolution(rosen, x0, lb, ub)
        assert_almost_equal(result.x, rosenx, decimal=2)
        x0 = [0.999, 0.5, 0.1]
        lb = [3., 0., 0.]
        ub = [2., 2., 2.]
        rosenx = [x0[0], 1., 1.]
        with warns(UserWarning):
            result = shuffled_complex_evolution(rosen, x0, lb, ub)
        assert_almost_equal(result.x, rosenx, decimal=2)

        # seed
        result = shuffled_complex_evolution(rosen, self.x0, self.lb, self.ub,
                                            polish=False, seed=1)
        result2 = shuffled_complex_evolution(rosen, self.x0, self.lb,
                                             self.ub, polish=False, seed=1)
        assert_equal(result.x, result2.x)
        assert_equal(result.nfev, result2.nfev)

        # restart
        restartfile1 = 'sce.restart.npz'
        restartfile2 = restartfile1 + '.txt'
        result = shuffled_complex_evolution(rosen, self.x0, self.lb,
                                            self.ub, polish=False, seed=1)
        _ = shuffled_complex_evolution(rosen, self.x0, self.lb, self.ub,
                                       polish=False, seed=1,
                                       restart=False,
                                       restartfile1=restartfile1, maxn=10)
        result2 = shuffled_complex_evolution(rosen, self.x0, self.lb,
                                             self.ub, polish=False, seed=1,
                                             restart=True)
        assert_equal(result.x, result2.x)
        assert_equal(result.nfev, result2.nfev)

        toremove = [restartfile1, restartfile2]
        for ff in toremove:
            if os.path.exists(ff):
                os.remove(ff)

        # printit
        result = shuffled_complex_evolution(rosen, self.x0, self.lb,
                                            self.ub, printit=0)
        assert_equal(result.fun, rosen(result.x))
        result = shuffled_complex_evolution(rosen, self.x0, self.lb,
                                            self.ub, printit=1)
        assert_equal(result.fun, rosen(result.x))
        result = shuffled_complex_evolution(rosen, self.x0, self.lb,
                                            self.ub, printit=2)
        assert_equal(result.fun, rosen(result.x))
        result = shuffled_complex_evolution(rosen, self.x0, self.lb,
                                            self.ub, printit=1, maxn=10)
        assert_equal(result.fun, rosen(result.x))
        result = shuffled_complex_evolution(rosen, self.x0, self.lb,
                                            self.ub, printit=2, maxn=10)
        assert_equal(result.fun, rosen(result.x))
        result = shuffled_complex_evolution(
            rosen, [0.5] * 5, [0.] * 5, [2.] * 5, printit=1,
            pcento=10)
        print(result.message)
        assert_equal(result.fun, rosen(result.x))

        # sampling
        result = shuffled_complex_evolution(
            rosen, self.x0, self.lb, self.ub,
            sampling=['half-open'] * len(self.lb))
        assert_almost_equal(result.x, self.rosenx, decimal=3)
        smpls = ['half-open', 'right-half-open', 'left-half-open', 'open',
                 'log']
        for sm in smpls:
            result = shuffled_complex_evolution(rosen, self.x0, self.lb,
                                                self.ub, sampling=sm)
            assert_almost_equal(result.x, self.rosenx, decimal=3)
        lb = [-1, -1.]
        result = shuffled_complex_evolution(rosen, self.x0, lb, self.ub,
                                            sampling='log')
        assert_almost_equal(result.x, self.rosenx, decimal=3)
        lb = [0., 0.]
        result = shuffled_complex_evolution(rosen, self.x0, lb, self.ub,
                                            sampling='log')
        assert_almost_equal(result.x, self.rosenx, decimal=3)

        # maxit
        result = shuffled_complex_evolution(self.negative_rosen, self.x0,
                                            self.lb, self.ub, maxn=100,
                                            polish=True, maxit=True)
        assert_equal(result.fun, self.negative_rosen(result.x))

    def test_errors(self):
        # test that the bounds checking works
        func = rosen
        x0 = [1.]
        bounds = [(-3)]
        assert_raises(TypeError, shuffled_complex_evolution, func, x0, bounds)
        bounds = [(-1, 1), (-1, 1)]
        assert_raises(ValueError, shuffled_complex_evolution, func, x0,
                      bounds, sampling='unknown')
        # test correct bool string
        assert_raises(ValueError, _strtobool, 'Ja')
        # test no initial population found
        func = deb03
        x0 = [-0.5, -0.5]
        bounds = [(-1, 0), (-1, 0.)]  # should be (0, 1) to work
        assert_raises(ValueError, shuffled_complex_evolution, func, x0,
                      bounds, printit=1)

    def test_shuffled_complex_evolution_args(self):
        # args
        a = 0.5
        result = shuffled_complex_evolution(rosenbrock_args, self.x0,
                                            self.lb, self.ub,
                                            args=(a,))
        assert_equal(result.fun, rosenbrock_args(result.x, a))
        assert_almost_equal(result.x, (a, a**2), decimal=3)
        # args and kwargs
        a = 0.5
        b = 200
        result = shuffled_complex_evolution(rosenbrock_args, self.x0,
                                            self.lb, self.ub,
                                            args=(a,), kwargs={'b': b})
        assert_equal(result.fun, rosenbrock_args(result.x, a, b=b))
        assert_almost_equal(result.x, (a, a**2), decimal=3)
        # kwargs
        a = 0.5
        b = 200
        result = shuffled_complex_evolution(rosenbrock_kwargs, self.x0,
                                            self.lb, self.ub,
                                            kwargs={'a': a, 'b': b})
        assert_equal(result.fun, rosenbrock_kwargs(result.x, a=a, b=b))
        assert_almost_equal(result.x, (a, a**2), decimal=3)


if __name__ == "__main__":
    unittest.main()
