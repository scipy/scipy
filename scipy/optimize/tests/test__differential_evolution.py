"""
Unit tests for the differential global minimization algorithm.
"""
import unittest
import scipy.optimize._differentialevolution as _differentialevolution
import numpy as np
from scipy.optimize import rosen
import numpy.testing as npt

SEED = 1    

class TestDifferentialEvolutionSolver(unittest.TestCase):

    def setUp(self):
        np.seterr(invalid='raise')

    def test_differential_evolution(self):
        '''
            test that the Jmin of DifferentialEvolutionSolver is the same as the function
            evaluation
        '''
        func = lambda x: np.cos(14.5 * x - 0.3) + (x + 0.2) * x
        bounds = [(-3, 3)]
        result = _differentialevolution.differential_evolution(func, bounds, xtol=1e-10)

        print result.x, result.fun
        npt.assert_almost_equal(result.fun, func(result.x))

    def test_bounds_checking(self):
        '''
            test that the bounds checking works
        '''
        f = lambda x: np.cos(14.5 * x - 0.3) + (x + 0.2) * x
        bounds = [(-3, None)]
        self.assertRaises(_differentialevolution.BoundsError,
                          _differentialevolution.differential_evolution, f, bounds)
        bounds = [(-3)]
        self.assertRaises(_differentialevolution.BoundsError,
                          _differentialevolution.differential_evolution, f, bounds)
        bounds = [(-3, 3), (3,4,5)]
        self.assertRaises(_differentialevolution.BoundsError,
                          _differentialevolution.differential_evolution, f, bounds)

    def test_select_samples(self):
        '''
            select_samples should return 5 separate random numbers.
        '''

        limits = np.arange(12.).reshape(2, 6)
        solver = _differentialevolution.DifferentialEvolutionSolver(None, limits, popsize=1)
        candidate = 0
        r1, r2, r3, r4, r5 = solver.select_samples(candidate, 1, 1, 1, 1, 1)
        assert len(np.unique(np.array([candidate, r1, r2, r3, r4, r5]))) == 6

    def test_rosen(self):
        '''
            test the Rosenbrock function
        '''
        limits = np.array([[0.,  0.,  0.,  0.,  0.],
                           [2.,  2.,  2.,  2.,  2.]])
        solver = _differentialevolution.DifferentialEvolutionSolver(rosen, limits)
        result = solver.solve()

    def test_rosen_from_diff_ev(self):
        '''
            test the Rosenbrock function from differential_evolution
        '''
        bounds = [(0, 2), (0, 2), (0, 2), (0, 2), (0, 2)]
        result = _differentialevolution.differential_evolution(rosen, bounds)

if __name__ == '__main__':
    unittest.main()
