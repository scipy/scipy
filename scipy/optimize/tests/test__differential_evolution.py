"""
Unit tests for the differential global minimization algorithm.
"""
import scipy.optimize._differentialevolution as _differentialevolution
import numpy as np
from scipy.optimize import rosen
import numpy.testing as npt

SEED = 1


class TestDifferentialEvolutionSolver(npt.TestCase):

    def setUp(self):
        np.seterr(invalid='raise')
        self.limits = np.array([[0., 0.],
                                [2., 2.]])
        self.bounds = [(0, 2), (0, 2)]

        dummy_limits = np.array([[0.], [100]])
        self.dummy_solver = _differentialevolution.DifferentialEvolutionSolver(
            self.dummy_function, dummy_limits)

    def dummy_function(self):
        pass

    def test__can_init_with_dithering(self):
        mutation = (0.5, 1)
        solver = _differentialevolution.DifferentialEvolutionSolver(
            self.dummy_function, self.limits, mutation=mutation)
        if hasattr(solver, 'dither') is False:
            raise AttributeError

    def test__scale_parameters(self):
        trial = np.array([0.3])
        npt.assert_equal(30, self.dummy_solver._scale_parameters(trial))

        # it should also work with the limits reversed
        self.dummy_solver.limits = np.array([[100], [0.]])
        npt.assert_equal(30, self.dummy_solver._scale_parameters(trial))

    def test__unscale_parameters(self):
        trial = np.array([30])
        npt.assert_equal(0.3, self.dummy_solver._unscale_parameters(trial))

        # it should also work with the limits reversed
        self.dummy_solver.limits = np.array([[100], [0.]])
        npt.assert_equal(0.3, self.dummy_solver._unscale_parameters(trial))

    def test__ensure_constraint(self):
        trial = np.array([1.1, -100, 2., 300., -0.00001])
        trial = self.dummy_solver._ensure_constraint(trial)
        npt.assert_equal(np.all(trial <= 1), True)

    def test_differential_evolution(self):
        # test that the Jmin of DifferentialEvolutionSolver
        # is the same as the function evaluation
        solver = _differentialevolution.DifferentialEvolutionSolver(
            rosen, self.limits)
        result = solver.solve()
        npt.assert_almost_equal(result.fun, rosen(result.x))

    def test_callback_terminates(self):
        # test that if the callback returns true, then the minimization halts
        bounds = [(0, 2), (0, 2)]

        def callback(param, convergence=0.):
            return True

        result = _differentialevolution.differential_evolution(
            rosen, bounds, callback=callback)

        npt.assert_string_equal(result.message,
                                'callback function requested stop early '
                                'by returning True')

    def test_bounds_checking(self):
        # test that the bounds checking works
        f = rosen
        bounds = [(-3, None)]
        self.assertRaises(_differentialevolution.BoundsError,
                          _differentialevolution.differential_evolution, f, bounds)
        bounds = [(-3)]
        self.assertRaises(_differentialevolution.BoundsError,
                          _differentialevolution.differential_evolution, f, bounds)
        bounds = [(-3, 3), (3, 4, 5)]
        self.assertRaises(_differentialevolution.BoundsError,
                          _differentialevolution.differential_evolution, f, bounds)

    def test_select_samples(self):
        # select_samples should return 5 separate random numbers.
        limits = np.arange(12.).reshape(2, 6)
        solver = _differentialevolution.DifferentialEvolutionSolver(
            None, limits, popsize=1)
        candidate = 0
        r1, r2, r3, r4, r5 = solver._select_samples(candidate, 5)
        npt.assert_equal(
            len(np.unique(np.array([candidate, r1, r2, r3, r4, r5]))), 6)

    @npt.dec.slow
    def test_rosen(self):
        # test the Rosenbrock function from object
        solver = _differentialevolution.DifferentialEvolutionSolver(
            rosen, self.limits)
        result = solver.solve()

    @npt.dec.slow
    def test_rosen_from_diff_ev(self):
        # test the Rosenbrock function from differential_evolution function

        result = _differentialevolution.differential_evolution(
            rosen, self.bounds)

    def test__make_random_gen(self):
#     If seed is None, return the RandomState singleton used by np.random.
#     If seed is an int, return a new RandomState instance seeded with seed.
#     If seed is already a RandomState instance, return it.
#     Otherwise raise ValueError.
        rsi = _differentialevolution._make_random_gen(1)
        npt.assert_equal(type(rsi), np.random.RandomState)
        rsi = _differentialevolution._make_random_gen(rsi)
        npt.assert_equal(type(rsi), np.random.RandomState)
        rsi = _differentialevolution._make_random_gen(None)
        npt.assert_equal(type(rsi), np.random.RandomState)
        self.assertRaises(
            ValueError, _differentialevolution._make_random_gen, 'a')


if __name__ == '__main__':
    unittest.main()
