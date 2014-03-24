"""
Unit tests for the differential global minimization algorithm.
"""
import scipy.optimize._differentialevolution as _differentialevolution
from scipy.optimize._differentialevolution import DifferentialEvolutionSolver
from scipy.optimize import differential_evolution
import numpy as np
from scipy.optimize import rosen
import numpy.testing as npt

SEED = 1


class TestDifferentialEvolutionSolver(npt.TestCase):

    def setUp(self):
        np.seterr(invalid='raise')
        self.limits = np.array([[0., 0.],
                                [2., 2.]])
        self.bounds = [(0., 2.), (0., 2.)]

        self.dummy_solver = DifferentialEvolutionSolver(self.dummy_function,
                                                        [(0, 100)])
            
        #dummy_solver2 will be used to test mutation strategies
        self.dummy_solver2 = DifferentialEvolutionSolver(self.dummy_function,
                                                         [(0, 1)],
                                                         popsize=7,
                                                         mutation=0.5)
        #create a population that's only 7 members long
        #[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
        population = np.atleast_2d(np.arange(0.1, 0.8, 0.1)).T
        self.dummy_solver2.population = population
        
    def dummy_function(self):
        pass

    """
    test all the mutation strategies
    """
    def test__strategy_resolves(self):
        #test that the correct mutation function is resolved by
        #different requested strategy arguments
        solver = DifferentialEvolutionSolver(rosen,
                                             self.bounds,
                                             strategy='best1exp')
        npt.assert_equal(solver.strategy, 'best1exp')
        npt.assert_equal(solver.mutation_func.__name__, '_best1')

        solver = DifferentialEvolutionSolver(rosen,
                                             self.bounds,
                                             strategy='best1bin')
        npt.assert_equal(solver.strategy, 'best1bin')
        npt.assert_equal(solver.mutation_func.__name__, '_best1')

        solver = DifferentialEvolutionSolver(rosen,
                                             self.bounds,
                                             strategy='rand1bin')
        npt.assert_equal(solver.strategy, 'rand1bin')
        npt.assert_equal(solver.mutation_func.__name__, '_rand1')

        solver = DifferentialEvolutionSolver(rosen,
                                             self.bounds,
                                             strategy='rand1exp')
        npt.assert_equal(solver.strategy, 'rand1exp')
        npt.assert_equal(solver.mutation_func.__name__, '_rand1')
        
        solver = DifferentialEvolutionSolver(rosen,
                                             self.bounds,
                                             strategy='rand2exp')
        npt.assert_equal(solver.strategy, 'rand2exp')
        npt.assert_equal(solver.mutation_func.__name__, '_rand2')

        solver = DifferentialEvolutionSolver(rosen,
                                             self.bounds,
                                             strategy='best2bin')
        npt.assert_equal(solver.strategy, 'best2bin')
        npt.assert_equal(solver.mutation_func.__name__, '_best2')

        solver = DifferentialEvolutionSolver(rosen,
                                             self.bounds,
                                             strategy='rand2bin')
        npt.assert_equal(solver.strategy, 'rand2bin')
        npt.assert_equal(solver.mutation_func.__name__, '_rand2')

        solver = DifferentialEvolutionSolver(rosen,
                                             self.bounds,
                                             strategy='rand2exp')
        npt.assert_equal(solver.strategy, 'rand2exp')
        npt.assert_equal(solver.mutation_func.__name__, '_rand2')

        solver = DifferentialEvolutionSolver(rosen,
                                             self.bounds,
                                             strategy='randtobest1bin')
        npt.assert_equal(solver.strategy, 'randtobest1bin')
        npt.assert_equal(solver.mutation_func.__name__, '_randtobest1')

        solver = DifferentialEvolutionSolver(rosen,
                                             self.bounds,
                                             strategy='randtobest1exp')
        npt.assert_equal(solver.strategy, 'randtobest1exp')
        npt.assert_equal(solver.mutation_func.__name__, '_randtobest1')

    def test__mutate1(self):
        #strategies */1/*, i.e. rand/1/bin, best/1/exp, etc.
        result = np.array([0.05])
        trial = self.dummy_solver2._best1(1, (2, 3, 4, 5, 6))
        npt.assert_allclose(trial, result)

        result = np.array([0.25])
        trial = self.dummy_solver2._rand1(1, (2, 3, 4, 5, 6))
        npt.assert_allclose(trial, result)

    def test__mutate2(self):
        #strategies */2/*, i.e. rand/2/bin, best/2/exp, etc.
        #[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]

        result = np.array([-0.1])
        trial = self.dummy_solver2._best2(1, (2, 3, 4, 5, 6))
        npt.assert_allclose(trial, result)

        result = np.array([0.1])
        trial = self.dummy_solver2._rand2(1, (2, 3, 4, 5, 6))
        npt.assert_allclose(trial, result)

    def test__randtobest1(self):
        #strategies randtobest/1/*
        result = np.array([0.1])
        trial = self.dummy_solver2._randtobest1(1, (2, 3, 4, 5, 6))
        npt.assert_allclose(trial, result)

    def test_can_init_with_dithering(self):
        mutation = (0.5, 1)
        solver = DifferentialEvolutionSolver(self.dummy_function,
                                             self.bounds,
                                             mutation=mutation)
        if solver.dither is False:
            raise AttributeError

        if solver.dither is None:
            raise AttributeError

    def test_invalid_mutation_values_arent_accepted(self):
        f = rosen
        mutation = (0.5, 3)
        self.assertRaises(ValueError,
                          DifferentialEvolutionSolver,
                          f,
                          self.bounds,
                          mutation=mutation)

        mutation = (-1, 3)
        self.assertRaises(ValueError,
                          DifferentialEvolutionSolver,
                          f,
                          self.bounds,
                          mutation=mutation)

        mutation = (-1, 'a')
        try:
            DifferentialEvolutionSolver(f,
                                        self.bounds,
                                        mutation=mutation)
            raise Exception("DifferentialEvolutionSolver should've raised some "
                            "error with mutation=(-1, 'a')")
        except Exception:
            #np.isfinite will either raise TypeError or NotImplementedError
            #for a string value
            pass                                  

        mutation = (0.1, np.nan)
        self.assertRaises(ValueError,
                          DifferentialEvolutionSolver,
                          f,
                          self.bounds,
                          mutation=mutation)

        mutation = (0.5)
        solver = DifferentialEvolutionSolver(f,
                                             self.bounds,
                                             mutation=mutation)
        npt.assert_equal(0.5, solver.scale)
        npt.assert_equal(False, solver.dither)

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
        self.dummy_solver._ensure_constraint(trial)
        npt.assert_equal(np.all(trial <= 1), True)

    def test_differential_evolution(self):
        # test that the Jmin of DifferentialEvolutionSolver
        # is the same as the function evaluation
        solver = DifferentialEvolutionSolver(rosen, self.bounds)
        result = solver.solve()
        npt.assert_almost_equal(result.fun, rosen(result.x))

    def test_callback_terminates(self):
        # test that if the callback returns true, then the minimization halts
        bounds = [(0, 2), (0, 2)]

        def callback(param, convergence=0.):
            return True

        result = differential_evolution(rosen, bounds, callback=callback)

        npt.assert_string_equal(result.message,
                                'callback function requested stop early '
                                'by returning True')

    def test_init_with_invalid_strategy(self):
        #test that passing an invalid strategy raises ValueError
        f = rosen
        bounds = [(-3, 3)]
        self.assertRaises(ValueError,
                          differential_evolution,
                          f,
                          bounds,
                          strategy='abc')
                          
    def test_bounds_checking(self):
        # test that the bounds checking works
        f = rosen
        bounds = [(-3, None)]
        self.assertRaises(ValueError,
                          differential_evolution,
                          f,
                          bounds)
        bounds = [(-3)]
        self.assertRaises(ValueError,
                          differential_evolution,
                          f,
                          bounds)
        bounds = [(-3, 3), (3, 4, 5)]
        self.assertRaises(ValueError,
                          differential_evolution,
                          f,
                          bounds)

    def test_select_samples(self):
        # select_samples should return 5 separate random numbers.
        limits = np.arange(12., dtype='float64').reshape(2, 6)
        bounds = list(zip(limits[0, :], limits[1, :]))
        solver = DifferentialEvolutionSolver(None, bounds, popsize=1)
        candidate = 0
        r1, r2, r3, r4, r5 = solver._select_samples(candidate, 5)
        npt.assert_equal(
            len(np.unique(np.array([candidate, r1, r2, r3, r4, r5]))), 6)

    def test_maxiter_stops_solve(self):
        #test that if the maximum number of iterations is exceeded
        #the solver stops
        solver = DifferentialEvolutionSolver(rosen, self.bounds, maxiter=1)
        result = solver.solve()
        npt.assert_equal(result.success, False)
        npt.assert_equal(result.message,
                        'Maximum number of iterations has been exceeded.')

    def test_maxfun_stops_solve(self):
        #test that if the maximum number of function evaluations is exceeded
        #the solver stops
        solver = DifferentialEvolutionSolver(rosen, self.bounds, maxfun=1)
        result = solver.solve()

        npt.assert_equal(result.nfev, 2)
        npt.assert_equal(result.success, False)
        npt.assert_equal(result.message,
                         'Maximum number of function evaluations has '
                              'been exceeded.')

    @npt.dec.slow
    def test_rosen(self):
        # test the Rosenbrock function from object
        solver = DifferentialEvolutionSolver(rosen, self.bounds)
        result = solver.solve()

    @npt.dec.slow
    def test_rosen_from_diff_ev(self):
        # test the Rosenbrock function from differential_evolution function

        result = differential_evolution(rosen, self.bounds)

    def test_exp_runs(self):
        # test whether exponential mutation loop runs
        solver = DifferentialEvolutionSolver(rosen,
                                             self.bounds,
                                             strategy='best1exp',
                                             maxiter=1)

        solver.solve()

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
    npt.run_module_suite()
