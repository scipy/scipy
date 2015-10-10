"""
differential_evolution: The differential evolution global optimization algorithm
Added by Andrew Nelson 2014
"""
from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.optimize import OptimizeResult, minimize
from scipy.optimize.optimize import _status_message
import numbers
import itertools

__all__ = ['differential_evolution']

_MACHEPS = np.finfo(np.float64).eps


def differential_evolution(func, bounds, args=(), strategy='best1bin',
                           maxiter=None, popsize=15, tol=0.01,
                           mutation=(0.5, 1), recombination=0.7, seed=None,
                           callback=None, disp=False, polish=True,
                           init='latinhypercube', workers=1):
    """Finds the global minimum of a multivariate function.
    Differential Evolution is stochastic in nature (does not use gradient
    methods) to find the minimium, and can search large areas of candidate
    space, but often requires larger numbers of function evaluations than
    conventional gradient based techniques.

    The algorithm is due to Storn and Price [1]_.

    Parameters
    ----------
    func : callable
        The objective function to be minimized.  Must be in the form
        ``f(x, *args)``, where ``x`` is the argument in the form of a 1-D array
        and ``args`` is a  tuple of any additional fixed parameters needed to
        completely specify the function. If you are using parallelisation with
        several `workers`, then this function must be picklable.
    bounds : sequence
        Bounds for variables.  ``(min, max)`` pairs for each element in ``x``,
        defining the lower and upper bounds for the optimizing argument of
        `func`. It is required to have ``len(bounds) == len(x)``.
        ``len(bounds)`` is used to determine the number of parameters in ``x``.
    args : tuple, optional
        Any additional fixed parameters needed to completely specify the
        objective function.
    strategy : str, optional
        The differential evolution strategy to use. Should be one of:

            - 'best1bin'
            - 'best1exp'
            - 'rand1exp'
            - 'randtobest1exp'
            - 'best2exp'
            - 'rand2exp'
            - 'randtobest1bin'
            - 'best2bin'
            - 'rand2bin'
            - 'rand1bin'

        The default is 'best1bin'.
    maxiter : int, optional
        The maximum number of generations over which the entire population is
        evolved. The maximum number of function evaluations (with no polishing)
        is: ``(maxiter + 1) * popsize * len(x)``
    popsize : int, optional
        A multiplier for setting the total population size.  The population has
        ``popsize * len(x)`` individuals.
    tol : float, optional
        When the mean of the population energies, multiplied by tol,
        divided by the standard deviation of the population energies
        is greater than 1 the solving process terminates:
        ``convergence = mean(pop) * tol / stdev(pop) > 1``
    mutation : float or tuple(float, float), optional
        The mutation constant. In the literature this is also known as
        differential weight, being denoted by F.
        If specified as a float it should be in the range [0, 2].
        If specified as a tuple ``(min, max)`` dithering is employed. Dithering
        randomly changes the mutation constant on a generation by generation
        basis. The mutation constant for that generation is taken from
        ``U[min, max)``. Dithering can help speed convergence significantly.
        Increasing the mutation constant increases the search radius, but will
        slow down convergence.
    recombination : float, optional
        The recombination constant, should be in the range [0, 1]. In the
        literature this is also known as the crossover probability, being
        denoted by CR. Increasing this value allows a larger number of mutants
        to progress into the next generation, but at the risk of population
        stability.
    seed : int or `np.random.RandomState`, optional
        If `seed` is not specified the `np.RandomState` singleton is used.
        If `seed` is an int, a new `np.random.RandomState` instance is used,
        seeded with seed.
        If `seed` is already a `np.random.RandomState instance`, then that
        `np.random.RandomState` instance is used.
        Specify `seed` for repeatable minimizations.
    disp : bool, optional
        Display status messages
    callback : callable, `callback(xk, convergence=val)`, optional
        A function to follow the progress of the minimization. ``xk`` is
        the current value of ``x0``. ``val`` represents the fractional
        value of the population convergence.  When ``val`` is greater than one
        the function halts. If callback returns `True`, then the minimization
        is halted (any polishing is still carried out).
    polish : bool, optional
        If True (default), then `scipy.optimize.minimize` with the `L-BFGS-B`
        method is used to polish the best population member at the end, which
        can improve the minimization slightly.
    init : string, optional
        Specify how the population initialization is performed. Should be
        one of:

            - 'latinhypercube'
            - 'random'

        The default is 'latinhypercube'. Latin Hypercube sampling tries to
        maximize coverage of the available parameter space. 'random' initializes
        the population randomly - this has the drawback that clustering can
        occur, preventing the whole of parameter space being covered.
    workers : pool object or int, optional
        Optional iterable object which is used for parallelization. It can be
        any object with a map method that follows the same calling sequence as
        the built-in map function, and with poolsize method.
        If int is given as the argument, then a multiprocessing-based pool is
        spawned internally with the corresponding number of parallel processes.
        'mpi4py'-based paralelization and 'joblib'-based parallelization 
        alternative pools can be also used here.

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a `OptimizeResult` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes.  If `polish`
        was employed, and a lower minimum was obtained by the polishing, then
        OptimizeResult also contains the ``jac`` attribute.

    Notes
    -----
    Differential evolution is a stochastic population based method that is
    useful for global optimization problems. At each pass through the population
    the algorithm mutates each candidate solution by mixing with other candidate
    solutions to create a trial candidate. There are several strategies [2]_ for
    creating trial candidates, which suit some problems more than others. The
    'best1bin' strategy is a good starting point for many systems. In this
    strategy two members of the population are randomly chosen. Their difference
    is used to mutate the best member (the `best` in `best1bin`), :math:`b_0`,
    so far:

    .. math::

        b' = b_0 + mutation * (population[rand0] - population[rand1])

    A trial vector is then constructed. Starting with a randomly chosen 'i'th
    parameter the trial is sequentially filled (in modulo) with parameters from
    `b'` or the original candidate. The choice of whether to use `b'` or the
    original candidate is made with a binomial distribution (the 'bin' in
    'best1bin') - a random number in [0, 1) is generated.  If this number is
    less than the `recombination` constant then the parameter is loaded from
    `b'`, otherwise it is loaded from the original candidate.  The final
    parameter is always loaded from `b'`.  Once the trial candidate is built
    its fitness is assessed. If the trial is better than the original candidate
    then it takes its place. If it is also better than the best overall
    candidate it also replaces that.
    To improve your chances of finding a global minimum use higher `popsize`
    values, with higher `mutation` and (dithering), but lower `recombination`
    values. This has the effect of widening the search radius, but slowing
    convergence.

    .. versionadded:: 0.15.0

    Examples
    --------
    Let us consider the problem of minimizing the Rosenbrock function. This
    function is implemented in `rosen` in `scipy.optimize`.

    >>> from scipy.optimize import rosen, differential_evolution
    >>> bounds = [(0,2), (0, 2), (0, 2), (0, 2), (0, 2)]
    >>> result = differential_evolution(rosen, bounds)
    >>> result.x, result.fun
    (array([1., 1., 1., 1., 1.]), 1.9216496320061384e-19)

    Next find the minimum of the Ackley function
    (http://en.wikipedia.org/wiki/Test_functions_for_optimization).

    >>> from scipy.optimize import differential_evolution
    >>> import numpy as np
    >>> def ackley(x):
    ...     arg1 = -0.2 * np.sqrt(0.5 * (x[0] ** 2 + x[1] ** 2))
    ...     arg2 = 0.5 * (np.cos(2. * np.pi * x[0]) + np.cos(2. * np.pi * x[1]))
    ...     return -20. * np.exp(arg1) - np.exp(arg2) + 20. + np.e
    >>> bounds = [(-5, 5), (-5, 5)]
    >>> result = differential_evolution(ackley, bounds)
    >>> result.x, result.fun
    (array([ 0.,  0.]), 4.4408920985006262e-16)

    This example demonstrates the usage of parallelization capabilities.

    >>> def example():
    ...     import numpy as np
    ...     from scipy.optimize import differential_evolution
    ...     from scipy.misc import PPool
    ...     import time
    ...     def ackley(x):
    ...         arg1 = -0.2 * np.sqrt(0.5 * (x[0] ** 2 + x[1] ** 2))
    ...         arg2 = 0.5 * (np.cos(2. * np.pi * x[0]) + np.cos(2. * np.pi * x[1]))
    ...         return -20. * np.exp(arg1) - np.exp(arg2) + 20. + np.e
    ...     
    ...     def objfuncheavy(params):
    ...         for it in range(100000):
    ...             it**2
    ...         return ackley(params)
    ...    
    ...     def objfunclight(params):
    ...         return ackley(params)
    ...     
    ...     bounds = [(-2,2), (-2, 2)]
    ...     p = PPool(n_jobs=10)
    ...     start_time = time.time()
    ...     result = differential_evolution(objfuncheavy,bounds, polish=False, workers=p)
    ...     print("Parallel heavy function: %s seconds ---" % (time.time() - start_time))
    ...     print(result)
    ...     start_time = time.time()
    ...     result = differential_evolution(objfuncheavy, bounds, polish=False)
    ...     print("Serial heavy function: %s seconds ---" % (time.time() - start_time))
    ...     print(result)
    ...     start_time = time.time()
    ...     result = differential_evolution(objfunclight, bounds, polish=False, workers=10)
    ...     print("Parallel light function: %s seconds ---" % (time.time() - start_time))
    ...     start_time = time.time()
    ...     result = differential_evolution(objfunclight, bounds, polish=False)
    ...     print("Serial light function: %s seconds ---" % (time.time() - start_time))
    >>> example() # doctest: +SKIP
    Parallel heavy function: 8.27454996109 seconds ---
        nfev: 3390
     success: True
         fun: 4.4408920985006262e-16
           x: array([ -2.22044605e-16,  -2.22044605e-16])
     message: 'Optimization terminated successfully.'
         nit: 112
    Serial heavy function: 29.4678740501 seconds ---
        nfev: 3150
     success: True
         fun: 4.4408920985006262e-16
           x: array([ 0.,  0.])
     message: 'Optimization terminated successfully.'
         nit: 104
    Parallel light function: 0.465750932693 seconds ---
    Serial light function: 0.209820985794 seconds ---

    
    Results show significant speedup in case of a heavy objective function. 
    The number of required iterations in parallel case can be higher, as the
    mutation in the parallel case is quasi-aggressive (the best solution is
    updated from the best available after evaluation of a subpopulation,
    whereas in the aggressive case the best solution is updated after
    evaluation of every population member). Nevertheless, the gained speedup
    can be still beneficial.
    In case when the objective function is computationally inexpensive, the
    computational overhand due to parallel execution can deteriorate the

    performance. In such a case, the usage of serial version of the algorithm
    is preferred.

    
    References
    ----------
    .. [1] Storn, R and Price, K, Differential Evolution - a Simple and
           Efficient Heuristic for Global Optimization over Continuous Spaces,
           Journal of Global Optimization, 1997, 11, 341 - 359.
    .. [2] http://www1.icsi.berkeley.edu/~storn/code.html
    .. [3] http://en.wikipedia.org/wiki/Differential_evolution
    """

    solver = DifferentialEvolutionSolver(func, bounds, args=args,
                                         strategy=strategy, maxiter=maxiter,
                                         popsize=popsize, tol=tol,
                                         mutation=mutation,
                                         recombination=recombination,
                                         seed=seed, polish=polish,
                                         callback=callback,
                                         disp=disp,
                                         init=init, workers=workers)
    return solver.solve()


class DifferentialEvolutionSolver(object):

    """This class implements the differential evolution solver

    Parameters
    ----------
    func : callable
        The objective function to be minimized.  Must be in the form
        ``f(x, *args)``, where ``x`` is the argument in the form of a 1-D array
        and ``args`` is a  tuple of any additional fixed parameters needed to
        completely specify the function. If you are using parallelisation with
        several `workers`, then this function must be picklable.
    bounds : sequence
        Bounds for variables.  ``(min, max)`` pairs for each element in ``x``,
        defining the lower and upper bounds for the optimizing argument of
        `func`. It is required to have ``len(bounds) == len(x)``.
        ``len(bounds)`` is used to determine the number of parameters in ``x``.
    args : tuple, optional
        Any additional fixed parameters needed to completely specify the
        objective function.
    strategy : str, optional
        The differential evolution strategy to use. Should be one of:

            - 'best1bin'
            - 'best1exp'
            - 'rand1exp'
            - 'randtobest1exp'
            - 'best2exp'
            - 'rand2exp'
            - 'randtobest1bin'
            - 'best2bin'
            - 'rand2bin'
            - 'rand1bin'

        The default is 'best1bin'

    maxiter : int, optional
        The maximum number of generations over which the entire population is
        evolved. The maximum number of function evaluations (with no polishing)
        is: ``(maxiter + 1) * popsize * len(x)``
    popsize : int, optional
        A multiplier for setting the total population size.  The population has
        ``popsize * len(x)`` individuals.
    tol : float, optional
        When the mean of the population energies, multiplied by tol,
        divided by the standard deviation of the population energies
        is greater than 1 the solving process terminates:
        ``convergence = mean(pop) * tol / stdev(pop) > 1``
    mutation : float or tuple(float, float), optional
        The mutation constant. In the literature this is also known as
        differential weight, being denoted by F.
        If specified as a float it should be in the range [0, 2].
        If specified as a tuple ``(min, max)`` dithering is employed. Dithering
        randomly changes the mutation constant on a generation by generation
        basis. The mutation constant for that generation is taken from
        U[min, max). Dithering can help speed convergence significantly.
        Increasing the mutation constant increases the search radius, but will
        slow down convergence.
    recombination : float, optional
        The recombination constant, should be in the range [0, 1]. In the
        literature this is also known as the crossover probability, being
        denoted by CR. Increasing this value allows a larger number of mutants
        to progress into the next generation, but at the risk of population
        stability.
    seed : int or `np.random.RandomState`, optional
        If `seed` is not specified the `np.random.RandomState` singleton is
        used.
        If `seed` is an int, a new `np.random.RandomState` instance is used,
        seeded with `seed`.
        If `seed` is already a `np.random.RandomState` instance, then that
        `np.random.RandomState` instance is used.
        Specify `seed` for repeatable minimizations.
    disp : bool, optional
        Display status messages
    callback : callable, `callback(xk, convergence=val)`, optional
        A function to follow the progress of the minimization. ``xk`` is
        the current value of ``x0``. ``val`` represents the fractional
        value of the population convergence.  When ``val`` is greater than one
        the function halts. If callback returns `True`, then the minimization
        is halted (any polishing is still carried out).
    polish : bool, optional
        If True, then `scipy.optimize.minimize` with the `L-BFGS-B` method
        is used to polish the best population member at the end. This requires
        a few more function evaluations.
    maxfun : int, optional
        Set the maximum number of function evaluations. However, it probably
        makes more sense to set `maxiter` instead.
    init : string, optional
        Specify which type of population initialization is performed. Should be
        one of:

            - 'latinhypercube'
            - 'random'
    workers : pool object or int, optional
        Optional iterable object which is used for parallelization. It can be
        any object with a map method that follows the same calling sequence as
        the built-in map function, and with poolsize method.
        If int is given as the argument, then a multiprocessing-based pool is
        spawned internally with the corresponding number of parallel processes.
        'mpi4py'-based paralelization and 'joblib'-based parallelization 
        alternative pools can be also used here.
    """

    # Dispatch of mutation strategy method (binomial or exponential).
    _binomial = {'best1bin': '_best1',
                 'randtobest1bin': '_randtobest1',
                 'best2bin': '_best2',
                 'rand2bin': '_rand2',
                 'rand1bin': '_rand1'}
    _exponential = {'best1exp': '_best1',
                    'rand1exp': '_rand1',
                    'randtobest1exp': '_randtobest1',
                    'best2exp': '_best2',
                    'rand2exp': '_rand2'}

    def __init__(self, func, bounds, args=(),
                 strategy='best1bin', maxiter=None, popsize=15,
                 tol=0.01, mutation=(0.5, 1), recombination=0.7, seed=None,
                 maxfun=None, callback=None, disp=False, polish=True,
                 init='latinhypercube', workers=1):
        
        self.disp = disp

        # default is serial
        self.__pool = None
        self.__created_pool = False
        self.pool_map = map
        self.poolsize = 1
        self.workers = workers

        if type(workers) is int:
            # We'll deal with creating our own pool in solve. This is because
            # we have to terminate them ourselves when the solve has finished.
            pass
        else:
            try:
                self.__pool = workers
                self.pool_map = workers.map
                self.poolsize = workers.poolsize()
            except:
                raise ValueError('Workers keyword expected a pool object with '
                                 'map and poolsize methods')
        
        if strategy in self._binomial:
            self.mutation_func = getattr(self, self._binomial[strategy])
        elif strategy in self._exponential:
            self.mutation_func = getattr(self, self._exponential[strategy])
        else:
            raise ValueError("Please select a valid mutation strategy")
        self.strategy = strategy

        self.callback = callback
        self.polish = polish
        self.tol = tol

        # Mutation constant should be in [0, 2). If specified as a sequence
        # then dithering is performed.
        self.scale = mutation
        if (not np.all(np.isfinite(mutation)) or
                np.any(np.array(mutation) >= 2) or
                np.any(np.array(mutation) < 0)):
            raise ValueError('The mutation constant must be a float in '
                             'U[0, 2), or specified as a tuple(min, max)'
                             ' where min < max and min, max are in U[0, 2).')

        self.dither = None
        if hasattr(mutation, '__iter__') and len(mutation) > 1:
            self.dither = [mutation[0], mutation[1]]
            self.dither.sort()

        self.cross_over_probability = recombination

        self.func = func
        self.args = args

        # convert tuple of lower and upper bounds to limits
        # [(low_0, high_0), ..., (low_n, high_n]
        #     -> [[low_0, ..., low_n], [high_0, ..., high_n]]
        self.limits = np.array(bounds, dtype='float').T
        if (np.size(self.limits, 0) != 2 or 
                not np.all(np.isfinite(self.limits))):
            raise ValueError('bounds should be a sequence containing '
                             'real valued (min, max) pairs for each value'
                             ' in x')

        self.maxiter = maxiter or 1000
        self.maxfun = (maxfun or ((self.maxiter + 1) * popsize *
                                  np.size(self.limits, 1)))

        # population is scaled to between [0, 1].
        # We have to scale between parameter <-> population
        # save these arguments for _scale_parameter and
        # _unscale_parameter. This is an optimization
        self.__scale_arg1 = 0.5 * (self.limits[0] + self.limits[1])
        self.__scale_arg2 = np.fabs(self.limits[0] - self.limits[1])

        parameter_count = np.size(self.limits, 1)
        self.random_number_generator = _make_random_gen(seed)

        # default initialization is a latin hypercube design, but there
        # are other population initializations possible.
        self.population = np.zeros((popsize * parameter_count,
                                    parameter_count))
        if init == 'latinhypercube':
            self.init_population_lhs()
        elif init == 'random':
            self.init_population_random()
        else:
            raise ValueError("The population initialization method must be one"
                             "of 'latinhypercube' or 'random'")

        self.population_energies = np.ones(
            popsize * parameter_count) * np.inf

    def init_population_lhs(self):
        """
        Initializes the population with Latin Hypercube Sampling
        Latin Hypercube Sampling ensures that the sampling of parameter space
        is maximised.
        """
        samples = np.size(self.population, 0)
        N = np.size(self.population, 1)
        rng = self.random_number_generator

        # Generate the intervals
        segsize = 1.0 / samples

        # Fill points uniformly in each interval
        rdrange = rng.rand(samples, N) * segsize
        rdrange += np.atleast_2d(
            np.linspace(0., 1., samples, endpoint=False)).T

        # Make the random pairings
        self.population = np.zeros_like(rdrange)

        for j in range(N):
            order = rng.permutation(range(samples))
            self.population[:, j] = rdrange[order, j]

    def init_population_random(self):
        """
        Initialises the population at random.  This type of initialization
        can possess clustering, Latin Hypercube sampling is generally better.
        """
        rng = self.random_number_generator
        self.population = rng.random_sample(self.population.shape)

    @property
    def x(self):
        """
        The best solution from the solver

        Returns
        -------
        x - ndarray
            The best solution from the solver.
        """
        return self._scale_parameters(self.population[0])

    def solve(self):
        """
        Runs the DifferentialEvolutionSolver.

        Returns
        -------
        res : OptimizeResult
            The optimization result represented as a ``OptimizeResult`` object.
            Important attributes are: ``x`` the solution array, ``success`` a
            Boolean flag indicating if the optimizer exited successfully and
            ``message`` which describes the cause of the termination. See
            `OptimizeResult` for a description of other attributes.  If `polish`
            was employed, and a lower minimum was obtained by the polishing,
            then OptimizeResult also contains the ``jac`` attribute.
        """

        # if int argument is supplied, then use either serial or parallel pool
        self.__created_pool = False
        if type(self.workers) is int:
            # parallel case
            if self.workers != 1:
                from scipy.misc import PPool
                if (self.workers > 1) and (self.workers < 100):
                    pool = PPool(n_jobs=self.workers)
                # if number of processes is negative or insane
                else:
                    pool = PPool()

                self.pool_map = pool.map
                self.poolsize = pool.poolsize()
                self.__pool = pool
                self.__created_pool = True

                if self.disp:
                    print("Starting %g workers in parallel." % self.poolsize)

        nfev, nit, warning_flag = 0, 0, False
        status_message = _status_message['success']
        
        # calculate energies to start with for the whole population
        parameters = []
        for candidate in self.population:
            # incapsulate additional fixed arguments to parameters
            parameters.append(self._scale_parameters(candidate))

        # calculate starting energies for the whole population
        parameters = self._scale_parameters(self.population)
        energies = self.pool_map(_wrapper,
                                 zip(itertools.repeat(self.func),
                                 parameters,
                                 itertools.repeat(self.args)))

        # the squeeze is necessary because some objective functions return
        # arrays instead of floats.
        self.population_energies = np.fromiter(energies, float).squeeze()

        nfev += len(self.population)
        
        # check the number of evaluation of the functions
        # (this, perhaps, should be deprecated as unnecessary)
        if nfev > self.maxfun:
            warning_flag = True
            status_message = _status_message['maxfev']

        minval = np.argmin(self.population_energies)

        # put the lowest energy into the best solution position

        self.population_energies[[minval, 0]] = (
            self.population_energies[[0, minval]])

        # and exchange places of previous and new best solutions
        self.population[[0, minval], :] = self.population[[minval, 0], :]

        if warning_flag:
            return OptimizeResult(
                x=self.x,
                fun=self.population_energies[0],
                nfev=nfev,
                nit=nit,
                message=status_message,
                success=(warning_flag is not True))

        # divide populations into subpopulations depending on poolsize
        # determine number of sub-populations 'nsp', or number of
        # parallel evaluations, which is 'nsp + 1'
        # and the length of the reminder 'lenrem'
        nsp, lenrem = divmod(np.size(self.population, 0), self.poolsize)
        sp_sizes = [self.poolsize] * nsp
        if lenrem:
            sp_sizes += [lenrem]

        # do the optimisation, evolving over several generations
        for nit in range(1, self.maxiter + 1):
        
            if self.dither is not None:
                self.scale = self.random_number_generator.rand(
                ) * (self.dither[1] - self.dither[0]) + self.dither[0]

            # within a generation iterate among sub-populations
            # when self.poolsize == 1 the mutation is aggressive
            # when self.poolsize > 1 the mutation is quasi-aggressive
            for it, sp_size in enumerate(sp_sizes):
                # range of candidates we're evaluating in the subpopulation
                cnd_range = range(self.poolsize*it, self.poolsize*it + sp_size)

                # create trial vectors through mutation
                trials = np.array([self._mutate(cand) for cand in cnd_range])

                # ensure parameters are within the limits
                trials[trials < 0] = (
                    np.random.random(np.count_nonzero(trials < 0)))
                trials[trials > 1] = (
                    np.random.random(np.count_nonzero(trials > 1)))

                # scale trials to parameter values
                parameters = self._scale_parameters(trials)
                nfev += sp_size

                # calculate the energies of the trials
                spenergies = self.pool_map(_wrapper,
                                           zip(itertools.repeat(self.func),
                                           parameters,
                                           itertools.repeat(self.args)))
                spenergies = np.fromiter(spenergies, float).squeeze()

                # check the number of evaluation of the functions
                # (this, perhaps, should be deprecated as unnecessary)
                if nfev > self.maxfun:
                    warning_flag = True
                    status_message = _status_message['maxfev']
                    break

                # find out which trial candidates have have lower energy than
                # the existing population and replace the original population
                # members
                improved = spenergies < self.population_energies[cnd_range]

                new_pop = np.where(improved[:, np.newaxis],
                                   trials,
                                   self.population[cnd_range])
                self.population[cnd_range] = new_pop

                # also replace the energies if they got lower
                new_energy = np.where(improved,
                                      spenergies,
                                      self.population_energies[cnd_range])
                self.population_energies[cnd_range] = new_energy

                # the overall best solution may have changed. If so, replace it
                minval = np.argmin(self.population_energies)
                if minval:
                    self.population_energies[[minval, 0]] = (
                        self.population_energies[[0, minval]])
                    self.population[[0, minval], :] = (
                        self.population[[minval, 0], :])
                    if self.disp:
                        print(" Best updated: f(x)= %g" %
                              self.population_energies[0])
                        print(self._scale_parameters(self.population[0]))

            # report on the results of the current generation
            if self.disp:
                print("differential_evolution step %d: f(x)= %g"
                      % (nit,
                         self.population_energies[0]))

            # stop when the fractional s.d. of the population is less than tol
            # of the mean energy
            convergence = (np.std(self.population_energies) /
                           np.abs(np.mean(self.population_energies) + _MACHEPS))

            if (self.callback and
                    self.callback(self._scale_parameters(self.population[0]),
                                  convergence=self.tol / convergence) is True):

                warning_flag = True
                status_message = ('callback function requested stop early '
                                  'by returning True')
                break

            if convergence < self.tol or warning_flag:
                break

        else:
            status_message = _status_message['maxiter']
            warning_flag = True

        DE_result = OptimizeResult(
            x=self.x,
            fun=self.population_energies[0],
            nfev=nfev,
            nit=nit,
            message=status_message,
            success=(warning_flag is not True))

        if self.polish:
            result = minimize(self.func,
                              np.copy(DE_result.x),
                              method='L-BFGS-B',
                              bounds=self.limits.T,
                              args=self.args)

            nfev += result.nfev
            DE_result.nfev = nfev

            if result.fun < DE_result.fun:
                DE_result.fun = result.fun
                DE_result.x = result.x
                DE_result.jac = result.jac
                # to keep internal state consistent
                self.population_energies[0] = result.fun
                self.population[0] = self._unscale_parameters(result.x)

        if self.__created_pool:
            self.__pool.terminate()

        return DE_result

    def _scale_parameters(self, trial):
        """
        scale from a number between 0 and 1 to parameters
        """
        return self.__scale_arg1 + (trial - 0.5) * self.__scale_arg2

    def _unscale_parameters(self, parameters):
        """
        scale from parameters to a number between 0 and 1.
        """
        return (parameters - self.__scale_arg1) / self.__scale_arg2 + 0.5

    def _ensure_constraint(self, trial):
        """
        make sure the parameters lie between the limits
        """
        for index, param in enumerate(trial):
            if param > 1 or param < 0:
                trial[index] = self.random_number_generator.rand()

    def _mutate(self, candidate):
        """
        create a trial vector based on a mutation strategy
        """
        trial = np.copy(self.population[candidate])
        parameter_count = np.size(trial, 0)

        fill_point = self.random_number_generator.randint(0, parameter_count)

        if (self.strategy == 'randtobest1exp' or
                self.strategy == 'randtobest1bin'):
            bprime = self.mutation_func(candidate,
                                        self._select_samples(candidate, 5))
        else:
            bprime = self.mutation_func(self._select_samples(candidate, 5))

        if self.strategy in self._binomial:
            crossovers = self.random_number_generator.rand(parameter_count)
            crossovers = crossovers < self.cross_over_probability
            # the last one is always from the bprime vector for binomial
            # If you fill in modulo with a loop you have to set the last one to
            # true. If you don't use a loop then you can have any random entry
            # be True.
            crossovers[fill_point] = True
            trial = np.where(crossovers, bprime, trial)
            return trial

        elif self.strategy in self._exponential:
            i = 0
            while (i < parameter_count and
                   self.random_number_generator.rand() <
                   self.cross_over_probability):

                trial[fill_point] = bprime[fill_point]
                fill_point = (fill_point + 1) % parameter_count
                i += 1

            return trial

    def _best1(self, samples):
        """
        best1bin, best1exp
        """
        r0, r1 = samples[:2]
        return (self.population[0] + self.scale *
                (self.population[r0] - self.population[r1]))

    def _rand1(self, samples):
        """
        rand1bin, rand1exp
        """
        r0, r1, r2 = samples[:3]
        return (self.population[r0] + self.scale *
                (self.population[r1] - self.population[r2]))

    def _randtobest1(self, candidate, samples):
        """
        randtobest1bin, randtobest1exp
        """
        r0, r1 = samples[:2]
        bprime = np.copy(self.population[candidate])
        bprime += self.scale * (self.population[0] - bprime)
        bprime += self.scale * (self.population[r0] -
                                self.population[r1])
        return bprime

    def _best2(self, samples):
        """
        best2bin, best2exp
        """
        r0, r1, r2, r3 = samples[:4]
        bprime = (self.population[0] + self.scale *
                  (self.population[r0] + self.population[r1] -
                   self.population[r2] - self.population[r3]))

        return bprime

    def _rand2(self, samples):
        """
        rand2bin, rand2exp
        """
        r0, r1, r2, r3, r4 = samples
        bprime = (self.population[r0] + self.scale *
                  (self.population[r1] + self.population[r2] -
                   self.population[r3] - self.population[r4]))

        return bprime

    def _select_samples(self, candidate, number_samples):
        """
        obtain random integers from range(np.size(self.population, 0)),
        without replacement.  You can't have the original candidate either.
        """
        idxs = list(range(np.size(self.population, 0)))
        idxs.remove(candidate)
        self.random_number_generator.shuffle(idxs)
        idxs = idxs[:number_samples]
        return idxs


def _wrapper(xarg):
    """
    Wrapper to allow the objective function to be called in a multiprocessing
    context.

    Parameters
    ----------
    xarg: tuple
        Tuple containing the objective function and all its parameters.
        xarg[0]: callable
            The objective function.
        xarg[1]: np.ndarray
            Solution array.
        xarg[2]: tuple
            The extra parameters required to fully specify the objective
            function.
    Returns
    -------
    Evaluated function value
    """
    func = xarg[0]
    return func(xarg[1], *xarg[2])


def _make_random_gen(seed):
    """Turn seed into a np.random.RandomState instance

    If seed is None, return the RandomState singleton used by np.random.
    If seed is an int, return a new RandomState instance seeded with seed.
    If seed is already a RandomState instance, return it.
    Otherwise raise ValueError.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (numbers.Integral, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)
