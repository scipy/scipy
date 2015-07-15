"""
differential_evolution: The differential evolution global optimization algorithm
Added by Andrew Nelson 2014
Parallel version by Pavel Ponomarev 2015
"""
from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.optimize import OptimizeResult, minimize
from scipy.optimize.optimize import _status_message
from scipy.optimize.pools.Spool import Spool
import numbers

__all__ = ['differential_evolution']

_MACHEPS = np.finfo(np.float64).eps


def differential_evolution(func, bounds, args=(), strategy='best1bin',
                           maxiter=None, popmul=15, tol=0.01,
                           mutation=(0.5, 1), recombination=0.7, seed=None,
                           callback=None, disp=False, polish=True,
                           init='latinhypercube', pool=None):
    """Finds the global minimum of a multivariate function.
    Differential Evolution is stochastic in nature (does not use gradient
    methods) to find the minimium, and can search large areas of candidate
    space, but often requires larger numbers of function evaluations than
    conventional gradient based techniques. Suitable to find the minimum of 
    non-differentiable, non-linear and noizy functions.

    The algorithm is due to Storn and Price [1]_.

    Parameters
    ----------
    func : callable
        The objective function to be minimized.  Must be in the form
        ``f(x, *args)``, where ``x`` is the argument in the form of a 1-D array
        and ``args`` is a  tuple of any additional fixed parameters needed to
        completely specify the function.
    bounds : sequence
        Bounds for variables.  ``(min, max)`` pairs for each element in ``x``,
        defining the lower and upper bounds for the optimizing argument of
        `func`. It is required to have ``len(bounds) == len(x)``.
        ``len(bounds)`` is used to determine the number of parameters in ``x``.
    args : tuple, optional
        Any additional fixed parameters needed to
        completely specify the objective function.
    strategy : str, optional
        The mutation strategy to use. Should be one of:

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
        The maximum number of times the entire population is evolved.
        The maximum number of function evaluations (with no polishing) is:
        ``(maxiter + 1) * popmul * len(x)``
        Default is 1000.
    popmul : int, optional
        A multiplier for setting the total population size. The population has
        ``popmul * len(x)`` individuals. Default is `popmul = 15`.
    tol : float, optional
        When the mean of the population energies, multiplied by tol,
        divided by the standard deviation of the population energies
        is greater than 1 the solving process terminates:
        ``convergence = mean(pop) * tol / stdev(pop) > 1``
        Default value is 'tol = 0.01'
    mutation : float or tuple(float, float), optional
        The mutation constant.
        If specified as a float it should be in the range [0, 2].
        If specified as a tuple ``(min, max)`` dithering is employed. Dithering
        randomly changes the mutation constant on a generation by generation
        basis. The mutation constant for that generation is taken from
        ``U[min, max)``. Dithering can help speed convergence significantly.
        Increasing the mutation constant increases the search radius, but will
        slow down convergence.
    pool : iterable, optional
        Optional iterable object which is used for parallelization. It can be
        any object with a map method that follows the same calling sequence as
        the built-in map function. There are two helper classes -- MPIpool and 
        JLpool that provide 'mpi4py'-based paralelization and 'joblib'-based
        parallelization alternatives.
    recombination : float, optional
        The recombination constant, should be in the range [0, 1]. Increasing
        this value allows a larger number of mutants to progress into the next
        generation, but at the risk of population stability.
    seed : int or `np.random.RandomState`, optional
        If `seed` is not specified the `np.RandomState` singleton is used.
        If `seed` is an int, a new `np.random.RandomState` instance is used,
        seeded with seed.
        If `seed` is already a `np.random.RandomState instance`, then that
        `np.random.RandomState` instance is used.
        Specify `seed` for repeatable minimizations.
    disp : bool, optional
        Display status messages.
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

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a `OptimizeResult` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes. If `polish`
        was employed, then OptimizeResult also contains the `jac` attribute.

    Notes
    -----
    Differential evolution is a stochastic population based method that is
    useful for constrained, parameter bound, nonlinear, disconinuous
    global optimization problems. At each pass through the population
    the algorithm mutates each candidate by mixing its parameters with other
    candidates to create a trial candidate. There are several strategies [2]_ for
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
    parameter is always loaded from `b'`. Once the trial candidate is built
    its fitness is assessed. If the trial is better than the original candidate
    then it takes its place. If it is also better than the best overall
    candidate it also replaces that.
    To improve your chances of finding a global minimum use higher `popmul`
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
                                         popmul=popmul, tol=tol,
                                         mutation=mutation,
                                         recombination=recombination,
                                         seed=seed, polish=polish,
                                         callback=callback,
                                         disp=disp,
                                         init=init, pool=pool)
    return solver.solve()


class DifferentialEvolutionSolver(object):

    """This class implements the differential evolution solver

    Parameters are equal to he parameters of the module.

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
                 strategy='best1bin', maxiter=None, popmul=15,
                 tol=0.01, mutation=(0.5, 1), recombination=0.7, seed=None,
                 callback=None, disp=False, polish=True,
                 init='latinhypercube', pool=None):
        
        # use aggressive mutation strategy (or quasi-aggressive in case of parallel)
        # original algorythm was aggressive
        self.aggressive = True
        
        if pool is None:
            self.pool = Spool()
        else:
            self.pool = pool
            print("Starting %g workers in parallel." % self.pool.poolsize())
        
        
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

        #Mutation constant should be in [0, 2). If specified as a sequence
        #then dithering is performed.
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
        if (np.size(self.limits, 0) != 2
                or not np.all(np.isfinite(self.limits))):
            raise ValueError('bounds should be a sequence containing '
                             'real valued (min, max) pairs for each value'
                             ' in x')

        self.maxiter = maxiter or 1000


        # population is scaled to between [0, 1].
        # We have to scale between parameter <-> population
        # save these arguments for _scale_parameter and
        # _unscale_parameter. This is an optimization
        self.__scale_arg1 = 0.5 * (self.limits[0] + self.limits[1])
        self.__scale_arg2 = np.fabs(self.limits[0] - self.limits[1])

        parameter_count = np.size(self.limits, 1)
        self.random_number_generator = _make_random_gen(seed)

        #default initialization is a latin hypercube design, but there
        #are other population initializations possible.
        self.population = np.zeros((popmul * parameter_count,
                                    parameter_count))
        if init == 'latinhypercube':
            self.init_population_lhs()
        elif init == 'random':
            self.init_population_random()
        else:
            raise ValueError("The population initialization method must be one"
                             "of 'latinhypercube' or 'random'")

        self.population_energies = np.ones(
            popmul * parameter_count) * np.inf

        self.disp = disp

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
            `OptimizeResult` for a description of other attributes. If polish
            was employed, then OptimizeResult also contains the ``hess_inv`` and
            ``jac`` attributes.
        """

        nfev, nit, warning_flag = 0, 0, False
        status_message = _status_message['success']

        # calculate energies to start with
        params = []
        for candidate in self.population:
            params.append(self._scale_parameters(candidate), *(self.args))
            nfev += 1

        self.population_energies = self.pool.map(self.func, params)
        
        minval = np.argmin(self.population_energies)

        # put the lowest energy into the best solution position.
        lowest_energy = self.population_energies[minval]
        self.population_energies[minval] = self.population_energies[0]
        self.population_energies[0] = lowest_energy
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

        # do the optimisation.
        for nit in range(1, self.maxiter + 1):
        
            if self.dither is not None:
                self.scale = self.random_number_generator.rand(
                    ) * (self.dither[1] - self.dither[0]) + self.dither[0]

            # brake the population to subpopulation depending on the size of the pool
            # determine number of sub-populations 'nsp', or number of 
            # parallel evaluations, which is 'nsp + 1'
            # and the length of the reminder 'lenrem'
            nsp, lenrem = divmod(np.size(self.population, 0), self.pool.poolsize())
            
            itsp = 0
            # iterate among sub-populations
            # when self.pool.poolsize() == 1 the mutation is aggressive
            # when self.pool.poolsize() > 1 the mutation is quasi-aggressive
            while itsp <= nsp:
                # determine the length of the subpopulation
                lensp = self.pool.poolsize() if (itsp < nsp) else lenrem

                # initialize list for all members (parameters) of the current subpopulation
                spparams = []
                # mutate parameters for the sub-population
                for candidate in xrange(lensp):
                    trial = self._mutate(candidate + self.pool.poolsize()*itsp)
                    self._ensure_constraint(trial)
                    spparams.append(self._scale_parameters(trial), *(self.args))
                    nfev += 1
                    
                # in parallel case the self.func must return a list of energies
                # for the whole subpopulation
                spenergies = self.pool.map(self.func, spparams)

                # update population and their energies if subpopulation members are 
                # better by iteration among all members (or jobs) of the subpopulation
                for itjob in xrange(lensp):
                    energy = spenergies[itjob]
                    if energy < self.population_energies[itjob + self.pool.poolsize()*itsp]:
                        self.population[itjob + self.pool.poolsize()*itsp] =\
                            self._unscale_parameters(spparams[itjob])
                        self.population_energies[itjob + self.pool.poolsize()*itsp] = energy
                        
                        # update global best if there is a better in the current sub-population
                        # and strategy is aggressive
                        if self.aggressive:
                            if energy < self.population_energies[0]:
                                self.population_energies[0] = energy
                                self.population[0] =\
                                    self._unscale_parameters(spparams[itjob])
                                if self.disp:
                                    print(" Best updated: f(x)= %g"
                                          % (self.population_energies[0]))
                                    print(self._scale_parameters(self.population[0]))

                                # exchange places between old and new global best    
                                self.population[[0, itjob + self.pool.poolsize()*itsp], :] =\
                                    self.population[[itjob + self.pool.poolsize()*itsp, 0], :]
                itsp += 1
            


                                  
            # report on the results of the current generation
            if self.disp:
                print("differential_evolution step %d: f(x)= %g"
                      % (nit,
                         self.population_energies[0]))
                print(self._scale_parameters(self.population[0]))

            if (self.callback and
                    self.callback(self._scale_parameters(self.population[0]),
                                  convergence=self.tol / convergence) is True):

                warning_flag = True
                status_message = ('callback function requested stop early '
                                  'by returning True')
                break

            # stop when the fractional s.d. of the population is less than tol
            # of the mean energy
            convergence = (np.std(self.population_energies) /
                           np.abs(np.mean(self.population_energies) +
                                  _MACHEPS))

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

        if (self.strategy == 'randtobest1exp'
                or self.strategy == 'randtobest1bin'):
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
                            (self.population[r0] + self.population[r1]
                           - self.population[r2] - self.population[r3]))

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
