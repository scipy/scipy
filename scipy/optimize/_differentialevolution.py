from __future__ import division
import numpy as np
from scipy.optimize import OptimizeResult, minimize
import numbers
import numpy.testing as npt

__all__ = ['differential_evolution']

MACHEPS = 2.22e-16

# standard status messages of optimizers
_status_message = {'success': 'Optimization terminated successfully.',
                   'maxfev': 'Maximum number of function evaluations has '
                              'been exceeded.',
                   'maxiter': 'Maximum number of iterations has been '
                              'exceeded.',
                   'pr_loss': 'Desired error not necessarily achieved due '
                              'to precision loss.',
                   'terminated': 'callback function requested stop early by '
                              'returning True'}


class BoundsError(Exception):

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


def bounds_to_limits(bounds):
    """
        convert tuple of lower and upper bounds to limits
        [(low_0, high_0), ..., (low_n, high_n] 
            -> [[low_0, ..., low_n], [high_0, ..., high_n]]
    """
    return np.array(bounds, float).T


def limits_to_bounds(limits):
    """
        convert limits to tuple of lower and upper bounds
        [[low_0, ..., low_n], [high_0, ..., high_n]] -->
            [(low_0, high_0), ..., (low_n, high_n] 
    """
    return [(limits[0, idx], limits[1, idx])
            for idx in range(np.size(limits, 1))]


def differential_evolution(func, bounds, args=(), strategy='best1bin',
                           maxiter=None, popsize=15, tol=0.01,
                           mutation=(0.5, 1), recombination=0.7, seed=None,
                           callback=None, disp=False, polish=False):
    """
    Find the global minimum of a multivariate function using the differential
    evolution algorithm. It is stochastic in nature (does not use gradient
    methods) to find the minimium, and can search large areas of candidate
    space, but often requires large numbers of function evaluations.

    The algorithm is due to Storn and Price [1]_.

    .. versionadded:: 0.15.0

    Parameters
    ----------
    func : callable
        The objective function to be minimized.  Must be in the form
        ``f(x, *args)``, where `x` is the argument in the form of a 1-D array
        and ``args`` is a  tuple of any additional fixed parameters needed to
        completely specify the function.
    bounds: sequence
        Bounds for variables.  ``(min, max)`` pairs for each element in ``x``,
        defining the lower and upper bounds for the optimizing argument of
        ``func``. It is required to have ``len(bounds) == len(x)``.
        ``len(bounds)`` is used to determine the number of parameters in ``x``.
    args : tuple, optional
        Any additional fixed parameters needed to
        completely specify the objective function.
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

    maxiter: int, optional
        The maximum number of times the entire population is evolved.
        The maximum number of function evaluations is:
        ``maxiter * popsize * len(x)``
    popsize : int, optional
        A multiplier for setting the total population size.  The population has
        ``popsize * len(x)`` individuals.
    tol : float, optional:
        When the mean of the population energies, multiplied by tol,
        divided by the standard deviation of the population energies
        is greater than 1 the solving process terminates:
        ``convergence = mean(pop) * tol / stdev(pop) > 1``
    mutation : float or tuple(float, float), optional:
        The mutation constant.
        If specified as a float it should be in the range [0, 2].
        If specified as a tuple ``(min, max)`` dithering is employed. Dithering
        randomly changes the mutation constant on a generation by generation
        basis. The mutation constant for that generation is taken from
        U[min, max). Dithering can help speed convergence significantly.
        Increasing the mutation constant increases the search radius, but will
        slow down convergence.
    recombination : float, optional:
        The recombination constant, should be in the range [0, 1]. Increasing
        this value allows a larger number of mutants to progress into the next
        generation, but at the risk of population stability.
    seed : int or np.RandomState, optional:
        If seed is not specified the np.RandomState singleton is used.
        If seed is an int, a new np.RandomState instance is used,
        seeded with seed.
        If seed is already a np.RandomState instance, then that
        np.RandomState instance is used.
        Specify seed for repeatable minimizations.
    disp : bool, optional:
        Display status messages
    callback : callable, ``callback(xk, convergence=val)``, optional:
        A function to follow the progress of the minimization. ``xk`` is
        the current value of ``x0``. ``val`` represents the fractional
        value of the population convergence.  When ``val`` is greater than one
        the function halts. If callback returns ``True``, then the minimization
        is halted (any polishing is still carried out).
    polish : bool, optional
        If true, then scipy.optimize.minimize with the `L-BFGS-B` method
        is used to polish the best population member at the end. This requires
        a few more function evaluations.

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes. If ``polish``
        was employed, then OptimizeResult also contains the ``jac`` attribute.
    
    Notes
    -----
    Differential evolution is a stochastic population based method that is
    useful for global optimization problems. At each pass through the population
    the algorithm mutates each candidate solution by mixing with other candidate
    solutions to create a trial candidate. There are several strategies [2]_ for
    creating trial candidates, which suit some problems more than others. The
    'best1bin' strategy is a good starting point for many systems. In this
    strategy two members of the population are randomly chosen. Their difference
    is used to mutate the best member (the `best` in `best1bin`), b_0, so far:
    .. math:: 
        b' = b_0 + mutation * (population[rand0] - population[rand1])
    
    A trial vector is then constructed. Starting with a randomly chosen ``i``'th
    parameter the trial is sequentially filled (in modulo) with parameters from
    ``b'`` or the original candidate. The choice of whether to use ``b'`` or the
    original candidate is made with a binomial distribution (the `bin` in
    `best1bin`) - a random number in ``[0, 1)`` is generated.  If this number is
    less than the ``recombination`` constant then the parameter is loaded from
    ``b'``, otherwise it is loaded from the original candidate.  The final
    parameter is always loaded from ``b'``.  Once the trial candidate is built
    its fitness is assessed. If the trial is better than the original candidate
    then it takes its place. If it is also better than the best overall
    candidate it also replaces that.
    To improve your chances of finding a global minimum use higher ``popsize``
    values, with higher ``mutation``, but lower ``recombination``, values. This
    has the effect of widening the search radius, but slowing convergence.
    
    Examples
    --------
    Let us consider the problem of minimizing the Rosenbrock function. This
    function is implemented in `rosen` in `scipy.optimize`.

    >>> from scipy.optimize import rosen, differential_evolution
    >>> bounds = [(0,2), (0, 2), (0, 2), (0, 2), (0, 2)]
    >>> result = differential_evolution(rosen, bounds)
    >>> print result
    
    Next find the minimum of the Ackley function
    (http://en.wikipedia.org/wiki/Test_functions_for_optimization).
    
    >>> from scipy.optimize import differential_evolution
    >>> from numpy import exp, sqrt, cos, pi, e
    >>> def ackley(x):
    ...     arg1 = -0.2 * sqrt(0.5 * (x[0] ** 2 + x[1] ** 2))
    ...     arg2 = 0.5 * (cos(2. * pi * x[0]) + cos(2. * pi * x[1]))
    ...     return -20. * exp(arg1) - exp(arg2) + 20. + e
    >>> bounds = [(-5, 5), (-5, 5)]
    >>> result = differential_evolution(ackley, bounds)
    >>> print result

    References
    ----------
    .. [1] Storn, R and Price, K, Differential Evolution - a Simple and
           Efficient Heuristic for Global Optimization over Continuous Spaces,
           Journal of Global Optimization, 1997, 11, 341 - 359.
    .. [2] http://www1.icsi.berkeley.edu/~storn/code.html
    .. [3] http://en.wikipedia.org/wiki/Differential_evolution
    """

    # assemble the bounds into the limits
    try:
        limits = bounds_to_limits(bounds)
        npt.assert_equal(np.size(limits, 0), 2)
    except (ValueError, AssertionError):
        # it is required to have (min, max) pairs for each value in x
        raise BoundsError('Bounds should be a sequence containing '
                          'real valued (min, max) pairs for each value'
                          ' in x')

    if type(mutation) is tuple:
        mutation_err_message = 'The mutation constant must be a float in '
        'U[0, 2), or specified as a tuple(min, max) where min < max and min, '
        'max are in U[0, 2.'
        if len(mutation) < 2:
            raise ValueError(mutation_err_message)
        if mutation[0] > mutation[1]:
            mutation = (mutation[1], mutation[0])

    solver = DifferentialEvolutionSolver(func, limits, args=args,
                                         strategy=strategy, maxiter=maxiter,
                                         popsize=popsize, tol=tol,
                                         mutation=mutation,
                                         recombination=recombination,
                                         seed=seed, polish=polish,
                                         callback=callback,
                                         disp=disp)

    result = solver.solve()
    return result


class DifferentialEvolutionSolver(object):

    """
    Parameters
    ----------
    func : callable
        The objective function to be minimized.  Must be in the form
        `f(x, *args)`, where `x` is the argument in the form of a 1-D array
        and `args` is a  tuple of any additional fixed parameters needed to
        completely specify the function.
    limits: 2-D ndarray
        lower and upper limits for the optimizing argument of func. Must have
        shape (2, len(x))
    args : tuple, optional
        Any additional fixed parameters needed to completely
        specify the objective function.
    strategy : str, optional
        The differential evolution strategy to use.
    maxiter: int, optional
        The maximum number of times the entire population is evolved. The
        maximum number of function evaluations is: maxiter * popsize * len(x)
    popsize : int, optional
        A multiplier for setting the total population size.  The population has
        popsize * len(x) individuals.
    tol : float, optional:
        When the mean of the population energies, multiplied by tol,
        divided by the standard deviation of the population energies is greater
        than 1 the solving process terminates.
        i.e. mean(pop) * tol / stdev(pop) > 1
    mutation : float, optional:
        The mutation constant, should be in the range [0, 2].
    recombination : float, optional:
        The recombination constant, should be in the range [0, 1].
    seed : int or np.RandomState, optional:
        If seed is not specified the np.RandomState singleton is used.
        If seed is an int, a new np.RandomState instance is used,
        seeded with seed.
        If seed is already a np.RandomState instance, then that
        np.RandomState instance is used.
        Specify seed for repeatable minimizations.
    callback : callable, optional:
        A function to follow the progress of the minimization.
        Called as ``callback(xk, convergence=val)``, where ``xk`` is the
        current value of ``x0``. ``val`` represents the fractional value
        of the population convergence.  When ``val`` is greater than one
        the function halts. If callback returns True, then the minimization
        is halted (any polishing is still carried out).
    disp : bool, optional
    polish : bool, optional
        If true, then scipy.optimize.minimize with the `L-BFGS-B` method
        is used to polish the best population member at the end.
    """

    def __init__(self, func, limits, args=(),
                 strategy=None, maxiter=None, popsize=15,
                 tol=0.01, mutation=(0.5, 1), recombination=0.7, seed=None,
                 maxfun=None, callback=None, disp=False, polish=False):

        if strategy is not None:
            self.strategy = getattr(
                DifferentialEvolutionSolver, '_' + strategy.lower())
        else:
            self.strategy = getattr(DifferentialEvolutionSolver, '_best1bin')

        self.callback = callback
        self.polish = polish

        self.maxiter = 1000
        if maxiter is not None:
            self.maxiter = maxiter

        self.maxfun = (self.maxiter + 1) * popsize * np.size(limits, 1)
        if maxfun is not None:
            self.maxfun = maxfun

        self.tol = tol
        self.scale = mutation
        self.dither = None
        if type(mutation) is tuple:
            self.dither = mutation

        self.cross_over_probability = recombination

        self.func = func
        self.args = args
        self.limits = limits
        if np.any(np.isnan(limits)):
            raise BoundsError('Bounds should be a sequence'
                              ' containing real valued '
                              '(min, max) pairs for each value in x')

        self.bounds = limits_to_bounds(self.limits)
        # population is scaled to between [0, 1].
        # We have to scale between parameter <-> population
        # save these arguments for _scale_parameter and
        #_unscale_parameter. This is an optimization
        self.__scale_arg1 = 0.5 * (self.limits[0] + self.limits[1])
        self.__scale_arg2 = np.fabs(self.limits[0] - self.limits[1])

        self.nfev = 0
        self.nit = 0

        self.parameter_count = np.size(self.limits, 1)
        self.population_size = popsize * self.parameter_count

        self.random_number_generator = _make_random_gen(seed)

        self.population = self.random_number_generator.rand(
            popsize *
            self.parameter_count,
            self.parameter_count)

        self.population_energies = np.ones(
            popsize * self.parameter_count) * 1.e300

        self.disp = disp

    def solve(self):
        """
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

        self.nfev = 0
        self.nit = 0
        status_message = _status_message['success']
        warning_flag = False

        # population is scaled to between [0, 1].
        # We have to scale between parameter <-> population
        # save these arguments for _scale_parameter and
        #_unscale_parameter. This is an optimization
        self.__scale_arg1 = 0.5 * (self.limits[0] + self.limits[1])
        self.__scale_arg2 = np.fabs(self.limits[0] - self.limits[1])

        # calculate energies to start with
        for index, candidate in enumerate(self.population):
            parameters = self._scale_parameters(candidate)
            self.population_energies[
                index] = self.func(
                parameters,
                *self.args)
            self.nfev += 1

            if self.nfev > self.maxfun:
                warning_flag = True
                status_message = _status_message['maxfev']

        minval = np.argmin(self.population_energies)

        # put the lowest energy into the best solution position.
        lowest_energy = self.population_energies[minval]
        self.population_energies[minval] = self.population_energies[0]
        self.population_energies[0] = lowest_energy

        self.population[[0, minval], :] = self.population[[minval, 0], :]

        # do the optimisation.
        for iteration in xrange(self.maxiter):
            if self.dither is not None:
                self.scale = self.random_number_generator.rand(
                ) * (self.dither[1] - self.dither[0]) + self.dither[0]
            for candidate in xrange(self.population_size):
                if self.nfev >= self.maxfun:
                    warning_flag = True
                    status_message = _status_message['maxfev']
                    break

                trial = self.strategy(self, candidate)
                self._ensure_constraint(trial)
                parameters = self._scale_parameters(trial)

                energy = self.func(parameters, *self.args)
                self.nfev += 1

                if energy < self.population_energies[candidate]:
                    self.population[candidate] = trial
                    self.population_energies[candidate] = energy

                    if energy < self.population_energies[0]:
                        self.population_energies[0] = energy
                        self.population[0] = trial

            # stop when the fractional s.d. of the population is less than tol
            # of the mean energy
            self.convergence = np.std(self.population_energies) / \
                np.abs(np.mean(self.population_energies) + MACHEPS)

            self.nit = iteration + 1

            if self.disp:
                print("differential_evolution step %d: f(x)= %g"
                      % (self.nit,
                         self.population_energies[0]))

            if self.callback:
                if self.callback(
                        self._scale_parameters(self.population[0]),
                        convergence=self.tol / self.convergence) \
                        is True:

                    warning_flag = True
                    status_message = _status_message['terminated']
                    break

            if self.convergence < self.tol:
                break

            if warning_flag:
                break

        if self.nit == self.maxiter:
            status_message = _status_message['maxiter']
            warning_flag = True

        DE_result = OptimizeResult(
            x=self._scale_parameters(self.population[0]),
            fun=self.population_energies[0],
            nfev=self.nfev,
            nit=self.nit,
            message=status_message,
            success=(warning_flag != True))

        if self.polish:
            result = minimize(self.func,
                              np.copy(DE_result.x),
                              method='L-BFGS-B',
                              bounds=self.bounds,
                              args=self.args)

            self.nfev += result.nfev
            DE_result.nfev = self.nfev

            if result.fun < DE_result.fun:
                DE_result.fun = result.fun
                DE_result.x = result.x
#                 DE_result.hess_inv = result.hess_inv
                DE_result.jac = result.jac
                # to keep internal state consistent
                self.population_energies[0] = result.fun
                self.population[0] = self._unscale_parameters(result.x)

        return DE_result

    def _scale_parameters(self, trial):
        # scale from a number between 0 and 1 to parameters
        return self.__scale_arg1 + (trial - 0.5) * self.__scale_arg2

    def _unscale_parameters(self, parameters):
        # scale from parameters to a number between 0 and 1.
        return (parameters - self.__scale_arg1) / self.__scale_arg2 + 0.5

    def _ensure_constraint(self, trial):
        for index, param in enumerate(trial):
            if param > 1 or param < 0:
                trial[index] = self.random_number_generator.rand()

    def _best1bin(self, candidate):
        r0, r1 = self._select_samples(candidate, 2)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        crossovers = self.random_number_generator.rand(self.parameter_count)
        crossovers = crossovers < self.cross_over_probability

        while i < self.parameter_count:
            if crossovers[i] or i == self.parameter_count - 1:
                trial[n] = self.population[0, n] + self.scale * \
                    (self.population[r0, n] - self.population[r1, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _rand1bin(self, candidate):
        r0, r1, r2 = self._select_samples(candidate, 3)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        crossovers = self.random_number_generator.rand(self.parameter_count)
        crossovers = crossovers < self.cross_over_probability

        while i < self.parameter_count:
            if crossovers[i] or i == self.parameter_count - 1:
                trial[n] = self.population[r0, n] + self.scale * \
                    (self.population[r1, n]
                     - self.population[r2, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _best1exp(self, candidate):
        r0, r1 = self._select_samples(candidate, 2)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0
        while i < self.parameter_count and self.random_number_generator.rand() < self.cross_over_probability:
            trial[n] = self.population[0, n] + self.scale * \
                (self.population[r0, n] - self.population[r1, n])
            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _rand1exp(self, candidate):
        r0, r1, r2 = self._select_samples(candidate, 3)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count and self.random_number_generator.rand() < self.cross_over_probability:
            trial[n] = self.population[r0, n] + self.scale * \
                (self.population[r1, n] - self.population[r2, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _randtobest1exp(self, candidate):
        r0, r1 = self._select_samples(candidate, 2)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count and self.random_number_generator.rand() < self.cross_over_probability:
            trial[n] += self.scale * (self.population[0, n] - trial[n]) + \
                self.scale * \
                (self.population[r0, n]
                 - self.population[r1, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _best2exp(self, candidate):
        r0, r1, r2, r3 = self._select_samples(candidate, 4)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count and self.random_number_generator.rand() < self.cross_over_probability:
            trial[n] = self.population[0, n]
            + self.scale * (self.population[r0, n]
                            + self.population[r1, n]
                            - self.population[r2, n]
                            - self.population[r3, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _rand2exp(self, candidate):
        r0, r1, r2, r3, r4 = self._select_samples(candidate, 5)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count and self.random_number_generator.rand() < self.cross_over_probability:
            trial[n] = self.population[r0, n] + self.scale * \
                (self.population[r1, n] + self.population[r2, n]
                 - self.population[r3, n] - self.population[r4, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _randtobest1bin(self, candidate):
        r0, r1 = self._select_samples(candidate, 2)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        crossovers = self.random_number_generator.rand(self.parameter_count)
        crossovers = crossovers < self.cross_over_probability

        while i < self.parameter_count:
            if crossovers[i] or i == self.parameter_count - 1:
                trial[n] += self.scale * \
                    (self.population[0, n] - trial[n]) + \
                    self.scale * \
                    (self.population[r0, n] - self.population[r1, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _best2bin(self, candidate):
        r0, r1, r2, r3 = self._select_samples(candidate, 4)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        crossovers = self.random_number_generator.rand(self.parameter_count)
        crossovers = crossovers < self.cross_over_probability

        while i < self.parameter_count:
            if crossovers[i] or i == self.parameter_count - 1:
                trial[n] = self.population[0, n]    + self.scale * \
                    (self.population[r0, n] + self.population[r1, n]
                     - self.population[r2, n] - self.population[r3, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _rand2bin(self, candidate):
        r0, r1, r2, r3, r4 = self._select_samples(candidate, 5)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        crossovers = self.random_number_generator.rand(self.parameter_count)
        crossovers = crossovers < self.cross_over_probability

        while i < self.parameter_count:
            if crossovers[i] or i == self.parameter_count - 1:
                trial[n] = self.population[r0, n] + self.scale * \
                    (self.population[r1, n] + self.population[r2, n] -
                     self.population[r3, n] - self.population[r4, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _select_samples(self, candidate, number_samples):
        """
            obtain random integers from range(self.population_size), without
            replacement.  You can't have the original candidate either
        """
        idxs = range(self.population_size)
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


if __name__ == "__main__":
    def test():
        from scipy.optimize import rosen
        from numpy import exp, sqrt, cos, pi, e
        # minimum expected at [1, 1, 1, 1, 1]
        bounds = [(0, 2), (0, 2), (0, 2), (0, 2), (0, 2)]
        result = differential_evolution(rosen,
                                        bounds,
                                        seed=1,
                                        polish=True,
                                        disp=False)
        print (result)

        # now do Ackley function
        def ackley(x):
            arg1 = -0.2 * sqrt(0.5 * (x[0] ** 2 + x[1] ** 2))
            arg2 = 0.5 * (cos(2. * pi * x[0]) + cos(2. * pi * x[1]))
            return -20. * exp(arg1) - exp(arg2) + 20. + e
        bounds = [(-5, 5), (-5, 5)]
        result = differential_evolution(ackley,
                                        bounds,
                                        disp=False,
                                        polish=True)
        print (result)

#     import cProfile
#     cProfile.run('test()')
    test()
