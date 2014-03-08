from __future__ import division
import numpy as np
from scipy.optimize import OptimizeResult, minimize
import numbers

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

def differential_evolution(func, bounds, args=(), DEstrategy=None,
                           maxiter=None, popsize=15, tol=0.01,
                           mutation=0.8, recombination=0.8, seed=None,
                           callback=None, disp=False, polish=False, **options):
    """
    Find the global minimum of a multivariate function using the differential
    evolution algorithm. It is stochastic in nature (does not use gradient
    methods) to find the minimium, and can search large areas of candidate
    space, but often requires large numbers of function evaluations.

    The algorithm is originally due to Storn and Price [1]_.

    .. versionadded:: 0.15.0

    Parameters
    ----------
    func : callable
        The objective function to be minimized.  Must be in the form
        `f(x, *args)`, where `x` is the argument in the form of a 1-D array
        and `args` is a  tuple of any additional fixed parameters needed to
        completely specify the function.
    bounds: sequence
        Bounds for variables.  ``(min, max)`` pairs for each element in ``x``,
        defining the lower and upper bounds for the optimizing argument of func.
        It is required to have len(bounds) == len(x)
    args : tuple, optional
        Any additional fixed parameters needed to
        completely specify the objective function.
    DEstrategy : str, optional
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
        maxiter * popsize * len(x)
    popsize : int, optional
        A multiplier for setting the total population size.  The population has
        popsize * len(x) individuals.
    tol : float, optional:
        When the mean of the population energies, multiplied by tol,
        divided by the standard deviation of the population energies
        is greater than 1 the solving process terminates.
        i.e. convergence = mean(pop) * tol / stdev(pop) > 1
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
    disp : bool, optional:
        Display status messages
    callback : callable, optional:
        A function to follow the progress of the minimization.
        Called as ``callback(xk, convergence=val)``, where ``xk`` is
        the current value of ``x0``. ``val`` represents the fractional
        value of the population convergence.  When ``val`` is greater than one
        the function halts.
        If callback returns True, then the minimization is halted (any
        polishing is still carried out).
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
        `OptimizeResult` for a description of other attributes. If polish was
        employed, then OptimizeResult also contains the ``jac`` attribute.

    References
    ----------
    .. [1] http://www1.icsi.berkeley.edu/~storn/code.html
    .. [2] http://en.wikipedia.org/wiki/Differential_evolution

    Examples
    --------
    Let us consider the problem of minimizing the Rosenbrock function. This
    function is implemented in `rosen` in `scipy.optimize`.

    >>> from scipy.optimize import rosen, differential_evolution
    >>> bounds = [(0,2), (0, 2), (0, 2), (0, 2), (0, 2)]
    >>> result = differential_evolution(rosen, bounds)
    >>> print result

    """

    # assemble the bounds into the limits
    try:
        limits = bounds_to_limits(bounds)
        assert np.size(limits, 0) == 2
    except (ValueError, AssertionError):
        # it is required to have (min, max) pairs for each value in x
        raise BoundsError('Bounds should be a sequence containing '
                          'real valued (min, max) pairs for each value'
                          ' in x')

    solver = DifferentialEvolutionSolver(func, limits, args=args,
                                         DEstrategy=DEstrategy, maxiter=maxiter,
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
    DEstrategy : str, optional
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
                 DEstrategy=None, maxiter=None, popsize=15,
                 tol=0.01, mutation=0.8, recombination=0.8, seed=None,
                 maxfun=None, callback=None, disp=False, polish=False,
                 **options):

        if DEstrategy is not None:
            self.DEstrategy = getattr(
                DifferentialEvolutionSolver, '_' + DEstrategy.lower())
        else:
            self.DEstrategy = getattr(DifferentialEvolutionSolver, '_best1bin')

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
        self.cross_over_probability = recombination

        self.func = func
        self.args = args
        self.limits = limits
        if np.any(np.isnan(limits)):
            raise BoundsError('Bounds should be a sequence'
                              ' containing real valued '
                              '(min, max) pairs for each value in x')

        self.bounds = limits_to_bounds(self.limits)

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
            for candidate in xrange(self.population_size):
                if self.nfev >= self.maxfun:
                    warning_flag = True
                    status_message = _status_message['maxfev']
                    break

                trial = self.DEstrategy(self, candidate)
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
        return (
            0.5 * (self.limits[0] + self.limits[1]) +
            (trial - 0.5) * np.fabs(self.limits[0] - self.limits[1])
        )

    def _unscale_parameters(self, parameters):
        return (parameters - 0.5 * (self.limits[0] + self.limits[1])) / \
            np.fabs(self.limits[0] - self.limits[1]) + 0.5

    def _ensure_constraint(self, trial):
        for index, param in enumerate(trial):
            if param > 1 or param < 0:
                trial[index] = self.random_number_generator.rand()

    def _best1bin(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 0, 0, 0)
        n = self.random_number_generator.randint(0, self.parameter_count)
        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count:
            if self.random_number_generator.rand() < self.cross_over_probability or i == self.parameter_count - 1:
                trial[n] = self.population[0, n] + self.scale * \
                    (self.population[r1, n] - self.population[r2, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _best1exp(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 0, 0, 0)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0
        while i < self.parameter_count and self.random_number_generator.rand() < self.cross_over_probability:
            trial[n] = self.population[0, n] + self.scale * \
                (self.population[r1, n] - self.population[r2, n])
            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _rand1exp(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 1, 0, 0)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count and self.random_number_generator.rand() < self.cross_over_probability:
            trial[n] = self.population[r1, n] + self.scale * \
                (self.population[r2, n] - self.population[r3, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _randtobest1exp(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 0, 0, 0)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count and self.random_number_generator.rand() < self.cross_over_probability:
            trial[n] += self.scale * (self.population[0, n] - trial[n]) + \
                self.scale * \
                (self.population[r1, n]
                 - self.population[r2, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _best2exp(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 1, 1, 0)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count and self.random_number_generator.rand() < self.cross_over_probability:
            trial[n] = self.population[0, n]
            + self.scale * (self.population[r1, n]
                            + self.population[r2, n]
                            - self.population[r3, n]
                            - self.population[r4, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _rand2exp(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 1, 1, 1)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count and self.random_number_generator.rand() < self.cross_over_probability:
            trial[n] = self.population[r1, n]
            + self.scale * (self.population[r2, n]
                            + self.population[r3, n]
                            - self.population[r4, n]
                            - self.population[r5, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _randtobest1bin(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 0, 0, 0)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count:
            if self.random_number_generator.rand() < self.cross_over_probability or i == self.parameter_count - 1:
                trial[n] += self.scale * (self.population[0, n] - trial[n])
                + self.scale * \
                    (self.population[r1, n] - self.population[r2, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _best2bin(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 1, 1, 0)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count:
            if self.random_number_generator.rand() < self.cross_over_probability or i == self.parameter_count - 1:
                trial[n] = self.population[0, n]
                + self.scale * (self.population[r1, n]
                                + self.population[r2, n]
                                - self.population[r3, n]
                                - self.population[r4, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _rand2bin(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 1, 1, 1)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count:
            if self.random_number_generator.rand() < self.cross_over_probability or i == self.parameter_count - 1:
                trial[n] = self.population[r1, n]
                + self.scale * (self.population[r2, n]
                                + self.population[r3, n]
                                - self.population[r4, n]
                                - self.population[r5, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def _rand1bin(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 1, 0, 0)

        n = self.random_number_generator.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0
        while i < self.parameter_count:
            if self.random_number_generator.rand() < self.cross_over_probability or i == self.parameter_count - 1:
                trial[n] = self.population[r1, n]
                + self.scale * (self.population[r2, n]
                                - self.population[r3, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def select_samples(self, candidate, r1, r2, r3, r4, r5):
        if r1:
            while True:
                r1 = self.random_number_generator.randint(
                    0, self.population_size)
                if r1 != candidate:
                    break
        if r2:
            while True:
                r2 = self.random_number_generator.randint(
                    0, self.population_size)
                if r2 != candidate and r1 != r2:
                    break
        if r3:
            while True:
                r3 = self.random_number_generator.randint(
                    0, self.population_size)
                if r3 != candidate and r3 != r2 and r3 != r1:
                    break
        if r4:
            while True:
                r4 = self.random_number_generator.randint(
                    0, self.population_size)
                if r4 != candidate and r4 != r3 and r4 != r2 and r4 != r1:
                    break
        if r5:
            while True:
                r5 = self.random_number_generator.randint(
                    0, self.population_size)
                if r5 != candidate and r5 != r4 and r5 != r3 and r5 != r2 and r5 != r1:
                    break

        return r1, r2, r3, r4, r5


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
    from scipy.optimize import rosen
    # minimum expected at [1, 1, 1, 1, 1]
    bounds = [(0, 2), (0, 2), (0, 2), (0, 2), (0, 2)]
    result = differential_evolution(rosen,
                                    bounds,
                                    mutation=0.7,
                                    recomb=0.9,
                                    polish=True,
                                    disp=True)
    print (result)
