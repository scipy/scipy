from __future__ import division
import numpy as np
import numpy.random as npr
import scipy.optimize

__all__ = ['differential_evolution']

# standard status messages of optimizers
_status_message = {'success': 'Optimization terminated successfully.',
                   'maxfev': 'Maximum number of function evaluations has '
                              'been exceeded.',
                   'maxiter': 'Maximum number of iterations has been '
                              'exceeded.',
                   'pr_loss': 'Desired error not necessarily achieved due '
                              'to precision loss.',
                   'aborted': 'Minimization aborted by callback function'}


class BoundsException(Exception):

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


def differential_evolution(func, bounds, args=(), DEstrategy=None,
                           maxiter=None, popsize=20, tol=0.01,
                           km=0.7, recomb=0.5, seed=None, callback=None,
                           disp=False, **options):
    """
    Find the global minimum of a multivariate function using the differential
    evolution algorithm. It is stochastic in nature (does not use gradient
    methods) to find the minimium, and can search large areas of candidate
    space, but often requires large numbers of function evaluations.

    The algorithm is originally due to Storn and Price [1]_.

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

            - 'Best1Bin'
            - 'Best1Exp'
            - 'Rand1Exp'
            - 'RandToBest1Exp'
            - 'Best2Exp'
            - 'Rand2Exp'
            - 'Best2Bin'
            - 'Rand2Bin'
            - 'Rand1Bin'

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
        i.e. mean(pop) * tol / stdev(pop) > 1
    km : float, optional:
        The mutation constant, should be in the range [0, 1].
    recomb : float, optional:
        The recombination constant, should be in the range [0, 1].
    seed : float, optional:
        Seeds the random number generator for repeatable minimizations.
    disp : bool, optional:
        Display status messages
    callback : callable, optional:
        A function to follow the progress of the minimization.
        Called as ``callback(xk, convergence=val)``, where ``xk`` is
        the current value of ``x0``. ``val`` represents the fractional
        value of the population convergence.  When ``val`` is greater than one
        the function halts.
        If this function returns False, then the minimization is halted.

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes.

    References
    ----------
    .. [1] http://www1.icsi.berkeley.edu/~storn/code.html
    .. [2] http://en.wikipedia.org/wiki/Differential_evolution

    Examples
    --------
    Let us consider the problem of minimizing the Rosenbrock function. This
    function is implemented in `rosen` in `scipy.optimize`.

    >>> from scipy.optimize import rosen, differential_evolution
    >>> func = lambda x: np.cos(14.5 * x - 0.3) + (x + 0.2) * x
    >>> bounds = [(-3, 3)]
    >>> result = differential_evolution(func, bounds,
    ... tol=1e-2, popsize=40, km=0.6, recomb=0.9, DEstrategy='Best1Bin')
    >>> print result

    """

    # assemble the bounds into the limits
    try:
        limits = np.array(bounds, float).T
        assert np.size(limits, 0) == 2
    except (ValueError, AssertionError) as e:
        # it is required to have (min, max) pairs for each value in x
        raise BoundsException('Bounds should be a sequence containing '
                              'real valued (min, max) pairs for each value'
                              ' in x')

    solver = DEsolver(func, limits, args=args, DEstrategy=DEstrategy,
                      maxiter=maxiter, popsize=popsize,
                      tol=tol, km=km, recomb=recomb, seed=seed,
                      callback=callback)

    result = solver.solve()
    return result


class DEsolver(object):

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
    km : float, optional:
        The mutation constant, should be in the range [0, 1].
    recomb : float, optional:
        The recombination constant, should be in the range [0, 1].
    seed : int, optional:
        Seed initializing the pseudo-random number generator. Can be an integer,
        an array (or other sequence) of integers of any length, or None (the
        default). If you use the seed you will get repeatable minimizations.
    callback : callable, optional:
        A function to follow the progress of the minimization.
        Called as ``callback(xk)``, where ``xk`` is the current value of ``x0``.
    disp : bool, optional
    """

    def __init__(self, func, limits, args=(),
                 DEstrategy=None, maxiter=None, popsize=20,
                 tol=0.01, km=0.7, recomb=0.5, seed=None, maxfun=None,
                 callback=None, disp=False, **options):

        if DEstrategy is not None:
            self.DEstrategy = getattr(DEsolver, DEstrategy)
        else:
            self.DEstrategy = getattr(DEsolver, 'Best1Bin')

        self.callback = callback

        self.maxiter = 1000
        if maxiter is not None:
            self.maxiter = maxiter

        self.maxfun = (self.maxiter + 1) * popsize * np.size(limits, 1)
        if maxfun is not None:
            self.maxfun = maxfun

        self.tol = tol
        self.scale = km
        self.crossOverProbability = recomb

        self.func = func
        self.args = args
        self.limits = limits
        if np.any(np.isnan(limits)):
            raise BoundsException('Bounds should be a sequence'
                                  ' containing real valued '
                                  '(min, max) pairs for each value in x')

        self.nfev = 0
        self.nit = 0

        self.parameter_count = np.size(self.limits, 1)
        self.population_size = popsize * self.parameter_count

        self.RNG = npr.RandomState()
        self.RNG.seed(seed)

        self.population = self.RNG.rand(
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
            `OptimizeResult` for a description of other attributes.
        """
        status_message = _status_message['success']
        warning_flag = False

        # calculate energies to start with
        for index, candidate in enumerate(self.population):
            params = self.__scale_parameters(candidate)
            self.population_energies[
                index] = self.func(
                params,
                *self.args)
            self.nfev += 1

            if self.nfev == self.maxfun:
                warning_flag = True
                status_message = status_message['maxfev']

        minval = np.argmin(self.population_energies)

        # put the lowest energy into the best solution position.
        lowest_energy = self.population_energies[minval]
        self.population_energies[minval] = self.population_energies[0]
        self.population_energies[0] = lowest_energy

        self.population[[0, minval], :] = self.population[[minval, 0], :]

        # do the optimisation.
        for iteration in xrange(self.maxiter):
            for candidate in xrange(self.population_size):
                trial = self.DEstrategy(self, candidate)
                self.__ensure_constraint(trial)
                params = self.__scale_parameters(trial)

                energy = self.func(params, *self.args)
                self.nfev += 1
                if self.nfev == self.maxfun:
                    warning_flag == True
                    status_message = _status_message['maxfev']
                    break

                if energy < self.population_energies[candidate]:
                    self.population[candidate] = trial
                    self.population_energies[candidate] = energy

                    if energy < self.population_energies[0]:
                        self.population_energies[0] = energy
                        self.population[0] = trial
                        if self.disp:
                            print("differential_evolution step %d: f %g" \
                                 % (self.iteration,
                                    self.population_energies[0]))

            # stop when the fractional s.d. of the population is less than tol
            # of the mean energy
            self.convergence = np.std(self.population_energies) / \
                np.mean(self.population_energies)

            if self.callback:
                should_continue = self.callback(
                    self.__scale_parameters(self.population[0]),
                    convergence=self.tol / self.convergence)
                if should_continue is False:
                    warning_flag = True
                    status_message = _status_message['aborted']
                    break

            if self.convergence < self.tol:
                break

            if warning_flag:
                break

            self.nit += 1

        if self.nit == self.maxiter:
            status_message = _status_message['maxiter']
            warning_flag = True

        result = scipy.optimize.OptimizeResult(
            x=self.__scale_parameters(self.population[0]),
            fun=self.population_energies[0],
            nfev=self.nfev,
            nit=self.nit,
            message=status_message,
            success=(warning_flag != True))

        return result

    def __scale_parameters(self, trial):
        return (
            0.5 * (self.limits[0] + self.limits[1]) +
            (trial - 0.5) * np.fabs(self.limits[0] - self.limits[1])
        )

    def __ensure_constraint(self, trial):
        for index, param in enumerate(trial):
            if param > 1 or param < 0:
                trial[index] = self.RNG.rand()

    def Best1Bin(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 0, 0, 0)
        n = self.RNG.randint(0, self.parameter_count)
        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count:
            if self.RNG.rand() < self.crossOverProbability or i == self.parameter_count - 1:
                trial[n] = self.population[0, n] + self.scale * \
                    (self.population[r1, n] - self.population[r2, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def Best1Exp(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 0, 0, 0)

        n = self.RNG.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0
        while i < self.parameter_count and self.RNG.rand() < self.crossOverProbability:
            trial[n] = self.population[0, n] + self.scale * \
                (self.population[r1, n] - self.population[r2, n])
            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def Rand1Exp(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 1, 0, 0)

        n = self.RNG.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count and self.RNG.rand() < self.crossOverProbability:
            trial[n] = self.population[r1, n] + self.scale * \
                (self.population[r2, n] - self.population[r3, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def RandToBest1Exp(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 0, 0, 0)

        n = self.RNG.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count and self.RNG.rand() < self.crossOverProbability:
            trial[n] += self.scale * (self.population[0, n] - trial[n]) + \
                self.scale * \
                (self.population[r1, n]
                 - self.population[r2, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def Best2Exp(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 1, 1, 0)

        n = self.RNG.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count and self.RNG.rand() < self.crossOverProbability:
            trial[n] = self.population[0, n]
            + self.scale * (self.population[r1, n]
                            + self.population[r2, n]
                            - self.population[r3, n]
                            - self.population[r4, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def Rand2Exp(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 1, 1, 1)

        n = self.RNG.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count and self.RNG.rand() < self.crossOverProbability:
            trial[n] = self.population[r1, n]
            + self.scale * (self.population[r2, n]
                            + self.population[r3, n]
                            - self.population[r4, n]
                            - self.population[r5, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def RandToBest1Bin(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 0, 0, 0)

        n = self.RNG.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count:
            if self.RNG.rand() < self.crossOverProbability or i == self.parameter_count - 1:
                trial[n] += self.scale * (self.population[0, n] - trial[n])
                + self.scale * \
                    (self.population[r1, n] - self.population[r2, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def Best2Bin(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 1, 1, 0)

        n = self.RNG.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count:
            if self.RNG.rand() < self.crossOverProbability or i == self.parameter_count - 1:
                trial[n] = self.population[0, n]
                + self.scale * (self.population[r1, n]
                                + self.population[r2, n]
                                - self.population[r3, n]
                                - self.population[r4, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def Rand2Bin(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 1, 1, 1)

        n = self.RNG.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0

        while i < self.parameter_count:
            if self.RNG.rand() < self.crossOverProbability or i == self.parameter_count - 1:
                trial[n] = self.population[r1, n]
                + self.scale * (self.population[r2, n]
                                + self.population[r3, n]
                                - self.population[r4, n]
                                - self.population[r5, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def Rand1Bin(self, candidate):
        r1, r2, r3, r4, r5 = self.select_samples(candidate, 1, 1, 1, 0, 0)

        n = self.RNG.randint(0, self.parameter_count)

        trial = np.copy(self.population[candidate])
        i = 0
        while i < self.parameter_count:
            if self.RNG.rand() < self.crossOverProbability or i == self.parameter_count - 1:
                trial[n] = self.population[r1, n]
                + self.scale * (self.population[r2, n]
                                - self.population[r3, n])

            n = (n + 1) % self.parameter_count
            i += 1

        return trial

    def select_samples(self, candidate, r1, r2, r3, r4, r5):
        if r1:
            while True:
                r1 = self.RNG.randint(0, self.population_size)
                if r1 != candidate:
                    break
        if r2:
            while True:
                r2 = self.RNG.randint(0, self.population_size)
                if r2 != candidate and r1 != r2:
                    break
        if r3:
            while True:
                r3 = self.RNG.randint(0, self.population_size)
                if r3 != candidate and r3 != r2 and r3 != r1:
                    break
        if r4:
            while True:
                r4 = self.RNG.randint(0, self.population_size)
                if r4 != candidate and r4 != r3 and r4 != r2 and r4 != r1:
                    break
        if r5:
            while True:
                r5 = self.RNG.randint(0, self.population_size)
                if r5 != candidate and r5 != r4 and r5 != r3 and r5 != r2 and r5 != r1:
                    break

        return r1, r2, r3, r4, r5

if __name__ == "__main__":
    # minimum expected at ~-0.195
    func = lambda x: np.cos(14.5 * x - 0.3) + (x + 0.2) * x
    bounds = [(-3, 3)]
    result = differential_evolution(func,
                                    bounds,
                                    tol=1e-2,
                                    popsize=40,
                                    km=0.6,
                                    recomb=0.9,
                                    DEstrategy='Best1Bin')
    print result
