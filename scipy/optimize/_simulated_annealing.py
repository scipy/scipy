# Simulated Dual Annealing implementation.
# Copyright (c) 2018 Sylvain Gubian <sylvain.gubian@pmi.com>,
# Yang Xiang <yang.xiang@pmi.com>
# Author: Sylvain Gubian, PMP S.A.
"""
SDA: A Simulated Dual Annealing global optimization algorithm
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.optimize import OptimizeResult
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy._lib._util import check_random_state

__all__ = ['simulated_annealing']
BIG_VALUE = 1e16


class VisitingDistribution(object):
    """
    Class used to generate new coordinates based on the distorted
    Cauchy-Lorentz distribution. Depending on the steps within the Markov
    chain, the class implements the strategy for generating new location
    changes.
    """
    TAIL_LIMIT = 1.e8
    MIN_VISIT_BOUND = 1.e-10

    def __init__(self, lb, ub, visiting_param, rand_state):
        self.visiting_param = visiting_param
        self.rand_state = rand_state
        self.lower = lb
        self.upper = ub
        self.bound_range = ub - lb
        self.x_gaussian_sample = None
        self.square_gaussian = 0
        self.square_root_gaussian = None

    def visiting(self, x, step, temperature):
        dim = x.size
        if step < dim:
            # Changing all coordinates with a new visting value
            visits = np.array([self.visit_fn(
                temperature) for _ in range(dim)])
            upper_sample = self.rand_state.random_sample()
            lower_sample = self.rand_state.random_sample()
            visits[visits > self.TAIL_LIMIT] = self.TAIL_LIMIT * upper_sample
            visits[visits < -self.TAIL_LIMIT] = -self.TAIL_LIMIT * lower_sample
            x_visit = visits + x
            a = x_visit - self.lower
            b = np.fmod(a, self.bound_range) + self.bound_range
            x_visit = np.fmod(b, self.bound_range) + self.lower
            x_visit[np.fabs(
                x_visit - self.lower) < self.MIN_VISIT_BOUND] += 1.e-10
        else:
            # Changing only one coordinate at a time based on Markov chain step
            x_visit = np.copy(x)
            visit = self.visit_fn(temperature)
            if visit > self.TAIL_LIMIT:
                visit = self.TAIL_LIMIT * self.rand_state.random_sample()
            elif visit < -self.TAIL_LIMIT:
                visit = -self.TAIL_LIMIT * self.rand_state.random_sample()
            index = step - dim
            x_visit[index] = visit + x[index]
            a = x_visit[index] - self.lower[index]
            b = np.fmod(a, self.bound_range[index]) + self.bound_range[index]
            x_visit[index] = np.fmod(b, self.bound_range[
                index]) + self.lower[index]
            if np.fabs(x_visit[index] - self.lower[
                    index]) < self.MIN_VISIT_BOUND:
                x_visit[index] += self.MIN_VISIT_BOUND
        return x_visit

    def visit_fn(self, temperature):
        # Formula Visita from p. 405 of reference [2]
        factor1 = np.exp(np.log(temperature) / (self.visiting_param - 1.0))
        factor2 = np.exp((4.0 - self.visiting_param) * np.log(
            self.visiting_param - 1.0))
        factor3 = np.exp((2.0 - self.visiting_param) * np.log(2.0) / (
            self.visiting_param - 1.0))
        factor4 = np.sqrt(np.pi) * factor1 * factor2 / (factor3 * (
            3.0 - self.visiting_param))
        factor5 = 1.0 / (self.visiting_param - 1.0) - 0.5
        d1 = 2.0 - factor5
        factor6 = np.pi * (1.0 - factor5) / np.sin(
            np.pi * (1.0 - factor5)) / np.exp(gammaln(d1))
        sigmax = np.exp(-(self.visiting_param - 1.0) * np.log(
            factor6 / factor4) / (3.0 - self.visiting_param))
        x = sigmax * self.gaussian_fn(1)
        y = self.gaussian_fn(0)
        den = np.exp(
            (self.visiting_param - 1.0) * np.log((np.fabs(y))) / (
                3.0 - self.visiting_param))
        return x / den

    def gaussian_fn(self, axis):
        if axis == 1:
            enter = True
            while enter or (self.square_gaussian <= 0 or
                            self.square_gaussian >= 1):
                enter = False
                sample1 = self.rand_state.random_sample()
                self.x_gaussian_sample = sample1 * 2.0 - 1.0
                sample2 = self.rand_state.random_sample()
                y_gaussian_sample = sample2 * 2.0 - 1.0
                self.square_gaussian = self.x_gaussian_sample ** 2 + \
                    y_gaussian_sample ** 2
            self.square_root_gaussian = np.sqrt(-2.0 / self.square_gaussian * 
                np.log(self.square_gaussian))
            return self.square_root_gaussian * y_gaussian_sample
        else:
            return self.square_root_gaussian * self.x_gaussian_sample


class EnergyState():
    """
    Class used to record the energy state. At any time, it knows what is the
    currently used coordinates and the most recent best location
    """
    # Maximimum number of trials for generating a valid starting point
    MAX_REINIT_COUNT = 1000

    def __init__(self, lower, upper):
        self.ebest = None
        self.current_energy = None
        self.current_location = None
        self.xbest = None
        self.lower = lower
        self.upper = upper

    def __str__(self):
        return 'Current: {0}@{1} Best: {2}@{3}'.format(
            self.current_energy, self.current_location, self._ebest,
            self.xbest,
        )

    def __repr__(self):
        return self.__str__

    def reset(self, obj_fun_wrapper, rand_state, x0=None):
        if x0 is None:
            self.current_location = self.lower + rand_state.random_sample(
                len(self.lower)) * (self.upper - self.lower)
        else:
            self.current_location = np.copy(x0)
        init_error = True
        reinit_counter = 0
        while init_error:
            self.current_energy = obj_fun_wrapper.func(self.current_location,
                *obj_fun_wrapper.fun_args)
            obj_fun_wrapper.nb_fun_call += 1
            if self.current_energy is None:
                raise ValueError('Objective function is returning None')
            if (self.current_energy >= BIG_VALUE or
                    np.isnan(self.current_energy)):
                if reinit_counter >= EnergyState.MAX_REINIT_COUNT:
                    init_error = False
                    message = (
                        'Stopping algorithm because function '
                        'create NaN or (+/-) inifinity values even with '
                        'trying new random parameters'
                    )
                    raise ValueError(message)
                self.current_location = self.lower + rand_state.random_sample(
                    self.lower.size) * (self.upper - self.lower)
                reinit_counter += 1
            else:
                init_error = False
            # If first time reset, initialize ebest and xbest
            if self.ebest is None and self.xbest is None:
                self.ebest = self.current_energy
                self.xbest = np.copy(self.current_location)
            # Otherwise, keep them in case of reannealing reset


class MarkovChain(object):
    """
    Class used for the Markov chain and related strategy for local search
    decision
    """
    def __init__(self, acceptance_param, visit_dist, obj_fun_wrapper,
                 rand_state, state):
        # Local markov chain minimum energy and location
        self.emin = state.current_energy
        self.xmin = np.array(state.current_location)
        # Global optimizer state
        self.state = state
        # Acceptance parameter
        self.acceptance_param = acceptance_param
        # Visiting distribution instance
        self.visit_dist = visit_dist
        # Wrapper to objective function and related local minimizer
        self.obj_fun_wrapper = obj_fun_wrapper
        self.not_improved_idx = 0
        self.not_improved_max_idx = 1000
        self._rand_state = rand_state
        self.temperature_step = 0
        self.K = 100 * len(state.current_location)

    def run(self, step, temperature):
        self.temperature_step = temperature / float(step + 1)
        self.not_improved_idx += 1
        for j in range(self.state.current_location.size * 2):
            if j == 0:
                self.state_improved = False
            if step == 0 and j == 0:
                self.state_improved = True
            x_visit = self.visit_dist.visiting(
                self.state.current_location, j, temperature)
            # Calling the objective function
            e = self.obj_fun_wrapper.func(x_visit,
                *self.obj_fun_wrapper.fun_args)
            self.obj_fun_wrapper.nb_fun_call += 1
            if e < self.state.current_energy:
                #  print('Better energy: {0}'.format(e))
                # We have got a better ernergy value
                self.state.current_energy = e
                self.state.current_location = np.copy(x_visit)
                if e < self.state.ebest:
                    self.state.ebest = e
                    self.state.xbest = np.copy(x_visit)
                    self.state_improved = True
                    self.not_improved_idx = 0
            else:
                # We have not improved but do we accept the new location?
                r = self._rand_state.random_sample()
                pqv_temp = (self.acceptance_param - 1.0) * (
                    e - self.state.current_energy) / self.temperature_step + 1.
                if pqv_temp < 0.:
                    pqv = 0.
                else:
                    pqv = np.exp(np.log(pqv_temp) / (
                        1. - self.acceptance_param))
                if r <= pqv:
                    # We accept the new location and update state
                    self.state.current_energy = e
                    self.state.current_location = np.copy(x_visit)
                    self.xmin = np.copy(self.state.current_location)

                # No improvement since long time
                if self.not_improved_idx >= self.not_improved_max_idx:
                    if j == 0 or self.state.current_energy < self.emin:
                        self.emin = self.state.current_energy
                        self.xmin = np.copy(self.state.current_location)
        # End of MarkovChain loop

    def local_search(self):
        # Decision making for performing a local search
        # based on Markov chain results
        # If energy has been improved or no improvement since too long,
        # performing a local search with the best Markov chain location
        if self.state_improved:
            # Global energy has improved, let's see if LS improved further
            e, x = self.obj_fun_wrapper.local_search(self.state.xbest)
            if e < self.state.ebest:
                self.not_improved_idx = 0
                self.state.ebest = e
                self.state.xbest = np.copy(x)
                self.state.current_energy = e
                self.state.current_location = np.copy(x)
                return
        # Check probability of a need to perform a LS even if no improvment
        # (Dual annealing principle)
        do_ls = False
        if self.K < 90 * len(self.state.current_location):
            pls = np.exp(self.K * (self.state.ebest - self.state.current_energy
                                   ) / self.temperature_step)
            if pls >= self._rand_state.random_sample():
                do_ls = True
        # Global energy not improved, let's see what LS gives
        # on the best Markov chain location
        if self.not_improved_idx >= self.not_improved_max_idx:
            do_ls = True
        if do_ls:
            e, x = self.obj_fun_wrapper.local_search(self.xmin)
            self.xmin = np.copy(x)
            self.emin = e
            self.not_improved_idx = 0
            self.not_improved_max_idx = self.state.current_location.size
            if e < self.state.ebest:
                self.state.ebest = self.emin
                self.state.xbest = np.copy(self.xmin)
                self.state.current_energy = e
                self.state.current_location = np.copy(x)


class ObjectiveFunWrapper(object):
    """
    Class used to wrap around the objective function in order to apply local
    search.
    Default local minimizer is SciPy minimizer L-BFGS-B
    """
    
    MAXITER_RATIO = 6
    MAXITER_MIN = 100
    MAXITER_MAX = 1000

    def __init__(self, bounds, func, **kwargs):
        self.func = func
        self.nb_fun_call = 0
        self.kwargs = kwargs
        self.minimizer = minimize
        self.fun_args = ()
        bounds_list = list(zip(*bounds))
        self.lower = np.array(bounds_list[0])
        self.upper = np.array(bounds_list[1])
        if 'args' in self.kwargs:
            self.fun_args = self.kwargs.get('args')
        
    def set_default_minimizer(self):
        # By default, use SciPy minimize with 'L-BFGS-B' method
        n = len(self.lower)
        ls_max_iter = min(max(n * self.MAXITER_RATIO, self.MAXITER_MIN),
                          self.MAXITER_MAX)
        self.kwargs['method'] = 'L-BFGS-B'
        self.kwargs['options'] = {
            'maxiter': ls_max_iter,
        }
        self.kwargs['bounds'] = list(zip(self.lower, self.upper))

    def local_search(self, x):
        x_tmp = np.copy(x)
        fun_temp = self.func(x, *self.fun_args)
        mres = self.minimizer(self.func, x, *self.fun_args, **self.kwargs)
        self.nb_fun_call += mres.nfev 
        # Check if is valid value
        is_finite = np.all(np.isfinite(mres.x)) and np.isfinite(mres.fun)
        in_bounds = np.all(mres.x >= self.lower) and np.all(
            mres.x <= self.upper)
        is_valid = is_finite and in_bounds

        # Use the new point only if it is valid and return a better results
        if is_valid and mres.fun < fun_temp:
            return mres.fun, mres.x
        else:
            return fun_temp, x_tmp 


def simulated_annealing(func, x0, bounds, maxiter=1000, local_search_options={},
        initial_temp=5230., visit=2.62, accept=-5.0, maxfun=1e7, seed=None,
        no_local_search=False):
    """
    Find the global minimum of a function using the Simulated Dual Annealing
    algorithm.

    Parameters
    ----------
    func : callable
        The objective function to be minimized.  Must be in the form
        ``f(x, *args)``, where ``x`` is the argument in the form of a 1-D array
        and ``args`` is a  tuple of any additional fixed parameters needed to
        completely specify the function.
    x0 : ndarray, shape(n,)
        A single initial starting point coordinates. If ``None`` is provided,
        initial coordinates are automatically generated.
    bounds : sequence, shape (n, 2) 
        Bounds for variables.  ``(min, max)`` pairs for each element in ``x``,
        defining bounds for the objective function parameter.
    maxiter : int, optional
        The maximum number of global search iterations. Default value is 1000.
    local_search_options : dict, optional
        Extra keyword arguments to be passed to the local minimizer
            ``scipy.optimize.minimize()`` Some important options could be:
            method : str
                The minimization method (e.g. ``"L-BFGS-B"``)
            args : tuple
                Extra arguments passed to the objective function (``func``) and
                its derivatives (Jacobian, Hessian).
    initial_temp : float, optional
        The initial temperature, use higher values to facilitates a wider
        search of the energy landscape, allowing simulated_annealing to escape
        local minima that it is trapped in. Default value is 5230.
    visit : float, optional
        Parameter for visiting distribution. Higher values give the visiting
        distribution a heavier tail, this makes the algorithm jump to a more
        distant region. The value range is (0, 3]. Default value is 2.62. 
    accept : float, optional
        Parameter for acceptance distribution. It is used to control the
        probability of acceptance. The lower the acceptance parameter, the
        smaller the probability of acceptance. Default value is -5.0.
    maxfun : int, optional
        Soft limit for the number of objective function calls. If the
        algorithm is in the middle of a local search, this number will be
        exceeded, the algorithm will stop just after the local search is
        done.
    seed : {int, `np.random.RandomState`}, optional
        If `seed` is not specified the `np.RandomState` singleton is used.
        If `seed` is an int, a new `np.random.RandomState` instance is used,
        seeded with seed.
        If `seed` is already a `np.random.RandomState instance`, then that
        `np.random.RandomState` instance is used.
        Specify `seed` for repeatable minimizations. The random numbers
        generated with this seed only affect the visiting distribution
        function and new coordinates generation.
    no_local_search : bool, optional
        If `no_local_search` is set to `True`, a traditional Generalized
        Simulated Annealing will be performed with no local search
        strategy applied.

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``fun`` the value
        of the function at the solution, and ``message`` which describes the
        cause of the termination.
        See `OptimizeResult` for a description of other attributes.

    Notes
    -----
    SDA is an implementation of the Simulated Dual Annealing. This stochastic
    approach [3]_ implements the generalization of CSA (Classical Simulated
    Annealing) and FSA (Fast Simulated Annealing) [1]_ [2]_ coupled to a
    strategy for applying a local search on accepted locations [4]_.
    An alternative implementation of this same algorithm is described in [5]_
    and benchmarks are presented in [6]_. SDA introduces an advanced dual
    annealing method to refine the solution found by the generalized annealing
    process. This algorithm uses a distorted Cauchy-Lorentz visiting
    distribution, with its shape controlled by the parameter :math:`q_{v}`

    .. math::

        g_{q_{v}}(\\Delta x(t)) \\propto \\frac{ \\
        \\left[T_{q_{v}}(t) \\right]^{-\\frac{D}{3-q_{v}}}}{ \\
        \\left[{1+(q_{v}-1)\\frac{(\Delta x(t))^{2}} { \\
        \\left[T_{q_{v}}(t)\\right]^{\\frac{2}{3-q_{v}}}}}\\right]^{ \\
        \\frac{1}{q_{v}-1}+\\frac{D-1}{2}}}

    Where :math:`t` is the artificial time. This visiting distribution is used
    to generate a trial jump distance :math:`\Delta x(t)` of variable
    :math:`x(t)` under artificial temperature :math:`T_{q_{v}}(t)`.

    From the starting point, after calling the visiting distribution
    function, the acceptance probability is computed as follows:

    .. math::

        p_{q_{a}} = \min{\{1,\\left[1-(1-q_{a}) \\beta \\Delta E \\right]^{ \\
        \\frac{1}{1-q_{a}}}\\}}

    Where :math:`q_{a}` is a acceptance parameter. For :math:`q_{a}<1`, zero
    acceptance probability is assigned to the cases where

    .. math::

        [1-(1-q_{a}) \\beta \\Delta E] < 0

    The artificial temperature :math:`T_{q_{v}}(t)` is decreased according to

    .. math::

        T_{q_{v}}(t) = T_{q_{v}}(1) \\frac{2^{q_{v}-1}-1}{\\left( \\
        1 + t\\right)^{q_{v}-1}-1}

    Where :math:`q_{v}` is the visiting parameter.

    .. versionadded:: 1.1.0

    References
    ----------
    .. [1] Tsallis C. Possible generalization of Boltzmann-Gibbs
        statistics. Journal of Statistical Physics, 52, 479-487 (1998).
    .. [2] Tsallis C, Stariolo DA. Generalized Simulated Annealing.
        Physica A, 233, 395-406 (1996).
    .. [3] Xiang Y, Sun DY, Fan W, Gong XG. Generalized Simulated
        Annealing Algorithm and Its Application to the Thomson Model.
        Physics Letters A, 233, 216-220 (1997).
    .. [4] Xiang Y, Gong XG. Efficiency of Generalized Simulated
        Annealing. Physical Review E, 62, 4473 (2000).
    .. [5] Xiang Y, Gubian S, Suomela B, Hoeng J. Generalized
        Simulated Annealing for Efficient Global Optimization: the GenSA
        Package for R. The R Journal, Volume 5/1 (2013).
    .. [6] Mullen, K. Continuous Global Optimization in R. Journal of
        Statistical Software, 60(6), 1 - 45, (2014). DOI:10.18637/jss.v060.i06

    Examples
    --------
    The following example is a 10-dimensional problem, with many local minima.
    The function involved is called Rastrigin
    (https://en.wikipedia.org/wiki/Rastrigin_function)

    >>> from scipy.optimize import simulated_annealing
    >>> func = lambda x: np.sum(x * x - 10 * np.cos(
    ...    2 * np.pi * x)) + 10 * np.size(x)
    >>> lw = [-5.12] * 10
    >>> up = [5.12] * 10
    >>> ret = simulated_annealing(func, None, bounds=list(zip(lw, up)))
    >>> print("global minimum: xmin = {0}, f(xmin) = {1:.6f}".format(
    ...    ret.x, ret.fun))

    global minimum: xmin = [ -5.30926309e-09  -9.47607022e-09  -4.09044159e-09
    -8.91916554e-09 -7.58266242e-09  -7.51256718e-09  -5.29363641e-09
    -5.77883504e-09 -4.36341509e-09  -5.22583018e-09], f(xmin) = 0.000000
    """

    if x0 is not None and not len(x0) == len(bounds):
        raise ValueError('Bounds size does not match x0')
    lu = list(zip(*bounds))
    lower = np.array(lu[0])
    upper = np.array(lu[1])
    # Checking that bounds are consistent
    if not np.all(lower < upper):
        raise ValueError('Bounds are note consistent min < max')
    # Wrapper for the objective function and minimizer
    obj_fun_wrapper = ObjectiveFunWrapper(bounds, func, **local_search_options)
    if not local_search_options:
        obj_fun_wrapper.set_default_minimizer()
    # Initialization of RandomState for reproducible runs if seed provided
    rand_state = check_random_state(seed)
    # Initialization of the energy state
    energy_state = EnergyState(lower, upper)
    energy_state.reset(obj_fun_wrapper, rand_state, x0)
    # Minimum value of annealing temperature reached to perform
    # re-annealing
    temperature_restart = 0.1
    # VisitingDistribution instance
    visit_dist = VisitingDistribution(lower, upper, visit, rand_state)
    # Markov chain instance
    markov_chain = MarkovChain(accept, visit_dist, obj_fun_wrapper, rand_state,
                               energy_state)

    # Run the search lopp 
    max_steps_reached = False
    iteration = 0
    t1 = np.exp((visit - 1) * np.log(2.0)) - 1.0
    while(not max_steps_reached):
        for i in range(maxiter):
            # Compute temperature for this step
            s = float(i) + 2.0
            t2 = np.exp((visit - 1) * np.log(s)) - 1.0
            temperature = initial_temp * t1 / t2
            iteration += 1
            if iteration == maxiter:
                max_steps_reached = True
                break
            # Need a re-annealing process?
            if temperature < temperature_restart:
                energy_state.reset(obj_fun_wrapper, rand_state)
                break
            # starting Markov chain
            markov_chain.run(i, temperature)
            if obj_fun_wrapper.nb_fun_call >= maxfun:
                break
            if not no_local_search:
                markov_chain.local_search()
                if obj_fun_wrapper.nb_fun_call >= maxfun:
                    break

    # Return the OptimizeResult
    res = OptimizeResult()
    res.x = energy_state.xbest
    res.fun = energy_state.ebest
    res.nit = iteration
    res.ncall = obj_fun_wrapper.nb_fun_call
    return res
