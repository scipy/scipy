# Simulated Dual Annealing implementation.
# Copyright (c) 2016 Sylvain Gubian <sylvain.gubian@pmi.com>,
# Yang Xiang <yang.xiang@pmi.com>
# Author: Sylvain Gubian, PMP S.A.
"""
SDA: A Simulated Dual Annealing global optimization algorithm
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.optimize import OptimizeResult
from scipy.optimize import minimize
from scipy.optimize._numdiff import approx_derivative
from scipy.special import gammaln
from scipy._lib._util import check_random_state

__all__ = ['sda']
BIG_VALUE = 1e16


class VisitingDistribution(object):
    """
    Class used to generate new coordinates based on the distorted
    Cauchy-Lorentz distribution. Depending on the steps within the Markov
    chain, the class implements the strategy for generating new location
    changes.
    """
    tail_limit = 1.e8
    min_visit_bound = 1.e-10

    def __init__(self, lb, ub, qv, rs):
        self.qv = qv
        self.rs = rs
        self.lower = lb
        self.upper = ub
        self.b_range = ub - lb
        self.x_gauss = None
        self.s_gauss = 0
        self.root_gauss = None

    def visiting(self, x, step, temperature):
        dim = x.size
        if step < dim:
            # Changing all coordinates with a new visting value
            visits = np.array([self.visit_fn(
                temperature) for _ in range(dim)])
            upper_sample = self.rs.random_sample()
            lower_sample = self.rs.random_sample()
            visits[visits > self.tail_limit] = self.tail_limit * upper_sample
            visits[visits < -self.tail_limit] = -self.tail_limit * lower_sample
            x_visit = visits + x
            a = x_visit - self.lower
            b = np.fmod(a, self.b_range) + self.b_range
            x_visit = np.fmod(b, self.b_range) + self.lower
            x_visit[np.fabs(
                x_visit - self.lower) < self.min_visit_bound] += 1.e-10
        else:
            # Changing only one coordinate at a time based on Markov chain step
            x_visit = np.copy(x)
            visit = self.visit_fn(temperature)
            if visit > self.tail_limit:
                visit = self.tail_limit * self.rs.random_sample()
            elif visit < -self.tail_limit:
                visit = -self.tail_limit * self.rs.random_sample()
            index = step - dim
            x_visit[index] = visit + x[index]
            a = x_visit[index] - self.lower[index]
            b = np.fmod(a, self.b_range[index]) + self.b_range[index]
            x_visit[index] = np.fmod(b, self.b_range[
                index]) + self.lower[index]
            if np.fabs(x_visit[index] - self.lower[
                    index]) < self.min_visit_bound:
                x_visit[index] += self.min_visit_bound
        return x_visit

    def visit_fn(self, temperature):
        factor1 = np.exp(np.log(temperature) / (self.qv - 1.0))
        factor2 = np.exp((4.0 - self.qv) * np.log(self.qv - 1.0))
        factor3 = np.exp((2.0 - self.qv) * np.log(2.0) / (self.qv - 1.0))
        factor4 = np.sqrt(np.pi) * factor1 * factor2 / (factor3 * (
            3.0 - self.qv))
        factor5 = 1.0 / (self.qv - 1.0) - 0.5
        d1 = 2.0 - factor5
        factor6 = np.pi * (1.0 - factor5) / np.sin(
            np.pi * (1.0 - factor5)) / np.exp(gammaln(d1))
        sigmax = np.exp(-(self.qv - 1.0) * np.log(
            factor6 / factor4) / (3.0 - self.qv))
        x = sigmax * self.gaussian_fn(1)
        y = self.gaussian_fn(0)
        den = np.exp(
            (self.qv - 1.0) * np.log((np.fabs(y))) / (3.0 - self.qv))
        return x / den

    def gaussian_fn(self, axis):
        if axis == 1:
            enter = True
            while enter or (self.s_gauss <= 0 or self.s_gauss >= 1):
                enter = False
                sample1 = self.rs.random_sample()
                self.x_gauss = sample1 * 2.0 - 1.0
                sample2 = self.rs.random_sample()
                y_gauss = sample2 * 2.0 - 1.0
                self.s_gauss = self.x_gauss ** 2 + y_gauss ** 2
            self.root_gauss = np.sqrt(-2.0 / self.s_gauss * np.log(
                self.s_gauss))
            return self.root_gauss * y_gauss
        else:
            return self.root_gauss * self.x_gauss


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

    def reset(self, owf, rs, x0=None):
        if x0 is None:
            self.current_location = self.lower + rs.random_sample(
                len(self.lower)) * (self.upper - self.lower)
        else:
            self.current_location = np.copy(x0)
        init_error = True
        reinit_counter = 0
        while init_error:
            self.current_energy = owf.func(self.current_location)
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
                self.current_location = self.lower + rs.random_sample(
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
    def __init__(self, qa, vd, ofw, rs, state):
        # Local markov chain minimum energy and location
        self.emin = state.current_energy
        self.xmin = np.array(state.current_location)
        # Global optimizer state
        self.state = state
        # Acceptance parameter
        self.qa = qa
        # Visiting distribution instance
        self.vd = vd
        # Wrapper to objective function and related local minimizer
        self.ofw = ofw
        self.not_improved_idx = 0
        self.not_improved_max_idx = 1000
        self._rs = rs
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
            x_visit = self.vd.visiting(
                self.state.current_location, j, temperature)
            # Calling the objective function
            e = self.ofw.func_wrapper(x_visit)
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
                r = self._rs.random_sample()
                pqa_temp = (self.qa - 1.0) * (
                    e - self.state.current_energy) / self.temperature_step + 1.
                if pqa_temp < 0.:
                    pqa = 0.
                else:
                    pqa = np.exp(np.log(pqa_temp) / (1. - self.qa))
                if r <= pqa:
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
            e, x = self.ofw.local_search(self.state.xbest)
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
            if pls >= self._rs.random_sample():
                do_ls = True
        # Global energy not improved, let's see what LS gives
        # on the best Markov chain location
        if self.not_improved_idx >= self.not_improved_max_idx:
            do_ls = True
        if do_ls:
            e, x = self.ofw.local_search(self.xmin)
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
    search and default gradient computation.
    Default local minimizer is L-BFGS-B
    """
    def __init__(self, bounds, func, **kwargs):
        self.func = func
        self.nb_fun_call = 0
        self.kwargs = kwargs
        self.minimizer = minimize
        self.fun_args = None
        lu = list(zip(*bounds))
        self.lower = np.array(lu[0])
        self.upper = np.array(lu[1])
        self.ls_max_iter = self.lower.size * 6

        if self.ls_max_iter < 100:
            self.ls_max_iter = 100
        if self.ls_max_iter > 1000:
            self.ls_max_iter = 1000

        # By default, scipy L-BFGS-B is used with a custom 3 points gradient
        # computation
        else:
            self.fun_args = ()
        if not self.kwargs or 'method' not in self.kwargs:
            self.kwargs['method'] = 'L-BFGS-B'
            self.kwargs['options'] = {
                'disp': None, 'maxls': 100, 'iprint': -1, 'gtol': 1e-06,
                'eps': 1e-06,
                'maxiter': self.ls_max_iter,
                'maxcor': 10, 'maxfun': 15000
            }
            if 'jac' not in self.kwargs:
                self.kwargs['jac'] = self.gradient
            if 'bounds' not in self.kwargs:
                self.kwargs['bounds'] = bounds
        if 'args' in self.kwargs:
            self.fun_args = self.kwargs.get('args')

    def func_wrapper(self, x):
        self.nb_fun_call += 1
        return self.func(x, *self.fun_args)

    def gradient(self, x):
        g = approx_derivative(self.func_wrapper, x)
        return g

    def local_search(self, x, maxlsiter=None):
        mres = self.minimizer(self.func_wrapper, x, **self.kwargs)
        if not mres.success:
            return BIG_VALUE, None
        return (mres.fun, mres.x)


class SDARunner(object):
    MAX_REINIT_COUNT = 1000

    def __init__(self, fun, x0, bounds, seed=None, local_search_options=None,
                 temperature_start=5230, qv=2.62, qa=-5.0,
                 maxfun=1e7, max_steps=500, no_local_search=False):
        if x0 is not None and not len(x0) == len(bounds):
            raise ValueError('Bounds size does not match x0')
        lu = list(zip(*bounds))
        lower = np.array(lu[0])
        upper = np.array(lu[1])
        # Checking that bounds are consistent
        if not np.all(lower < upper):
            raise ValueError('Bounds are note consistent min < max')
        # Wrapper for the objective function and minimizer
        if local_search_options is None:
            local_search_options = dict()
        self.owf = ObjectiveFunWrapper(bounds, fun, **local_search_options)
        # Initialization of RandomState for reproducible runs if seed provided
        self.rs = check_random_state(seed)
        # Initialization of the energy state
        self.es = EnergyState(lower, upper)
        self.es.reset(self.owf, self.rs, x0)
        # Maximum number of function call that can be used a stopping criterion
        self.maxfun = maxfun
        # Maximum number of step (main iteration)  that can be used as
        # stopping criterion
        self.max_steps = max_steps
        # Minimum value of annealing temperature reached to perform
        # re-annealing
        self.temperature_start = temperature_start
        self.temperature_restart = 0.1
        # VisitingDistribution instance
        vd = VisitingDistribution(lower, upper, qv, self.rs)
        # Markov chain instance
        self.mc = MarkovChain(qa, vd, self.owf, self.rs, self.es)
        self.qv = qv
        self.no_local_search = no_local_search

    def search(self):
        max_steps_reached = False
        self._iter = 0
        t1 = np.exp((self.qv - 1) * np.log(2.0)) - 1.0
        while(not max_steps_reached):
            for i in range(self.max_steps):
                # Compute temperature for this step
                s = float(i) + 2.0
                t2 = np.exp((self.qv - 1) * np.log(s)) - 1.0
                temperature = self.temperature_start * t1 / t2
                self._iter += 1
                if self._iter == self.max_steps:
                    max_steps_reached = True
                    break
                # Need a re-annealing process?
                if temperature < self.temperature_restart:
                    self.es.reset(self.owf, self.rs)
                    break
                # starting Markov chain
                self.mc.run(i, temperature)
                if self.owf.nb_fun_call >= self.maxfun:
                    break
                if not self.no_local_search:
                    self.mc.local_search()
                    if self.owf.nb_fun_call >= self.maxfun:
                        break

    @property
    def result(self):
        """ The OptimizeResult """
        res = OptimizeResult()
        res.x = self.es.xbest
        res.fun = self.es.ebest
        res.nit = self._iter
        res.ncall = self.owf.nb_fun_call
        return res


def sda(func, x0, bounds, maxiter=1000, local_search_options=None,
        initial_temp=5230., visit=2.62, accept=-5.0, maxfun=1e7, seed=None,
        no_local_search=False):
    """
    Find the global minimum of a function using the Simulated Dual Annealing
    algorithm

    Parameters
    ----------
    func : callable
        The objective function to be minimized.  Must be in the form
        ``f(x, *args)``, where ``x`` is the argument in the form of a 1-D array
        and ``args`` is a  tuple of any additional fixed parameters needed to
        completely specify the function.
    x0 : ndarray
        The starting coordinates. If ``None`` is provided, initial
        coordinates are automatically generated.
    bounds : sequence
        Bounds for variables.  ``(min, max)`` pairs for each element in ``x``,
        defining the lower and upper bounds for the optimizing argument of
        `func`. It is required to have ``len(bounds) == len(x)``.
        ``len(bounds)`` is used to determine the number of parameters in ``x``.
    maxiter : int, optional
        The maximum number of sda iterations. Increase this value if the
        objective function is very complicated with high dimensions.
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
        search of the energy landscape, allowing sda to escape local minima
        that it is trapped in.
    visit : float, optional
        Parameter for visiting distribution. Higher values give the visiting
        distribution a heavier tail, this makes the algorithm jump to a more
        distant region. The value range is (0, 3]
    accept : float, optional
        Parameter for acceptance distribution. It is used to control the
        probability of acceptance. The lower the acceptance parameter, the
        smaller the probability of acceptance. It has to be any negative value.
    maxfun : int, optional
        Soft limit for the number of objective function calls. If the
        algorithm is in the middle of a local search, this number will be
        exceeded, the algorithm will stop just after the local search is
        done.
    seed : int or `np.random.RandomState`, optional
        If `seed` is not specified the `np.RandomState` singleton is used.
        If `seed` is an int, a new `np.random.RandomState` instance is used,
        seeded with seed.
        If `seed` is already a `np.random.RandomState instance`, then that
        `np.random.RandomState` instance is used.
        Specify `seed` for repeatable minimizations. The random numbers
        generated with this seed only affect the visiting distribution
        function and new coordinates generation.
    no_local_search : boolean, optional
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
    strategy for applying a local search on accepted locations [4]_. A first
    implementation was done in C++ for the R language [5]_ and has been
    benchmarked [6]_. SDA introduces an advanced dual annealing method to
    refine the solution found by the generalized annealing process.
    This algorithm uses a distorted Cauchy-Lorentz visiting distribution, with
    its shape controlled by the parameter :math:`q_{v}`

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
    .. [1] Tsallis C (1988). "Possible generalization of Boltzmann-Gibbs
        statistics." Journal of Statistical Physics, 52, 479-487.
    .. [2] Tsallis C, Stariolo DA (1996). "Generalized Simulated Annealing."
        Physica A, 233, 395-406.
    .. [3] Xiang Y, Sun DY, Fan W, Gong XG (1997). "Generalized Simulated
        Annealing Algorithm and Its Application to the Thomson Model."
        Physics Letters A, 233, 216-220.
    .. [4] Xiang Y, Gong XG (2000a). "Efficiency of Generalized Simulated
        Annealing." PHYSICAL REVIEW E, 62, 4473.
    .. [5] Xiang Y, Gubian S, Suomela B, Hoeng (2013). "Generalized Simulated
        Annealing for Efficient Global Optimization: the GenSA Package for
        R". The R Journal, Volume 5/1, June 2013.
        http://journal.r-project.org/.
    .. [6] Mullen, K. (2014). Continuous Global Optimization in R. Journal of
        Statistical Software, 60(6), 1 - 45.
        http://dx.doi.org/10.18637/jss.v060.i06

    Examples
    --------
    The following example is a 10-dimensional problem, with many local minima.
    The function involved is called Rastrigin
    (https://en.wikipedia.org/wiki/Rastrigin_function)

    >>> from scipy.optimize import sda
    >>> func = lambda x: np.sum(x * x - 10 * np.cos(
    ...    2 * np.pi * x)) + 10 * np.size(x)
    >>> lw = [-5.12] * 10
    >>> up = [5.12] * 10
    >>> ret = sda(func, None, bounds=list(zip(lw, up)))
    >>> print("global minimum: xmin = {0}, f(xmin) = {1}".format(
    ...    ret.x, ret.fun))
    """
    gr = SDARunner(func, x0, bounds, seed, local_search_options,
                   temperature_start=initial_temp, qv=visit, qa=accept,
                   maxfun=maxfun, max_steps=maxiter,
                   no_local_search=no_local_search)
    gr.search()
    return gr.result
