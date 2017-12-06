# Generalized simulated annealing implementation.
# Copyright (c) 2016 Sylvain Gubian <sylvain.gubian@pmi.com>,
# Yang Xiang <yang.xiang@pmi.com>
# Author: Sylvain Gubian, PMP S.A.
"""
gensa: A generalized simulated annealing global optimization algorithm
"""
from __future__ import division, print_function, absolute_import

import numpy as np
import time
from scipy.optimize import OptimizeResult
from scipy.optimize import _lbfgsb
from scipy.special import gammaln
from scipy._lib._util import check_random_state

__all__ = ['gensa']


class GenSARunner(object):
    """This class implements the core of the gensa algorithm.

    fun : callable
        The objective function
    x0 : ndarray
        The starting coordinates.
    bounds : sequence
        Bounds for variables.  ``(min, max)`` pairs for each element in ``x``,
        defining the lower and upper bounds for the optimizing argument of
        `func`. It is required to have ``len(bounds) == len(x)``.
        ``len(bounds)`` is used to determine the number of parameters in ``x``.
    args : tuple, optional
        Any additional fixed parameters needed to
        completely specify the objective function.
    seed : int or `np.random.RandomState`, optional
        If `seed` is not specified the `np.RandomState` singleton is used.
        If `seed` is an int, a new `np.random.RandomState` instance is used,
        seeded with seed.
        If `seed` is already a `np.random.RandomState instance`, then that
        `np.random.RandomState` instance is used.
        Specify `seed` for repeatable minimizations. The random numbers
        generated with this seed only affect the visiting distribution
        function and new coordinates generation.
    temp_start : float, optional
        The initial temperature, use higher values to facilitates a wider
        search of the energy landscape, allowing gensa to escape local minima
        that it is trapped in.
    qv : float, optional
        Parameter for visiting distribution. Higher values give the visiting
        distribution a heavier tail, this makes the algorithm jump to a more
        distant region. The value range is (0, 3]
    qa : float, optional
        Parameter for acceptance distribution. It is used to control the
        probability of acceptance. The lower the acceptance parameter, the
        smaller the probability of acceptance. It has to be any negative value.
    maxfun : int, optional
        Soft limit for the number of objective function calls. If the
        algorithm is in the middle of a local search, this number will be
        exceeded, the algorithm will stop just after the local search is
        done.
    maxsteps : int, optional
        The maximum number of gensa iterations will perform.
    """

    KSPRING = 1.e8
    BIG_VALUE = 1.e13
    MAX_REINIT_COUNT = 100

    def __init__(self, fun, x0, bounds, args=(), seed=None,
            temperature_start=5230, qv=2.62, qa=-5.0, maxfun=1e7, maxsteps=500,
            pure_sa=False):
        self.fun = fun
        self.args = args
        self.pure_sa = pure_sa
        if x0 is not None and not len(x0) == len(bounds):
            raise ValueError('Bounds size does not match x0')
        lu = list(zip(*bounds))
        self._lower = np.array(lu[0])
        self._upper = np.array(lu[1])
        # Checking that bounds are consistent
        if not np.all(self._lower < self._upper):
            raise ValueError('Bounds are note consistent min < max')
        # Initialization of RandomState for reproducible runs if seed provided
        self._random_state = check_random_state(seed)
        if x0 is None:
            x0 = self._lower + self._random_state.random_sample(
                    len(bounds)) * (self._upper - self._lower)
        self.seed = seed
        # Number of maximum sof iteration for local search
        self.itsoftmax = x0.size * 6
        # Size of the markov chain. Twice the dimension problem is a
        # recommended value.
        self.markov_length = x0.size * 2
        # In case the real value of the global minimum is known
        # it can be used as stopping criterion
        self.know_real = False
        self.real_threshold = -np.inf
        # Maximum duration time of execution as a stopping criterion
        # Default is unlimited duration
        self.maxtime = np.inf
        # Maximum number of function call that can be used a stopping criterion
        self.maxfuncall = maxfun
        # Maximum number of step (main iteration)  that ca be used as
        # stopping criterion
        self.maxsteps = maxsteps
        # Minimum value of annealing temperature reached to perform
        # re-annealing
        self.temperature_restart = 0.1
        # Visiting distribution parameter
        self.qv = qv
        # Acceptance parameter value
        self.qa = qa
        # Initial temperature value for annealing
        self.temperature_start = temperature_start
        # Not yet implemented contraint function that would be used in the
        # future
        self.has_constraint = False
        self.judge_constraint = None
        self.factr = 1000
        self.pgtol = 1.e-6
        self.reps = 1.e-6

        self._x = np.array(x0)
        self._xrange = self._upper - self._lower
        self._xbackup = np.array(self._x)
        self._xmin = np.array(self._x)
        self._xbuffer = np.zeros(self._x.size)
        self._temperature = self.temperature_start
        self._usey = 1
        self._xgas = 0.0
        self._ygas = 0.0
        self._ranbyx = 0.0
        self._ranbyy = 0.0
        self._sgas = -1.0
        self._step_record = 0
        self._nbfuncall = 0
        self._emin_unchanged = True
        self._index_no_emin_update = 0
        self._temperature_qa = 0
        self._initialize()

    def _initialize(self):
        """
        Random coordinate generation in case given initial coordinates
        are giving invalid objective value
        """
        self._nbfuncall = 0
        if self.markov_length % self._x.size != 0:
            raise ValueError('Incorrect markov length.')
        in_constraint = True
        init_error = True
        reinit_counter = 0
        while(init_error):
            self._energy(self._x)
            if self._etot >= self.BIG_VALUE:
                if reinit_counter >= self.MAX_REINIT_COUNT:
                    init_error = False
                    self._message = [(
                        'Stopping algorithm because function '
                        'create NaN or (+/-) inifinity values even with '
                        'trying new random parameters')]
                    raise ValueError(self._message[0])
                self._x = self._lower + self._random_state.random_sample(
                        self._x.size) * (self._upper - self._lower)
                reinit_counter += 1
            else:
                init_error = False

    def _visiting(self, step):
        """
        Assignement of components values based on visiting distribution.
        The way of exploring space depends on the Markov chain stepping
        """
        # It it is the first part of the markov chain
        # Changing all components at the same time
        if step < self._x.size:
            visits = np.array([self._visita() for _ in range(self._x.size)])
            visits[visits > 1.e8] = 1.e8 * self._random_state.random_sample()
            visits[visits < -1e8] = -1.e8 * self._random_state.random_sample()
            self._x = visits + self._xbackup
            a = self._x - self._lower
            b = np.fmod(a, self._xrange) + self._xrange
            self._x = np.fmod(b, self._xrange) + self._lower
            self._x[np.fabs(self._x - self._lower) < 1.e-10] += 1.e-10
        else:
            # Second part of the markov chain
            # Now change only one component at a time
            visit = self._visita()
            if visit > 1.e8:
                visit = 1.e8 * self._random_state.random_sample()
            elif visit < -1e8:
                visit = -1.e8 * self._random_state.random_sample()
            index = step - self._x.size
            self._x[index] = visit + self._xbackup[index]
            a = self._x[index] - self._lower[index]
            b = np.fmod(
                a, self._xrange[index]) + self._xrange[index]
            self._x[index] = np.fmod(
                b, self._xrange[index]) + self._lower[index]
            if np.fabs(self._x[index] - self._lower[
                    index]) < 1.e-10:
                self._x[index] += 1.e-10

    def _run_markov_chain(self, step):
        for j in range(self.markov_length):
            if j == 0:
                self._emin_unchanged = True
            if step == 0 and j == 0:
                self._emin_unchanged = False
            # Keeping old location before visiting new place
            self._xbackup = np.array(self._x)
            self._visiting(j)
            self._energy(self._x)
            if self._etot < self._etot0:
                # We get a better energy value
                self._etot0 = np.array(self._etot)
                if self._etot < self._emin:
                    self._emin = np.array(self._etot)
                    self._xmin = np.array(self._x)
                    self._emin_unchanged = False
                    self._index_no_emin_update = 0
            else:
                # We do not have improvement but do we accept location
                r = self._random_state.random_sample()
                pqa1 = (self.qa - 1.0) * (
                    self._etot - self._etot0) / self._temperature_qa + 1.0
                if pqa1 < 0.0:
                    pqa = 0.0
                else:
                    pqa = np.exp(np.log(pqa1) / (1.0 - self.qa))
                if r > pqa:
                    # We reject the new visiting location
                    self._x = self._xbackup
                else:
                    # The new visiting location is accepted
                    self._etot0 = self._etot
            if self._check_stopping():
                self._stop_search()
                return 0
            if self._index_no_emin_update >= self._index_tol_emin_update - 1:
                if j == 0:
                    self._emin_markov = np.array(self._etot0)
                    self._xmin_markov = np.array(self._x)
                else:
                    if self._etot0 < self._emini_markov:
                        self._emin_markov = np.array(self._etot0)
                        self._xmin_markov = np.array(self._x)

    def start_search(self):
        """ Start annealing process with eventually re-annealing
        """
        in_constraint = True
        self._emin_unchanged = True
        self._emin_markov = 0.0
        self._xmin_markov = np.zeros(self._x.size)
        self._index_no_emin_update = 0
        self._index_tol_emin_update = 1000
        self._starttime = time.time()
        self._emin = np.array(self._etot)
        self._xmin = np.array(self._x)
        self._etot0 = np.array(self._etot)
        if self._etot < self._emin:
            self._emin = np.array(self._etot)
            self._xmin = np.array(self._x)
        self._etot0 = self._etot
        if self._check_stopping():
            self._stop_search()
        self._step_record = 0
        self._temperature = self.temperature_start
        max_steps_not_exceeded = True
        while(max_steps_not_exceeded):
            for i in range(self.maxsteps):
                # Evaluating iteration artificial temperature
                s = float(i) + 2.0
                t1 = np.exp((self.qv - 1) * np.log(2.0)) - 1.0
                t2 = np.exp((self.qv - 1) * np.log(s)) - 1.0
                self._temperature = self.temperature_start * t1 / t2
                self._step_record += 1
                self._temperature_qa = self._temperature / float(i + 1)
                self._index_no_emin_update += 1
                # break out of both for-loop and while loop because, annealing
                # and eventually re-anneling reached maximum number of
                # iteration set.
                if self._step_record == self.maxsteps:
                    max_steps_not_exceeded = False
                    break
                # Need a re-annealing process? - Restarting main loop
                if self._temperature < self.temperature_restart:
                    break
                # Starting Markov chain
                self._run_markov_chain(i)

                # Decision making for performing a local search
                # based on the markov chain results
                if not self._emin_unchanged and not self.pure_sa:
                    temp = np.array(self._xmin)
                    etemp = self._ls_energy(temp)
                    temp = self._xbuffer
                    if etemp < self._emin:
                        self._xmin = np.array(temp)
                        self._emin = np.array(etemp)
                        self._index_no_emin_update = 0
                if self._index_no_emin_update >= (
                        self._index_tol_emin_update - 1) and not self.pure_sa:
                    self._emin_markov = np.array(self._ls_energy(
                        self._xmin_markov))
                    self._index_no_emin_update = 0
                    self._index_tol_emin_update = self._x.size
                    if self._emin_markov < self._emin:
                        self._xmin = np.array(self._xmin_markov)
                        self._emin = np.array(self._emin_markov)
                if self._check_stopping():
                    self._stop_search()
                    return 0
        self._stop_search()
        self._message = ["Number of iteration reached"]
        return 0

    def _check_stopping(self):
        """Check if the search has to be stopped
        """
        if self.know_real:
            if self._emin <= self.real_threshold:
                self._message = ["Known value for minimum reached"]
                return True
        self._endtime = time.time()
        delta = self._endtime - self._starttime
        if delta >= self.maxtime:
            self._message = ["Time limit reached"]
            return True
        if self._nbfuncall >= self.maxfuncall:
            self._message = ["Number of function call reached"]
            return True

    def _stop_search(self):
        """Record time stamp when stop searching
        """
        if self.pure_sa:
            # In case of pure sa approach, doing ad LS at the end
            temp = np.array(self._xmin)
            etemp = self._ls_energy(temp)
            temp = self._xbuffer
            self._xmin = np.array(temp)
            self._emin = np.array(etemp)
        self._endtime = time.time()

    def _coordin(self):
        """Random generation of new coordinates
        """
        self._x = self._random_state.random_sample(
            self._x.size) * self._xrange + self._lower

    def _energy(self, x):
        """Calling objective function and adding elasticity if needed
        """
        delta_energy = 0
        if self.has_constraint:
            in_constraint = self.judge_constraint()
            if not in_constraint:
                self._etot = self.BIG_VALUE
                return 0
        if np.all(np.logical_and(
                x >= self._lower, x <= self._upper)):
            delta_energy = 0
        else:
            lcomp = x < self._lower
            ucomp = x > self._upper
            delta_energy_l = np.fabs(
                x[lcomp] - self._lower[lcomp]) * self.KSPRING
            delta_energy_u = np.fabs(
                x[ucomp] - self._upper[ucomp]) * self.KSPRING
            delta_energy = np.sum(delta_energy_l) + np.sum(
                delta_energy_u)
        self._etot = self.fun(x, *self.args)
        self._nbfuncall += 1
        self._etot = self._etot + delta_energy
        if np.isinf(self._etot) or np.isnan(self._etot):
            self._etot = self.BIG_VALUE

    def _ls_energy(self, x):
        """Performing a local search on the current location
        """
        self._xbuffer = np.array(x)
        self._smooth_search()
        self._x = np.array(self._xbuffer)
        return self._fvalue

    def _visita(self):
        """Visiting distribution function
        """
        pi = np.arcsin(1.0) * 2.0
        fator1 = np.exp(np.log(self._temperature) / (self.qv - 1.0))
        fator2 = np.exp((4.0 - self.qv) * np.log(self.qv - 1.0))
        fator3 = np.exp((2.0 - self.qv) * np.log(2.0) / (self.qv - 1.0))
        fator4 = np.sqrt(pi) * fator1 * fator2 / (fator3 * (3.0 - self.qv))
        fator5 = 1.0 / (self.qv - 1.0) - 0.5
        d1 = 2.0 - fator5
        fator6 = pi * (1.0 - fator5) / \
            np.sin(pi * (1.0 - fator5)) / np.exp(gammaln(d1))
        sigmax = np.exp(-(self.qv - 1.0) *
                        np.log(fator6 / fator4) / (3.0 - self.qv))
        x = sigmax * self._yygas()
        y = self._yygas()
        den = np.exp(
            (self.qv - 1.0) * np.log((np.fabs(y))) / (3.0 - self.qv))
        return x / den

    def _fobjective(self, x):
        """Wrapper to the objective function
        """
        self._x = np.array(x)
        self._energy(self._x)
        return self._etot

    def _smooth_search(self):
        """Local search implementation using setulb core function
        """
        m = 5
        iteration = 0
        f = 0
        n = self._x.size
        ndb = np.zeros(self._x.size)
        ndb.fill(2)
        iprint = -1
        x = np.array(self._xbuffer, np.float64)
        f = np.array(0.0, np.float64)
        l = np.array(self._lower, np.float64)
        u = np.array(self._upper, np.float64)
        g = np.zeros(n, np.float64)
        wa = np.zeros(2 * m * n + 5 * n + 11 * m * m + 8 * m, np.float64)
        iwa = np.zeros(3 * n, np.int32)
        task = np.zeros(1, 'S60')
        csave = np.zeros(1, 'S60')
        lsave = np.zeros(4, np.int32)
        isave = np.zeros(44, np.int32)
        dsave = np.zeros(29, np.float64)
        task[:] = 'START'
        if self.itsoftmax < 100:
            self.itsoftmax = 100
        elif self.itsoftmax > 1000:
            self.itsoftmax = 1000
        while True:
            if iteration >= self.itsoftmax:
                self._xbuffer = np.array(x)
                self._fvalue = f
                return
            _lbfgsb.setulb(m, x, l, u, ndb,
                           f, g, self.factr, self.pgtol, wa, iwa, task,
                           iprint, csave, lsave, isave, dsave, 100)
            iteration += 1
            task_str = task.tostring()
            if task_str.startswith(b'FG'):
                self._xbuffer = np.array(x)
                f = self._fobjective(self._xbuffer)
                if self.know_real:
                    if f <= self.real_threshold:
                        self._xbuffer = np.array(x)
                        self._fvalue = f
                        return
                g = np.array(self._compute_gradient(x), np.float64)
            elif task_str.startswith(b'NEW_X'):
                pass
            else:
                self._fvalue = f
                self._xbuffer = np.array(x)
                return

    def _compute_gradient(self, x):
        """Computation of numerical derivatives for local search
        """
        g = np.zeros(self._x.size, np.float64)
        for i in range(self._x.size):
            x1 = np.array(x)
            x2 = np.array(x)
            respl = self.reps
            respr = self.reps
            x1[i] = x[i] + respr
            if x1[i] > self._upper[i]:
                x1[i] = self._upper[i]
                respr = x1[i] - x[i]
            x2[i] = x[i] - respl
            if x2[i] < self._lower[i]:
                x2[i] = self._lower[i]
                respl = x[i] - x2[i]
            f1 = self._fobjective(x1)
            f2 = self._fobjective(x2)
            g[i] = ((f1 - f2)) / (respl + respr)
        idx = np.logical_or(np.isnan(g), np.isinf(g))
        g[idx] = 101.0
        return g

    def _yygas(self):
        if self._usey == 1:
            enter = True
            while(enter or (self._sgas <= 0 or self._sgas >= 1)):
                enter = False
                self._xgas = self._random_state.random_sample() * 2.0 - 1.0
                self._ygas = self._random_state.random_sample() * 2.0 - 1.0
                self._sgas = self._xgas * self._xgas + self._ygas * self._ygas
            root = np.sqrt(-2.0 / self._sgas * np.log(self._sgas))
            self._ranbyx = self._xgas * root
            self._ranbyy = self._ygas * root
            retval = self._ranbyy
            self._usey = 0
        else:
            retval = self._ranbyx
            self._usey = 1
        return retval

    @property
    def result(self):
        """ The OptimizeResult """
        res = OptimizeResult()
        res.x = self._xmin
        res.fun = self._fvalue
        res.message = self._message
        res.nit = self._step_record
        return res


def gensa(func, x0, bounds, maxiter=500, initial_temp=5230., visit=2.62,
        accept=-5.0, maxfun=1e7, args=(), seed=None, pure_sa=False):
    """
    Find the global minimum of a function using the Generalized Simulated
    Annealing algorithm

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
        The maximum number of gensa iterations
    initial_temp : float, optional
        The initial temperature, use higher values to facilitates a wider
        search of the energy landscape, allowing gensa to escape local minima
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
    args : tuple, optional
        Any additional fixed parameters needed to
        completely specify the objective function.
    seed : int or `np.random.RandomState`, optional
        If `seed` is not specified the `np.RandomState` singleton is used.
        If `seed` is an int, a new `np.random.RandomState` instance is used,
        seeded with seed.
        If `seed` is already a `np.random.RandomState instance`, then that
        `np.random.RandomState` instance is used.
        Specify `seed` for repeatable minimizations. The random numbers
        generated with this seed only affect the visiting distribution
        function and new coordinates generation.

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
    GenSA is an implementation of the General Simulated Annealing algorithm
    (GSA [2]_). This stochastic approach generalizes CSA (Classical Simulated
    Annealing) and FSA (Fast Simulated Annealing) to find the neighborhood of
    minima, then calls a local method (lbfgsb) to find their exact value.
    GenSA can process complicated and high dimension non-linear objective
    functions with a large number of local minima as described by [6]_.

    GSA uses a distorted Cauchy-Lorentz visiting distribution, with its shape
    controlled by the parameter :math:`q_{v}`

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

    .. versionadded:: 0.19.0

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

    >>> from scipy.optimize import gensa
    >>> func = lambda x: np.sum(x * x - 10 * np.cos(
    ...    2 * np.pi * x)) + 10 * np.size(x)
    >>> lw = [-5.12] * 10
    >>> up = [5.12] * 10
    >>> ret = gensa(func, None, bounds=(zip(lw, up)))
    """
    gr = GenSARunner(func, x0, bounds, args, seed,
            temperature_start=initial_temp, qv = visit, qa = accept,
            maxfun=maxfun, maxsteps=maxiter,pure_sa=pure_sa)
    gr.start_search()
    return gr.result

