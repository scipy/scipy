"""
gensa: A generalized simulated annealing global optimization algorithm
"""
from __future__ import division, print_function, absolute_import

import numpy as np
import time
from scipy.optimize import OptimizeResult
from scipy.optimize import _lbfgsb
from scipy.special import gammaln

__all__ = ['gensa']

KSPRING = 1.e8
BIG_VALUE = 1.e13
MAX_REINIT_COUNT = 100


class GenSARunner(object):
    """This class implements the core of the gensa algorithm.

    x0 : ndarray
        The starting coordinates.
    fun : callable
        The objective function
    lower: ndarray
        The lower bounds of the search domain
    upper: ndarray
        The upper bounds of the search domain
    fargs: dict
        Additional arguments to pass to the objective function

    """

    def __init__(self, fun, x0, lower, upper, fargs=()):
        self._lower = np.array(lower)
        self._upper = np.array(upper)
        if x0 is not None:
            if self._lower.size != x0.size:
                raise ValueError("Lower bounds size does not match x0")
            if self._upper.size != x0.size:
                raise ValueError("Upper bounds size does not match x0")
        if x0 is None:
            x0 = self._lower + np.random.rand(len(self._lower)) * (
                    self._upper - self._lower)
        else:
            x0 = np.array(x0)
        self._fun = fun
        self._fargs = fargs
        self._x = np.array(x0)
        self._dim, = self._x.shape
        self._xrange = self._upper - self._lower
        self._xbackup = np.array(self._x)
        self._xmini = np.array(self._x)
        self._xbuffer = np.zeros(self._dim)
        self._itsoftmax = self._dim * 6
        self._factr = 1000
        self._pgtol = 1.e-6
        self._reps = 1.e-6
        self._markovlength = self._dim * 2
        self._hasconstraint = False
        self._tempsta = 5230
        self._temp = self._tempsta
        self._knowreal = False
        self._maxtime = 3600
        self._maxfuncall = 1e7
        self._maxsteps = 500
        self._temprestart = 0.1
        self._qv = 2.62
        self._usey = 1
        self._xgas = 0.0
        self._ygas = 0.0
        self._ranbyx = 0.0
        self._ranbyy = 0.0
        self._sgas = -1.0
        self._qa = -5.0  # Acceptance parameter
        self.nblocal = 0
        self._step_record = 0
        self._nbfuncall = 0
        self._extensive = False

    def _judge_constraint(self):
        # Not implemented yet
        return True

    def initialize(self):
        """
        Random coordinate generation in case given initial coordinates
        are giving invalid objective value
        """
        self._nbfuncall = 0
        self._indTrace = 0
        if self._markovlength % self._dim != 0:
            raise ValueError('Incorrect markov length.')
        inconstraint = True
        initerror = True
        reinit_counter = 0
        while(initerror):
            if inconstraint:
                if self._hasconstraint:
                    inconstraint = self._judge_constraint()
                    while(not inconstraint):
                        self._coordin()
                        inconstraint = self._judge_constraint()
            self._energy(self._x)
            if self._etot >= BIG_VALUE:
                if reinit_counter >= MAX_REINIT_COUNT:
                    initerror = False
                    self._message = [(
                        'Stopping algorithm because function '
                        'create NaN or (+/-) inifinity values even with '
                        'trying new random parameters')]
                    raise ValueError(self._message)
                for i in range(self._dim):
                    self._x[i] = self._lower[i] + np.random.ranf() * (
                        self._upper[i] - self._lower[i])
                reinit_counter += 1
            else:
                initerror = False

    def start_search(self):
        """ Start annealing search process """
        self.nblocal += 1
        inconstraint = True
        itnew = 0
        emini_unchanged = True
        emini_markov = 0.0
        index_no_emini_update = 0
        index_tol_emini_update = 1000
        xminimarkov = np.zeros(self._dim)
        etot_ls_ex = BIG_VALUE
        # Here check that the function is smooth
        self._starttime = time.time()
        self._emini = np.array(self._etot)
        self._xmini = np.array(self._x)
        self._etot0 = np.array(self._etot)
        if self._etot < self._emini:
            self._emini = np.array(self._etot)
            self._xmini = np.array(self._x)
        self._etot0 = self._etot
        if self._check_stopping():
            self._stop_search()
        self._step_record = 0
        self._temp = self._tempsta
        tolowtemp = True
        while(tolowtemp):
            for i in range(self._maxsteps):
                itnew = i + 1
                s1 = float(itnew)
                s = s1 + 1.0
                t1 = np.exp((self._qv - 1) * np.log(2.0)) - 1.0
                t2 = np.exp((self._qv - 1) * np.log(s)) - 1.0
                self._temp = self._tempsta * t1 / t2
                self._step_record += 1
                # Maximum number of iterations reached
                if self._step_record == self._maxsteps:
                    tolowtemp = False
                    break
                # Need a reannealing process
                if self._temp < self._temprestart:
                    break
                temp_qa = self._temp / float(itnew)
                index_no_emini_update += 1
                # Starting markov chain loop
                for j in range(self._markovlength):
                    if j == 0:
                        emini_unchanged = True
                    if i == 0 and j == 0:
                        emini_unchanged = False
                    self._xbackup = np.array(self._x)
                    # inconstraint = False
                    # while (not inconstraint):
                    if j < self._dim:
                        for k in range(self._dim):
                            visit = self._visita()
                            if visit > 1.e8:
                                visit = 1.e8 * np.random.ranf()
                            elif visit < -1e8:
                                visit = -1e8 * np.random.ranf()
                            self._x[k] = visit + self._xbackup[k]
                            a = self._x[k] - self._lower[k]
                            b = np.fmod(a, self._xrange[k]) + self._xrange[k]
                            self._x[k] = np.fmod(
                                b, self._xrange[k]) + self._lower[k]
                            if np.fabs(self._x[k] - self._lower[k] < 1.e-10):
                                self._x[k] += 1.e-10
                        # end coordinates loop
                    else:
                        # Now change only one component at a time
                        visit = self._visita()
                        if visit > 1.e8:
                            visit = 1.e8 * np.random.ranf()
                        elif visit < -1e8:
                            visit = -1.e8 * np.random.ranf()
                        index = j - self._dim
                        self._x[index] = visit + self._xbackup[index]
                        a = self._x[index] - self._lower[index]
                        b = np.fmod(
                            a, self._xrange[index]) + self._xrange[index]
                        self._x[index] = np.fmod(
                            b, self._xrange[index]) + self._lower[index]
                        if np.fabs(self._x[index] - self._lower[
                                index]) < 1.e-10:
                            self._x[index] += 1.e-10
                    if self._hasconstraint:
                        inconstraint = self._judge_constraint()
                    else:
                        inconstraint = True
                    if inconstraint:
                        self._energy(self._x)
                        if self._etot < self._etot0:
                            self._etot0 = np.array(self._etot)
                            if self._etot < self._emini:
                                self._emini = np.array(self._etot)
                                self._xmini = np.array(self._x)
                                emini_unchanged = False
                                index_no_emini_update = 0
                        else:
                            r = np.random.ranf()
                            pqa1 = (self._qa - 1.0) * (
                                self._etot - self._etot0) / temp_qa + 1.0
                            if pqa1 < 0.0:
                                pqa = 0.0
                            else:
                                pqa = np.exp(np.log(pqa1) / (1.0 - self._qa))
                            if r > pqa:
                                self._x = self._xbackup
                            else:
                                self._etot0 = self._etot
                        if self._check_stopping():
                            self._stop_search()
                            return 0
                    # End of while not inconstraint
                    if index_no_emini_update >= index_tol_emini_update - 1:
                        if j == 0:
                            emini_markov = np.array(self._etot0)
                            xminimarkov = np.array(self._x)
                        else:
                            if self._etot0 < emini_markov:
                                emini_markov = np.array(self._etot0)
                                xminimarkov = np.array(self._x)
                # End of markov chain loop
                if self._extensive:
                    if j == 0:
                        #print("-> First extensive")
                        temp = np.array(self._x)
                        etot_ls_ex = self._ls_energy(temp)
                        index_no_emini_update = 0
                    else:
                        if np.sum((self._etot - etot_ls_ex)**2) > 1e-10:
                            #print("->extensive")
                            temp = np.array(self._x)
                            etot_ls_ex = self._ls_energy(temp)
                            index_no_emini_update = 0
                    if etot_ls_ex < self._emini:
                        #print("-> better VALUE in extensive")
                        self._xmini = np.array(temp)
                        self._emini = np.array(etot_ls_ex)
                else:
                    if not emini_unchanged:
                        temp = np.array(self._xmini)
                        etemp = self._ls_energy(temp)
                        temp = self._xbuffer
                        if etemp < self._emini:
                            self._xmini = np.array(temp)
                            self._emini = np.array(etemp)
                            index_no_emini_update = 0
                if index_no_emini_update >= index_tol_emini_update - 1:
                    emini_markov = np.array(self._ls_energy(xminimarkov))
                    index_no_emini_update = 0
                    index_tol_emini_update = self._dim
                    if emini_markov < self._emini:
                        self._xmini = np.array(xminimarkov)
                        self._emini = np.array(emini_markov)
                if self._check_stopping():
                    self._stop_search()
                    return 0
            # End main loop
            #tolowtemp = False
        self._stop_search()
        self._message = ["Number of iteration reached"]
        return 0

    def _check_stopping(self):
        """Check if the search has to be stopped"""
        if self._knowreal:
            if self._emini <= self._realthreshold:
                self._message = ["Known value for minimum reached"]
                return True
        self._endtime = time.time()
        delta = self._endtime - self._starttime
        if delta >= self._maxtime:
            self._message = ["Time limit reached"]
            return True
        if self._nbfuncall >= self._maxfuncall:
            self._message = ["Number of function call reached"]
            return True
        #TODO:
        # Check if no more improvments from tracer
        # Compare tracer value minernergy with self._emini
        # if lower than 1e-10, return true

    def _stop_search(self):
        """Record time stamp when stop searching"""
        self._endtime = time.time()

    def _coordin(self):
        """Random generation of new coordinates"""
        for i in range(self._dim):
            self._x[i] = np.random.ranf() * self._xrange[i] + self._lower[i]

    def _energy(self, x):
        """Calling objective function and adding elasticity if needed"""
        delta_energy = 0
        if self._hasconstraint:
            inconstraint = self._judge_constraint()
            if not inconstraint:
                self._etot = BIG_VALUE
                return 0
        if np.all(np.logical_and(
                x >= self._lower, x <= self._upper)):
            delta_energy = 0
        else:
            lcomp = x < self._lower
            ucomp = x > self._upper
            delta_energy_l = np.fabs(
                x[lcomp] - self._lower[lcomp]) * KSPRING
            delta_energy_u = np.fabs(
                x[ucomp] - self._upper[ucomp]) * KSPRING
            delta_energy = np.sum(delta_energy_l) + np.sum(
                delta_energy_u)
        self._etot = self._fun(x, *self._fargs)
        self._nbfuncall += 1
        self._etot = self._etot + delta_energy
        if np.isinf(self._etot) or np.isnan(self._etot):
            self._etot = BIG_VALUE

    def _ls_energy(self, x):
        """Performing a local search on the current location"""
        self._xbuffer = np.array(x)
        self._smooth_search()
        self._x = np.array(self._xbuffer)
        return self._fvalue

    def _visita(self):
        """Visiting function"""
        pi = np.arcsin(1.0) * 2.0
        fator1 = np.exp(np.log(self._temp) / (self._qv - 1.0))
        fator2 = np.exp((4.0 - self._qv) * np.log(self._qv - 1.0))
        fator3 = np.exp((2.0 - self._qv) * np.log(2.0) / (self._qv - 1.0))
        fator4 = np.sqrt(pi) * fator1 * fator2 / (fator3 * (3.0 - self._qv))
        fator5 = 1.0 / (self._qv - 1.0) - 0.5
        d1 = 2.0 - fator5
        fator6 = pi * (1.0 - fator5) / \
            np.sin(pi * (1.0 - fator5)) / np.exp(gammaln(d1))
        sigmax = np.exp(-(self._qv - 1.0) *
                        np.log(fator6 / fator4) / (3.0 - self._qv))
        x = sigmax * self._yygas()
        y = self._yygas()
        den = np.exp(
            (self._qv - 1.0) * np.log((np.fabs(y))) / (3.0 - self._qv))
        return x / den

    def _fobjective(self, x):
        """Wrapper to the objective function"""
        self._x = np.array(x)
        self._energy(self._x)
        return self._etot

    def _smooth_search(self):
        """Local search implementation using setulb core function"""
        m = 5
        iteration = 0
        f = 0
        n = self._dim
        ndb = np.zeros(self._dim)
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
        if self._itsoftmax < 100:
            self._itsoftmax = 100
        elif self._itsoftmax > 1000:
            self._itsoftmax = 1000
        while True:
            if iteration >= self._itsoftmax:
                self._xbuffer = np.array(x)
                self._fvalue = f
                return
            _lbfgsb.setulb(m, x, l, u, ndb,
                           f, g, self._factr, self._pgtol, wa, iwa, task,
                           iprint, csave, lsave, isave, dsave, 100)
            iteration += 1
            task_str = task.tostring()
            if task_str.startswith(b'FG'):
                self._xbuffer = np.array(x)
                f = self._fobjective(self._xbuffer)
                if self._knowreal:
                    if f <= self._realthreshold:
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
        """Computation of numerical derivatives for local search"""
        # Can be done later here:
        # Check here if the expression of the derivatives is provided
        # For that, add an argument for jacobian in the main function
        # and use MemoizeJac from optimize.py
        g = np.zeros(self._dim, np.float64)
        for i in range(self._dim):
            x1 = np.array(x)
            x2 = np.array(x)
            respl = self._reps
            respr = self._reps
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
                self._xgas = np.random.ranf() * 2.0 - 1.0
                self._ygas = np.random.ranf() * 2.0 - 1.0
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
        res.x = self._xmini
        res.fun = self._fvalue
        res.message = self._message
        res.nit = self._step_record
        return res

def gensa(func, x0, bounds, niter=500, T=5230., visitparam=2.62,
        acceptparam=-5.0, maxtime=3600, maxcall=1e7, args=()):
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
    bounds : sequence
        Bounds for variables.  ``(min, max)`` pairs for each element in ``x``,
        defining the lower and upper bounds for the optimizing argument of
        `func`. It is required to have ``len(bounds) == len(x)``.
        ``len(bounds)`` is used to determine the number of parameters in ``x``.
    bounds: sequence
        Bounds for variable. (min, max) pairs for each element in x,
        defining the bounds on that parameter.
    niter: integer, optional
        The number of gensa iterations
    T: float, optional
        Initial value for temperature. The higher the initial temperature, the
        the higher the probability of hill-climbing for the initial iterations.
    visitparam: float, optional
        Parameter for visiting distribution. Higher value means the
        visiting distribution has a heavier tail which makes the
        algorithm jump to a farer region. The value range is ]0; 3]
    acceptparam: float, optional
        Parameter for acceptance distribution. It is used to control the
        probability of acceptance. The lower the acceptance parameter, the
        smaller the probability of acceptance. It has to be any negative value.
    maxtime: integer, optional
        Time limit for running the algorithm in seconds
    maxcall: integer, optional
        Soft limit for the number of objective function calls. If the
        algorithm is in the middle of a local search, this number will be
        exceeded, the algorithm will stop just after the local search is
        done.
    args : tuple, optional
        Any additional fixed parameters needed to
        completely specify the objective function.


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
    minima, then calls a local method (lbfgs) to find their exact value.
    GenSA can process complicated and high dimension non-linear objective
    functions with a large number of local minima as described by Muller
    paper [6]_.
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
    >>> func = lambda x: np.sum(x * x - 10 * np.cos(2 * np.pi * x))\
            + 10 * np.size(x)
    >>> lw = [-5.12] * 10
    >>> up = [5.12] * 10
    >>> ret = gensa(func, None, bounds=(zip(lw, up)))
    >>> print("global minimum: xmin = {0}, f(xmin) = {1}".format(
    ...    ret.x, ret.fun))
    """
    lower, upper = zip(*bounds)
    gr = GenSARunner(func, x0, np.array(lower), np.array(upper), fargs=args)
    gr._tempsta = T
    gr._qv = visitparam
    gr._qa = acceptparam
    gr._maxtime = maxtime
    gr._maxfuncall = maxcall
    gr._maxsteps = niter
    gr._extensive = False
    gr.initialize()
    gr.start_search()
    return gr.result


def main():
    func = lambda x: np.sum(x * x - 10 * np.cos(2 * np.pi * x))\
        + 10 * np.size(x)
    dim = 5
    lw = np.array([-5.12] * dim)
    up = np.array([5.12] * dim)
    xinit = lw + np.random.rand(dim) * (up - lw)
    ret = gensa(func, xinit, (zip(lw, up)))
    print("global minimum: xmin = {0}, f(xmin) = {1}".format(
        ret.x, ret.fun))

if __name__ == '__main__':
    main()

