##############################################################################
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# Author: Sylvain Gubian, PMP SA
##############################################################################
import numpy as np
from numpy.testing import assert_allclose
from numpy.random import RandomState
import sys
import time
from scipy.optimize import _lbfgsb
from scipy.special import gammaln

KSPRING = 1.e8
BIG_VALUE = 1.e13

# Here put explanation for why we have testing functions here
def rastrigin(x, y=None):
    # 2D implementation (used for plotting)
    if y is not None:
        return 20 + x ** 2 - 10 * np.cos(
            2 * np.pi * x) + y ** 2 - 10 * np.cos(2 * np.pi * y)
    # Generalized version
    return np.sum(x * x - 10 * np.cos(2 * np.pi * x)) + 10 * np.size(x)

def test_rastrigin():
    x = np.random.ranf(10) / 1.e+13
    assert(rastrigin(x) == 0.0)



class GenSARunnerException(Exception):
    pass


class Tracer():
    pass


class GenSARunner(object):

    def __init__(self, x0, fun, lower, upper, fargs={}):
        self._fun = fun
        self._fargs = fargs
        self._x = np.array(x0)
        self._lower = np.array(lower)
        self._upper = np.array(upper)
        self._dim, = self._x.shape
        self._tracer = Tracer()
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
        self._knowreal = False
        self._maxtime = 3600 # Using hours units
        self._maxfuncall = 1e7
        self._maxsteps = 5000
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

    def _judge_constraint(self):
        return True

    def initialize(self):
        self._nbfuncall = 0
        self._indTrace = 0
        if self._markovlength % self._dim != 0:
            raise GenSARunnerException('Incorrect markov length.')
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
                    raise GenSARunnerException(
                        ('Stopping algorithm because function create NaN',
                         'or (+/-) inifinity values even with trying new',
                         'random parameters'))
                rd = 0.0
                for i in range(self._dim):
                    self._x[i] = self._lower[i] + np.random.ranf() * (
                        self._upper[i] - self._lower[i])
                reinit_counter += 1
            else:
                initerror = False

    def start_search(self):
        self.nblocal += 1
        inconstraint = True
        itnew = 0
        emini_unchanged = True
        emini_markov = 0.0
        index_no_emini_update = 0
        index_tol_emini_update = 1000
        xminimarkov = np.zeros(self._dim)
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
        step_record = 0
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
                step_record += 1
                if step_record == self._maxsteps:
                    tolowtemp = False
                    break
                if self._temp < self._temprestart:
                    break
                temp_qa = self._temp / float(itnew)
                index_no_emini_update += 1
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
                        inconstraint = judge_constraint()
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
                        # Add data to tracer
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
                if not emini_unchanged:
                    temp = np.array(self._xmini)
                    etemp = self._ls_energy(temp)
                    temp = self._xbuffer
                    if etemp < self._emini:
                        self._xmini = np.array(temp)
                        self._emini = np.array(etemp)
                        index_no_emini_update = 0
                        # Add data to tracer
                if index_no_emini_update >= index_tol_emini_update - 1:
                    emini_markov = np.array(self._ls_energy(xminimarkov))
                    index_no_emini_update = 0
                    index_tol_emini_update = self._dim
                    if emini_markov < self._emini:
                        self._xmini = np.array(xminimarkov)
                        self._emini = np.array(emini_markov)
                        # Add data to tracer
                        if self._check_stopping():
                            self._stop_search()
                            return 0
            # End main loop
            #tolowtemp = False
        self._stop_search()
        return 0

    def _check_stopping(self):
        if self._knowreal:
            if self._emini <= self._realthreshold:
                return True
        self._endtime = time.time()
        delta = self._endtime - self._starttime
        if delta >= self._maxtime:
            return True
        if self._nbfuncall >= self._maxfuncall:
            return True
        # Check if no more improvments from tracer
        # Compare tracer value minernergy with self._emini
        # if lower than 1e-10, return true

    def _stop_search(self):
        self._endtime = time.time()

    def _coordin(self):
        for i in range(self._dim):
            self._x[i] = np.random.ranf() * self._xrange[i] + self._lower[i]

    def _energy(self, x):
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
        self._etot = self._fun(x, **self._fargs)
        self._nbfuncall += 1
        self._etot = self._etot + delta_energy
        if np.isinf(self._etot) or np.isnan(self._etot):
            self._etot = BIG_VALUE

    def _ls_energy(self, x):
        self._xbuffer = np.array(x)
        self._smooth_search()
        self._x = np.array(self._xbuffer)
        return self._fvalue

    # @profile
    def _visita(self):
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
        self._x = np.array(x)
        self._energy(self._x)
        return self._etot

    def _smooth_search(self):
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
                           iprint, csave, lsave, isave, dsave)
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
        # The below lines can be vectorized using np.vectorized on fobjective
        # to prevent doing a loop
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

    # @profile
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
    def x(self):
        """ The solution array """
        return self._x

    @property
    def f(self):
        """ The objective function value """
        return self._fvalue


def test_optimizer():
    np.random.seed(123)
    for fname in TEST_CONFIG:
        yield optim_func, fname


def optim_func(fname):
    dim = TEST_CONFIG[fname]['dim']
    lower = np.array(TEST_CONFIG[fname]['lower'])
    upper = np.array(TEST_CONFIG[fname]['upper'])
    func = TEST_CONFIG[fname]['func']
    if 'args' in TEST_CONFIG[fname]:
        fargs = TEST_CONFIG[fname]['args']
    else:
        fargs = {}
    xinit = lower + np.random.rand(dim) * (upper - lower)
    sr = GenSARunner(xinit, func, lower, upper, fargs=fargs)
    sr._knowreal = True
    sr._realthreshold = TEST_CONFIG[fname]['fvalue']
    sr.initialize()
    sr.start_search()
    print('\nValue found for {0}: {1}'.format(
        fname, sr.f))
    assert(sr.f <= TEST_CONFIG[fname]['fvalue'])


def main():
    dim = 2
    lower = np.array([-5] * dim)
    upper = np.array([5] * dim)
    xinit = lower + np.random.rand(dim) * (upper - lower)
    results = []
    for i in range(100):
        np.random.seed(1234 + i)
        sr = GenSARunner(xinit, rosenbrock, lower, upper)
        sr._knowreal = True
        sr._realthreshold = 1e-8
        sr.initialize()
        sr.start_search()
        print('Solution is: {0}'.format(sr.x))
        print('f value is: {0}'.format(sr.f))
        print('Nb local search: {0}'.format(np.mean(results[1])))
        print('Nb fn call average: {0}'.format(np.mean(results[0])))
        results.append((sr._nbfuncall, sr.nblocal))

    # sr._maxsteps = 15000

if __name__ == '__main__':
    main()
