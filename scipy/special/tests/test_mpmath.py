"""
Test Scipy functions versus mpmath, if available.

"""
from __future__ import division, print_function, absolute_import

import sys
import os
import time

from distutils.version import LooseVersion

import numpy as np
from numpy.testing import dec, run_module_suite
from numpy import pi

import scipy.special as sc
from scipy.lib.six import reraise, with_metaclass
from scipy.special._testutils import FuncData, assert_func_equal

try:
    import mpmath
except ImportError:
    try:
        import sympy.mpmath as mpmath
    except ImportError:
        mpmath = None


def mpmath_check(min_ver):
    if mpmath is None:
        return dec.skipif(True, "mpmath is not installed")
    return dec.skipif(LooseVersion(mpmath.__version__) < LooseVersion(min_ver),
                      "mpmath version >= %s required" % min_ver)


#------------------------------------------------------------------------------
# expi
#------------------------------------------------------------------------------

@mpmath_check('0.10')
def test_expi_complex():
    dataset = []
    for r in np.logspace(-99, 2, 10):
        for p in np.linspace(0, 2*np.pi, 30):
            z = r*np.exp(1j*p)
            dataset.append((z, complex(mpmath.ei(z))))
    dataset = np.array(dataset, dtype=np.complex_)

    FuncData(sc.expi, dataset, 0, 1).check()


#------------------------------------------------------------------------------
# hyp2f1
#------------------------------------------------------------------------------

@mpmath_check('0.14')
def test_hyp2f1_strange_points():
    pts = [
        (2,-1,-1,0.7),
        (2,-2,-2,0.7),
    ]
    kw = dict(eliminate=True)
    dataset = [p + (float(mpmath.hyp2f1(*p, **kw)),) for p in pts]
    dataset = np.array(dataset, dtype=np.float_)

    FuncData(sc.hyp2f1, dataset, (0,1,2,3), 4, rtol=1e-10).check()


@mpmath_check('0.13')
def test_hyp2f1_real_some_points():
    pts = [
        (1,2,3,0),
        (1./3, 2./3, 5./6, 27./32),
        (1./4, 1./2, 3./4, 80./81),
        (2,-2,-3,3),
        (2,-3,-2,3),
        (2,-1.5,-1.5,3),
        (1,2,3,0),
        (0.7235, -1, -5, 0.3),
        (0.25, 1./3, 2, 0.999),
        (0.25, 1./3, 2, -1),
        (2,3,5,0.99),
        (3./2,-0.5,3,0.99),
        (2,2.5,-3.25,0.999),
        (-8, 18.016500331508873, 10.805295997850628, 0.90875647507000001),
        (-10,900,-10.5,0.99),
        (-10,900,10.5,0.99),
        (-1,2,1,1.0),
        (-1,2,1,-1.0),
        (-3,13,5,1.0),
        (-3,13,5,-1.0),
        (0.5, 1 - 270.5, 1.5, 0.999**2),  # from issue 1561
    ]
    dataset = [p + (float(mpmath.hyp2f1(*p)),) for p in pts]
    dataset = np.array(dataset, dtype=np.float_)

    olderr = np.seterr(invalid='ignore')
    try:
        FuncData(sc.hyp2f1, dataset, (0,1,2,3), 4, rtol=1e-10).check()
    finally:
        np.seterr(**olderr)


@mpmath_check('0.14')
def test_hyp2f1_some_points_2():
    # Taken from mpmath unit tests -- this point failed for mpmath 0.13 but
    # was fixed in their SVN since then
    pts = [
        (112, (51,10), (-9,10), -0.99999),
        (10,-900,10.5,0.99),
        (10,-900,-10.5,0.99),
    ]

    def fev(x):
        if isinstance(x, tuple):
            return float(x[0]) / x[1]
        else:
            return x

    dataset = [tuple(map(fev, p)) + (float(mpmath.hyp2f1(*p)),) for p in pts]
    dataset = np.array(dataset, dtype=np.float_)

    FuncData(sc.hyp2f1, dataset, (0,1,2,3), 4, rtol=1e-10).check()


@mpmath_check('0.13')
def test_hyp2f1_real_some():
    dataset = []
    for a in [-10, -5, -1.8, 1.8, 5, 10]:
        for b in [-2.5, -1, 1, 7.4]:
            for c in [-9, -1.8, 5, 20.4]:
                for z in [-10, -1.01, -0.99, 0, 0.6, 0.95, 1.5, 10]:
                    try:
                        v = float(mpmath.hyp2f1(a, b, c, z))
                    except:
                        continue
                    dataset.append((a, b, c, z, v))
    dataset = np.array(dataset, dtype=np.float_)

    olderr = np.seterr(invalid='ignore')
    try:
        FuncData(sc.hyp2f1, dataset, (0,1,2,3), 4, rtol=1e-9,
                 ignore_inf_sign=True).check()
    finally:
        np.seterr(**olderr)


@mpmath_check('0.12')
@dec.slow
def test_hyp2f1_real_random():
    dataset = []

    npoints = 500
    dataset = np.zeros((npoints, 5), np.float_)

    np.random.seed(1234)
    dataset[:,0] = np.random.pareto(1.5, npoints)
    dataset[:,1] = np.random.pareto(1.5, npoints)
    dataset[:,2] = np.random.pareto(1.5, npoints)
    dataset[:,3] = 2*np.random.rand(npoints) - 1

    dataset[:,0] *= (-1)**np.random.randint(2, npoints)
    dataset[:,1] *= (-1)**np.random.randint(2, npoints)
    dataset[:,2] *= (-1)**np.random.randint(2, npoints)

    for ds in dataset:
        if mpmath.__version__ < '0.14':
            # mpmath < 0.14 fails for c too much smaller than a, b
            if abs(ds[:2]).max() > abs(ds[2]):
                ds[2] = abs(ds[:2]).max()
        ds[4] = float(mpmath.hyp2f1(*tuple(ds[:4])))

    FuncData(sc.hyp2f1, dataset, (0,1,2,3), 4, rtol=1e-9).check()

#------------------------------------------------------------------------------
# erf (complex)
#------------------------------------------------------------------------------


@mpmath_check('0.14')
def test_erf_complex():
    # need to increase mpmath precision for this test
    old_dps, old_prec = mpmath.mp.dps, mpmath.mp.prec
    try:
        mpmath.mp.dps = 70
        x1, y1 = np.meshgrid(np.linspace(-10, 1, 31), np.linspace(-10, 1, 11))
        x2, y2 = np.meshgrid(np.logspace(-80, .8, 31), np.logspace(-80, .8, 11))
        points = np.r_[x1.ravel(),x2.ravel()] + 1j*np.r_[y1.ravel(),y2.ravel()]

        assert_func_equal(sc.erf, lambda x: complex(mpmath.erf(x)), points,
                          vectorized=False, rtol=1e-13)
        assert_func_equal(sc.erfc, lambda x: complex(mpmath.erfc(x)), points,
                          vectorized=False, rtol=1e-13)
    finally:
        mpmath.mp.dps, mpmath.mp.prec = old_dps, old_prec


#------------------------------------------------------------------------------
# lpmv
#------------------------------------------------------------------------------

@mpmath_check('0.15')
def test_lpmv():
    pts = []
    for x in [-0.99, -0.557, 1e-6, 0.132, 1]:
        pts.extend([
            (1, 1, x),
            (1, -1, x),
            (-1, 1, x),
            (-1, -2, x),
            (1, 1.7, x),
            (1, -1.7, x),
            (-1, 1.7, x),
            (-1, -2.7, x),
            (1, 10, x),
            (1, 11, x),
            (3, 8, x),
            (5, 11, x),
            (-3, 8, x),
            (-5, 11, x),
            (3, -8, x),
            (5, -11, x),
            (-3, -8, x),
            (-5, -11, x),
            (3, 8.3, x),
            (5, 11.3, x),
            (-3, 8.3, x),
            (-5, 11.3, x),
            (3, -8.3, x),
            (5, -11.3, x),
            (-3, -8.3, x),
            (-5, -11.3, x),
        ])

    def mplegenp(nu, mu, x):
        if mu == int(mu) and x == 1:
            # mpmath 0.17 gets this wrong
            if mu == 0:
                return 1
            else:
                return 0
        return mpmath.legenp(nu, mu, x)

    dataset = [p + (mplegenp(p[1], p[0], p[2]),) for p in pts]
    dataset = np.array(dataset, dtype=np.float_)

    def evf(mu, nu, x):
        return sc.lpmv(mu.astype(int), nu, x)

    olderr = np.seterr(invalid='ignore')
    try:
        FuncData(evf, dataset, (0,1,2), 3, rtol=1e-10, atol=1e-14).check()
    finally:
        np.seterr(**olderr)


#------------------------------------------------------------------------------
# beta
#------------------------------------------------------------------------------

@mpmath_check('0.15')
def test_beta():
    np.random.seed(1234)

    b = np.r_[np.logspace(-200, 200, 4),
              np.logspace(-10, 10, 4),
              np.logspace(-1, 1, 4),
              np.arange(-10, 11, 1),
              np.arange(-10, 11, 1) + 0.5,
              -1, -2.3, -3, -100.3, -10003.4]
    a = b

    ab = np.array(np.broadcast_arrays(a[:,None], b[None,:])).reshape(2, -1).T

    old_dps, old_prec = mpmath.mp.dps, mpmath.mp.prec
    try:
        mpmath.mp.dps = 400

        assert_func_equal(sc.beta,
                          lambda a, b: float(mpmath.beta(a, b)),
                          ab,
                          vectorized=False,
                          rtol=1e-10,
                          ignore_inf_sign=True)

        assert_func_equal(
            sc.betaln,
            lambda a, b: float(mpmath.log(abs(mpmath.beta(a, b)))),
            ab,
            vectorized=False,
            rtol=1e-10)
    finally:
        mpmath.mp.dps, mpmath.mp.prec = old_dps, old_prec


#------------------------------------------------------------------------------
# Machinery for systematic tests
#------------------------------------------------------------------------------

class Arg(object):
    """
    Generate a set of numbers on the real axis, concentrating on
    'interesting' regions and covering all orders of magnitude.
    """

    def __init__(self, a=-np.inf, b=np.inf, inclusive_a=True, inclusive_b=True):
        self.a = a
        self.b = b
        self.inclusive_a = inclusive_a
        self.inclusive_b = inclusive_b
        if self.a == -np.inf:
            self.a = -np.finfo(float).max/2
        if self.b == np.inf:
            self.b = np.finfo(float).max/2

    def values(self, n):
        """Return an array containing approximatively `n` numbers."""
        n1 = max(2, int(0.3*n))
        n2 = max(2, int(0.2*n))
        n3 = max(8, n - n1 - n2)

        v1 = np.linspace(-1, 1, n1)
        v2 = np.r_[np.linspace(-10, 10, max(0, n2-4)),
                   -9, -5.5, 5.5, 9]
        if self.a >= 0 and self.b > 0:
            v3 = np.r_[
                np.logspace(-30, -1, 2 + n3//4),
                np.logspace(5, np.log10(self.b), 1 + n3//4),
                ]
            v4 = np.logspace(1, 5, 1 + n3//2)
        elif self.a < 0 and self.b > 0:
            v3 = np.r_[
                np.logspace(-30, -1, 2 + n3//8),
                np.logspace(5, np.log10(self.b), 1 + n3//8),
                -np.logspace(-30, -1, 2 + n3//8),
                -np.logspace(5, np.log10(-self.a), 1 + n3//8)
                ]
            v4 = np.r_[
                np.logspace(1, 5, 1 + n3//4),
                -np.logspace(1, 5, 1 + n3//4)
                ]
        elif self.b < 0:
            v3 = np.r_[
                -np.logspace(-30, -1, 2 + n3//4),
                -np.logspace(5, np.log10(-self.b), 1 + n3//4),
                ]
            v4 = -np.logspace(1, 5, 1 + n3//2)
        else:
            v3 = []
            v4 = []
        v = np.r_[v1, v2, v3, v4, 0]
        if self.inclusive_a:
            v = v[v >= self.a]
        else:
            v = v[v > self.a]
        if self.inclusive_b:
            v = v[v <= self.b]
        else:
            v = v[v < self.b]
        return np.unique(v)


class FixedArg(object):
    def __init__(self, values):
        self._values = np.asarray(values)

    def values(self, n):
        return self._values


class ComplexArg(object):
    def __init__(self, a=complex(-np.inf, -np.inf), b=complex(np.inf, np.inf)):
        self.real = Arg(a.real, b.real)
        self.imag = Arg(a.imag, b.imag)

    def values(self, n):
        m = max(2, int(np.sqrt(n)))
        x = self.real.values(m)
        y = self.imag.values(m)
        return (x[:,None] + 1j*y[None,:]).ravel()


class IntArg(object):
    def __init__(self, a=-1000, b=1000):
        self.a = a
        self.b = b

    def values(self, n):
        v1 = Arg(self.a, self.b).values(max(1 + n//2, n-5)).astype(int)
        v2 = np.arange(-5, 5)
        v = np.unique(np.r_[v1, v2])
        v = v[(v >= self.a) & (v < self.b)]
        return v


class MpmathData(object):
    def __init__(self, scipy_func, mpmath_func, arg_spec, name=None,
                 dps=None, prec=None, n=5000, rtol=1e-7, atol=1e-300,
                 ignore_inf_sign=False):
        self.scipy_func = scipy_func
        self.mpmath_func = mpmath_func
        self.arg_spec = arg_spec
        self.dps = dps
        self.prec = prec
        self.n = n
        self.rtol = rtol
        self.atol = atol
        self.ignore_inf_sign = ignore_inf_sign
        if isinstance(self.arg_spec, np.ndarray):
            self.is_complex = np.issubdtype(self.arg_spec.dtype, np.complexfloating)
        else:
            self.is_complex = any([isinstance(arg, ComplexArg) for arg in self.arg_spec])
        self.ignore_inf_sign = ignore_inf_sign
        if not name or name == '<lambda>':
            name = getattr(scipy_func, '__name__', None)
        if not name or name == '<lambda>':
            name = getattr(mpmath_func, '__name__', None)
        self.name = name

    def check(self):
        np.random.seed(1234)

        # Generate values for the arguments
        if isinstance(self.arg_spec, np.ndarray):
            argarr = self.arg_spec.copy()
        else:
            num_args = len(self.arg_spec)
            ms = np.asarray([1.5 if isinstance(arg, ComplexArg) else 1.0
                             for arg in self.arg_spec])
            ms = (self.n**(ms/sum(ms))).astype(int) + 1

            argvals = []
            for arg, m in zip(self.arg_spec, ms):
                argvals.append(arg.values(m))

            argarr = np.array(np.broadcast_arrays(*np.ix_(*argvals))).reshape(num_args, -1).T

        # Check
        old_dps, old_prec = mpmath.mp.dps, mpmath.mp.prec
        try:
            if self.dps is not None:
                dps_list = [self.dps]
            else:
                dps_list = [20]
            if self.prec is not None:
                mpmath.mp.prec = self.prec

            # Proper casting of mpmath input and output types. Using
            # native mpmath types as inputs gives improved precision
            # in some cases.
            if np.issubdtype(argarr.dtype, np.complexfloating):
                pytype = complex
                mptype = lambda x: mpmath.mpc(complex(x))
            else:
                mptype = lambda x: mpmath.mpf(float(x))

                def pytype(x):
                    if abs(x.imag) > 1e-16*(1 + abs(x.real)):
                        return np.nan
                    else:
                        return float(x.real)

            # Try out different dps until one (or none) works
            for j, dps in enumerate(dps_list):
                mpmath.mp.dps = dps

                try:
                    assert_func_equal(self.scipy_func,
                                      lambda *a: pytype(self.mpmath_func(*map(mptype, a))),
                                      argarr,
                                      vectorized=False,
                                      rtol=self.rtol, atol=self.atol,
                                      ignore_inf_sign=self.ignore_inf_sign,
                                      nan_ok=True)
                    break
                except AssertionError:
                    if j >= len(dps_list)-1:
                        reraise(*sys.exc_info())
        finally:
            mpmath.mp.dps, mpmath.mp.prec = old_dps, old_prec

    def __repr__(self):
        if self.is_complex:
            return "<MpmathData: %s (complex)>" % (self.name,)
        else:
            return "<MpmathData: %s>" % (self.name,)


def assert_mpmath_equal(*a, **kw):
    d = MpmathData(*a, **kw)
    d.check()


def nonfunctional_tooslow(func):
    return dec.skipif(True, "    Test not yet functional (too slow), needs more work.")(func)


def knownfailure_overridable(msg=None):
    if not msg:
        msg = "Undiagnosed issues (corner cases, wrong comparison values, or otherwise)"
    msg = msg + " [Set environment variable SCIPY_XFAIL=1 to run this test nevertheless.]"

    def deco(func):
        try:
            if bool(os.environ['SCIPY_XFAIL']):
                return func
        except (ValueError, KeyError):
            pass
        return dec.knownfailureif(True, msg)(func)
    return deco


class _SystematicMeta(type):
    """
    Metaclass which decorates all of the test_* methods with

    - @mpmath_check(...)
    - @dec.slow

    """

    mpmath_min_version = '0.17'

    def __new__(cls, cls_name, bases, dct):
        for name, item in list(dct.items()):
            if name.startswith('test_'):
                item = dec.slow(item)
                item = mpmath_check(cls.mpmath_min_version)(item)
                dct[name] = item
        return type.__new__(cls, cls_name, bases, dct)


#------------------------------------------------------------------------------
# Dealing with mpmath quirks
#------------------------------------------------------------------------------

def _trace_args(func):
    def tofloat(x):
        if isinstance(x, mpmath.mpc):
            return complex(x)
        else:
            return float(x)

    def wrap(*a, **kw):
        sys.stderr.write("%r: " % (tuple(map(tofloat, a)),))
        sys.stderr.flush()
        try:
            r = func(*a, **kw)
            sys.stderr.write("-> %r" % r)
        finally:
            sys.stderr.write("\n")
            sys.stderr.flush()
        return r
    return wrap

try:
    import posix
    import signal
    POSIX = ('setitimer' in dir(signal))
except ImportError:
    POSIX = False


class _TimeoutError(Exception):
    pass


def _time_limited(timeout=0.5, return_val=np.nan, use_sigalrm=True):
    """
    Decorator for setting a timeout for pure-Python functions.

    If the function does not return within `timeout` seconds, the
    value `return_val` is returned instead.

    On POSIX this uses SIGALRM by default. On non-POSIX, settrace is
    used. Do not use this with threads: the SIGALRM implementation
    does probably not work well. The settrace implementation only
    traces the current thread.

    The settrace implementation slows down execution speed. Slowdown
    by a factor around 10 is probably typical.
    """
    if POSIX and use_sigalrm:
        def sigalrm_handler(signum, frame):
            raise _TimeoutError()

        def deco(func):
            def wrap(*a, **kw):
                old_handler = signal.signal(signal.SIGALRM, sigalrm_handler)
                signal.setitimer(signal.ITIMER_REAL, timeout)
                try:
                    return func(*a, **kw)
                except _TimeoutError:
                    return return_val
                finally:
                    signal.setitimer(signal.ITIMER_REAL, 0)
                    signal.signal(signal.SIGALRM, old_handler)
            return wrap
    else:
        def deco(func):
            def wrap(*a, **kw):
                start_time = time.time()

                def trace(frame, event, arg):
                    if time.time() - start_time > timeout:
                        raise _TimeoutError()
                    return None  # turn off tracing except at function calls
                sys.settrace(trace)
                try:
                    return func(*a, **kw)
                except _TimeoutError:
                    sys.settrace(None)
                    return return_val
                finally:
                    sys.settrace(None)
            return wrap
    return deco


def _exception_to_nan(func):
    """Decorate function to return nan if it raises an exception"""
    def wrap(*a, **kw):
        try:
            return func(*a, **kw)
        except Exception:
            return np.nan
    return wrap


def _inf_to_nan(func):
    """Decorate function to return nan if it returns inf"""
    def wrap(*a, **kw):
        v = func(*a, **kw)
        if not np.isfinite(v):
            return np.nan
        return v
    return wrap


#------------------------------------------------------------------------------
# Systematic tests
#------------------------------------------------------------------------------

HYPERKW = dict(maxprec=200, maxterms=200)


class TestSystematic(with_metaclass(_SystematicMeta, object)):
    def test_airyai(self):
        # oscillating function, limit range
        assert_mpmath_equal(lambda z: sc.airy(z)[0],
                            mpmath.airyai,
                            [Arg(-1e8, 1e8)],
                            rtol=1e-5)
        assert_mpmath_equal(lambda z: sc.airy(z)[0],
                            mpmath.airyai,
                            [Arg(-1e3, 1e3)])

    def test_airyai_complex(self):
        assert_mpmath_equal(lambda z: sc.airy(z)[0],
                            mpmath.airyai,
                            [ComplexArg()])

    def test_airyai_prime(self):
        # oscillating function, limit range
        assert_mpmath_equal(lambda z: sc.airy(z)[1], lambda z:
                            mpmath.airyai(z, derivative=1),
                            [Arg(-1e8, 1e8)],
                            rtol=1e-5)
        assert_mpmath_equal(lambda z: sc.airy(z)[1], lambda z:
                            mpmath.airyai(z, derivative=1),
                            [Arg(-1e3, 1e3)])

    def test_airyai_prime_complex(self):
        assert_mpmath_equal(lambda z: sc.airy(z)[1], lambda z:
                            mpmath.airyai(z, derivative=1),
                            [ComplexArg()])

    def test_airybi(self):
        # oscillating function, limit range
        assert_mpmath_equal(lambda z: sc.airy(z)[2], lambda z:
                            mpmath.airybi(z),
                            [Arg(-1e8, 1e8)],
                            rtol=1e-5)
        assert_mpmath_equal(lambda z: sc.airy(z)[2], lambda z:
                            mpmath.airybi(z),
                            [Arg(-1e3, 1e3)])

    def test_airybi_complex(self):
        assert_mpmath_equal(lambda z: sc.airy(z)[2], lambda z:
                            mpmath.airybi(z),
                            [ComplexArg()])

    def test_airybi_prime(self):
        # oscillating function, limit range
        assert_mpmath_equal(lambda z: sc.airy(z)[3], lambda z:
                            mpmath.airybi(z, derivative=1),
                            [Arg(-1e8, 1e8)],
                            rtol=1e-5)
        assert_mpmath_equal(lambda z: sc.airy(z)[3], lambda z:
                            mpmath.airybi(z, derivative=1),
                            [Arg(-1e3, 1e3)])

    def test_airybi_prime_complex(self):
        assert_mpmath_equal(lambda z: sc.airy(z)[3], lambda z:
                            mpmath.airybi(z, derivative=1),
                            [ComplexArg()])

    def test_bei(self):
        assert_mpmath_equal(sc.bei,
                            _exception_to_nan(lambda z: mpmath.bei(0, z, **HYPERKW)),
                            [Arg(-1e3, 1e3)])

    def test_ber(self):
        assert_mpmath_equal(sc.ber,
                            _exception_to_nan(lambda z: mpmath.ber(0, z, **HYPERKW)),
                            [Arg(-1e3, 1e3)])

    def test_bernoulli(self):
        assert_mpmath_equal(lambda n: sc.bernoulli(int(n))[int(n)],
                            lambda n: float(mpmath.bernoulli(int(n))),
                            [IntArg(0, 13000)],
                            rtol=1e-9, n=13000)

    def test_besseli(self):
        assert_mpmath_equal(sc.iv,
                            _exception_to_nan(lambda v, z: mpmath.besseli(v, z, **HYPERKW)),
                            [Arg(-1e100, 1e100), Arg()],
                            atol=1e-270)

    def test_besseli_complex(self):
        assert_mpmath_equal(lambda v, z: sc.iv(v.real, z),
                            _exception_to_nan(lambda v, z: mpmath.besseli(v, z, **HYPERKW)),
                            [Arg(-1e100, 1e100), ComplexArg()])

    def test_besselj(self):
        assert_mpmath_equal(sc.jv,
                            _exception_to_nan(lambda v, z: mpmath.besselj(v, z, **HYPERKW)),
                            [Arg(-1e100, 1e100), Arg(-1e3, 1e3)],
                            ignore_inf_sign=True)

        # loss of precision at large arguments due to oscillation
        assert_mpmath_equal(sc.jv,
                            _exception_to_nan(lambda v, z: mpmath.besselj(v, z, **HYPERKW)),
                            [Arg(-1e100, 1e100), Arg(-1e8, 1e8)],
                            ignore_inf_sign=True,
                            rtol=1e-5)

    def test_besselj_complex(self):
        assert_mpmath_equal(lambda v, z: sc.jv(v.real, z),
                            _exception_to_nan(lambda v, z: mpmath.besselj(v, z, **HYPERKW)),
                            [Arg(), ComplexArg()])

    def test_besselk(self):
        def mpbesselk(v, x):
            r = float(mpmath.besselk(v, x, **HYPERKW))
            if abs(r) > 1e305:
                # overflowing to inf a bit earlier is OK
                r = np.inf * np.sign(r)
            if abs(v) == abs(x) and abs(r) == np.inf and abs(x) > 1:
                # wrong result (kv(x,x) -> 0 for x > 1),
                # try with higher dps
                old_dps = mpmath.mp.dps
                mpmath.mp.dps = 200
                try:
                    r = float(mpmath.besselk(v, x, **HYPERKW))
                finally:
                    mpmath.mp.dps = old_dps
            return r
        assert_mpmath_equal(sc.kv,
                            _exception_to_nan(mpbesselk),
                            [Arg(-1e100, 1e100), Arg()])

    def test_besselk_int(self):
        assert_mpmath_equal(sc.kn,
                            _exception_to_nan(lambda v, z: mpmath.besselk(v, z, **HYPERKW)),
                            [IntArg(-1000, 1000), Arg()])

    def test_besselk_complex(self):
        assert_mpmath_equal(lambda v, z: sc.kv(v.real, z),
                            _exception_to_nan(lambda v, z: mpmath.besselk(v, z, **HYPERKW)),
                            [Arg(-1e100, 1e100), ComplexArg()])

    def test_bessely(self):
        def mpbessely(v, x):
            r = float(mpmath.bessely(v, x, **HYPERKW))
            if abs(r) > 1e305:
                # overflowing to inf a bit earlier is OK
                r = np.inf * np.sign(r)
            if abs(r) == 0 and x == 0:
                # invalid result from mpmath, point x=0 is a divergence
                return np.nan
            return r
        assert_mpmath_equal(sc.yv,
                            _exception_to_nan(mpbessely),
                            [Arg(-1e100, 1e100), Arg(-1e8, 1e8)],
                            n=5000)

    def test_bessely_complex(self):
        def mpbessely(v, x):
            r = complex(mpmath.bessely(v, x, **HYPERKW))
            if abs(r) > 1e305:
                # overflowing to inf a bit earlier is OK
                olderr = np.seterr(invalid='ignore')
                try:
                    r = np.inf * np.sign(r)
                finally:
                    np.seterr(**olderr)
            return r
        assert_mpmath_equal(lambda v, z: sc.yv(v.real, z),
                            _exception_to_nan(mpbessely),
                            [Arg(), ComplexArg()],
                            n=15000)

    def test_bessely_int(self):
        def mpbessely(v, x):
            r = float(mpmath.bessely(v, x))
            if abs(r) == 0 and x == 0:
                # invalid result from mpmath, point x=0 is a divergence
                return np.nan
            return r
        assert_mpmath_equal(lambda v, z: sc.yn(int(v), z),
                            _exception_to_nan(mpbessely),
                            [IntArg(-1000, 1000), Arg(-1e8, 1e8)])

    def test_beta(self):
        bad_points = []

        def beta(a, b, nonzero=False):
            if a < -1e12 or b < -1e12:
                # Function is defined here only at integers, but due
                # to loss of precision this is numerically
                # ill-defined. Don't compare values here.
                return np.nan
            if (a < 0 or b < 0) and (abs(float(a + b)) % 1) == 0:
                # close to a zero of the function: mpmath and scipy
                # will not round here the same, so the test needs to be
                # run with an absolute tolerance
                if nonzero:
                    bad_points.append((float(a), float(b)))
                    return np.nan
            return mpmath.beta(a, b)

        assert_mpmath_equal(sc.beta,
                            lambda a, b: beta(a, b, nonzero=True),
                            [Arg(), Arg()],
                            dps=400,
                            ignore_inf_sign=True)

        assert_mpmath_equal(sc.beta,
                            beta,
                            np.array(bad_points),
                            dps=400,
                            ignore_inf_sign=True,
                            atol=1e-14)

    def test_betainc(self):
        assert_mpmath_equal(sc.betainc,
                            _time_limited()(_exception_to_nan(lambda a, b, x: mpmath.betainc(a, b, 0, x, regularized=True))),
                            [Arg(), Arg(), Arg()])

    def test_binom(self):
        bad_points = []

        def binomial(n, k, nonzero=False):
            if abs(k) > 1e8*(abs(n) + 1):
                # The binomial is rapidly oscillating in this region,
                # and the function is numerically ill-defined. Don't
                # compare values here.
                return np.nan
            if n < k and abs(float(n-k) - np.round(float(n-k))) < 1e-15:
                # close to a zero of the function: mpmath and scipy
                # will not round here the same, so the test needs to be
                # run with an absolute tolerance
                if nonzero:
                    bad_points.append((float(n), float(k)))
                    return np.nan
            return mpmath.binomial(n, k)

        assert_mpmath_equal(sc.binom,
                            lambda n, k: binomial(n, k, nonzero=True),
                            [Arg(), Arg()],
                            dps=400)

        assert_mpmath_equal(sc.binom,
                            binomial,
                            np.array(bad_points),
                            dps=400,
                            atol=1e-14)

    def test_chebyt_int(self):
        assert_mpmath_equal(lambda n, x: sc.eval_chebyt(int(n), x),
                            _exception_to_nan(lambda n, x: mpmath.chebyt(n, x, **HYPERKW)),
                            [IntArg(), Arg()], dps=50)

    @knownfailure_overridable("some cases in hyp2f1 not fully accurate")
    def test_chebyt(self):
        assert_mpmath_equal(sc.eval_chebyt,
                            lambda n, x: _time_limited()(_exception_to_nan(mpmath.chebyt))(n, x, **HYPERKW),
                            [Arg(-101, 101), Arg()], n=10000)

    def test_chebyu_int(self):
        assert_mpmath_equal(lambda n, x: sc.eval_chebyu(int(n), x),
                            _exception_to_nan(lambda n, x: mpmath.chebyu(n, x, **HYPERKW)),
                            [IntArg(), Arg()], dps=50)

    @knownfailure_overridable("some cases in hyp2f1 not fully accurate")
    def test_chebyu(self):
        assert_mpmath_equal(sc.eval_chebyu,
                            lambda n, x: _time_limited()(_exception_to_nan(mpmath.chebyu))(n, x, **HYPERKW),
                            [Arg(-101, 101), Arg()])

    def test_chi(self):
        def chi(x):
            return sc.shichi(x)[1]
        assert_mpmath_equal(chi, mpmath.chi, [Arg()])
        # check asymptotic series cross-over
        assert_mpmath_equal(chi, mpmath.chi, [FixedArg([88 - 1e-9, 88, 88 + 1e-9])])

    def test_ci(self):
        def ci(x):
            return sc.sici(x)[1]
        # oscillating function: limit range
        assert_mpmath_equal(ci,
                            mpmath.ci,
                            [Arg(-1e8, 1e8)])

    def test_digamma(self):
        assert_mpmath_equal(sc.digamma,
                            _exception_to_nan(mpmath.digamma),
                            [Arg()],
                            dps=50)

    @knownfailure_overridable()
    def test_digamma_complex(self):
        assert_mpmath_equal(sc.digamma,
                            _time_limited()(_exception_to_nan(mpmath.digamma)),
                            [ComplexArg()],
                            n=200)

    def test_e1(self):
        assert_mpmath_equal(sc.exp1,
                            mpmath.e1,
                            [Arg()])

    def test_e1_complex(self):
        # E_1 oscillates as Im[z] -> +- inf, so limit range
        assert_mpmath_equal(sc.exp1,
                            mpmath.e1,
                            [ComplexArg(complex(-np.inf, -1e8), complex(np.inf, 1e8))],
                            rtol=1e-11)

        # Check cross-over reqion
        assert_mpmath_equal(sc.exp1,
                            mpmath.e1,
                            (np.linspace(-50, 50, 171)[:,None]
                             + np.r_[0, np.logspace(-3, 2, 61),
                                       -np.logspace(-3, 2, 11)]*1j
                             ).ravel(),
                            rtol=1e-11)
        assert_mpmath_equal(sc.exp1,
                            mpmath.e1,
                            (np.linspace(-50, -35, 10000) + 0j),
                            rtol=1e-11)

    def test_ei(self):
        assert_mpmath_equal(sc.expi,
                            mpmath.ei,
                            [Arg()],
                            rtol=1e-11)

    def test_ei_complex(self):
        # Ei oscillates as Im[z] -> +- inf, so limit range
        assert_mpmath_equal(sc.expi,
                            mpmath.ei,
                            [ComplexArg(complex(-np.inf, -1e8), complex(np.inf, 1e8))],
                            rtol=1e-9)

    def test_ellipe(self):
        assert_mpmath_equal(sc.ellipe,
                            mpmath.ellipe,
                            [Arg()])

    @knownfailure_overridable("insufficient accuracy from Cephes at phi < 0, or for extremely large m")
    def test_ellipf(self):
        assert_mpmath_equal(sc.ellipkinc,
                            mpmath.ellipf,
                            [Arg(), Arg(b=1.0)])

    def test_ellipk(self):
        assert_mpmath_equal(sc.ellipk,
                            mpmath.ellipk,
                            [Arg(b=1.0)])
        assert_mpmath_equal(sc.ellipkm1,
                            lambda m: mpmath.ellipk(1 - m),
                            [Arg(a=0.0)],
                            dps=400)

    def test_ellipfun_sn(self):
        # Oscillating function --- limit range of first argument; the
        # loss of precision there is an expected numerical feature
        # rather than an actual bug
        assert_mpmath_equal(lambda u, m: sc.ellipj(u, m)[0],
                            lambda u, m: mpmath.ellipfun("sn", u=u, m=m),
                            [Arg(-1e6, 1e6), Arg(a=0, b=1)],
                            atol=1e-20)

    def test_ellipfun_cn(self):
        # see comment in ellipfun_sn
        assert_mpmath_equal(lambda u, m: sc.ellipj(u, m)[1],
                            lambda u, m: mpmath.ellipfun("cn", u=u, m=m),
                            [Arg(-1e6, 1e6), Arg(a=0, b=1)],
                            atol=1e-20)

    def test_ellipfun_dn(self):
        # see comment in ellipfun_sn
        assert_mpmath_equal(lambda u, m: sc.ellipj(u, m)[2],
                            lambda u, m: mpmath.ellipfun("dn", u=u, m=m),
                            [Arg(-1e6, 1e6), Arg(a=0, b=1)],
                            atol=1e-20)

    def test_erf(self):
        assert_mpmath_equal(sc.erf,
                            lambda z: mpmath.erf(z),
                            [Arg()])

    def test_erf_complex(self):
        assert_mpmath_equal(sc.erf,
                            lambda z: mpmath.erf(z),
                            [ComplexArg()], n=200)

    def test_erfc(self):
        assert_mpmath_equal(sc.erfc,
                            _exception_to_nan(lambda z: mpmath.erfc(z)),
                            [Arg()])

    def test_erfc_complex(self):
        assert_mpmath_equal(sc.erfc,
                            _exception_to_nan(lambda z: mpmath.erfc(z)),
                            [ComplexArg()], n=200)

    def test_erfi(self):
        assert_mpmath_equal(sc.erfi,
                            mpmath.erfi,
                            [Arg()], n=200)

    def test_erfi_complex(self):
        assert_mpmath_equal(sc.erfi,
                            mpmath.erfi,
                            [ComplexArg()], n=200)

    def test_eulernum(self):
        assert_mpmath_equal(lambda n: sc.euler(n)[-1],
                            mpmath.eulernum,
                            [IntArg(1, 10000)], n=10000)

    @knownfailure_overridable("spurious(?) inf for negative x")
    def test_expint(self):
        assert_mpmath_equal(sc.expn,
                            _exception_to_nan(mpmath.expint),
                            [IntArg(0, 100), Arg()])

    def test_fresnels(self):
        def fresnels(x):
            return sc.fresnel(x)[0]
        assert_mpmath_equal(fresnels,
                            mpmath.fresnels,
                            [Arg()])

    def test_fresnelc(self):
        def fresnelc(x):
            return sc.fresnel(x)[1]
        assert_mpmath_equal(fresnelc,
                            mpmath.fresnelc,
                            [Arg()])

    def test_gamma(self):
        assert_mpmath_equal(sc.gamma,
                            _exception_to_nan(mpmath.gamma),
                            [Arg()])

    @dec.knownfailureif(True, "BUG: special.gammainc(1e20, 1e20) never returns")
    def test_gammainc(self):
        assert_mpmath_equal(sc.gammainc,
                            _exception_to_nan(
                                lambda z, b: mpmath.gammainc(z, b=b)/mpmath.gamma(z)),
                            [Arg(a=0), Arg(a=0)])

    @knownfailure_overridable()
    def test_gegenbauer(self):
        assert_mpmath_equal(sc.eval_gegenbauer,
                            _exception_to_nan(mpmath.gegenbauer),
                            [Arg(-1e3, 1e3), Arg(), Arg()])

    def test_gegenbauer_int(self):
        # Redefine functions to deal with numerical + mpmath issues
        def gegenbauer(n, a, x):
            # Avoid overflow at large `a` (mpmath would need an even larger
            # dps to handle this correctly, so just skip this region)
            if abs(a) > 1e100:
                return np.nan

            # Deal with n=0, n=1 correctly; mpmath 0.17 doesn't do these
            # always correctly
            if n == 0:
                r = 1.0
            elif n == 1:
                r = 2*a*x
            else:
                r = mpmath.gegenbauer(n, a, x)

            # Mpmath 0.17 gives wrong results (spurious zero) in some cases, so
            # compute the value by perturbing the result
            if float(r) == 0 and n <= 1-a and a < -1 and float(a) == int(float(a)):
                r = mpmath.gegenbauer(n, a + mpmath.mpf('1e-50'), x)

            # Differing overflow thresholds in scipy vs. mpmath
            if abs(r) > 1e270:
                return np.inf
            return r

        def sc_gegenbauer(n, a, x):
            r = sc.eval_gegenbauer(int(n), a, x)
            # Differing overflow thresholds in scipy vs. mpmath
            if abs(r) > 1e270:
                return np.inf
            return r
        assert_mpmath_equal(sc_gegenbauer,
                            _exception_to_nan(gegenbauer),
                            [IntArg(0, 100), Arg(), Arg()],
                            n=40000, dps=100,
                            ignore_inf_sign=True)

        # Check the small-x expansion
        assert_mpmath_equal(sc_gegenbauer,
                            _exception_to_nan(gegenbauer),
                            [IntArg(0, 100), Arg(), FixedArg(np.logspace(-30, -4, 30))],
                            dps=100,
                            ignore_inf_sign=True)

    @knownfailure_overridable()
    def test_gegenbauer_complex(self):
        assert_mpmath_equal(lambda n, a, x: sc.eval_gegenbauer(int(n), a.real, x),
                            _exception_to_nan(mpmath.gegenbauer),
                            [IntArg(0, 100), Arg(), ComplexArg()])

    @nonfunctional_tooslow
    def test_gegenbauer_complex_general(self):
        assert_mpmath_equal(lambda n, a, x: sc.eval_gegenbauer(n.real, a.real, x),
                            _exception_to_nan(mpmath.gegenbauer),
                            [Arg(-1e3, 1e3), Arg(), ComplexArg()])

    def test_hankel1(self):
        assert_mpmath_equal(sc.hankel1,
                            _exception_to_nan(lambda v, x: mpmath.hankel1(v, x,
                                                                          **HYPERKW)),
                            [Arg(-1e20, 1e20), Arg()])

    def test_hankel2(self):
        assert_mpmath_equal(sc.hankel2,
                            _exception_to_nan(lambda v, x: mpmath.hankel2(v, x, **HYPERKW)),
                            [Arg(-1e20, 1e20), Arg()])

    @knownfailure_overridable("issues at intermediately large orders")
    def test_hermite(self):
        assert_mpmath_equal(lambda n, x: sc.eval_hermite(int(n), x),
                            _exception_to_nan(mpmath.hermite),
                            [IntArg(0, 10000), Arg()])

    # hurwitz: same as zeta

    @nonfunctional_tooslow
    def test_hyp0f1(self):
        assert_mpmath_equal(sc.hyp0f1,
                            _exception_to_nan(lambda a, x: mpmath.hyp0f1(a, x, **HYPERKW)),
                            [Arg(), Arg()])

    @nonfunctional_tooslow
    def test_hyp0f1_complex(self):
        assert_mpmath_equal(lambda a, z: sc.hyp0f1(a.real, z),
                            _exception_to_nan(lambda a, x: mpmath.hyp0f1(a, x, **HYPERKW)),
                            [Arg(), ComplexArg()])

    @knownfailure_overridable()
    def test_hyp1f1(self):
        assert_mpmath_equal(_inf_to_nan(sc.hyp1f1),
                            _exception_to_nan(lambda a, b, x: mpmath.hyp1f1(a, b, x, **HYPERKW)),
                            [Arg(-1e5, 1e5), Arg(-1e5, 1e5), Arg()],
                            n=2000)

    @knownfailure_overridable()
    def test_hyp1f1_complex(self):
        assert_mpmath_equal(_inf_to_nan(lambda a, b, x: sc.hyp1f1(a.real, b.real, x)),
                            _exception_to_nan(lambda a, b, x: mpmath.hyp1f1(a, b, x, **HYPERKW)),
                            [Arg(-1e3, 1e3), Arg(-1e3, 1e3), ComplexArg()],
                            n=2000)

    @knownfailure_overridable()
    def test_hyp1f2(self):
        def hyp1f2(a, b, c, x):
            v, err = sc.hyp1f2(a, b, c, x)
            if abs(err) > max(1, abs(v)) * 1e-7:
                return np.nan
            return v
        assert_mpmath_equal(hyp1f2,
                            _exception_to_nan(lambda a, b, c, x: mpmath.hyp1f2(a, b, c, x, **HYPERKW)),
                            [Arg(), Arg(), Arg(), Arg()],
                            n=20000)

    @knownfailure_overridable()
    def test_hyp2f0(self):
        def hyp2f0(a, b, x):
            v, err = sc.hyp2f0(a, b, x, 1)
            if abs(err) > max(1, abs(v)) * 1e-7:
                return np.nan
            return v
        assert_mpmath_equal(hyp2f0,
                            lambda a, b, x: _time_limited(0.1)(_exception_to_nan(_trace_args(mpmath.hyp2f0)))(
                                a, b, x, **HYPERKW),
                            [Arg(), Arg(), Arg()])

    @knownfailure_overridable("spurious inf (or inf with wrong sign) for some argument values")
    def test_hyp2f1(self):
        assert_mpmath_equal(sc.hyp2f1,
                            _exception_to_nan(lambda a, b, c, x: mpmath.hyp2f1(a, b, c, x, **HYPERKW)),
                            [Arg(), Arg(), Arg(), Arg()])

    @nonfunctional_tooslow
    def test_hyp2f1_complex(self):
        # Scipy's hyp2f1 seems to have performance and accuracy problems
        assert_mpmath_equal(lambda a, b, c, x: sc.hyp2f1(a.real, b.real, c.real, x),
                            _exception_to_nan(lambda a, b, c, x: mpmath.hyp2f1(a, b, c, x, **HYPERKW)),
                            [Arg(-1e2, 1e2), Arg(-1e2, 1e2), Arg(-1e2, 1e2), ComplexArg()],
                            n=10)

    @knownfailure_overridable()
    def test_hyperu(self):
        assert_mpmath_equal(sc.hyperu,
                            _exception_to_nan(lambda a, b, x: mpmath.hyperu(a, b, x, **HYPERKW)),
                            [Arg(), Arg(), Arg()])

    def test_j0(self):
        # The Bessel function at large arguments is j0(x) ~ cos(x + phi)/sqrt(x)
        # and at large arguments the phase of the cosine loses precision.
        #
        # This is numerically expected behavior, so we compare only up to
        # 1e8 = 1e15 * 1e-7
        assert_mpmath_equal(sc.j0,
                            mpmath.j0,
                            [Arg(-1e3, 1e3)])
        assert_mpmath_equal(sc.j0,
                            mpmath.j0,
                            [Arg(-1e8, 1e8)],
                            rtol=1e-5)

    def test_j1(self):
        # See comment in test_j0
        assert_mpmath_equal(sc.j1,
                            mpmath.j1,
                            [Arg(-1e3, 1e3)])
        assert_mpmath_equal(sc.j1,
                            mpmath.j1,
                            [Arg(-1e8, 1e8)],
                            rtol=1e-5)

    @knownfailure_overridable()
    def test_jacobi(self):
        assert_mpmath_equal(sc.eval_jacobi,
                            _exception_to_nan(lambda a, b, c, x: mpmath.jacobi(a, b, c, x, **HYPERKW)),
                            [Arg(), Arg(), Arg(), Arg()])
        assert_mpmath_equal(lambda n, b, c, x: sc.eval_jacobi(int(n), b, c, x),
                            _exception_to_nan(lambda a, b, c, x: mpmath.jacobi(a, b, c, x, **HYPERKW)),
                            [IntArg(), Arg(), Arg(), Arg()])

    def test_jacobi_int(self):
        # Redefine functions to deal with numerical + mpmath issues
        def jacobi(n, a, b, x):
            # Mpmath does not handle n=0 case always correctly
            if n == 0:
                return 1.0
            return mpmath.jacobi(n, a, b, x)
        assert_mpmath_equal(lambda n, a, b, x: sc.eval_jacobi(int(n), a, b, x),
                            lambda n, a, b, x: _exception_to_nan(jacobi)(n, a, b, x, **HYPERKW),
                            [IntArg(), Arg(), Arg(), Arg()],
                            n=20000, dps=50)

    def test_kei(self):
        def kei(x):
            if x == 0:
                # work around mpmath issue at x=0
                return -pi/4
            return _exception_to_nan(mpmath.kei)(0, x, **HYPERKW)
        assert_mpmath_equal(sc.kei,
                            kei,
                            [Arg(-1e30, 1e30)], n=1000)

    def test_ker(self):
        assert_mpmath_equal(sc.ker,
                            _exception_to_nan(lambda x: mpmath.ker(0, x, **HYPERKW)),
                            [Arg(-1e30, 1e30)], n=1000)

    @nonfunctional_tooslow
    def test_laguerre(self):
        assert_mpmath_equal(_trace_args(sc.eval_laguerre),
                            lambda n, x: _exception_to_nan(mpmath.laguerre)(n, x, **HYPERKW),
                            [Arg(), Arg()])

    def test_laguerre_int(self):
        assert_mpmath_equal(lambda n, x: sc.eval_laguerre(int(n), x),
                            lambda n, x: _exception_to_nan(mpmath.laguerre)(n, x, **HYPERKW),
                            [IntArg(), Arg()], n=20000)

    def test_lambertw(self):
        assert_mpmath_equal(lambda x, k: sc.lambertw(x, int(k)),
                            lambda x, k: mpmath.lambertw(x, int(k)),
                            [Arg(), IntArg(0, 10)])

    @nonfunctional_tooslow
    def test_legendre(self):
        assert_mpmath_equal(sc.eval_legendre,
                            mpmath.legendre,
                            [Arg(), Arg()])

    def test_legendre_int(self):
        assert_mpmath_equal(lambda n, x: sc.eval_legendre(int(n), x),
                            lambda n, x: _exception_to_nan(mpmath.legendre)(n, x, **HYPERKW),
                            [IntArg(), Arg()],
                            n=20000)

        # Check the small-x expansion
        assert_mpmath_equal(lambda n, x: sc.eval_legendre(int(n), x),
                            lambda n, x: _exception_to_nan(mpmath.legendre)(n, x, **HYPERKW),
                            [IntArg(), FixedArg(np.logspace(-30, -4, 20))])

    def test_legenp(self):
        def lpnm(n, m, z):
            try:
                v = sc.lpmn(m, n, z)[0][-1,-1]
            except ValueError:
                return np.nan
            if abs(v) > 1e306:
                # harmonize overflow to inf
                v = np.inf * np.sign(v.real)
            return v

        def lpnm_2(n, m, z):
            v = sc.lpmv(m, n, z)
            if abs(v) > 1e306:
                # harmonize overflow to inf
                v = np.inf * np.sign(v.real)
            return v

        def legenp(n, m, z):
            if (z == 1 or z == -1) and int(n) == n:
                # Special case (mpmath may give inf, we take the limit by
                # continuity)
                if m == 0:
                    if n < 0:
                        n = -n - 1
                    return mpmath.power(mpmath.sign(z), n)
                else:
                    return 0

            if abs(z) < 1e-15:
                # mpmath has bad performance here
                return np.nan

            typ = 2 if abs(z) < 1 else 3
            v = _exception_to_nan(mpmath.legenp)(n, m, z, type=typ)

            if abs(v) > 1e306:
                # harmonize overflow to inf
                v = mpmath.inf * mpmath.sign(v.real)

            return v

        assert_mpmath_equal(lpnm,
                            legenp,
                            [IntArg(-100, 100), IntArg(-100, 100), Arg()])

        assert_mpmath_equal(lpnm_2,
                            legenp,
                            [IntArg(-100, 100), Arg(-100, 100), Arg(-1, 1)])

    def test_legenp_complex_2(self):
        def clpnm(n, m, z):
            try:
                return sc.clpmn(m.real, n.real, z, type=2)[0][-1,-1]
            except ValueError:
                return np.nan

        def legenp(n, m, z):
            if abs(z) < 1e-15:
                # mpmath has bad performance here
                return np.nan
            return _exception_to_nan(mpmath.legenp)(int(n.real), int(m.real), z, type=2)

        # mpmath is quite slow here
        x = np.array([-2, -0.99, -0.5, 0, 1e-5, 0.5, 0.99, 20, 2e3])
        y = np.array([-1e3, -0.5, 0.5, 1.3])
        z = (x[:,None] + 1j*y[None,:]).ravel()

        assert_mpmath_equal(clpnm,
                            legenp,
                            [FixedArg([-2, -1, 0, 1, 2, 10]), FixedArg([-2, -1, 0, 1, 2, 10]), FixedArg(z)],
                            rtol=1e-6,
                            n=500)

    def test_legenp_complex_3(self):
        def clpnm(n, m, z):
            try:
                return sc.clpmn(m.real, n.real, z, type=3)[0][-1,-1]
            except ValueError:
                return np.nan

        def legenp(n, m, z):
            if abs(z) < 1e-15:
                # mpmath has bad performance here
                return np.nan
            return _exception_to_nan(mpmath.legenp)(int(n.real), int(m.real), z, type=3)

        # mpmath is quite slow here
        x = np.array([-2, -0.99, -0.5, 0, 1e-5, 0.5, 0.99, 20, 2e3])
        y = np.array([-1e3, -0.5, 0.5, 1.3])
        z = (x[:,None] + 1j*y[None,:]).ravel()

        assert_mpmath_equal(clpnm,
                            legenp,
                            [FixedArg([-2, -1, 0, 1, 2, 10]), FixedArg([-2, -1, 0, 1, 2, 10]), FixedArg(z)],
                            rtol=1e-6,
                            n=500)

    @knownfailure_overridable("apparently picks wrong function at |z| > 1")
    def test_legenq(self):
        def lqnm(n, m, z):
            return sc.lqmn(m, n, z)[0][-1,-1]

        def legenq(n, m, z):
            if abs(z) < 1e-15:
                # mpmath has bad performance here
                return np.nan
            return _exception_to_nan(mpmath.legenq)(n, m, z, type=2)

        assert_mpmath_equal(lqnm,
                            legenq,
                            [IntArg(0, 100), IntArg(0, 100), Arg()])

    @nonfunctional_tooslow
    def test_legenq_complex(self):
        def lqnm(n, m, z):
            return sc.lqmn(int(m.real), int(n.real), z)[0][-1,-1]

        def legenq(n, m, z):
            if abs(z) < 1e-15:
                # mpmath has bad performance here
                return np.nan
            return _exception_to_nan(mpmath.legenq)(int(n.real), int(m.real), z, type=2)

        assert_mpmath_equal(lqnm,
                            legenq,
                            [IntArg(0, 100), IntArg(0, 100), ComplexArg()],
                            n=100)

    @knownfailure_overridable()
    def test_pcfd(self):
        def pcfd(v, x):
            return sc.pbdv(v, x)[0]
        assert_mpmath_equal(pcfd,
                            _exception_to_nan(lambda v, x: mpmath.pcfd(v, x, **HYPERKW)),
                            [Arg(), Arg()])

    @knownfailure_overridable("it's not the same as the mpmath function --- maybe different definition?")
    def test_pcfv(self):
        def pcfv(v, x):
            return sc.pbvv(v, x)[0]
        assert_mpmath_equal(pcfv,
                            lambda v, x: _time_limited()(_exception_to_nan(mpmath.pcfv))(v, x, **HYPERKW),
                            [Arg(), Arg()], n=1000)

    @knownfailure_overridable()
    def test_pcfw(self):
        def pcfw(a, x):
            return sc.pbwa(a, x)[0]
        assert_mpmath_equal(pcfw,
                            lambda v, x: _time_limited()(_exception_to_nan(mpmath.pcfw))(v, x, **HYPERKW),
                            [Arg(), Arg()], dps=50, n=1000)

    @knownfailure_overridable("issues at large arguments (atol OK, rtol not) and <eps-close to z=0")
    def test_polygamma(self):
        assert_mpmath_equal(sc.polygamma,
                            _time_limited()(_exception_to_nan(mpmath.polygamma)),
                            [IntArg(0, 1000), Arg()])

    def test_rgamma(self):
        def rgamma(x):
            if x < -8000:
                return np.inf
            else:
                v = mpmath.rgamma(x)
            return v
        assert_mpmath_equal(sc.rgamma,
                            rgamma,
                            [Arg()],
                            ignore_inf_sign=True)

    def test_rf(self):
        def mppoch(a, m):
            # deal with cases where the result in double precision
            # hits exactly a non-positive integer, but the
            # corresponding extended-precision mpf floats don't
            if float(a + m) == int(a + m) and float(a + m) <= 0:
                a = mpmath.mpf(a)
                m = int(a + m) - a
            return mpmath.rf(a, m)

        assert_mpmath_equal(sc.poch,
                            mppoch,
                            [Arg(), Arg()],
                            dps=400)

    def test_shi(self):
        def shi(x):
            return sc.shichi(x)[0]
        assert_mpmath_equal(shi, mpmath.shi, [Arg()])
        # check asymptotic series cross-over
        assert_mpmath_equal(shi, mpmath.shi, [FixedArg([88 - 1e-9, 88, 88 + 1e-9])])

    def test_si(self):
        def si(x):
            return sc.sici(x)[0]
        assert_mpmath_equal(si, mpmath.si, [Arg()])

    def test_spherharm(self):
        def spherharm(l, m, theta, phi):
            if m > l:
                return np.nan
            return sc.sph_harm(m, l, phi, theta)
        assert_mpmath_equal(spherharm,
                            mpmath.spherharm,
                            [IntArg(0, 100), IntArg(0, 100),
                             Arg(a=0, b=pi), Arg(a=0, b=2*pi)],
                            atol=1e-8, n=6000,
                            dps=150)

    def test_struveh(self):
        assert_mpmath_equal(sc.struve,
                            _exception_to_nan(mpmath.struveh),
                            [Arg(-1e4, 1e4), Arg(0, 1e4)],
                            rtol=5e-10)

    def test_struvel(self):
        def mp_struvel(v, z):
            if v < 0 and z < -v and abs(v) > 1000:
                # larger DPS needed for correct results
                old_dps = mpmath.mp.dps
                try:
                    mpmath.mp.dps = 300
                    return mpmath.struvel(v, z)
                finally:
                    mpmath.mp.dps = old_dps
            return mpmath.struvel(v, z)

        assert_mpmath_equal(sc.modstruve,
                            _exception_to_nan(mp_struvel),
                            [Arg(-1e4, 1e4), Arg(0, 1e4)],
                            rtol=5e-10,
                            ignore_inf_sign=True)

    def test_zeta(self):
        assert_mpmath_equal(sc.zeta,
                            _exception_to_nan(mpmath.zeta),
                            [Arg(a=1, b=1e10, inclusive_a=False),
                             Arg(a=0, inclusive_a=False)])

    def test_boxcox(self):

        def mp_boxcox(x, lmbda):
            x = mpmath.mp.mpf(x)
            lmbda = mpmath.mp.mpf(lmbda)
            if lmbda == 0:
                return mpmath.mp.log(x)
            else:
                return mpmath.mp.powm1(x, lmbda) / lmbda

        assert_mpmath_equal(sc.boxcox,
                            _exception_to_nan(mp_boxcox),
                            [Arg(a=0, inclusive_a=False), Arg()],
                            n=200,
                            dps=60,
                            rtol=1e-13)

    def test_boxcox1p(self):

        def mp_boxcox1p(x, lmbda):
            x = mpmath.mp.mpf(x)
            lmbda = mpmath.mp.mpf(lmbda)
            one = mpmath.mp.mpf(1)
            if lmbda == 0:
                return mpmath.mp.log(one + x)
            else:
                return mpmath.mp.powm1(one + x, lmbda) / lmbda

        assert_mpmath_equal(sc.boxcox1p,
                            _exception_to_nan(mp_boxcox1p),
                            [Arg(a=-1, inclusive_a=False), Arg()],
                            n=200,
                            dps=60,
                            rtol=1e-13)


if __name__ == "__main__":
    run_module_suite()
