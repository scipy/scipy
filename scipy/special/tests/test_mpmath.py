"""
Test Scipy functions versus mpmath, if available.

"""
from __future__ import division, print_function, absolute_import

import sys
import os
import time

from distutils.version import LooseVersion

import nose
import numpy as np
from numpy.testing import dec
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
                          rtol=1e-10)

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
        n1 = max(3, int(0.4*n))
        n2 = max(3, int(0.3*n))
        n3 = max(2, n - n1 - n2)

        v1 = np.linspace(-1, 1, n1)
        v2 = np.linspace(-10, 10, n2)
        if self.a >= 0 and self.b > 0:
            v3 = np.logspace(-30, np.log10(self.b), n3//2)
            v4 = np.logspace(-30, 5, n3//2)
        elif self.a < 0 and self.b > 0:
            v3 = np.r_[
                np.logspace(-30, np.log10(self.b), n3//4),
                -np.logspace(-30, np.log10(-self.a), n3//4)
                ]
            v4 = np.r_[
                np.logspace(-30, 5, n3//4),
                -np.logspace(-30, 5, n3//4)
                ]
        elif self.b < 0:
            v3 = -np.logspace(-30, np.log10(-self.b), n3//2)
            v3 = -np.logspace(-30, 5, n3//2)
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
                 dps=None, prec=None, n=5000, rtol=1e-7, atol=1e-300):
        self.scipy_func = scipy_func
        self.mpmath_func = mpmath_func
        self.arg_spec = arg_spec
        self.dps = dps
        self.prec = prec
        self.n = n
        self.rtol = rtol
        self.atol = atol
        self.is_complex = any([isinstance(arg, ComplexArg) for arg in self.arg_spec])
        if not name or name == '<lambda>':
            name = getattr(scipy_func, '__name__', None)
        if not name or name == '<lambda>':
            name =  getattr(mpmath_func, '__name__', None)
        self.name = name

    def check(self):
        np.random.seed(1234)

        # Generate values for the arguments
        num_args = len(self.arg_spec)
        m = int(self.n**(1./num_args)) + 1

        argvals = []
        for arg in self.arg_spec:
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
    def wrap(*a, **kw):
        sys.stderr.write("%r %r: " % (a, kw))
        sys.stderr.flush()
        r = func(*a, **kw)
        sys.stderr.write("-> %r\n" % r)
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
                    return None # turn off tracing except at function calls
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
    @knownfailure_overridable("accuracy issues at large arguments")
    def test_airyai(self):
        assert_mpmath_equal(lambda z: sc.airy(z)[0],
                            mpmath.airyai,
                            [Arg()])

    def test_airyai_complex(self):
        assert_mpmath_equal(lambda z: sc.airy(z)[0],
                            mpmath.airyai,
                            [ComplexArg()])

    @knownfailure_overridable("accuracy issues at large arguments")
    def test_airyai_prime(self):
        assert_mpmath_equal(lambda z: sc.airy(z)[1], lambda z:
                            mpmath.airyai(z, derivative=1),
                            [Arg()])

    def test_airyai_prime_complex(self):
        assert_mpmath_equal(lambda z: sc.airy(z)[1], lambda z:
                            mpmath.airyai(z, derivative=1),
                            [ComplexArg()])

    @knownfailure_overridable("accuracy issues at large arguments")
    def test_airybi(self):
        assert_mpmath_equal(lambda z: sc.airy(z)[2], lambda z:
                            mpmath.airybi(z),
                            [Arg()])

    def test_airybi_complex(self):
        assert_mpmath_equal(lambda z: sc.airy(z)[2], lambda z:
                            mpmath.airybi(z),
                            [ComplexArg()])

    @knownfailure_overridable("accuracy issues at large arguments")
    def test_airybi_prime(self):
        assert_mpmath_equal(lambda z: sc.airy(z)[3], lambda z:
                            mpmath.airybi(z, derivative=1),
                            [Arg()])

    def test_airybi_prime_complex(self):
        assert_mpmath_equal(lambda z: sc.airy(z)[3], lambda z:
                            mpmath.airybi(z, derivative=1),
                            [ComplexArg()])

    def test_bei(self):
        assert_mpmath_equal(sc.bei,
                            _exception_to_nan(lambda z: mpmath.bei(0, z, **HYPERKW)),
                            [Arg(-1e30, 1e30)])

    def test_ber(self):
        assert_mpmath_equal(sc.ber,
                            _exception_to_nan(lambda z: mpmath.ber(0, z, **HYPERKW)),
                            [Arg(-1e30, 1e30)])

    def test_bernoulli(self):
        assert_mpmath_equal(lambda n: sc.bernoulli(int(n))[int(n)],
                            lambda n: float(mpmath.bernoulli(int(n))),
                            [IntArg(0, 13000)],
                            rtol=1e-9, n=13000)

    def test_besseli(self):
        assert_mpmath_equal(sc.iv,
                            _exception_to_nan(lambda v, z: mpmath.besseli(v, z, **HYPERKW)),
                            [Arg(-1e100, 1e100), Arg()],
                            n=1000)

    def test_besseli_complex(self):
        assert_mpmath_equal(lambda v, z: sc.iv(v.real, z),
                            _exception_to_nan(lambda v, z: mpmath.besseli(v, z, **HYPERKW)),
                            [Arg(-1e100, 1e100), ComplexArg()])

    @knownfailure_overridable("Cephes jv at large arguments (absolute tolerance OK, relative not, Trac #1740)")
    def test_besselj(self):
        assert_mpmath_equal(sc.jv,
                            _exception_to_nan(lambda v, z: mpmath.besselj(v, z, **HYPERKW)),
                            [Arg(-1e100, 1e100), Arg()],
                            n=1000)

    def test_besselj_complex(self):
        assert_mpmath_equal(lambda v, z: sc.jv(v.real, z),
                            lambda v, z: mpmath.besselj(v, z, **HYPERKW),
                            [Arg(), ComplexArg()],
                            n=2000)


    def test_besselk(self):
        assert_mpmath_equal(sc.kv,
                            _exception_to_nan(lambda v, z: mpmath.besselk(v, z, **HYPERKW)),
                            [Arg(-1e100, 1e100), Arg()],
                            n=1000)

    @knownfailure_overridable("spurious inf at large arguments")
    def test_besselk_int(self):
        assert_mpmath_equal(sc.kn,
                            _exception_to_nan(lambda v, z: mpmath.besselk(v, z, **HYPERKW)),
                            [IntArg(), Arg()],
                            n=1000)

    def test_besselk_complex(self):
        assert_mpmath_equal(lambda v, z: sc.kv(v.real, z),
                            _exception_to_nan(lambda v, z: mpmath.besselk(v, z, **HYPERKW)),
                            [Arg(-1e100, 1e100), ComplexArg()])

    @knownfailure_overridable()
    def test_bessely(self):
        assert_mpmath_equal(sc.yv,
                            _exception_to_nan(lambda v, z: mpmath.bessely(v, z, **HYPERKW)),
                            [Arg(-1e100, 1e100), Arg()],
                            n=1000)

    def test_bessely_complex(self):
        assert_mpmath_equal(lambda v, z: sc.yv(v.real, z),
                            lambda v, z: _exception_to_nan(mpmath.bessely)(v, z, **HYPERKW),
                            [Arg(), ComplexArg()],
                            n=2000)

    @knownfailure_overridable()
    def test_beta(self):
        assert_mpmath_equal(sc.beta,
                            mpmath.beta,
                            [Arg(), Arg()],
                            dps=400)

    def test_betainc(self):
        assert_mpmath_equal(sc.betainc,
                            _exception_to_nan(lambda a, b, x: mpmath.betainc(a, b, 0, x, regularized=True)),
                            [Arg(), Arg(), Arg()])

    @knownfailure_overridable()
    def test_binom(self):
        assert_mpmath_equal(sc.binom,
                            mpmath.binomial,
                            [Arg(), Arg()],
                            dps=400)

    @knownfailure_overridable("issues at negative orders")
    def test_chebyt_int(self):
        assert_mpmath_equal(lambda n, x: sc.eval_chebyt(int(n), x),
                            _exception_to_nan(lambda n, x: mpmath.chebyt(n, x, **HYPERKW)),
                            [IntArg(), Arg()],
                            n=2000)

    @dec.knownfailureif(True, "test causes a segmentation fault")
    def test_chebyt(self):
        assert_mpmath_equal(sc.eval_chebyt,
                            _exception_to_nan(lambda n, x: mpmath.chebyt(n, x, **HYPERKW)),
                            [Arg(), Arg()],
                            n=2000)

    @knownfailure_overridable("issues at negative orders and at eps-small arguments")
    def test_chebyu_int(self):
        assert_mpmath_equal(sc.eval_chebyu,
                            _exception_to_nan(lambda n, x: mpmath.chebyu(n, x, **HYPERKW)),
                            [IntArg(), Arg()], n=2000)

    @dec.knownfailureif(True, "test causes a segmentation fault")
    def test_chebyu(self):
        assert_mpmath_equal(sc.eval_chebyu,
                            _exception_to_nan(lambda n, x: mpmath.chebyu(n, x, **HYPERKW)),
                            [Arg(), Arg()], n=2000)

    @knownfailure_overridable("overflows to inf somewhat too early at large arguments")
    def test_chi(self):
        def chi(x):
            return sc.shichi(x)[1]
        assert_mpmath_equal(chi,
                            mpmath.chi,
                            [Arg()])

    @knownfailure_overridable("Issues at extremely large arguments (absolute tolerance OK, relative not)")
    def test_ci(self):
        def ci(x):
            return sc.sici(x)[1]
        assert_mpmath_equal(ci,
                            mpmath.ci,
                            [Arg()])

    @knownfailure_overridable("Bad relative accuracy <eps-close to the pole at z=0")
    def test_digamma(self):
        assert_mpmath_equal(sc.digamma,
                            _exception_to_nan(mpmath.digamma),
                            [Arg()])

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

    @knownfailure_overridable("issues at very large arguments near imaginary axis (absolute tolerance OK, relative tolerance not)")
    def test_e1_complex(self):
        assert_mpmath_equal(sc.exp1,
                            mpmath.e1,
                            [ComplexArg()])

    def test_ei(self):
        assert_mpmath_equal(sc.expi,
                            mpmath.ei,
                            [Arg()])

    @knownfailure_overridable("picks the wrong branch of logarithm at real axis (Trac #1442)")
    def test_ei_complex(self):
        assert_mpmath_equal(sc.expi,
                            mpmath.ei,
                            [ComplexArg()])

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

    @knownfailure_overridable("issues at large arguments")
    def test_ellipfun_sn(self):
        assert_mpmath_equal(lambda u, m: sc.ellipj(u, m)[0],
                            lambda u, m: mpmath.ellipfun("sn", u=u, m=m),
                            [Arg(), Arg(a=0, b=1)])

    @knownfailure_overridable("issues at large arguments")
    def test_ellipfun_cn(self):
        assert_mpmath_equal(lambda u, m: sc.ellipj(u, m)[1],
                            lambda u, m: mpmath.ellipfun("cn", u=u, m=m),
                            [Arg(), Arg(a=0, b=1)])

    @knownfailure_overridable("issues at large arguments")
    def test_ellipfun_dn(self):
        assert_mpmath_equal(lambda u, m: sc.ellipj(u, m)[2],
                            lambda u, m: mpmath.ellipfun("dn", u=u, m=m),
                            [Arg(), Arg(a=0, b=1)])

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

    @knownfailure_overridable("behavior at large arguments (asymptotic expansion missing a sub-leading term)")
    def test_fresnels(self):
        def fresnels(x):
            return sc.fresnel(x)[0]
        assert_mpmath_equal(fresnels,
                            mpmath.fresnels,
                            [Arg()])

    @knownfailure_overridable("behavior at large arguments (asymptotic expansion missing a sub-leading term)")
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

    @knownfailure_overridable()
    def test_gegenbauer_int(self):
        assert_mpmath_equal(lambda n, a, x: sc.eval_gegenbauer(int(n), a, x),
                            _exception_to_nan(mpmath.gegenbauer),
                            [IntArg(0, 100), Arg(), Arg()])

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

    def test_hyp2f0(self):
        def hyp2f0(a, b, c, x):
            v, err = sc.hyp2f0(a, b, c, x)
            if abs(err) > max(1, abs(v)) * 1e-7:
                return np.nan
            return v
        assert_mpmath_equal(hyp2f0,
                            _exception_to_nan(lambda a, b, c, x: mpmath.hyp2f0(a, b, c, x, **HYPERKW)),
                            [Arg(), Arg(), Arg(), Arg()],
                            n=10000)

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

    @knownfailure_overridable("issues at large arguments")
    def test_j0(self):
        assert_mpmath_equal(sc.j0,
                            mpmath.j0,
                            [Arg()])

    @knownfailure_overridable("issues at large arguments")
    def test_j1(self):
        assert_mpmath_equal(sc.j1,
                            mpmath.j1,
                            [Arg()])

    @knownfailure_overridable()
    def test_jacobi(self):
        assert_mpmath_equal(sc.eval_jacobi,
                            _exception_to_nan(lambda a, b, c, x: mpmath.jacobi(a, b, c, x, **HYPERKW)),
                            [Arg(), Arg(), Arg(), Arg()])
        assert_mpmath_equal(lambda n, b, c, x: sc.eval_jacobi(int(n), b, c, x),
                            _exception_to_nan(lambda a, b, c, x: mpmath.jacobi(a, b, c, x, **HYPERKW)),
                            [IntArg(), Arg(), Arg(), Arg()])

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
        assert_mpmath_equal(sc.eval_laguerre,
                            mpmath.laguerre,
                            [Arg(), Arg(), Arg()])
        assert_mpmath_equal(sc.eval_laguerre,
                            mpmath.laguerre,
                            [IntArg(), Arg(), Arg()])

    def test_lambertw(self):
        assert_mpmath_equal(lambda x, k: sc.lambertw(x, int(k)),
                            lambda x, k: mpmath.lambertw(x, int(k)),
                            [Arg(), IntArg(0, 10)])

    @nonfunctional_tooslow
    def test_legendre(self):
        assert_mpmath_equal(sc.eval_legendre,
                            mpmath.legendre,
                            [Arg(), Arg()])
        assert_mpmath_equal(sc.eval_legendre,
                            mpmath.legendre,
                            [IntArg(), Arg()])

    @knownfailure_overridable("wrong function at |z| > 1 (Trac #1877)")
    def test_legenp(self):
        def lpnm(n, m, z):
            if m > n:
                return 0.0
            return sc.lpmn(m, n, z)[0][-1,-1]

        def lpnm_2(n, m, z):
            if m > n:
                return 0.0
            return sc.lpmv(m, n, z)

        def legenp(n, m, z):
            if abs(z) < 1e-15:
                # mpmath has bad performance here
                return np.nan
            return _exception_to_nan(mpmath.legenp)(n, m, z, **HYPERKW)

        assert_mpmath_equal(lpnm,
                            legenp,
                            [IntArg(0, 100), IntArg(0, 100), Arg()])

        assert_mpmath_equal(lpnm_2,
                            legenp,
                            [IntArg(0, 100), IntArg(0, 100), Arg(-1, 1)])

    @knownfailure_overridable("wrong function at |z| > 1 (Trac #1877)")
    def test_legenp_complex(self):
        def lpnm(n, m, z):
            if m > n:
                return 0.0
            return sc.lpmn(m, n, z)[0][-1,-1]

        def legenp(n, m, z):
            if abs(z) < 1e-15:
                # mpmath has bad performance here
                return np.nan
            return _exception_to_nan(mpmath.legenp)(int(n.real), int(m.real), **HYPERKW)

        assert_mpmath_equal(lpnm,
                            legenp,
                            [IntArg(0, 100), IntArg(0, 100), ComplexArg()])

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

    @knownfailure_overridable("all large negative arguments hit a pole --- the function itself is however numerically badly defined in this region, but maybe should return 0 instead?")
    def test_rgamma(self):
        assert_mpmath_equal(sc.rgamma,
                            mpmath.rgamma,
                            [Arg()])

    @knownfailure_overridable("invalid inf at very large negative arguments, accuracy issues at negative arguments eps-close to poles (some loss of precision in cancellation)")
    def test_rf(self):
        assert_mpmath_equal(sc.poch,
                            mpmath.rf,
                            [Arg(), Arg()], dps=100)

    @knownfailure_overridable("overflows to inf earlier than necessary")
    def test_shi(self):
        def shi(x):
            return sc.shichi(x)[0]
        assert_mpmath_equal(shi,
                            mpmath.shi,
                            [Arg()],
                            n=3000)

    def test_si(self):
        def si(x):
            return sc.sici(x)[0]
        assert_mpmath_equal(si,
                            mpmath.si,
                            [Arg()])

    def test_spherharm(self):
        def spherharm(l, m, theta, phi):
            if m > l:
                return np.nan
            return sc.sph_harm(m, l, phi, theta)
        assert_mpmath_equal(spherharm,
                            mpmath.spherharm,
                            [IntArg(0, 100), IntArg(0, 100),
                             Arg(a=0, b=2*pi), Arg(a=0, b=pi)],
                            atol=1e-13, n=6000)

    @knownfailure_overridable("problems at extremely large arguments (absolute tolerance OK, relative not)")
    def test_struve(self):
        assert_mpmath_equal(sc.struve,
                            _exception_to_nan(lambda v, x: mpmath.struveh(v, x, **HYPERKW)),
                            [Arg(-1e30, 1e30), Arg()],
                            n=2000)

    @knownfailure_overridable("accuracy problems at extremely large arguments (absolute tolerance OK, relative not)")
    def test_zeta(self):
        assert_mpmath_equal(sc.zeta,
                            _exception_to_nan(mpmath.zeta),
                            [Arg(a=1, inclusive_a=False),
                             Arg(a=0, inclusive_a=False)], n=1000)
