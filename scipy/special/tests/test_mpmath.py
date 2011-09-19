"""
Test Scipy functions versus mpmath, if available.

"""
import re
import numpy as np
from numpy.testing import dec
import scipy.special as sc

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
        return dec.skipif(True, "mpmath library is not present")

    def try_int(v):
        try: return int(v)
        except ValueError: return v

    def get_version(v):
        return map(try_int, re.split('[^0-9]', v))

    return dec.skipif(get_version(min_ver) > get_version(mpmath.__version__),
                      "mpmath %s required" % min_ver)


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
        FuncData(sc.hyp2f1, dataset, (0,1,2,3), 4, rtol=1e-9).check()
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
        x1, y1 = np.meshgrid(np.linspace(-10, 1, 11), np.linspace(-10, 1, 11))
        x2, y2 = np.meshgrid(np.logspace(-10, .8, 11), np.logspace(-10, .8, 11))
        points = np.r_[x1.ravel(),x2.ravel()] + 1j*np.r_[y1.ravel(),y2.ravel()]

        # note that the global accuracy of our complex erf algorithm is limited
        # roughly to 2e-8
        assert_func_equal(sc.erf, lambda x: complex(mpmath.erf(x)), points,
                          vectorized=False, rtol=2e-8)
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

    dataset = [p + (mpmath.legenp(p[1], p[0], p[2]),) for p in pts]
    dataset = np.array(dataset, dtype=np.float_)

    evf = lambda mu,nu,x: sc.lpmv(mu.astype(int), nu, x)

    olderr = np.seterr(invalid='ignore')
    try:
        FuncData(evf, dataset, (0,1,2), 3, rtol=1e-10, atol=1e-14).check()
    finally:
        np.seterr(**olderr)
