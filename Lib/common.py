# Functions which are common and require SciPy Base and Level 1 SciPy
# (stats, special, linalg)

# Needs eigenvalues
import Numeric
import types
from scipy import special, stats

__all__ = ['factorial','comb','rand','randn']
    
def factorial(n,exact=0):
    """n! = special.gamma(n+1)

    If exact==0, then floating point precision is used, otherwise
    exact long integer is computed."""
    if n < 0:
        raise ValueError, "n must be >= 0"
    if exact:
        n = int(n)
        val = 1L
        for k in xrange(1,n+1):
            val = val*k
        return val
    else:
        return special.gamma(n+1)


def comb(N,k,exact=0):
    """Combinations of N things taken k at a time.

    If exact==0, then floating point precision is used, otherwise
    exact long integer is computed.
    """
    if (k > N) or (N < 0) or (k < 0):
        raise ValueError, "N and k must be non-negative and k <= N"
    if exact:
        N,k = map(int,(N,k))
        val = 1L
        for n in xrange(N-k+1,N+1):
            val = val*n
        for n in xrange(1,k+1):
            val = val / n
        return val
    else:
        lgam = special.gammaln
        return fastumath.exp(lgam(N+1) - lgam(N-k+1) - lgam(k+1))


def rand(*args):
    """rand(d1,...,dn) returns a matrix of the given dimensions
    which is initialized to random numbers from a uniform distribution
    in the range [0,1).
    """
    return stats.random(args)

def randn(*args):
    """u = randn(d0,d1,...,dn) returns zero-mean, unit-variance Gaussian
    random numbers in an array of size (d0,d1,...,dn)."""
    return stats.stnorm(size=args)

def lena():
    import cPickle, os
    d,junk = os.path.split(os.path.abspath(scipy.__file__))
    fname = os.path.join(d,'plt','lena.dat')
    f = open(fname,'rb')
    lena = scipy.array(cPickle.load(f))
    f.close()
    return lena


def test(level=10):
    from scipy_test import module_test
    module_test(__name__,__file__,level=level)

def test_suite(level=1):
    from scipy_test import module_test_suite
    return module_test_suite(__name__,__file__,level=level)

