# Functions which are common and require SciPy Base and Level 1 SciPy
# (stats, special, linalg)

# Needs eigenvalues
import Numeric
import types, sys
from scipy import special, stats
from scipy_base import exp, amin, amax, ravel, asarray, cast, arange, \
     ones, NewAxis, transpose
import scipy_base.fastumath

__all__ = ['factorial','comb','rand','randn','disp','who','lena','lena','lena','lena','lena','lena','lena','lena']
    
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
        return exp(lgam(N+1) - lgam(N-k+1) - lgam(k+1))


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


#-----------------------------------------------------------------------------
# Matlab like functions for output and information on the variables used.
#-----------------------------------------------------------------------------
def disp(mesg, device=None, linefeed=1):
    """Display a message to device (default is sys.stdout) with(out) linefeed.
    """
    if device is None:
        device = sys.stdout
    if linefeed:
        device.write('%s\n' % mesg)
    else:
        device.write('%s' % mesg)
    device.flush()
    return

def who(vardict=None):
    """Print the Numeric arrays in the given dictionary (or globals() if None).
    """
    if vardict is None:
        print "Pass in a dictionary:  who(globals())"
        return
    sta = []
    cache = {}
    for name in vardict.keys():
        if isinstance(vardict[name],Numeric.ArrayType):
            var = vardict[name]
            idv = id(var)
            if idv in cache.keys():
                namestr = name + " (%s)" % cache[idv]
                original=0
            else:
                cache[idv] = name
                namestr = name
                original=1
            shapestr = " x ".join(map(str, var.shape))
            bytestr = str(var.itemsize()*Numeric.product(var.shape))
            sta.append([namestr, shapestr, bytestr, _namefromtype[var.typecode()], original])

    maxname = 0
    maxshape = 0
    maxbyte = 0
    totalbytes = 0
    for k in range(len(sta)):
        val = sta[k]
        if maxname < len(val[0]):
            maxname = len(val[0])
        if maxshape < len(val[1]):
            maxshape = len(val[1])
        if maxbyte < len(val[2]):
            maxbyte = len(val[2])
        if val[4]:
            totalbytes += int(val[2])

    max = Numeric.maximum
    if len(sta) > 0:
        sp1 = max(10,maxname)
        sp2 = max(10,maxshape)
        sp3 = max(10,maxbyte)
        prval = "Name %s Shape %s Bytes %s Type" % (sp1*' ', sp2*' ', sp3*' ')
        print prval + "\n" + "="*(len(prval)+5) + "\n"
        
    for k in range(len(sta)):
        val = sta[k]
        print "%s %s %s %s %s %s %s" % (val[0], ' '*(sp1-len(val[0])+4),
                                        val[1], ' '*(sp2-len(val[1])+5),
                                        val[2], ' '*(sp3-len(val[2])+5),
                                        val[3])
    print "\nUpper bound on total bytes  =       %d" % totalbytes
    return
    

#-----------------------------------------------------------------------------

def test(level=10):
    from scipy_base.testing import module_test
    module_test(__name__,__file__,level=level)

def test_suite(level=1):
    from scipy_base.testing import module_test_suite
    return module_test_suite(__name__,__file__,level=level)

