# Functions which are common and require SciPy Base and Level 1 SciPy
# (stats, special, linalg)

# Needs eigenvalues
import Numeric
import types, sys
from scipy import special, stats, linalg
from scipy_base import exp, amin, amax, ravel, asarray, cast, arange, \
     ones, NewAxis, transpose, hstack, product, array, typename
import scipy_base.fastumath

__all__ = ['factorial','comb','rand','randn','who','lena','central_diff_weights', 'derivative']
    
def factorial(n,exact=0):
    """n! = special.gamma(n+1)

    If exact==0, then floating point precision is used, otherwise
    exact long integer is computed."""
    if n < 0:
        raise ValueError, "n must be >= 0"
    if exact:
        n = long(n)
        val = 1L
        k = 1L
        while (k < n+1L):
            val = val*k
            k += 1
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
        N,k = map(long,(N,k))
        top = N
        val = 1L
        while (top > (N-k)):
            val *= top
            top -= 1
        n = 1L
        while (n < k+1L):
            val /= n
            n += 1
        return val
    else:
        lgam = special.gammaln
        return exp(lgam(N+1) - lgam(N-k+1) - lgam(k+1))

def central_diff_weights(Np,ndiv=1):
    """Return weights for an Np-point central derivative of order ndiv
       assuming equally-spaced function points.

       If weights are in the vector w, then 
       derivative is w[0] * f(x-ho*dx) + ... + w[-1] * f(x+h0*dx)
       
       Can be inaccurate for large number of points.
    """
    assert (Np >= ndiv+1), "Number of points must be at least the derivative order + 1."
    assert (Np % 2 == 1), "Odd-number of points only."
    ho = Np >> 1
    x = arange(-ho,ho+1.0)
    x = x[:,NewAxis]
    X = x**0.0
    for k in range(1,Np):
        X = hstack([X,x**k])
    w = product(arange(1,ndiv+1))*linalg.inv(X)[ndiv]
    return w

def derivative(func,x0,dx=1.0,n=1,args=(),order=3):
    """Given a function, use an N-point central differenece
       formula with spacing dx to compute the nth derivative at
       x0, where N is the value of order and must be odd.

       Warning: Decreasing the step size too small can result in
       round-off error.
    """
    assert (order >= n+1), "Number of points must be at least the derivative order + 1."
    assert (order % 2 == 1), "Odd number of points only."
    # pre-computed for n=1 and 2 and low-order for speed.
    if n==1:
        if order == 3:
            weights = array([-1,0,1])/2.0
        elif order == 5:
            weights = array([1,-8,0,8,-1])/12.0
        elif order == 7:
            weights = array([-1,9,-45,0,45,-9,1])/60.0
        elif order == 9:
            weights = array([3,-32,168,-672,0,672,-168,32,-3])/840.0
        else:
            weights = central_diff_weights(order,1)
    elif n==2:
        if order == 3:
            weights = array([1,-2.0,1])
        elif order == 5:
            weights = array([-1,16,-30,16,-1])/12.0
        elif order == 7:
            weights = array([2,-27,270,-490,270,-27,2])/180.0
        elif order == 9:
            weights = array([-9,128,-1008,8064,-14350,8064,-1008,128,-9])/5040.0
        else:
            weights = central_diff_weights(order,2)
    else:
        weights = central_diff_weights(order, n)
    val = 0.0
    ho = order >> 1
    for k in range(order):
        val += weights[k]*func(x0+(k-ho)*dx,*args)
    return val / product((dx,)*n)



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
            sta.append([namestr, shapestr, bytestr, typename(var.typecode()), original])

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

