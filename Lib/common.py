# Functions which are common and require SciPy Base and Level 1 SciPy
# (stats, special, linalg)

# Needs eigenvalues
import sys
import types
from scipy import special, stats, linalg
import Numeric

from scipy_base import exp, amin, amax, ravel, asarray, cast, arange, \
     ones, NewAxis, transpose, hstack, product, array, typename, where, \
     zeros, extract, insert, pi, sqrt, eye, poly1d, dot, r_

__all__ = ['factorial','factorial2','factorialk','comb','rand','randn','who',
           'lena','central_diff_weights', 'derivative', 'pade']
    
def factorial(n,exact=0):
    """n! = special.gamma(n+1)

    If exact==0, then floating point precision is used, otherwise
    exact long integer is computed.

    Notes:    
      - Array argument accepted only for exact=0 case.
      - If n<0, the return value is 0.
    """
    if exact:
        if n < 0:
            return 0L
        n = long(n)
        val = 1L
        k = 1L
        while (k < n+1L):
            val = val*k
            k += 1
        return val
    else:
        n = asarray(n)
        sv = special.errprint(0)
        vals = special.gamma(n+1)
        sv = special.errprint(sv)
        return where(n>=0,vals,0)


def factorial2(n,exact=0):
    """n!! = special.gamma(n/2+1)*2**((m+1)/2)/sqrt(pi)  n odd
           = 2**(n) * n!                                 n even

    If exact==0, then floating point precision is used, otherwise
    exact long integer is computed.

    Notes:    
      - Array argument accepted only for exact=0 case.
      - If n<0, the return value is 0.
    """
    if exact:
        if n < -1:
            return 0L
        if n <= 0:
            return 1L
        n = long(n)
        val = 1L
        k = n
        while (k > 0):
            val = val*k
            k -= 2
        return val
    else:
        n = asarray(n)
        vals = zeros(n.shape,'d')
        cond1 = (n % 2) & (n >= -1)
        cond2 = (1-(n % 2)) & (n >= -1)
        oddn = extract(cond1,n)
        evenn = extract(cond2,n)
        nd2o = oddn / 2.0
        nd2e = evenn / 2.0
        insert(vals,cond1,special.gamma(nd2o+1)/sqrt(pi)*pow(2.0,nd2o+0.5))
        insert(vals,cond2,special.gamma(nd2e+1) * pow(2.0,nd2e))
        return vals

def factorialk(n,k,exact=1):
    """n(!!...!)  = multifactorial of order k
        k times
    """
    if exact:
        if n < 1-k:
            return 0L
        if n<=0:
            return 1L
        n = long(n)
        val = 1L
        j = n
        while (j > 0):
            val = val*j
            j -= k
        return val
    else:
        raise NotImplementedError
        

def comb(N,k,exact=0):
    """Combinations of N things taken k at a time.

    If exact==0, then floating point precision is used, otherwise
    exact long integer is computed.

    Notes:    
      - Array arguments accepted only for exact=0 case.
      - If k > N, N < 0, or k < 0, then a 0 is returned.
    """
    if exact:
        if (k > N) or (N < 0) or (k < 0):
            return 0L
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
        k,N = asarray(k), asarray(N)
        lgam = special.gammaln
        cond = (k <= N) & (N >= 0) & (k >= 0)
        sv = special.errprint(0)
        vals = exp(lgam(N+1) - lgam(N-k+1) - lgam(k+1))
        sv = special.errprint(sv)
        return where(cond, vals, 0.0)

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

def pade(an, m):
    """Given Taylor series coefficients in an, return a Pade approximation to
    the function as the ratio of two polynomials p / q  where the order of q is m.
    """
    an = asarray(an)
    N = len(an) - 1
    n = N-m
    if (n < 0):
        raise ValueError, \
              "Order of q <m> must be smaller than len(an)-1."
    Akj = eye(N+1,n+1)
    Bkj = zeros((N+1,m),'d')
    for row in range(1,m+1):
        Bkj[row,:row] = -(an[:row])[::-1]
    for row in range(m+1,N+1):
        Bkj[row,:] = -(an[row-m:row])[::-1]
    C = hstack((Akj,Bkj))
    pq = dot(linalg.inv(C),an)
    p = pq[:n+1]
    q = r_[1.0,pq[n+1:]]
    return poly1d(p[::-1]), poly1d(q[::-1])




def rand(*args):
    """rand(d1,...,dn) returns a matrix of the given dimensions
    which is initialized to random numbers from a uniform distribution
    in the range [0,1).
    """
    return stats.random(args)

def randn(*args):
    """u = randn(d0,d1,...,dn) returns zero-mean, unit-variance Gaussian
    random numbers in an array of size (d0,d1,...,dn)."""
    return stats.norm.rvs(size=args)

def lena():
    import cPickle, os
    fname = os.path.join(os.path.dirname(__file__),'plt','lena.dat')
    f = open(fname,'rb')
    lena = array(cPickle.load(f))
    f.close()
    return lena


#-----------------------------------------------------------------------------
# Matlab like functions for output and information on the variables used.
#-----------------------------------------------------------------------------

def who(vardict=None):
    """Print the Numeric arrays in the given dictionary (or globals() if None).
    """
    if vardict is None:
        frame = sys._getframe().f_back
        vardict = frame.f_globals
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
    from scipy_test.testing import module_test
    module_test(__name__,__file__,level=level)

def test_suite(level=1):
    from scipy_test.testing import module_test_suite
    return module_test_suite(__name__,__file__,level=level)

