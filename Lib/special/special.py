# 
# Author:  Travis Oliphant, 2002
#

from cephes import *
from Numeric import *
import types
from scipy_base.fastumath import *
from scipy_base import squeeze, isscalar
import specfun
    
class general_function:
    """
 general_function(somefunction)  Genearlized Function class.

  Description:
 
    Define a generalized function which takes nested sequence
    objects or Numeric arrays as inputs and returns a
    Numeric array as output, evaluating the function over successive
    tuples of the input arrays like the python map function except it uses
    the broadcasting rules of Numeric Python.

  Input:

    somefunction -- a Python function or method

  Example:

    def myfunc(a,b):
        if a > b:
            return a-b
        else
            return a+b

    gfunc = general_function(myfunc)

    >>> gfunc([1,2,3,4],2)
    array([3,4,1,2])

    """
    def __init__(self,pyfunc,otypes=None,doc=None):
        if not callable(pyfunc) or type(pyfunc) is types.ClassType:
            raise TypeError, "Object is not a callable Python object."
        self.thefunc = pyfunc
        if doc is None:
            self.__doc__ = pyfunc.__doc__
        else:
            self.__doc__ = doc
        if otypes is None:
            self.otypes=''
        else:
            if isinstance(otypes,types.StringType):
                self.otypes=otypes
            else:
                raise ValueError, "Output types must be a string."

    def __call__(self,*args):
        return squeeze(arraymap(self.thefunc,args,self.otypes))


def djv(v,z,n=1):
    """Return the nth derivative of Jv(z) with respect to z.
    """
    if n == 0:
        return jv(v,z)
    else:
        return (djv(v-1,z,n-1) - djv(v+1,z,n-1))/2.0

def dyv(v,z,n=1):
    """Return the nth derivative of Yv(z) with respect to z.
    """
    if n == 0:
        return yv(v,z)
    else:
        return (dyv(v-1,z,n-1) - dyv(v+1,z,n-1))/2.0

def dkv(v,z,n=1):
    """Return the nth derivative of Kv(z) with respect to z.
    """
    if n == 0:
        return kv(v,z)
    else:
        return (dkv(v-1,z,n-1) - dkv(v+1,z,n-1))/2.0

def div(v,z,n=1):
    """Return the nth derivative of Iv(z) with respect to z.
    """
    if n <= 0:
        return iv(v,z)
    else:
        return (div(v-1,z,n-1) - div(v+1,z,n-1))/2.0

def dh1v(v,z,n=1):
    """Return the nth derivative of H1v(z) with respect to z.
    """
    if n <= 0:
        return hankel1(v,z)
    else:
        return (dh1v(v-1,z,n-1) - dh1v(v+1,z,n-1))/2.0

def dh2v(v,z,n=1):
    """Return the nth derivative of H2v(z) with respect to z.
    """
    if n <= 0:
        return hankel2(v,z)
    else:
        return (dh2v(v-1,z,n-1) - dh2v(v+1,z,n-1))/2.0

def erfinv(y):
    return ndtri((y+1)/2.0)/sqrt(2)

def erfcinv(y):
    return ndtri((2-y)/2.0)/sqrt(2)

def gammaincinv(a,y):
    """returns the inverse of the incomplete gamma integral in that it
    finds x such that gammainc(a,x)=y
    """
    return gammainccinv(a,1-y)

def hyp0f1(v,z):
    """Confluent hypergeometric limit function 0F1.
    Limit as q->infinity of 1F1(q;a;z/q)
    """
    z = asarray(z)
    if z.typecode() in ['F', 'D']:
        arg = 2*sqrt(abs(z))
        num = where(z>=0, iv(v-1,arg), jv(v-1,arg))
        den = abs(z)**((v-1.0)/2)
    else:
        num = iv(v-1,2*sqrt(z))
        den = z**((v-1.0)/2.0)
    num *= gamma(v)
    return where(z==0,1.0,num/ asarray(den))

def assoc_laguerre(x,n,k=0.0):
    gam = gamma
    fac = gam(k+1+n)/gam(k+1)/gam(n+1)
    return fac*hyp1f1(-n,k+1,x)

digamma = psi

def polygamma(n, x):
    """Polygamma function which is the nth derivative of the digamma (psi)
    function."""
    n, x = asarray(n), asarray(x)
    cond = (n==0)
    fac2 = (-1.0)**(n+1) * gamma(n+1.0) * zeta(n+1,x)
    if sometrue(cond):
        return where(cond, psi(x), fac2)
    return fac2

def mathieu_A(m,q):
    """Compute expansion coefficients for even Mathieu functions and
    modified Mathieu functions.
    """
    if not (isscalar(m) and isscalar(q)):
        raise ValueError, "m and q must be scalars."
    if (q < 0):
        raise ValueError, "q >=0"
    if (m != floor(m)) or (m<0):
        raise ValueError, "m must be an integer >=0."

    if (q <= 1):
        qm = 7.5+56.1*sqrt(Q)-134.7*Q+90.7*sqrt(Q)*Q
    else:
        qm=17.0+3.1*sqrt(Q)-.126*Q+.0037*sqrt(Q)*Q
    km = int(qm+0.5*m)
    if km > 251:
        print "Warning, too many predicted coefficients."
    kd = 1
    m = int(floor(m))
    if m % 2:
        kd = 2

    a = mathieu_a(m,q)
    fc = specfunc.fcoef(kd,m,q,a)
    return fc[:km]

def mathieu_B(m,q):
    """Compute expansion coefficients for even Mathieu functions and
    modified Mathieu functions.
    """
    if not (isscalar(m) and isscalar(q)):
        raise ValueError, "m and q must be scalars."
    if (q < 0):
        raise ValueError, "q >=0"
    if (m != floor(m)) or (m<=0):
        raise ValueError, "m must be an integer > 0"

    if (q <= 1):
        qm = 7.5+56.1*sqrt(Q)-134.7*Q+90.7*sqrt(Q)*Q
    else:
        qm=17.0+3.1*sqrt(Q)-.126*Q+.0037*sqrt(Q)*Q
    km = int(qm+0.5*m)
    if km > 251:
        print "Warning, too many predicted coefficients."
    kd = 4
    m = int(floor(m))
    if m % 2:
        kd = 3

    b = mathieu_b(m,q)
    fc = specfunc.fcoef(kd,m,q,b)
    return fc[:km]

