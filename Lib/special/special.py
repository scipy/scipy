from cephes import *
from Numeric import *
import types
from fastumath import *
from scipy.basic import squeeze
    
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
    def __init__(self,pyfunc,otypes=None):
        if not callable(pyfunc) or type(pyfunc) is types.ClassType:
            raise TypeError, "Object is not a callable Python object."
        self.thefunc = pyfunc
        self.__doc__ = pyfunc.__doc__
        if otypes is None:
            self.otypes=''
        else:
            if isinstance(otypes,types.StringType):
                self.otypes=otypes
            else:
                raise ValueError, "Output types must be a string."

    def __call__(self,*args):
        return squeeze(arraymap(self.thefunc,args,self.otypes))


def jv_prime(v,z):
    return (jv(v-1,z) - jv(v+1,z))/2.0

def yv_prime(v,z):
    return (yv(v-1,z) - yv(v+1,z))/2.0

def kv_prime(v,z):
    return (kv(v-1,z) - yv(v+1,z))/2.0

def iv_prime(v,z):
    return (iv(v-1,z) - iv(v+1,z))/2.0

def H1v_prime(v,z):
    return (hankel1(v-1,z)-hankel1(v+1,z))/2.0

def H2v_prime(v,z):
    return (hankel2(v-1,z)-hankel2(v+1,z))/2.0

def erfinv(y):
    return ndtri((y+1)/2.0)/sqrt(2)

def erfcinv(y):
    return ndtri((2-y)/2.0)/sqrt(2)

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
