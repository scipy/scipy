# Functions to implement several important functions for 
#   various Continous and Discrete Probability Distributions
#
# Author:  Travis Oliphant  2002-2003
# 

from __future__ import nested_scopes
import scipy
import scipy.special as special
import scipy.optimize as optimize
import Numeric as Num
import inspect
from Numeric import alltrue, where, arange, put, putmask, nonzero, \
     ravel, compress, take, ones, sum, shape, product, repeat, reshape, \
     zeros
from scipy_base.fastumath import *
from scipy_base import atleast_1d, polyval, angle, ceil, insert, extract, \
     any, argsort, argmax, argmin, vectorize, r_, asarray, nan, inf, select
import scipy_base

errp = special.errprint
arr = asarray
gam = special.gamma

import types, math
import stats as st
import rand
import rv


all = alltrue
sgf = vectorize
import new

def seed(x=0,y=0):
    """seed(x, y), set the seed using the integers x, y; 
    Set a random one from clock if  y == 0
    """
    if type (x) != types.IntType or type (y) != types.IntType :
        raise ArgumentError, "seed requires integer arguments."
    if y == 0:
        import random
        y = int(rv.initial_seed())
        x = random.randint(1,int(2l**31-2)) # Python 2.1 and 2.2 require this
                                            # bit of hackery to escape some 
                                            # wierd int/longint conversion
                                            # behavior
    rand.set_seeds(x,y)

seed()

def get_seed():
    "Return the current seed pair"
    return rand.get_seeds()

def _build_random_array(fun, args, size=None):
# Build an array by applying function fun to
# the arguments in args, creating an array with
# the specified shape.
# Allows an integer shape n as a shorthand for (n,).
    if isinstance(size, types.IntType): 
        size = [size]
    if size is not None and len(size) != 0:
        n = Num.multiply.reduce(size)
        s = apply(fun, args + (n,))
        s.shape = size
        return s
    else:
        n = 1
        s = apply(fun, args + (n,))
        return s[0]    

def random(size=None):
    "Returns array of random numbers between 0 and 1"
    return _build_random_array(rand.sample, (), size)

def random_integers(max, min=1, size=None):
    """random integers in range min to max inclusive"""
    return randint.rvs(min, max+1, size) 
     
def permutation(arg):
    """If arg is an integer, a permutation of indices arange(n), otherwise
    a permuation of the sequence"""
    if isinstance(arg,types.IntType):
        arg = Num.arange(arg)
    return rand.permutation(arg)


## Internal class to compute a ppf given a distribution.
##  (needs cdf function) and uses brentq from scipy.optimize
##  to compute ppf from cdf.
class general_cont_ppf:
    def __init__(self, dist, xa=-10.0, xb=10.0, xtol=1e-14):
        self.dist = dist
        self.cdf = eval('%scdf'%dist)
        self.xa = xa
        self.xb = xb
        self.xtol = xtol
        self.vecfunc = sgf(self._single_call,otypes='d')
    def _tosolve(self, x, q, *args):
        return apply(self.cdf, (x, )+args) - q
    def _single_call(self, q, *args):
        return scipy.optimize.brentq(self._tosolve, self.xa, self.xb, args=(q,)+args, xtol=self.xtol)
    def __call__(self, q, *args):
        return self.vecfunc(q, *args)


# Frozen RV class 
class rv_frozen:
    def __init__(self, dist, *args, **kwds):
        self.args = args
        self.kwds = kwds
        self.dist = dist
    def pdf(self,x):
        return self.dist.pdf(x,*self.args,**self.kwds)
    def cdf(self,x):
        return self.dist.cdf(x,*self.args,**self.kwds)
    def ppf(self,q):
        return self.dist.ppf(q,*self.args,**self.kwds)
    def isf(self,q):
        return self.dist.isf(q,*self.args,**self.kwds)
    def rvs(self, size=None):
        kwds = self.kwds
        kwds.update({'size':size})
        return self.dist.rvs(*self.args,**kwds)
    def sf(self,x):
        return self.dist.sf(x,*self.args,**self.kwds)
    def stats(self):
        return self.dist.stats(*self.args,**self.kwds)

    
            
##  NANs are returned for unsupported parameters.
##    location and scale parameters are optional for each distribution.
##    The shape parameters are generally required
##
##    The loc and scale parameters must be given as keyword parameters.
##    These are related to the common symbols in the .lyx file

##  skew is third central moment / variance**(1.5)
##  kurtosis is fourth central moment / variance**2 - 3


## References::

##  Documentation for ranlib, rv2, cdflib and
## 
##  Eric Wesstein's world of mathematics http://mathworld.wolfram.com/
##      http://mathworld.wolfram.com/topics/StatisticalDistributions.html
##
##  Documentation to Regress+ by Michael McLaughlin
##
##  Engineering and Statistics Handbook (NIST)
##      http://www.itl.nist.gov/div898/handbook/index.htm
##
##  Documentation for DATAPLOT from NIST
##      http://www.itl.nist.gov/div898/software/dataplot/distribu.htm
##
##  Norman Johnson, Samuel Kotz, and N. Balakrishnan "Continuous
##      Univariate Distributions", second edition,
##      Volumes I and II, Wiley & Sons, 1994.


## Each continuous random variable as the following methods
##
## rvs -- Random Variates (alternatively calling the class could produce these)
## pdf -- PDF
## cdf -- CDF
## sf  -- Survival Function (1-CDF)
## ppf -- Percent Point Function (Inverse of CDF)
## isf -- Inverse Survival Function (Inverse of SF)
## stats -- Return mean, variance, (Fisher's) skew, or (Fisher's) kurtosis
## nnlf  -- negative log likelihood function (to minimize)
## fit   -- Model-fitting
##
##  Maybe Later
##
##  hf  --- Hazard Function (PDF / SF)
##  chf  --- Cumulative hazard function (-log(SF))
##  psf --- Probability sparsity function (reciprocal of the pdf) in
##                units of percent-point-function (as a function of q).
##                Also, the derivative of the percent-point function.

## To define a new random variable you subclass the rv_continuous class
##   and re-define the 
##
##   _pdf method which will be given clean arguments (in between a and b)
##        and passing the argument check method
##
##      If postive argument checking is not correct for your RV
##      then you will also need to re-define
##   _argcheck

##   Correct, but potentially slow defaults exist for the remaining
##       methods but for speed and/or accuracy you can over-ride
##
##     _cdf, _ppf, _rvs, _isf, _sf
##
##   Rarely would you override _isf  and _sf but you could.
##
##   Statistics are computed using numerical integration by default.
##     For speed you can redefine this using
##
##    _stats  --- take shape parameters and return mu, mu2, g1, g2
##            --- If you can't compute one of these return it as None
##
##            --- Can also be defined with a keyword argument moments=<str>
##                  where <str> is a string composed of 'm', 'v', 's',
##                  and/or 'k'.  Only the components appearing in string
##                 should be computed and returned in the order 'm', 'v',
##                  's', or 'k'  with missing values returned as None
##
##    OR
##
##  You can override
##  
##    _munp    -- takes n and shape parameters and returns
##             --  then nth non-central moment of the distribution.
##       

def valarray(shape,value=nan,typecode=None):
    """Return an array of all value.
    """
    out = reshape(repeat([value],product(shape)),shape)
    if typecode is None:
        return out
    else:
        return out.astype(typecode)

def argsreduce(cond, *args):
    """Return a sequence of arguments converted to the dimensions of cond
    """
    newargs = list(args)
    expand_arr = (cond==cond)
    for k in range(len(args)):
        newargs[k] = extract(cond,arr(args[k])*expand_arr)
    return newargs    

class rv_continuous:
    """A Generic continuous random variable.

    Continuous random variables are defined from a standard form chosen
    for simplicity of representation.  The standard form may require
    some shape parameters to complete its specification.  The distributions
    also take optional location and scale parameters using loc= and scale=
    keywords (defaults: loc=0, scale=1)

    These shape, scale, and location parameters can be passed to any of the
    methods of the RV object such as the following:

    generic.rvs(<shape(s)>,loc=0,scale=1)
        - random variates 

    generic.pdf(x,<shape(s)>,loc=0,scale=1)
        - probability density function

    generic.cdf(x,<shape(s)>,loc=0,scale=1)
        - cumulative density function

    generic.sf(x,<shape(s)>,loc=0,scale=1)
        - survival function (1-cdf --- sometimes more accurate)

    generic.ppf(q,<shape(s)>,loc=0,scale=1)
        - percent point function (inverse of cdf --- percentiles)

    generic.isf(q,<shape(s)>,loc=0,scale=1)
        - inverse survival function (inverse of sf)

    generic.stats(<shape(s)>,loc=0,scale=1,moments='mv')
        - mean('m'), variance('v'), skew('s'), and/or kurtosis('k')

    generic.entropy(<shape(s)>,loc=0,scale=1)
        - (differential) entropy of the RV.

    Alternatively, the object may be called (as a function) to fix
       the shape, location, and scale parameters returning a
       "frozen" continuous RV object:

    myrv = generic(<shape(s)>,loc=0,scale=1)
        - frozen RV object with the same methods but holding the
            given shape, location, and scale fixed
    """
    def __init__(self, momtype=1, a=None, b=None, xa=-10.0, xb=10.0, xtol=1e-14, badvalue=None, name=None, longname=None, shapes=None, extradoc=None):
        if badvalue is None:
            badvalue = nan
        self.badvalue = badvalue
        self.name = name
        self.a = a
        self.b = b
        if a is None:
            self.a = -inf
        if b is None:
            self.b = inf
        self.xa = xa
        self.xb = xb
        self.xtol = xtol
        self._size = 1
        self.m = 0.0
        self.moment_type = momtype
        self.vecfunc = sgf(self._ppf_single_call,otypes='d')
        self.vecentropy = sgf(self._entropy,otypes='d')
        self.veccdf = sgf(self._cdf_single_call,otypes='d')
        self.expandarr = 1
        if momtype == 0:
            self.generic_moment = sgf(self._mom0_sc,otypes='d')
        else:
            self.generic_moment = sgf(self._mom1_sc,otypes='d')
        cdf_signature = inspect.getargspec(self._cdf.im_func)
        numargs1 = len(cdf_signature[0]) - 2
        pdf_signature = inspect.getargspec(self._pdf.im_func)
        numargs2 = len(pdf_signature[0]) - 2
        self.numargs = max(numargs1, numargs2)

        if longname is None:
            if name[0] in ['aeiouAEIOU']: hstr = "An "
            else: hstr = "A "
            longname = hstr + name
        if self.__doc__ is None:
            self.__doc__ = rv_continuous.__doc__
        if self.__doc__ is not None:
            self.__doc__ = self.__doc__.replace("A Generic",longname)
            if name is not None:
                self.__doc__ = self.__doc__.replace("generic",name)
            if shapes is None:
                self.__doc__ = self.__doc__.replace("<shape(s)>,","")
            else:
                self.__doc__ = self.__doc__.replace("<shape(s)>",shapes)
            if extradoc is not None:
                self.__doc__ = self.__doc__ + extradoc

    def _ppf_to_solve(self, x, q,*args):
        return apply(self.cdf, (x, )+args)-q
        
    def _ppf_single_call(self, q, *args):
        return scipy.optimize.brentq(self._ppf_to_solve, self.xa, self.xb, args=(q,)+args, xtol=self.xtol)

    # moment from definition
    def _mom_integ0(self, x,m,*args):
        return x**m * self.pdf(x,*args)
    def _mom0_sc(self, m,*args):
        return scipy.integrate.quad(self._mom_integ0, self.a,
                                    self.b, args=(m,)+args)[0]
    # moment calculated using ppf
    def _mom_integ1(self, q,m,*args):
        return (self.ppf(q,*args))**m
    def _mom1_sc(self, m,*args):
        return scipy.integrate.quad(self._mom_integ1, 0, 1,args=(m,)+args)[0]

    ## These are the methods you must define (standard form functions)
    def _argcheck(self, *args):
        # Default check for correct values on args and keywords.
        # Returns condition array of 1's where arguments are correct and
        #  0's where they are not.
        cond = 1
        for arg in args:
            cond = logical_and(cond,(arr(arg) > 0))
        return cond

    def _pdf(self,x,*args):
        return scipy.derivative(self._cdf,x,dx=1e-5,args=args,order=5)    

    ## Could also define any of these (return 1-d using self._size to get number)
    def _rvs(self, *args):
        ## Use basic inverse cdf algorithm for RV generation as default.
        U = rand.sample(self._size)
        Y = self._ppf(U,*args)
        return Y

    def _cdf_single_call(self, x, *args):
        return scipy.integrate.quad(self._pdf, self.a, x, args=args)[0]

    def _cdf(self, x, *args):
        return self.veccdf(x,*args)

    def _sf(self, x, *args):
        return 1.0-self._cdf(x,*args)

    def _ppf(self, q, *args):
        return self.vecfunc(q,*args)

    def _isf(self, q, *args):
        return self.vecfunc(1.0-q,*args)

    # The actual cacluation functions (no basic checking need be done)
    #  If these are defined, the others won't be looked at.
    #  Otherwise, the other set can be defined.
    def _stats(self,*args, **kwds):
        moments = kwds.get('moments')
        return None, None, None, None

    #  Central moments
    def _munp(self,n,*args):
        return self.generic_moment(n,*args)

    def __fix_loc_scale(self, args, loc, scale):
        N = len(args)
        if N > self.numargs:
            if N == self.numargs + 1 and loc is None:  # loc is given without keyword
                loc = args[-1]
            if N == self.numargs + 2 and scale is None: # loc and scale given without keyword
                loc, scale = args[-2:]
            args = args[:self.numargs]
        if scale is None:
            scale = 1.0
        if loc is None:
            loc = 0.0            
        return args, loc, scale

    # These are actually called, but should not
    #  be overwritten if you want to keep
    #  the error checking. 
    def rvs(self,*args,**kwds):
        """Random variates of given type.

        *args
        =====
        The shape parameter(s) for the distribution (see docstring of the
           instance object for more information)
        
        **kwds
        ======
        size  - number of random variates (default=1)
        loc   - location parameter (default=0)
        scale - scale parameter (default=1)
        """
        loc,scale,size=map(kwds.get,['loc','scale','size'])
        args, loc, scale = self.__fix_loc_scale(args, loc, scale)
        cond = logical_and(self._argcheck(*args),(scale >= 0))
        if not all(cond):
            raise ValueError, "Domain error in arguments."

        if size is None:
            size = 1
        else:
            self._size = product(size)
        if scipy.isscalar(size):
            self._size = size
            size = (size,)

        vals = reshape(self._rvs(*args),size)
        if scale == 0:
            return loc*ones(size,'d')
        else:
            return vals * scale + loc
        
    def pdf(self,x,*args,**kwds):
        """Probability density function at x of the given RV.

        *args
        =====
        The shape parameter(s) for the distribution (see docstring of the
           instance object for more information)
        
        **kwds
        ======
        loc   - location parameter (default=0)
        scale - scale parameter (default=1)
        """
        loc,scale=map(kwds.get,['loc','scale'])
        args, loc, scale = self.__fix_loc_scale(args, loc, scale)
        x,loc,scale = map(arr,(x,loc,scale))
        args = tuple(map(arr,args))
        x = arr((x-loc)*1.0/scale)
        cond0 = self._argcheck(*args) & (scale > 0)
        cond1 = (scale > 0) & (x > self.a) & (x < self.b)
        cond = cond0 & cond1
        output = zeros(shape(cond),'d')
        insert(output,(1-cond0)*(cond1==cond1),self.badvalue)
        goodargs = argsreduce(cond, *((x,)+args+(scale,)))
        scale, goodargs = goodargs[-1], goodargs[:-1]
        insert(output,cond,self._pdf(*goodargs) / scale)
        return output

    def cdf(self,x,*args,**kwds):
        """Cumulative distribution function at x of the given RV.

        *args
        =====
        The shape parameter(s) for the distribution (see docstring of the
           instance object for more information)

        **kwds
        ======
        loc   - location parameter (default=0)
        scale - scale parameter (default=1)
        """        
        loc,scale=map(kwds.get,['loc','scale'])
        args, loc, scale = self.__fix_loc_scale(args, loc, scale)
        x,loc,scale = map(arr,(x,loc,scale))
        args = tuple(map(arr,args))
        x = (x-loc)*1.0/scale
        cond0 = self._argcheck(*args) & (scale > 0)
        cond1 = (scale > 0) & (x > self.a) & (x < self.b)
        cond2 = (x >= self.b) & cond0
        cond = cond0 & cond1
        output = zeros(shape(cond),'d')
        insert(output,(1-cond0)*(cond1==cond1),self.badvalue)
        insert(output,cond2,1.0)
        goodargs = argsreduce(cond, *((x,)+args))
        insert(output,cond,self._cdf(*goodargs))
        return output

    def sf(self,x,*args,**kwds):
        """Survival function (1-cdf) at x of the given RV.

        *args
        =====
        The shape parameter(s) for the distribution (see docstring of the
           instance object for more information)

        **kwds
        ======
        loc   - location parameter (default=0)
        scale - scale parameter (default=1)
        """         
        loc,scale=map(kwds.get,['loc','scale'])
        args, loc, scale = self.__fix_loc_scale(args, loc, scale)
        x,loc,scale = map(arr,(x,loc,scale))
        args = tuple(map(arr,args))
        x = (x-loc)*1.0/scale
        cond0 = self._argcheck(*args) & (scale > 0)
        cond1 = (scale > 0) & (x > self.a) & (x < self.b)
        cond2 = cond0 & (x <= self.a)
        cond = cond0 & cond1
        output = zeros(shape(cond),'d')
        insert(output,(1-cond0)*(cond1==cond1),self.badvalue)
        insert(output,cond2,1.0)
        goodargs = argsreduce(cond, *((x,)+args))
        insert(output,cond,self._sf(*goodargs))
        return output

    def ppf(self,q,*args,**kwds):
        """Percent point function (inverse of cdf) at q of the given RV.

        *args
        =====
        The shape parameter(s) for the distribution (see docstring of the
           instance object for more information)

        **kwds
        ======
        loc   - location parameter (default=0)
        scale - scale parameter (default=1)
        """                 
        loc,scale=map(kwds.get,['loc','scale'])
        args, loc, scale = self.__fix_loc_scale(args, loc, scale)
        q,loc,scale = map(arr,(q,loc,scale))
        args = tuple(map(arr,args))
        cond0 = self._argcheck(*args) & (scale > 0) & (loc==loc)
        cond1 = (q > 0) & (q < 1)
        cond2 = (q==1) & cond0
        cond = cond0 & cond1
        output = valarray(shape(cond),value=self.a*scale + loc)
        insert(output,(1-cond0)+(1-cond1)*(q!=0.0), self.badvalue)
        insert(output,cond2,self.b*scale + loc)
        goodargs = argsreduce(cond, *((q,)+args+(scale,loc)))
        scale, loc, goodargs = goodargs[-2], goodargs[-1], goodargs[:-2]
        insert(output,cond,self._ppf(*goodargs)*scale + loc)
        return output
        
    def isf(self,q,*args,**kwds):
        """Inverse survival function at q of the given RV.

        *args
        =====
        The shape parameter(s) for the distribution (see docstring of the
           instance object for more information)

        **kwds
        ======
        loc   - location parameter (default=0)
        scale - scale parameter (default=1)
        """                         
        loc,scale=map(kwds.get,['loc','scale'])
        args, loc, scale = self.__fix_loc_scale(args, loc, scale)
        q,loc,scale = map(arr,(q,loc,scale))
        args = tuple(map(arr,args))
        cond0 = self._argcheck(*args) & (scale > 0) & (loc==loc)
        cond1 = (q > 0) & (q < 1)
        cond2 = (q==1) & cond0
        cond = cond0 & cond1
        output = valarray(shape(cond),value=self.b)
        insert(output,(1-cond0)*(cond1==cond1), self.badvalue)
        insert(output,cond2,self.a)
        goodargs = argsreduce(cond, *((1.0-q,)+args+(scale,loc)))
        scale, loc, goodargs = goodargs[-2], goodargs[-1], goodargs[:-2]
        insert(output,cond,self._ppf(*goodargs)*scale + loc)
        return output

    def stats(self,*args,**kwds):
        """Some statistics of the given RV

        *args
        =====
        The shape parameter(s) for the distribution (see docstring of the
           instance object for more information)

        **kwds
        ======
        loc     - location parameter (default=0)
        scale   - scale parameter (default=1)
        moments - a string composed of letters ['mvsk'] specifying
                   which moments to compute (default='mv')
                   'm' = mean,
                   'v' = variance,
                   's' = (Fisher's) skew,
                   'k' = (Fisher's) kurtosis.
        """                         
        loc,scale,moments=map(kwds.get,['loc','scale','moments'])

        N = len(args)
        if N > self.numargs:
            if N == self.numargs + 1 and loc is None:  # loc is given without keyword
                loc = args[-1]
            if N == self.numargs + 2 and scale is None: # loc and scale given without keyword
                loc, scale = args[-2:]
            if N == self.numargs + 3 and moments is None: # loc, scale, and moments
                loc, scale, moments = args[-3:]
            args = args[:self.numargs]
        if scale is None: scale = 1.0
        if loc is None: loc = 0.0
        if moments is None: moments = 'mv'
                        
        loc,scale = map(arr,(loc,scale))
        args = tuple(map(arr,args))
        cond = self._argcheck(*args) & (scale > 0) & (loc==loc)

        signature = inspect.getargspec(self._stats.im_func)
        if (signature[2] is not None) or ('moments' in signature[0]):
            mu, mu2, g1, g2 = self._stats(*args,**{'moments':moments})
        else:
            mu, mu2, g1, g2 = self._stats(*args)
        if g1 is None:
            mu3 = None
        else:
            mu3 = g1*(mu2**1.5)
        default = valarray(shape(cond), self.badvalue)
        output = []

        # Use only entries that are valid in calculation
        goodargs = argsreduce(cond, *(args+(scale,loc)))
        scale, loc, goodargs = goodargs[-2], goodargs[-1], goodargs[:-2]

        if 'm' in moments:
            if mu is None:
                mu = self._munp(1.0,*goodargs)
            out0 = default.copy()
            insert(out0,cond,mu*scale+loc)
            output.append(out0)
            
        if 'v' in moments:
            if mu2 is None:
                mu2p = self._munp(2.0,*goodargs)
                if mu is None:
                    mu = self._munp(1.0,*goodargs)
                mu2 = mu2p - mu*mu
            out0 = default.copy()
            insert(out0,cond,mu2*scale*scale)
            output.append(out0)
            
        if 's' in moments:
            if g1 is None:
                mu3p = self._munp(3.0,*goodargs)
                if mu is None:
                    mu = self._munp(1.0,*goodargs)                    
                if mu2 is None:
                    mu2p = self._munp(2.0,*goodargs)
                    mu2 = mu2p - mu*mu
                mu3 = mu3p - 3*mu*mu2 - mu**3
                g1 = mu3 / mu2**1.5
            out0 = default.copy()
            insert(out0,cond,g1)
            output.append(out0)
                
        if 'k' in moments:
            if g2 is None:
                mu4p = self._munp(4.0,*goodargs)
                if mu is None:
                    mu = self._munp(1.0,*goodargs)                    
                if mu2 is None:
                    mu2p = self._munp(2.0,*goodargs)
                    mu2 = mu2p - mu*mu
                if mu3 is None:
                    mu3p = self._munp(3.0,*goodargs)
                    mu3 = mu3p - 3*mu*mu2 - mu**3 
                mu4 = mu4p - 4*mu*mu3 - 6*mu*mu*mu2 - mu**4
                g2 = mu4 / mu2**2.0 - 3.0
            out0 = default.copy()
            insert(out0,cond,g2)
            output.append(out0)

        if len(output) == 1:
            return output[0]
        else:
            return tuple(output)

    def moment(self, n, *args):
        if (floor(n) != n):
            raise ValueError, "Moment must be an integer."
        if (n < 0): raise ValueError, "Moment must be positive."
        if (n == 0): return 1.0
        if (n > 0) and (n < 5):
            signature = inspect.getargspec(self._stats.im_func)
            if (signature[2] is not None) or ('moments' in signature[0]):
                dict = {'moments':{1:'m',2:'v',3:'vs',4:'vk'}[n]}
            else:
                dict = {}
            mu, mu2, g1, g2 = self._stats(*args,**dict)
            if (n==1):
                if mu is None: return self._munp(1,*args)
                else: return mu
            elif (n==2):
                if mu2 is None: return self._munp(2,*args)
                else: return mu
            elif (n==3):
                if g1 is None or mu2 is None: return self._munp(3,*args)
                else: return g1*(mu2**1.5)
            else: # (n==4)
                if g2 is None or mu2 is None: return self._munp(4,*args)
                else: return (g2+3.0)*(mu2**2.0)
        else:
            return self._munp(n,*args)

    def _nnlf(self, x, *args):
        return -sum(log(self._pdf(x, *args)))

    def nnlf(self, *args):
        # - sum (log pdf(x, theta))
        #   where theta are the parameters (including loc and scale)
        #
        try:
            x = args[-1]
            loc = args[-2]
            scale = args[-3]
            args = args[:-3]
        except IndexError:
            raise ValueError, "Not enough input arguments."
        if not self._argcheck(*args) or scale <= 0:
            return inf
        x = arr((x-loc) / scale)
        cond0 = (x <= self.a) | (x >= self.b)
        if (any(cond0)):
            return inf
        else:
            N = len(x)
            return self._nnlf(self, x, *args) + N*log(scale)

    def fit(self, data, *args, **kwds):
        loc0, scale0 = map(kwds.get, ['loc', 'scale'],[0.0, 1.0])
        Narg = len(args)
        if Narg != self.numargs:
            if Narg > self.numargs:
                raise ValueError, "Too many input arguments."
            else:
                args += (1.0,)*(self.numargs-Narg)
        # location and scale are at the end                
        x0 = args + (loc0, scale0)
        return optimize.fmin(self.nnlf,x0,args=(ravel(data),),disp=0)

    def est_loc_scale(self, data, *args):
        mu, mu2, g1, g2 = self.stats(*args,**{'moments':'mv'})
        muhat = stats.nanmean(data)
        mu2hat = stats.nanstd(data)
        Shat = sqrt(mu2hat / mu2)
        Lhat = muhat - Shat*mu
        return Lhat, Shat

    def freeze(self,*args,**kwds):
        return rv_frozen(self,*args,**kwds)
                
    def __call__(self, *args, **kwds):
        return self.freeze(*args, **kwds)

    def _entropy(self, *args):
        def integ(x):
            val = self._pdf(x, *args)
            return val*log(val)
        return -scipy.integrate.quad(integ,self.a,self.b)[0]

    def entropy(self, *args, **kwds):
        loc,scale=map(kwds.get,['loc','scale'])
        args, loc, scale = self.__fix_loc_scale(args, loc, scale)
        args = map(arr,args)
        cond0 = self._argcheck(*args) & (scale > 0) & (loc==loc)
        output = zeros(shape(cond0),'d')
        insert(output,(1-cond0),self.badvalue)
        goodargs = argsreduce(cond0, *args)
        insert(output,cond0,self.vecentropy(*goodargs)+log(scale))
        return output                

_EULER = 0.577215664901532860606512090082402431042  # -special.psi(1)
_ZETA3 = 1.202056903159594285399738161511449990765  # special.zeta(3,1)  Apery's constant

## Kolmogorov-Smirnov one-sided and two-sided test statistics

class ksone_gen(rv_continuous):
    def _cdf(self,x,n):
        return 1.0-special.smirnov(n,x)
    def _ppf(self,q,n):
        return special.smirnovi(n,1.0-q)
ksone = ksone_gen(a=0.0,name='ksone', longname="Kolmogorov-Smirnov "\
                  "A one-sided test statistic.", shapes="n",
                  extradoc="""

General Kolmogorov-Smirnov one-sided test.
"""
                  )

class kstwobign_gen(rv_continuous):
    def _cdf(self,x):
        return 1.0-special.kolmogorov(x)
    def _ppf(self,q):
        return special.kolmogi(1.0-q)
kstwobign = kstwobign_gen(a=0.0,name='kstwobign', longname='Kolmogorov-Smirnov two-sided (for large N)', extradoc="""

Kolmogorov-Smirnov two-sided test for large N
"""
                          )


## Normal distribution

# loc = mu, scale = std
class norm_gen(rv_continuous):
    def _rvs(self):
        return rand.standard_normal(self._size)
    def _pdf(self,x):
        return 1.0/sqrt(2*pi)*exp(-x**2/2.0)
    def _cdf(self,x):
        return special.ndtr(x)
    def _ppf(self,q):
        return special.ndtri(q)
    def _stats(self):
        return 0.0, 1.0, 0.0, 0.0
    def _entropy(self):
        return 0.5*(log(2*pi)+1)
norm = norm_gen(name='norm',longname='A normal',extradoc="""

Normal distribution

The location (loc) keyword specifies the mean.
The scale (scale) keyword specifies the standard deviation.
""")

## Alpha distribution
##
class alpha_gen(rv_continuous):
    def _pdf(self, x, a):
        return 1.0/arr(x**2)/special.ndtr(a)*norm.pdf(a-1.0/x)
    def _cdf(self, x, a):
        return special.ndtr(a-1.0/x) / special.ndtr(a)
    def _ppf(self, q, a):
        return 1.0/arr(a-special.ndtri(q*special.ndtr(a)))
    def _stats(self):
        return [inf]*2 + [nan]*2
alpha = alpha_gen(a=0.0,name='alpha',shapes='a',extradoc="""

Alpha distribution
""")

## Anglit distribution
##
class anglit_gen(rv_continuous):
    def _pdf(self, x):
        return cos(2*x)
    def _cdf(self, x):
        return sin(x+pi/4)**2.0
    def _ppf(self, q):
        return (arcsin(sqrt(q))-pi/4)
    def _stats(self):
        return 0.0, pi*pi/16-0.5, 0.0, -2*(pi**4 - 96)/(pi*pi-8)**2
    def _entropy(self):
        return 1-log(2)
anglit = anglit_gen(a=-pi/4,b=pi/4,name='anglit', extradoc="""

Anglit distribution

""")


## Arcsine distribution
##
class arcsine_gen(rv_continuous):
    def _pdf(self, x):
        return 1.0/pi/sqrt(x*(1-x))
    def _cdf(self, x):
        return 2.0/pi*arcsin(sqrt(x))
    def _ppf(self, q):
        return sin(pi/2.0*q)**2.0
    def _stats(self):
        mup = 0.5, 3.0/8.0, 15.0/48.0, 35.0/128.0
        mu = 0.5
        mu2 = 1.0/8
        g1 = 0
        g2 = -3.0/2.0
        return mu, mu2, g1, g2
    def _entropy(self):
        return -0.24156447527049044468
arcsine = arcsine_gen(a=0.0,b=1.0,name='arcsine',extradoc="""

Arcsine distribution
""")


## Beta distribution
##
class beta_gen(rv_continuous):
    def _rvs(self, a, b):
        return rand.beta(a,b,self._size)
    def _pdf(self, x, a, b):
        Px = (1.0-x)**(b-1.0) * x**(a-1.0)
        Px /= special.beta(a,b)
        return Px
    def _cdf(self, x, a, b):
        return special.btdtr(a,b,x)
    def _ppf(self, q, a, b):
        return special.btdtri(a,b,q)
    def _stats(self, a, b):
        mn = a *1.0 / (a + b)
        var = (a*b*1.0)*(a+b+1.0)/(a+b)**2.0
        g1 = 2.0*(b-a)*sqrt((1.0+a+b)/(a*b)) / (2+a+b)
        g2 = 6.0*(a**3 + a**2*(1-2*b) + b**2*(1+b) - 2*a*b*(2+b))
        g2 /= a*b*(a+b+2)*(a+b+3)
        return mn, var, g1, g2  
beta = beta_gen(a=0.0, b=1.0, name='beta',shapes='a,b',extradoc="""

Beta distribution
""")

## Beta Prime
class betaprime_gen(rv_continuous):
    def _rvs(self, a, b):
        u1 = gamma.rvs(a,size=self._size)
        u2 = gamma.rvs(b,size=self._size)
        return (u1 / u2)
    def _pdf(self, x, a, b):
        return 1.0/special.beta(a,b)*x**(a-1.0)/(1+x)**(a+b)
    def _cdf(self, x, a, b):
        x = where(x==1.0, 1.0-1e-6,x)
        return pow(x,a)*special.hyp2f1(a+b,a,1+a,-x)/a/special.beta(a,b)
    def _munp(self, n, a, b):
        if (n == 1.0):
            return where(b > 1, a/(b-1.0), inf)
        elif (n == 2.0):
            return where(b > 2, a*(a+1.0)/((b-2.0)*(b-1.0)), inf)
        elif (n == 3.0):
            return where(b > 3, a*(a+1.0)*(a+2.0)/((b-3.0)*(b-2.0)*(b-1.0)),
                         inf)
        elif (n == 4.0):
            return where(b > 4,
                         a*(a+1.0)*(a+2.0)*(a+3.0)/((b-4.0)*(b-3.0) \
                                                    *(b-2.0)*(b-1.0)), inf)
        else:
            raise NotImplementedError
betaprime = betaprime_gen(a=0.0, b=500.0, name='betaprime', shapes='a,b',
                          extradoc="""

Beta prime distribution
""")

## Bradford
##

class bradford_gen(rv_continuous):
    def _pdf(self, x, c):
        return  c / (c*x + 1.0) / log(1.0+c)
    def _cdf(self, x, c):
        return log(1.0+c*x) / log(c+1.0)
    def _ppf(self, q, c):
        return ((1.0+c)**q-1)/c
    def _stats(self, c, moments='mv'):
        k = log(1.0+c)
        mu = (c-k)/(c*k)
        mu2 = ((c+2.0)*k-2.0*c)/(2*c*k*k)
        g1 = None
        g2 = None
        if 's' in moments:
            g1 = sqrt(2)*(12*c*c-9*c*k*(c+2)+2*k*k*(c*(c+3)+3))
            g1 /= sqrt(c*(c*(k-2)+2*k))*(3*c*(k-2)+6*k)
        if 'k' in moments:
            g2 = c**3*(k-3)*(k*(3*k-16)+24)+12*k*c*c*(k-4)*(k-3) \
                 + 6*c*k*k*(3*k-14) + 12*k**3
            g2 /= 3*c*(c*(k-2)+2*k)**2
        return mu, mu2, g1, g2
    def _entropy(self, c):
        k = log(1+c)
        return k/2.0 - log(c/k)
bradford = bradford_gen(a=0.0, b=1.0, name='bradford', longname="A Bradford",
                        shapes='c', extradoc="""

Bradford distribution
""")


## Burr

# burr with d=1 is called the fisk distribution
class burr_gen(rv_continuous):
    def _pdf(self, x, c, d):
        return c*d*(x**(-c-1.0))*((1+x**(-c*1.0))**(-d-1.0))
    def _cdf(self, x, c, d):
        return (1+x**(-c*1.0))**(-d**1.0)
    def _ppf(self, q, c, d):
        return (q**(-1.0/d)-1)**(-1.0/c)
    def _stats(self, c, d, moments='mv'):
        g2c, g2cd = gam(1-2.0/c), gam(2.0/c+d)
        g1c, g1cd = gam(1-1.0/c), gam(1.0/c+d)
        gd = gam(d)
        k = gd*g2c*g2cd - g1c**2 * g1cd**2
        mu = g1c*g1cd / gd
        mu2 = k / gd**2.0
        g1, g2 = None, None
        g3c, g3cd = None, None
        if 's' in moments:
            g3c, g3cd = gam(1-3.0/c), gam(3.0/c+d)
            g1 = 2*g1c**3 * g1cd**3 + gd*gd*g3c*g3cd - 3*gd*g2c*g1c*g1cd*g2cd
            g1 /= sqrt(k**3)
        if 'k' in moments:
            if g3c is None:
                g3c = gam(1-3.0/c)
            if g3cd is None:
                g3cd = gam(3.0/c+d)
            g4c, g4cd = gam(1-4.0/c), gam(4.0/c+d)
            g2 = 6*gd*g2c*g2cd * g1c**2 * g1cd**2 + gd**3 * g4c*g4cd
            g2 -= 3*g1c**4 * g1cd**4 -4*gd**2*g3c*g1c*g1cd*g3cd
        return mu, mu2, g1, g2
burr = burr_gen(a=0.0, name='burr', longname="Burr",
                shapes="c,d", extradoc="""

Burr distribution
""")
    
# Fisk distribution
# burr is a generalization

class fisk_gen(burr_gen):
    def _pdf(self, x, c):
        return burr_gen._pdf(self, x, c, 1.0)
    def _cdf(self, x, c):
        return burr_gen._cdf(self, x, c, 1.0)
    def _ppf(self, x, c):
        return burr_gen._ppf(self, x, c, 1.0)
    def _stats(self, c):
        return burr_gen._stats(self, c, 1.0)
    def _entropy(self, c):
        return 2 - log(c)
fisk = fisk_gen(a=0.0, name='fink', longname="A funk",
                shapes='c', extradoc="""

Fink distribution.
"""
                )

## Cauchy

# median = loc

class cauchy_gen(rv_continuous):
    def _pdf(self, x):
        return 1.0/pi/(1.0+x*x)
    def _cdf(self, x):
        return 0.5 + 1.0/pi*arctan(x)
    def _ppf(self, q):
        return tan(pi*q-pi/2.0)
    def _sf(self, x):
        return 0.5 - 1.0/pi*arctan(x)
    def _isf(self, q):
        return tan(pi/2.0-pi*q)
    def _stats(self):
        return inf, inf, nan, nan
    def _entropy(self):
        return log(4*pi)
cauchy = cauchy_gen(name='cauchy',longname='Cauchy',extradoc="""

Cauchy distribution
"""
                    )

## Chi
##   (positive square-root of chi-square)
##   chi(1, loc, scale) = halfnormal
##   chi(2, 0, scale) = Rayleigh
##   chi(3, 0, scale) = MaxWell

class chi_gen(rv_continuous):
    def _rvs(self, df):
        return sqrt(chi2.rvs(df,size=self._size))
    def _pdf(self, x, df):
        return x**(df-1.)*exp(-x*x*0.5)/(2.0)**(df*0.5-1)/gam(df*0.5)
    def _cdf(self, x, df):
        return special.gammainc(df*0.5,0.5*x*x)
    def _ppf(self, q, df):
        return sqrt(2*special.gammaincinv(df*0.5,q))
    def _stats(self, df):
        mu = sqrt(2)*special.gamma(df/2.0+0.5)/special.gamma(df/2.0)
        mu2 = df - mu*mu
        g1 = (2*mu**3.0 + mu*(1-2*df))/arr(mu2**1.5)
        g2 = 2*df*(1.0-df)-6*mu**4 + 4*mu**2 * (2*df-1)
        g2 /= arr(mu2**2.0)
        return mu, mu2, g1, g2
chi = chi_gen(a=0.0,name='chi',shapes='df',extradoc="""

Chi distribution
"""
              )
    

## Chi-squared (gamma-distributed with loc=0 and scale=2 and shape=df/2)
class chi2_gen(rv_continuous):
    def _rvs(self, df):
        return rand.chi2(df,self._size)
    def _pdf(self, x, df):
        Px = x**(df/2.0-1)*exp(-x/2.0)
        Px /= special.gamma(df/2.0)* 2**(df/2.0)
        return Px
    def _cdf(self, x, df):
        return special.chdtr(df, x)
    def _sf(self, x, df):
        return special.chdtrc(df, x)
    def _isf(self, p, df):
        return special.chdtri(df, p)
    def _ppf(self, p, df):
        return self._isf(1.0-p, df)
    def _stats(self, df):
        mu = df
        mu2 = 2*df
        g1 = 2*sqrt(2.0/df)
        g2 = 12.0/df
        return mu, mu2, g1, g2
chi2 = chi2_gen(a=0.0,name='chi2',longname='A chi-squared',shapes='df',
                extradoc="""

Chi-squared distribution
"""
                )

## Cosine (Approximation to the Normal)
class cosine_gen(rv_continuous):
    def _pdf(self, x):
        return 1.0/2/pi*(1+cos(x))
    def _cdf(self, x):
        return 1.0/2/pi*(pi + x + sin(x))
    def _stats(self):
        return 0.0, pi*pi/3.0-2.0, 0.0, -6.0*(pi**4-90)/(5.0*(pi*pi-6)**2)
    def _entropy(self):
        return log(4*pi)-1.0
cosine = cosine_gen(a=-pi,b=pi,name='cosine',extradoc="""

Cosine distribution (Approximation to the normal)
""")

## Double Gamma distribution
class dgamma_gen(rv_continuous):
    def _rvs(self, a):
        u = random(size=self._size)
        return (gamma.rvs(a,size=self._size)*Num.where(u>=0.5,1,-1))
    def _pdf(self, x, a):
        ax = abs(x)
        return 1.0/(2*special.gamma(a))*ax**(a-1.0) * exp(-ax)
    def _cdf(self, x, a):
        fac = 0.5*special.gammainc(a,abs(x))
        return where(x>0,0.5+fac,0.5-fac)
    def _sf(self, x, a):
        fac = 0.5*special.gammainc(a,abs(x))
        return where(x>0,0.5-0.5*fac,0.5+0.5*fac)        
    def _ppf(self, q, a):
        fac = special.gammainccinv(a,1-abs(2*q-1))
        return where(q>0.5, fac, -fac)
    def _stats(self, a):
        mu2 = a*(a+1.0)
        return 0.0, mu2, 0.0, (a+2.0)*(a+3.0)/mu2-3.0
dgamma = dgamma_gen(name='dgamma',longname="A double gamma",
                    shapes='a',extradoc="""

Double gamma distribution
"""
                    )

## Double Weibull distribution
##
class dweibull_gen(rv_continuous):
    def _rvs(self, c):
        u = random(size=self._size)
        return weibull_min.rvs(c, size=self._size)*(Num.where(u>=0.5,1,-1))    
    def _pdf(self, x, c):
        ax = abs(x)
        Px = c/2.0*ax**(c-1.0)*exp(-ax**c)
        return Px
    def _cdf(self, x, c):
        Cx1 = 0.5*exp(-abs(x)**c)
        return where(x > 0, 1-Cx1, Cx1)
    def _ppf(self, q, c):
        fac = where(q<=0.5,2*q,2*q-1)
        fac = pow(arr(log(1.0/fac)),1.0/c)
        return where(q>0.5,fac,-fac)
    def _stats(self, c):
        var = gam(1+2.0/c)
        return 0.0, var, 0.0, gam(1+4.0/c)/var
dweibull = dweibull_gen(name='dweibull',longname="A double Weibull",
                        shapes='c',extradoc="""

Double Weibull distribution
"""
                        )

## ERLANG
##
## Special case of the Gamma distribution with shape parameter an integer.
##
class erlang_gen(rv_continuous):
    def _rvs(self, n):
        return gamma.rvs(n,size=self._size)
    def _arg_check(self, n):
        return (n > 0) & (floor(n)==n)
    def _pdf(self, x, n):
        Px = (x)**(n-1.0)*exp(-x)/special.gamma(n)
        return Px
    def _cdf(self, x, n):
        return special.gdtr(1.0,n,x)
    def _sf(self, x, n):
        return special.gdtrc(1.0,n,x)
    def _ppf(self, q, n):
        return special.gdtrix(1.0, n, q)
    def _stats(self, n):
        n = n*1.0
        return n, n, 2/sqrt(n), 6/n
    def _entropy(self, n):
        return special.psi(n)*(1-n) + 1 + special.gammaln(n)        
erlang = erlang_gen(a=0.0,name='erlang',longname='An Erlang',
                    shapes='n',extradoc="""

Erlang distribution (Gamma with integer shape parameter)
"""
                    )
       
## Exponential (gamma distributed with a=1.0, loc=loc and scale=scale)
## scale == 1.0 / lambda

class expon_gen(rv_continuous):
    def _rvs(self):
        return rand.standard_exp(self._size)
    def _pdf(self, x):
        return exp(-x)
    def _cdf(self, x):
        return 1.0-exp(-x)
    def _ppf(self, q):
        return -log(1.0-q)
    def _stats(self):
        return 1.0, 1.0, 2.0, 6.0
    def _entropy(self):
        return 1.0
expon = expon_gen(a=0.0,name='expon',longname="An exponential",
                  extradoc="""

Exponential distribution
"""
                  )


## Exponentiated Weibull
class exponweib_gen(rv_continuous):
    def _pdf(self, x, a, c):
        exc = exp(-x**c)
        return a*c*(1-exc)**arr(a-1) * exc * x**arr(c-1)
    def _cdf(self, x, a, c):
        exc = exp(-x**c)
        return arr((1-exc)**a)
    def _ppf(self, q, a, c):
        return (-log(1-q**(1.0/a)))**arr(1.0/c)
exponweib = exponweib_gen(a=0.0,name='exponweib',
                          longname="An exponentiated Weibull",
                          shapes="a,c",extradoc="""

Exponentiated Weibull distribution                          
"""
                          )
                   
## Exponential Power

class exponpow_gen(rv_continuous):
    def _pdf(self, x, b):
        xbm1 = arr(x**(b-1.0))
        xb = xbm1 * x
        return exp(1)*b*xbm1 * exp(xb - exp(xb))
    def _cdf(self, x, b):
        xb = arr(x**b)
        return 1.0-exp(1-exp(xb))
    def _ppf(self, q, b):
        return pow(log(1.0-log(1.0-q)), 1.0/b)
exponpow = exponpow_gen(a=0.0,name='exponpow',longname="An exponential power",
                        shapes='b',extradoc="""

Exponential Power distribution
"""
                        )

## Faigue-Life (Birnbaum-Sanders)
class fatiguelife_gen(rv_continuous):
    def _rvs(self, c):
        z = norm.rvs(size=self._size)
        U = random(size=self._size)
        fac = 2 + c*c*z*z
        det = sqrt(fac*fac - 4)
        t1 = fac + det
        t2 = fac - det
        return t1*(U>0.5) + t2*(U<0.5)        
    def _pdf(self, x, c):
        return (x+1)/arr(2*c*sqrt(2*pi*x**3))*exp(-(x-1)**2/arr((2.0*x*c**2)))
    def _cdf(self, x, c):
        return special.ndtr(1.0/c*(sqrt(x)-1.0/arr(sqrt(x))))
    def _ppf(self, q, c):
        tmp = c*special.ndtri(q)
        return 0.25*(tmp + sqrt(tmp**2 + 4))**2
    def _stats(self, c):
        c2 = c*c
        mu = c2 / 2.0 + 1
        den = 5*c2 + 4
        mu2 = c2*den /4.0
        g1 = 4*c*sqrt(11*c2+6.0)/den**1.5
        g2 = 6*c2*(93*c2+41.0) / den**2.0
        return mu, mu2, g1, g2
fatiguelife = fatiguelife_gen(a=0.0,name='fatiguelife',
                              longname="A fatigue-life (Birnbaum-Sanders)",
                              shapes='c',extradoc="""

Fatigue-life (Birnbaum-Sanders) distribution
"""
                              )

## Folded Cauchy

class foldcauchy_gen(rv_continuous):
    def _rvs(self, c):
        return abs(cauchy.rvs(loc=c,size=self._size))
    def _pdf(self, x, c):
        return 1.0/pi*(1.0/(1+(x-c)**2) + 1.0/(1+(x+c)**2))
    def _cdf(self, x, c):
        return 1.0/pi*(arctan(x-c) + arctan(x+c))
    def _stats(self, x, c):
        return inf, inf, nan, nan
foldcauchy = foldcauchy_gen(a=0.0, name='foldcauchy',
                            longname = "A folded Cauchy",
                            shapes='c',extradoc="""

A folded Cauchy distributions
"""
                            )
        
## F

class f_gen(rv_continuous):
    def _rvs(self, dfn, dfd):
        return rand.f(dfn, dfd, self._size)
    def _pdf(self, x, dfn, dfd):
        n = arr(1.0*dfn)
        m = arr(1.0*dfd)
        Px = m**(m/2) * n**(n/2) * x**(n/2-1)
        Px /= (m+n*x)**((n+m)/2)*special.beta(n/2,m/2)
        return Px
    def _cdf(self, x, dfn, dfd):
        return special.fdtr(dfn, dfd, x)
    def _sf(self, x, dfn, dfd):
        return special.fdtrc(dfn, dfd, x)
    def _ppf(self, q, dfn, dfd):
        return special.fdtri(dfn, dfd, q)        
    def _stats(self, dfn, dfd):
        v2 = arr(dfd*1.0)
        v1 = arr(dfn*1.0)
        mu = where (v2 > 2, v2 / arr(v2 - 2), inf)
        mu2 = 2*v2*v2*(v2+v1-2)/(v1*(v2-2)**2 * (v2-4))
        mu2 = where(v2 > 4, mu2, inf)
        g1 = 2*(v2+2*v1-2)/(v2-6)*sqrt((2*v2-4)/(v1*(v2+v1-2)))
        g1 = where(v2 > 6, g1, nan)
        g2 = 3/(2*v2-16)*(8+g1*g1*(v2-6))
        g2 = where(v2 > 8, g2, nan)
        return mu, mu2, g1, g2
f = f_gen(a=0.0,name='f',longname='An F',shapes="dfn,dfd",
          extradoc="""

F distribution
"""
          )

## Folded Normal  
##   abs(Z) where (Z is normal with mu=L and std=S so that c=abs(L)/S)
##
##  note: regress docs have scale parameter correct, but first parameter
##    he gives is a shape parameter A = c * scale

##  Half-normal is folded normal with shape-parameter c=0.

class foldnorm_gen(rv_continuous):
    def _rvs(self, c):
        return abs(norm.rvs(loc=c,size=self._size))
    def _pdf(self, x, c):
        return sqrt(2.0/pi)*cosh(c*x)*exp(-(x*x+c*c)/2.0)
    def _cdf(self, x, c,):
        return special.ndtr(x-c) + special.ndtr(x+c) - 1.0
    def _stats(self, c):
        fac = special.erf(c/sqrt(2))
        mu = sqrt(2.0/pi)*exp(-0.5*c*c)+c*fac
        mu2 = c*c + 1 - mu*mu
        c2 = c*c
        g1 = sqrt(2/pi)*exp(-1.5*c2)*(4-pi*exp(c2)*(2*c2+1.0))
        g1 += 2*c*fac*(6*exp(-c2) + 3*sqrt(2*pi)*c*exp(-c2/2.0)*fac + \
                       pi*c*(fac*fac-1))
        g1 /= pi*mu2**1.5
    
        g2 = c2*c2+6*c2+3+6*(c2+1)*mu*mu - 3*mu**4
        g2 -= 4*exp(-c2/2.0)*mu*(sqrt(2.0/pi)*(c2+2)+c*(c2+3)*exp(c2/2.0)*fac)
        g2 /= mu2**2.0
        return mu, mu2, g1, g2
foldnorm = foldnorm_gen(a=0.0,name='foldnorm',longname='A folded normal',
                        shapes='c',extradoc="""

Folded normal distribution
"""
                        )

## Extreme Value Type II or Frechet
## (defined in Regress+ documentation as Extreme LB) as
##   a limiting value distribution.
##
class frechet_r_gen(rv_continuous):
    def _pdf(self, x, c):
        return c*pow(x,c-1)*exp(-pow(x,c))
    def _cdf(self, x, c):
        return 1-exp(-pow(x,c))
    def _ppf(self, q, c):
        return pow(-log(1-q),1.0/c)
    def _munp(self, n, c):
        return special.gamma(1.0+n*1.0/c)
    def _entropy(self, c):
        return -_EULER / c - log(c) + _EULER + 1
frechet_r = frechet_r_gen(a=0.0,name='frechet_r',longname="A Frechet right",
                          shapes='c',extradoc="""

A Frechet (right) distribution (also called Weibull minimum)
"""
                          )
weibull_min = frechet_r_gen(a=0.0,name='weibull_min',
                            longname="A Weibull minimum",
                            shapes='c',extradoc="""

A Weibull minimum distribution
"""
                            )

class frechet_l_gen(rv_continuous):
    def _pdf(self, x, c):
        return c*pow(-x,c-1)*exp(-pow(-x,c))
    def _cdf(self, x, c):
        return exp(-pow(-x,c))
    def _ppf(self, q, c):
        return -pow(-log(q),1.0/c)
    def _munp(self, n, c):
        val = special.gamma(1.0+n*1.0/c)
        if (int(n) % 2): sgn = -1
        else:            sgn = 1
        return sgn*val
    def _entropy(self, c):
        return -_EULER / c - log(c) + _EULER + 1    
frechet_l = frechet_l_gen(b=0.0,name='frechet_l',longname="A Frechet left",
                          shapes='c',extradoc="""

A Frechet (left) distribution (also called Weibull maximum)
"""
                          )
weibull_max = frechet_l_gen(b=0.0,name='weibull_max',
                            longname="A Weibull maximum",
                            shapes='c',extradoc="""

A Weibull maximum distribution
"""
                            )


## Generalized Logistic
##
class genlogistic_gen(rv_continuous):
    def _pdf(self, x, c):
        Px = c*exp(-x)/(1+exp(-x))**(c+1.0)
        return Px
    def _cdf(self, x, c):
        Cx = (1+exp(-x))**(-c)
        return Cx
    def _ppf(self, q, c):
        vals = -log(pow(q,-1.0/c)-1)
        return vals
    def _stats(self, c):
        zeta = special.zeta
        mu = _EULER + special.psi(c)
        mu2 = pi*pi/6.0 + zeta(2,c)
        g1 = -2*zeta(3,c) + 2*_ZETA3
        g1 /= mu2**1.5
        g2 = pi**4/15.0 + 6*zeta(4,c)
        g2 /= mu2**2.0
        return mu, mu2, g1, g2
genlogistic = genlogistic_gen(name='genlogistic',
                              longname="A generalized logistic",
                              shapes='c',extradoc="""

Generalized logistic distribution
"""
                              )

## Generalized Pareto
class genpareto_gen(rv_continuous):
    def _argcheck(self, c):
        c = arr(c)
        self.b = where(c < 0, 1.0/abs(c), inf)
        return where(c==0, 0, 1)
    def _pdf(self, x, c):
        Px = pow(1+c*x,arr(-1.0-1.0/c))
        return Px
    def _cdf(self, x, c):
        return 1.0 - pow(1+c*x,arr(-1.0/c))
    def _ppf(self, q, c):
        vals = 1.0/c * (pow(1-q, -c)-1)
        return vals
    def _munp(self, n, c):
        k = arange(0,n+1)
        val = (-1.0/c)**n * sum(scipy.comb(n,k)*(-1)**k / (1.0-c*k))
        return where(c*n < 1, val, inf)
    def _entropy(self, c):
        if (c > 0):
            return 1+c
        else:
            self.b = -1.0 / c
            return rv_continuous._entropy(self, c)
genpareto = genpareto_gen(a=0.0,name='genpareto',
                          longname="A generalized Pareto",
                          shapes='c',extradoc="""

Generalized Pareto distribution
"""                          
                          )

## Generalized Exponential

class genexpon_gen(rv_continuous):
    def _pdf(self, x, a, b, c):
        return (a+b*(1-exp(-c*x)))*exp((a-b)*x+b*(1-exp(-c*x))/c)
    def _cdf(self, x, a, b, c):
        return 1.0-exp((a-b)*x + b*(1-exp(-c*x))/c)
genexpon = genexpon_gen(a=0.0,name='genexpon',
                        longname='A generalized exponential',
                        shapes='a,b,c',extradoc="""

Generalized exponential distribution
"""
                        )

## Generalized Extreme Value
##  c=0 is just gumbel distribution.
##  This version does not accept c==0
##  Use gumbel_r for c==0

class genextreme_gen(rv_continuous):
    def _argcheck(self, c):
        self.b = where(c > 0, 1.0 / c, inf)
        self.a = where(c < 0, 1.0 / c, -inf)
        return (c!=0)
    def _pdf(self, x, c):
        ex2 = 1-c*x
        pex2 = pow(ex2,1.0/c)
        p2 = exp(-pex2)*pex2/ex2
        return p2
    def _cdf(self, x, c):
        return exp(-pow(1-c*x,1.0/c))
    def _ppf(self, q, c):
        return 1.0/c*(1-(-log(q))**c)
    def _munp(self, n, c):
        k = arange(0,n+1)
        vals = 1.0/c**n * sum(scipy.comb(n,k) * (-1)**k * special.gamma(c*k + 1))
        return where(c*n > -1, vals, inf)
genextreme = genextreme_gen(name='genextreme',
                            longname="A generalized extreme value",
                            shapes='c',extradoc="""

Generalized extreme value (see gumbel_r for c=0)
"""
                            )
        
## Gamma (Use MATLAB and MATHEMATICA (b=theta=scale, a=alpha=shape) definition)

## gamma(a, loc, scale)  with a an integer is the Erlang distribution
## gamma(1, loc, scale)  is the Exponential distribution
## gamma(df/2, 0, 2) is the chi2 distribution with df degrees of freedom.

class gamma_gen(rv_continuous):
    def _rvs(self, a):
        return rand.standard_gamma(a, self._size)
    def _pdf(self, x, a):
        return x**(a-1)*exp(-x)/special.gamma(a)
    def _cdf(self, x, a):
        return special.gammainc(a, x)
    def _ppf(self, q, a):
        return special.gammaincinv(a,q)
    def _stats(self, a):
        return a, a, 2.0/sqrt(a), 6.0/a
    def _entropy(self, a):
        return special.psi(a)*(1-a) + 1 + special.gammaln(a)
gamma = gamma_gen(a=0.0,name='gamma',longname='A gamma',
                  shapes='a',extradoc="""

Gamma distribution
"""
                  )

# Generalized Gamma
class gengamma_gen(rv_continuous):
    def _argcheck(self, a, c):
        return (a > 0) & (c != 0)
    def _pdf(self, x, a, c):
        return abs(c)* x**(c*a-1) / special.gamma(a) * exp(-x**c)
    def _cdf(self, x, a, c):
        val = special.gammainc(a,x**c)
        cond = c + 0*val
        return where(cond>0,val,1-val)
    def _ppf(self, q, a, c):
        val1 = special.gammaincinv(a,q)
        val2 = special.gammaincinv(a,1.0-q)
        ic = 1.0/c
        cond = c+0*val1
        return where(cond > 0,val1**ic,val2**ic)
    def _munp(self, n, a, c):
        return special.gamma(a+n*1.0/c) / special.gamma(a)
    def _entropy(a,c):
        val = special.psi(a)
        return a*(1-val) + 1.0/c*val + special.gammaln(a)-log(abs(c))
gengamma = gengamma_gen(a=0.0, name='gengamma',
                        longname='A generalized gamma',
                        shapes="a,c", extradoc="""

Generalized gamma distribution
"""
                        )

##  Generalized Half-Logistic
##

class genhalflogistic_gen(rv_continuous):
    def _argcheck(self, c):
        self.b = 1.0 / c
        return (c > 0)
    def _pdf(self, x, c):
        limit = 1.0/c
        tmp = arr(1-c*x)
        tmp0 = tmp**(limit-1)
        tmp2 = tmp0*tmp
        return 2*tmp0 / (1+tmp2)**2
    def _cdf(self, x, c):
        limit = 1.0/c
        tmp = arr(1-c*x)
        tmp2 = tmp**(limit)
        return (1.0-tmp2) / (1+tmp2)
    def _ppf(self, q, c):
        return 1.0/c*(1-((1.0-q)/(1.0+q))**c)
    def _entropy(self,c):
        return 2 - (2*c+1)*log(2)
genhalflogistic = genhalflogistic_gen(a=0.0, name='genhalflogistic',
                                      longname="A generalized half-logistic",
                                      shapes='c',extradoc="""

Generalized half-logistic
"""
                                      )

## Gompertz (Truncated Gumbel)
##  Defined for x>=0

class gompertz_gen(rv_continuous):
    def _pdf(self, x, c):
        ex = exp(x)
        return c*ex*exp(-c*(ex-1))
    def _cdf(self, x, c):
        return 1.0-exp(-c*(exp(x)-1))
    def _ppf(self, q, c):
        return log(1-1.0/c*log(1-q))
    def _entropy(self, c):
        return 1.0 - log(c) - exp(c)*special.expn(1,c)
gompertz = gompertz_gen(a=0.0, name='gompertz',
                        longname="A Gompertz (truncated Gumbel) distribution",
                        shapes='c',extradoc="""

Gompertz (truncated Gumbel) distribution
"""
                        )
    
## Gumbel, Log-Weibull, Fisher-Tippett, Gompertz
## The left-skewed gumbel distribution.
## and right-skewed are available as gumbel_l  and gumbel_r

class gumbel_r_gen(rv_continuous):
    def _pdf(self, x):
        ex = exp(-x)
        return ex*exp(-ex)
    def _cdf(self, x):
        return exp(-exp(-x))
    def _ppf(self, q):
        return -log(-log(q))
    def _stats(self):
        return _EULER, pi*pi/6.0, \
               12*sqrt(6)/pi**3 * _ZETA3, 12.0/5
    def _entropy(self):
        return 1.0608407169541684911
gumbel_r = gumbel_r_gen(name='gumbel_r',longname="A (right-skewed) Gumbel",
                        extradoc="""

Right-skewed Gumbel (Log-Weibull, Fisher-Tippett, Gompertz) distribution
"""
                        )
class gumbel_l_gen(rv_continuous):
    def _pdf(self, x):
        ex = exp(x)
        return ex*exp(-ex)
    def _cdf(self, x):
        return 1.0-exp(-exp(x))
    def _ppf(self, q):
        return log(-log(1-q))
    def _stats(self):
        return _EULER, pi*pi/6.0, \
               12*sqrt(6)/pi**3 * _ZETA3, 12.0/5
    def _entropy(self):
        return 1.0608407169541684911
gumbel_l = gumbel_l_gen(name='gumbel_l',longname="A left-skewed Gumbel",
                        extradoc="""

Left-skewed Gumbel distribution
"""
                        )

# Half-Cauchy

class halfcauchy_gen(rv_continuous):
    def _pdf(self, x):
        return 2.0/pi/(1.0+x*x)
    def _cdf(self, x):
        return 2.0/pi*arctan(x)
    def _ppf(self, q):
        return tan(pi/2*q)
    def _stats(self):
        return inf, inf, nan, nan
    def _entropy(self):
        return log(2*pi)
halfcauchy = halfcauchy_gen(a=0.0,name='halfcauchy',
                            longname="A Half-Cauchy",extradoc="""

Half-Cauchy distribution
"""
                            )


## Half-Logistic
##  

class halflogistic_gen(rv_continuous):
    def _pdf(self, x):
        return 0.5/(cosh(x/2.0))**2.0
    def _cdf(self, x):
        return tanh(x/2.0)
    def _ppf(self, q):
        return 2*arctanh(q)
    def _munp(self, n):
        if n==1: return 2*log(2)
        if n==2: return pi*pi/3.0
        if n==3: return 9*_ZETA3
        if n==4: return 7*pi**4 / 15.0
        return 2*(1-pow(2.0,1-n))*special.gamma(n+1)*special.zeta(n,1)
    def _entropy(self):
        return 2-log(2)
halflogistic = halflogistic_gen(a=0.0, name='halflogistic',
                                longname="A half-logistic",
                                extradoc="""

Half-logistic distribution                                
"""
                                )


## Half-normal = chi(1, loc, scale)

class halfnorm_gen(rv_continuous):
    def _rvs(self):
        return abs(norm.rvs(size=self._size))
    def _pdf(self, x):
        return sqrt(2.0/pi)*exp(-x*x/2.0)
    def _cdf(self, x):
        return special.ndtr(x)*2-1.0
    def _ppf(self, q):
        return special.ndtri((1+q)/2.0)
    def _stats(self):
        return sqrt(2.0/pi), 1-2.0/pi, sqrt(2)*(4-pi)/(pi-2)**1.5, \
               8*(pi-3)/(pi-2)**2
    def _entropy(self):
        return 0.5*log(pi/2.0)+0.5
halfnorm = halfnorm_gen(a=0.0, name='halfnorm',
                        longname="A half-normal",
                        extradoc="""
                        
Half-normal distribution
"""                        
                        )

## Hyperbolic Secant

class hypsecant_gen(rv_continuous):
    def _pdf(self, x):
        return 1.0/(pi*cosh(x))
    def _cdf(self, x):
        return 2.0/pi*arctan(exp(x))
    def _ppf(self, q):
        return log(tan(pi*q/2.0))
    def _stats(self):
        return 0, pi*pi/4, 0, 2
    def _entropy(self):
        return log(2*pi)
hypsecant = hypsecant_gen(name='hypsecant',longname="A hyperbolic secant",
                          extradoc="""

Hyperbolic secant distribution
"""
                          )

## Gauss Hypergeometric

class gausshyper_gen(rv_continuous):
    def _argcheck(self, a, b, c, z):
        return (a > 0) & (b > 0) & (c==c) & (z==z)
    def _pdf(self, x, a, b, c, z):
        Cinv = gam(a)*gam(b)/gam(a+b)*special.hyp2f1(c,a,a+b,-z)
        return 1.0/Cinv * x**(a-1.0) * (1.0-x)**(b-1.0) / (1.0+z*x)**c
    def _munp(self, n, a, b, c, z):
        fac = special.beta(n+a,b) / special.beta(a,b)
        num = special.hyp2f1(c,a+n,a+b+n,-z)
        den = special.hyp2f1(c,a,a+b,-z)
        return fac*num / den
gausshyper = gausshyper_gen(a=0.0, b=1.0, name='gausshyper',
                            longname="A Gauss hypergeometric",
                            shapes="a,b,c,z",
                            extradoc="""

Gauss hypergeometric distribution
"""
                            )

##  Inverted Gamma
#     special case of generalized gamma with c=-1
#

class invgamma_gen(rv_continuous):
    def _pdf(self, x, a):
        return x**(-a-1) / special.gamma(a) * exp(-1.0/x)
    def _cdf(self, x, a):
        return 1.0-special.gammainc(a, 1.0/x)
    def _ppf(self, q, a):
        return 1.0/special.gammaincinv(a,1-q)
    def _munp(self, n, a):
        return special.gamma(a-n) / special.gamma(a)
    def _entropy(self, a):
        return a - (a+1.0)*special.psi(a) + special.gammaln(a)
invgamma = invgamma_gen(a=0.0, name='invgamma',longname="An inverted gamma",
                        shapes='a',extradoc="""

Inverted gamma distribution
"""
                        )


## Inverse Normal Distribution
# scale is gamma from DATAPLOT and B from Regress

class invnorm_gen(rv_continuous):
    def _rvs(self, mu):
        return rv._inst._Wald(mu,size=(self._size,))
    def _pdf(self, x, mu):
        return 1.0/sqrt(2*pi*x**3.0)*exp(-1.0/(2*x)*((x-mu)/mu)**2)
    def _cdf(self, x, mu):
        fac = sqrt(1.0/x)
        C1 = norm.cdf(fac*(x-mu)/mu)
        C1 += exp(2.0/mu)*norm.cdf(-fac*(x+mu)/mu)
        return C1
    def _stats(self, mu):
        return mu, mu**3.0, 3*sqrt(mu), 15*mu
invnorm = invnorm_gen(a=0.0, name='invnorm', longname="An inverse normal",
                      shapes="mu",extradoc="""

Inverse normal distribution
"""
                      )

## Inverted Weibull

class invweibull_gen(rv_continuous):
    def _pdf(self, x, c):
        xc1 = x**(-c-1.0)
        xc2 = xc1*x
        return c*xc1*xc2
    def _cdf(self, x, c):
        xc1 = x**(-c)
        return exp(-xc1)
    def _ppf(self, q, c):
        return pow(-log(q),arr(-1.0/c))
    def _entropy(self, c):
        return 1+_EULER + _EULER / c - log(c)
invweibull = invweibull_gen(a=0,name='invweibull',
                            longname="An inverted Weibull",
                            shapes='c',extradoc="""

Inverted Weibull distribution
"""
                            )

## Johnson SB

class johnsonsb_gen(rv_continuous):
    def _argcheck(self, a, b):
        return (b > 0) & (a==a)
    def _pdf(self, x, a, b):
        trm = norm.pdf(a+b*log(x/(1.0-x)))
        return b*1.0/(x*(1-x))*trm
    def _cdf(self, x, a, b):
        return norm.cdf(a+b*log(x/(1.0-x)))
    def _ppf(self, q, a, b):
        return 1.0/(1+exp(-1.0/b*norm.ppf(q)-a))
johnsonsb = johnsonsb_gen(a=0.0,b=1.0,name='johnsonb',
                          longname="A Johnson SB",
                          shapes="a,b",extradoc="""

Johnson SB distribution
"""
                          )

## Johnson SU
class johnsonsu_gen(rv_continuous):
    def _argcheck(self, a, b):
        return (b > 0) & (a==a)
    def _pdf(self, x, a, b):
        x2 = x*x
        trm = norm.pdf(a+b*log(x+sqrt(x2+1)))
        return b*1.0/sqrt(x2+1.0)*trm
    def _cdf(self, x, a, b):
        return norm.cdf(a+b*log(x+sqrt(x*x+1)))
    def _ppf(self, q, a, b):
        return sinh((norm.ppf(q)-a)/b)
johnsonsu = johnsonsu_gen(name='johnsonsu',longname="A Johnson SU",
                          shapes="a,b", extradoc="""

Johnson SB distribution
"""
                          )


## Laplace Distribution

class laplace_gen(rv_continuous):
    def _pdf(self, x):
        return 0.5*exp(-abs(x))
    def _cdf(self, x):
        return where(x > 0, 1.0-0.5*exp(-x), 0.5*exp(x))
    def _ppf(self, q):
        return where(q > 0.5, -log(2*(1-q)), log(2*q))
    def _stats(self):
        return 0, 2, 0, 3
    def _entropy(self):
        return log(2)+1
laplace = laplace_gen(name='laplace', longname="A Laplace",
                      extradoc="""

Laplacian distribution
"""
                      )


## Levy Distribution

class levy_gen(rv_continuous):
    def _pdf(self, x):
        return 1/sqrt(2*x)/x*exp(-1/(2*x))
    def _cdf(self, x):
        return 2*(1-norm._cdf(1/sqrt(x)))
    def _ppf(self, q):
        val = norm._ppf(1-q/2.0)
        return 1.0/(val*val)
    def _stats(self):
        return inf, inf, nan, nan
levy = levy_gen(a=0.0,name="levy", longname = "A Levy", extradoc="""

Levy distribution
"""
                )

## Left-skewed Levy Distribution

class levy_l_gen(rv_continuous):
    def _pdf(self, x):
        ax = abs(x)
        return 1/sqrt(2*ax)/ax*exp(-1/(2*ax))
    def _cdf(self, x):
        ax = abs(x)
        return 2*norm._cdf(1/sqrt(ax))-1
    def _ppf(self, q):
        val = norm._ppf((q+1.0)/2)
        return -1.0/(val*val)
    def _stats(self):
        return inf, inf, nan, nan
levy_l = levy_l_gen(b=0.0,name="levy_l", longname = "A left-skewed Levy", extradoc="""

Left-skewed Levy distribution
"""
                )

## Levy-stable Distribution (only random variates)

class levy_stable_gen(rv_continuous):
    def _rvs(self, alpha, beta):
        sz = self._size
        TH = uniform.rvs(loc=-pi/2.0,scale=pi,size=sz)
        W = expon.rvs(size=sz)
        if alpha==1:
            return 2/pi*(pi/2+beta*TH)*tan(TH)-beta*log((pi/2*W*cos(TH))/(pi/2+beta*TH))
        # else
        ialpha = 1.0/alpha
        aTH = alpha*TH
        if beta==0:
            return W/(cos(TH)/tan(aTH)+sin(TH))*((cos(aTH)+sin(aTH)*tan(TH))/W)**ialpha
        # else
        val0 = beta*tan(pi*alpha/2)
        th0 = arctan(val0)/alpha
        val3 = W/(cos(TH)/tan(alpha*(th0+TH))+sin(TH))
        res3 = val3*((cos(aTH)+sin(aTH)*tan(TH)-val0*(sin(aTH)-cos(aTH)*tan(TH)))/W)**ialpha
        return res3
        
    def _argcheck(self, alpha, beta):        
        if beta == -1:
            self.b = 0.0
        elif beta == 1:
            self.a = 0.0
        return (alpha > 0) & (alpha <= 2) & (beta <= 1) & (beta >= -1)
    
    def _pdf(self, x, alpha, beta):
        raise NotImplementedError

levy_stable = levy_stable_gen(name='levy_stable', longname="A Levy-stable",
                    shapes="alpha, beta", extradoc="""

Levy-stable distribution (only random variates available -- ignore other docs)
"""
                    )
    

## Logistic (special case of generalized logistic with c=1)
## Sech-squared

class logistic_gen(rv_continuous):
    def _pdf(self, x):
        ex = exp(-x)
        return ex / (1+ex)**2.0
    def _cdf(self, x):
        return 1.0/(1+exp(-x))
    def _ppf(self, q):
        return -log(1.0/q-1)
    def _stats(self):
        return 0, pi*pi/3.0, 0, 6.0/5.0
    def _entropy(self):
        return 1.0
logistic = logistic_gen(name='logistic', longname="A logistic",
                        extradoc="""

Logistic distribution
"""
                        )


## Log Gamma
#
class loggamma_gen(rv_continuous):
    def _pdf(self, x, c):
        return exp(c*x-exp(x))/ special.gamma(c)
    def _cdf(self, x, c):
        return special.gammainc(c, exp(x))/ special.gamma(c)
    def _ppf(self, q, c):
        return log(special.gammaincinv(c,q*special.gamma(c)))
loggamma = loggamma_gen(name='loggamma', longname="A log gamma",
                        extradoc="""

Log gamma distribution
"""
                        )

## Log-Laplace  (Log Double Exponential)
##

class loglaplace_gen(rv_continuous):
    def _pdf(self, x, c):
        cd2 = c/2.0
        c = where(x < 1, c, -c)
        return cd2*x**(c-1)
    def _cdf(self, x, c):
        return where(x < 1, 0.5*x**c, 1-0.5*x**(-c))
    def _ppf(self, q, c):
        return where(q < 0.5, (2.0*q)**(1.0/c), (2*(1.0-q))**(-1.0/c))
    def _entropy(self, c):
        return log(2.0/c) + 1.0
loglaplace = loglaplace_gen(a=0.0, name='loglaplace',
                            longname="A log-Laplace",shapes='c',
                            extradoc="""

Log-Laplace distribution (Log Double Exponential)
"""
                            )

## Lognormal (Cobb-Douglass)
## std is a shape parameter and is the variance of the underlying
##    distribution.
## the mean of the underlying distribution is log(scale)

class lognorm_gen(rv_continuous):
    def _rvs(self, s):
        return exp(s * norm.rvs(size=self._size))
    def _pdf(self, x, s):
        Px = exp(-log(x)**2 / (2*s**2))
        return Px / (s*x*sqrt(2*pi))
    def _cdf(self, x, s):
        return norm.cdf(log(x)/s)
    def _ppf(self, q, s):
        return exp(s*norm._ppf(q))
    def _stats(self, s):
        p = exp(s*s)
        mu = sqrt(p)
        mu2 = p*(p-1)
        g1 = sqrt((p-1))*(2+p)
        g2 = scipy.polyval([1,2,3,0,-6.0],p)
        return mu, mu2, g1, g2
    def _entropy(self, s):
        return 0.5*(1+log(2*pi)+2*log(s))
lognorm = lognorm_gen(a=0.0, name='lognorm',
                      longname='A lognormal', shapes='s',
                      extradoc="""

Lognormal distribution
"""
                      )

# Gibrat's distribution is just lognormal with s=1

class gilbrat_gen(lognorm_gen):
    def _rvs(self):
        return lognorm_gen._rvs(self, 1.0)
    def _pdf(self, x):
        return lognorm_gen._pdf(self, x, 1.0)
    def _cdf(self, x):
        return lognorm_gen._cdf(self, x, 1.0)
    def _ppf(self, q):
        return lognorm_gen._ppf(self, q, 1.0)
    def _stats(self):
        return lognorm_gen._stats(self, 1.0)
    def _entropy(self):
        return 0.5*log(2*pi) + 0.5
gilbrat = gilbrat_gen(a=0.0, name='gilbrat', longname='A Gilbrat',
                      extradoc="""

Gilbrat distribution
"""
                      )


# MAXWELL
#  a special case of chi with df = 3, loc=0.0, and given scale = 1.0/sqrt(a)
#    where a is the parameter used in mathworld description

class maxwell_gen(rv_continuous):
    def _rvs(self):
        return chi.rvs(3.0,size=self._size)
    def _pdf(self, x):
        return sqrt(2.0/pi)*x*x*exp(-x*x/2.0)
    def _cdf(self, x):
        return special.gammainc(1.5,x*x/2.0)
    def _ppf(self, q):
        return sqrt(2*special.gammaincinv(1.5,q))
    def _stats(self):
        val = 3*pi-8
        return 2*sqrt(2.0/pi), 3-8/pi, sqrt(2)*(32-10*pi)/val**1.5, \
               (-12*pi*pi + 160*pi - 384) / val**2.0
    def _entropy(self):
        return _EULER + 0.5*log(2*pi)-0.5
maxwell = maxwell_gen(a=0.0, name='maxwell', longname="A Maxwell",
                      extradoc="""

Maxwell distribution
"""
                      )

# Mielke's Beta-Kappa

class mielke_gen(rv_continuous):
    def _pdf(self, x, k, s):
        return k*x**(k-1.0) / (1.0+x**s)**(1.0+k*1.0/s)
    def _cdf(self, x, k, s):
        return x**k / (1.0+x**s)**(k*1.0/s)
    def _ppf(self, q, k, s):
        qsk = pow(q,s*1.0/k)
        return pow(qsk/(1.0-qsk),1.0/s)
mielke = mielke_gen(a=0.0, name='mielke', longname="A Mielke's Beta-Kappa",
                    shapes="k,s", extradoc="""

Mielke's Beta-Kappa distribution
"""
                    )
     
# Nakagami (cf Chi)

class nakagami_gen(rv_continuous):
    def _pdf(self, x, nu):
        return 2*nu**nu/gam(nu)*(x**(2*nu-1.0))*exp(-nu*x*x)
    def _cdf(self, x, nu):
        return special.gammainc(nu,nu*x*x)
    def _ppf(self, q, nu):
        return sqrt(1.0/nu*special.gammaincinv(nu,q))
    def _stats(self, nu):
        mu = gam(nu+0.5)/gam(nu)/sqrt(nu)
        mu2 = 1.0-mu*mu
        g1 = mu*(1-4*nu*mu2)/2.0/nu/mu2**1.5
        g2 = -6*mu**4*nu + (8*nu-2)*mu**2-2*nu + 1
        g2 /= nu*mu2**2.0
        return mu, mu2, g1, g2
nakagami = nakagami_gen(a=0.0, name="nakagami", longname="A Nakagami",
                        shapes='nu', extradoc="""

Nakagami distribution
"""
                        )
    

# Non-central chi-squared
# nc is lambda of definition, df is nu

class ncx2_gen(rv_continuous):
    def _rvs(self, df, nc):
        return rand.noncentral_chi2(df,nc,self._size)        
    def _pdf(self, x, df, nc):
        a = arr(df/2.0)
        Px = exp(-nc/2.0)*special.hyp0f1(a,nc*x/4.0)
        Px *= exp(-x/2.0)*x**(a-1) / arr(2**a * special.gamma(a))
        return Px
    def _cdf(self, x, df, nc):
        return special.chndtr(x,df,nc)
    def _ppf(self, q, df, nc):
        return special.chndtrix(q,df,nc)
    def _stats(self, df, nc):
        val = df + 2.0*nc
        return df + nc, 2*val, sqrt(8)*(val+nc)/val**1.5, \
               12.0*(val+2*nc)/val**2.0
ncx2 = ncx2_gen(a=0.0, name='ncx2', longname="A non-central chi-squared",
                shapes="df,nc", extradoc="""

Non-central chi-squared distribution
"""
                )

# Non-central F

class ncf_gen(rv_continuous):
    def _rvs(self, dfn, dfd, nc):
        return rand.noncentral_f(dfn,dfd,nc,self._size)
    def _pdf(self, x, dfn, dfd, nc):
        n1,n2 = dfn, dfd
        Px = exp(-nc/2+nc*n1*x/(2*(n2+n1*x)))
        Px *= n1**(n1/2) * n2**(n2/2) * x**(n1/2-1)
        Px *= (n2+n1*x)**(-(n1+n2)/2)
        Px *= special.gamma(n1/2)*special.gamma(1+n2/2)
        Px *= special.assoc_laguerre(-nc*n1*x/(2.0*(n2+n1*x)),n2/2,n1/2-1)
        Px /= special.beta(n1/2,n2/2)*special.gamma((n1+n2)/2.0)
    def _cdf(self, x, dfn, dfd, nc):
        return special.ncfdtr(dfn,dfd,nc,x)
    def _ppf(self, q, dfn, dfd, nc):
        return special.ncfdtri(dfn, dfd, nc, q)
    def _munp(self, n, dfn, dfd, nc):
        val = (dfn *1.0/dfd)**n
        val *= gam(n+0.5*dfn)*gam(0.5*dfd-n) / gam(dfd*0.5)
        val *= exp(-nc / 2.0)
        val *= special.hyp1f1(n+0.5*dfn, 0.5*dfn, 0.5*nc)
        return val
    def _stats(self, dfn, dfd, nc):        
        mu = where(dfd <= 2, inf, dfd / (dfd-2.0)*(1+nc*1.0/dfn))
        mu2 = where(dfd <=4, inf, 2*(dfd*1.0/dfn)**2.0 * \
                    ((dfn+nc/2.0)**2.0 + (dfn+nc)*(dfd-2.0)) / \
                    ((dfd-2.0)**2.0 * (dfd-4.0)))
        return mu, mu2, None, None
ncf = ncf_gen(a=0.0, name='ncf', longname="A non-central F distribution",
              shapes="dfn,dfd,nc", extradoc="""

Non-central F distribution
"""
              )

## Student t distribution

class t_gen(rv_continuous):
    def _rvs(self, df):
        Y = f.rvs(df, df, size=self._size)
        sY = sqrt(Y)
        return 0.5*sqrt(df)*(sY-1.0/sY)
    def _pdf(self, x, df):
        r = arr(df*1.0)
        Px = exp(special.gammaln((r+1)/2)-special.gammaln(r/2))
        Px /= sqrt(r*pi)*(1+(x**2)/r)**((r+1)/2)
        return Px
    def _cdf(self, x, df):
        return special.stdtr(df, x)
    def _ppf(self, q, df):
        return special.stdtrit(df, q)
    def _stats(self, df):
        mu2 = where(df > 2, df / (df-2.0), inf)
        g1 = where(df > 3, 0.0, nan)
        g2 = where(df > 4, 6.0/(df-4.0), nan)
        return 0, mu2, g1, g2
t = t_gen(name='t',longname="Student's T",
          shapes="df", extradoc="""

Student's T distribution
"""
          )

## Non-central T distribution

class nct_gen(rv_continuous):
    def _rvs(self, df, nc):
        return norm.rvs(loc=nc,size=self._size)*sqrt(df) / sqrt(chi2.rvs(df,size=self._size))
    def _pdf(self, x, df, nc):
        n = df*1.0
        nc = nc*1.0
        x2 = x*x
        ncx2 = nc*nc*x2
        fac1 = n + x2
        Px = n**(n/2) * special.gamma(n+1)
        Px /= arr(2.0**n*exp(nc*nc/2)*fac1**(n/2)*special.gamma(n/2))
        valF = ncx2 / (2*fac1)
        trm1 = sqrt(2)*nc*x*special.hyp1f1(n/2+1,1.5,valF)
        trm1 /= arr(fac1*special.gamma((n+1)/2))
        trm2 = special.hyp1f1((n+1)/2,0.5,valF)
        trm2 /= arr(sqrt(fac1)*special.gamma(n/2+1))
        Px *= trm1+trm2
        return Px
    def _cdf(self, x, df, nc):
        return special.nctdtr(df, nc, x)
    def _ppf(self, q, df, nc):
        return special.nctdtrit(df, nc, q)
    def _stats(self, df, nc, moments='mv'):
        mu, mu2, g1, g2 = None, None, None, None
        val1 = gam((df-1.0)/2.0)
        val2 = gam(df/2.0)
        if 'm' in moments:
            mu = nc*sqrt(df/2.0)*val1/val2
        if 'v' in moments:
            var = (nc*nc+1.0)*df/(df-2.0)
            var -= nc*nc*df* val1**2 / 2.0 / val2**2
            mu2 = var
        if 's' in moments:
            g1n = 2*nc*sqrt(df)*val1*((nc*nc*(2*df-7)-3)*val2**2 \
                                      -nc*nc*(df-2)*(df-3)*val1**2)
            g1d = (df-3)*sqrt(2*df*(nc*nc+1)/(df-2) - \
                              nc*nc*df*(val1/val2)**2) * val2 * \
                              (nc*nc*(df-2)*val1**2 - \
                               2*(nc*nc+1)*val2**2)
            g1 = g1n/g1d
        if 'k' in moments:            
            g2n = 2*(-3*nc**4*(df-2)**2 *(df-3) *(df-4)*val1**4 + \
                     2**(6-2*df) * nc*nc*(df-2)*(df-4)* \
                     (nc*nc*(2*df-7)-3)*pi* gam(df+1)**2 - \
                     4*(nc**4*(df-5)-6*nc*nc-3)*(df-3)*val2**4)
            g2d = (df-3)*(df-4)*(nc*nc*(df-2)*val1**2 - \
                                 2*(nc*nc+1)*val2)**2
            g2 = g2n / g2d
        return mu, mu2, g1, g2
nct = nct_gen(name="nct", longname="A Noncentral T",
              shapes="df,nc", extradoc="""

Non-central Student T distribution
"""
              )

# Pareto

class pareto_gen(rv_continuous):
    def _pdf(self, x, b):
        return b * x**(-b-1)
    def _cdf(self, x, b):
        return 1 -  x**(-b)
    def _ppf(self, q, b):
        return pow(1-q, -1.0/b)
    def _stats(self, b, moments='mv'):
        mu, mu2, g1, g2 = None, None, None, None
        if 'm' in moments:
            mask = b > 1
            bt = extract(mask,b)
            mu = valarray(shape(b),value=inf)
            insert(mu, mask, bt / (bt-1.0))
        if 'v' in moments:
            mask = b > 2
            bt = extract( mask,b)
            mu2 = valarray(shape(b), value=inf)
            insert(mu2, mask, bt / (bt-2.0) / (bt-1.0)**2)
        if 's' in moments:
            mask = b > 3
            bt = extract( mask,b)
            g1 = valarray(shape(b), value=nan)
            vals = 2*(bt+1.0)*sqrt(b-2.0)/((b-3.0)*sqrt(b))
            insert(g1, mask, vals)
        if 'k' in moments:
            mask = b > 4
            bt = extract( mask,b)
            g2 = valarray(shape(b), value=nan)
            vals = 6.0*polyval([1.0,1.0,-6,-2],bt)/ \
                   polyval([1.0,-7.0,12.0,0.0],bt)
            insert(g2, mask, vals)
        return mu, mu2, g1, g2
    def _entropy(self, c):
        return 1 + 1.0/c - log(c)
pareto = pareto_gen(a=1.0, name="pareto", longname="A Pareto",
                    shapes="b", extradoc="""

Pareto distribution
"""
                    )

# LOMAX (Pareto of the second kind.)
#  Special case of Pareto of the first kind (location=-1.0)

class lomax_gen(rv_continuous):
    def _pdf(self, x, c):
        return c*1.0/(1.0+x)**(c+1.0)
    def _cdf(self, x, c):
        return 1.0-1.0/(1.0+x)**c
    def _ppf(self, q, c):
        return pow(1.0-q,-1.0/c)-1
    def _stats(self, c):
        mu, mu2, g1, g2 = pareto.stats(c, loc=-1.0, moments='mvsk')
        return mu, mu2, g1, g2
    def _entropy(self, c):
        return 1+1.0/c-log(c)
lomax = lomax_gen(a=0.0, name="lomax",
                  longname="A Lomax (Pareto of the second kind)",
                  shapes="c", extradoc="""

Lomax (Pareto of the second kind) distribution
"""
                  )
## Power-function distribution
##   Special case of beta dist. with d =1.0

class powerlaw_gen(rv_continuous):
    def _pdf(self, x, a):
        return a*x**(a-1.0)
    def _cdf(self, x, a):
        return x**(a*1.0)
    def _ppf(self, q, a):
        return pow(q, 1.0/a)
    def _stats(self, a):
        return a/(a+1.0), a*(a+2.0)/(a+1.0)**2, \
               2*(1.0-a)*sqrt((a+2.0)/(a*(a+3.0))), \
               6*polyval([1,-1,-6,2],a)/(a*(a+3.0)*(a+4))
    def _entropy(self, a):
        return 1 - 1.0/a - log(a)
powerlaw = powerlaw_gen(a=0.0, b=1.0, name="powerlaw",
                        longname="A power-function",
                        shapes="a", extradoc="""

Power-function distribution
"""
                        )

# Power log normal

class powerlognorm_gen(rv_continuous):
    def _pdf(self, x, c, s):
        return c/(x*s)*norm.pdf(log(x)/s)
    def _cdf(self, x, c, s):
        return 1.0 - pow(norm.cdf(-log(x)/s),c*1.0)
    def _ppf(self, q, c, s):
        return exp(-s*norm.ppf(pow(1.0-q,1.0/c)))
powerlognorm = powerlognorm_gen(a=0.0, name="powerlognorm",
                                longname="A power log-normal",
                                shapes="c,s", extradoc="""

Power log-normal distribution
"""
                                )

# Power Normal

class powernorm_gen(norm_gen):
    def _pdf(self, x, c):
        return c*norm_gen._pdf(self, x)* \
               (norm_gen._cdf(self, -x)**(c-1.0))
    def _cdf(self, x, c):
        return 1.0-norm_gen._cdf(self, -x)**(c*1.0)
    def _ppf(self, q, c):
        return -norm_gen._ppf(self, pow(1.0-q,1.0/c))
powernorm = powernorm_gen(name='powernorm', longname="A power normal",
                          shapes="c", extradoc="""

Power normal distribution
"""
                          )

# R-distribution ( a general-purpose distribution with a
#  variety of shapes.

class rdist_gen(rv_continuous):
    def _pdf(self, x, c):
        return pow((1.0-x*x),c/2.0-1) / special.beta(0.5,c/2.0)
    def _cdf(self, x, c):
        return 0.5 + x/special.beta(0.5,c/2.0)* \
               special.hyp2f1(0.5,1.0-c/2.0,1.5,x*x)
    def _munp(self, c):
        return (1-(n % 2))*special.beta((n+1.0)/2,c/2.0)
rdist = rdist_gen(a=-1.0,b=1.0, name="rdist", longname="An R-distributed",
                  shapes="c", extradoc="""

R-distribution
"""
                  )

# Rayleigh distribution (this is chi with df=2 and loc=0.0)
# scale is the mode.

class rayleigh_gen(rv_continuous):
    def _rvs(self):
        return chi.rvs(2,size=self._size)
    def _pdf(self, r):
        return r*exp(-r*r/2.0)
    def _cdf(self, r):
        return 1.0-exp(-r*r/2.0)
    def _ppf(self, q):
        return sqrt(-2*log(1-q))
    def _stats(self):
        val = 4-pi
        return pi/2, val/2, 2*(pi-3)*sqrt(pi)/val**1.5, \
               6*pi/val-16/val**2
    def _entropy(self):
        return _EULER/2.0 + 1 - 0.5*log(2)
rayleigh = rayleigh_gen(a=0.0, name="rayleigh",
                        longname="A Rayleigh",
                        extradoc="""

Rayleigh distribution
"""
                        )

# Reciprocal Distribution
class reciprocal_gen(rv_continuous):
    def _argcheck(self, a, b):
        self.a = a
        self.b = b
        self.d = log(b*1.0 / a)
        return (a > 0) & (b > 0) & (b > a)
    def _pdf(self, x, a, b):
        # argcheck should be called before _pdf
        return 1.0/(x*self.d)
    def _cdf(self, x, a, b):
        return (log(x)-log(a)) / self.d
    def _ppf(self, q, a, b):
        return a*pow(b*1.0/a,q)
    def _munp(self, n, a, b):
        return 1.0/self.d / n * (pow(b*1.0,n) - pow(a*1.0,n))
    def _entropy(self,a,b):
        return 0.5*log(a*b)+log(log(b/a))
reciprocal = reciprocal_gen(name="reciprocal",
                            longname="A reciprocal",
                            shapes="a,b", extradoc="""

Reciprocal distribution
"""
                            )

# Rice distribution

class rice_gen(rv_continuous):
    def _pdf(self, x, b):
        return x*exp(-(x*x+b*b)/2.0)*special.i0(x*b)
    def _munp(self, n, b):
        nd2 = n/2.0
        n1 = 1+nd2
        b2 = b*b/2.0
        return 2.0**(nd2)*exp(-b2)*special.gamma(n1) * \
               special.hyp1f1(n1,1,b2)
rice = rice_gen(a=0.0, name="rice", longname="A Rice",
                shapes="b", extradoc="""

Rician distribution
"""
                )

# Reciprocal Inverse Gaussian

class recipinvgauss_gen(rv_continuous):
    def _pdf(self, x, mu):
        return 1.0/sqrt(2*pi*x)*exp(-(1-mu*x)**2.0 / (2*x*mu**2.0))
    def _cdf(self, x, mu):
        trm1 = 1.0/mu - x
        trm2 = 1.0/mu + x
        isqx = 1.0/sqrt(x)
        return 1.0-norm.cdf(isqx*trm1)-exp(2.0/mu)*norm.cdf(-isqx*trm2)
recipinvgauss = recipinvgauss_gen(a=0.0, name='recipinvgauss',
                                  longname="A reciprocal inverse Gaussian",
                                  shapes="mu", extradoc="""

Reciprocal inverse Gaussian
"""
                                  )

# Semicircular

class semicircular_gen(rv_continuous):
    def _pdf(self, x):
        return 2.0/pi*sqrt(1-x*x)
    def _cdf(self, x):
        return 0.5+1.0/pi*(x*sqrt(1-x*x) + arcsin(x))
    def _stats(self):
        return 0, 0.25, 0, -1.0
    def _entropy(self):
        return 0.64472988584940017414        
semicircular = semicircular_gen(a=-1.0,b=1.0, name="semicircular",
                                longname="A semicircular",
                                extradoc="""

Semicircular distribution
"""
                                )

# Triangular
# up-sloping line from loc to (loc + c*scale) and then downsloping line from
#    loc + c*scale to loc + scale

# _trstr = "Left must be <= mode which must be <= right with left < right"
class triang_gen(rv_continuous):
    def _argcheck(self, c):
        return (c >= 0) & (c <= 1)
    def _pdf(self, x, c):
        return where(x < c, 2*x/c, 2*(1-x)/(1-c))
    def _cdf(self, x, c):
        return where(x < c, x*x/c, (x*x-2*x+c)/(c-1))
    def _ppf(self, q, c):
        return where(q < c, sqrt(c*q), 1-sqrt((1-c)*(1-q)))
    def _stats(self, c):
        return (c+1.0)/3.0, (1.0-c+c*c)/18, sqrt(2)*(2*c-1)*(c+1)*(c-2) / \
               (5*(1.0-c+c*c)**1.5), -3.0/5.0
    def _entropy(self,c):
        return 0.5-log(2)
triang = triang_gen(a=0.0, b=1.0, name="triang", longname="A Triangular",
                    shapes="c", extradoc="""

Triangular distribution

    up-sloping line from loc to (loc + c*scale) and then downsloping
    for (loc + c*scale) to (loc+scale).

    - standard form is in the range [0,1] with c the mode.
    - location parameter shifts the start to loc
    - scale changes the width from 1 to scale
"""
                    )

# Truncated Exponential

class truncexpon_gen(rv_continuous):
    def _argcheck(self, b):
        self.b = b
        return (b > 0)
    def _pdf(self, x, b):
        return exp(-x)/(1-exp(-b))
    def _cdf(self, x, b):
        return (1.0-exp(-x))/(1-exp(-b))
    def _ppf(self, q, b):
        return -log(1-q+q*exp(-b))
    def _munp(self, n, b):
        return gam(n+1)-special.gammainc(1+n,b)
    def _entropy(self, b):
        eB = exp(b)
        return log(eB-1)+(1+eB*(b-1.0))/(1.0-eB)
truncexpon = truncexpon_gen(a=0.0, name='truncexpon',
                            longname="A truncated exponential",
                            shapes="b", extradoc="""

Truncated exponential distribution
"""
                            )

# Truncated Normal

class truncnorm_gen(rv_continuous):
    def _argcheck(self, a, b):
        self.a = a
        self.b = b
        self.nb = norm._cdf(b)
        self.na = norm._cdf(a)
        return (a != b)
    def _pdf(self, x, a, b):
        return norm._pdf(x) / (self.nb - self.na)
    def _cdf(self, x, a, b):
        return (norm._cdf(x) - self.na) / (self.nb - self.na)
    def _ppf(self, q, a, b):
        return norm._ppf(q*self.nb + self.na*(1.0-q))
    def _stats(self, a, b):
        nA, nB = self.na, self.nb
        d = nB - nA
        pA, pB = norm._pdf(a), norm._pdf(b)
        mu = (pB - pA) / d
        mu2 = 1 + (a*pA - b*pB) / d - mu*mu
        return mu, mu2, None, None
truncnorm = truncnorm_gen(name='truncnorm', longname="A truncated normal",
                          shapes="a,b", extradoc="""

Truncated Normal distribution.

  The standard form of this distribution is a standard normal truncated to the
  range [a,b] --- notice that a and b are defined over the domain
  of the standard normal.  To convert clip values for a specific mean and
  standard deviation use a,b = (myclip_a-my_mean)/my_std, (myclip_b-my_mean)/my_std
"""
                          )

# Tukey-Lambda
# A flexible distribution ranging from Cauchy (lam=-1)
#   to logistic (lam=0.0)
#   to approx Normal (lam=0.14)
#   to u-shape (lam = 0.5)
#   to Uniform from -1 to 1 (lam = 1)

class tukeylambda_gen(rv_continuous):
    def _pdf(self, x, lam):
        Fx = arr(special.tklmbda(x,lam))
        Px = Fx**(lam-1.0) + (arr(1-Fx))**(lam-1.0)
        Px = 1.0/arr(Px)
        return where((lam > 0) & (abs(x) < 1.0/lam), Px, 0.0)
    def _cdf(self, x, lam):
        return special.tklmbda(x, lam)
    def _ppf(self, q, lam):
        q = q*1.0
        vals1 = (q**lam - (1-q)**lam)/lam
        vals2 = log(q/(1-q))
        return where((lam == 0)&(q==q), vals2, vals1)
    def _stats(self, lam):
        mu2 = 2*gam(lam+1.5)-lam*pow(4,-lam)*sqrt(pi)*gam(lam)*(1-2*lam)
        mu2 /= lam*lam*(1+2*lam)*gam(1+1.5)
        mu4 = 3*gam(lam)*gam(lam+0.5)*pow(2,-2*lam) / lam**3 / gam(2*lam+1.5)
        mu4 += 2.0/lam**4 / (1+4*lam)
        mu4 -= 2*sqrt(3)*gam(lam)*pow(2,-6*lam)*pow(3,3*lam) * \
               gam(lam+1.0/3)*gam(lam+2.0/3) / (lam**3.0 * gam(2*lam+1.5) * \
                                                gam(lam+0.5))
        g2 = mu4 / mu2 / mu2 - 3.0
                                               
        return 0, mu2, 0, g2
    def _entropy(self, lam):
        def integ(p):
            return log(pow(p,lam-1)+pow(1-p,lam-1))
        return scipy.integrate.quad(integ,0,1)[0]
tukeylambda = tukeylambda_gen(name='tukeylambda', longname="A Tukey-Lambda",
                              shapes="lam", extradoc="""

Tukey-Lambda distribution

   A flexible distribution ranging from Cauchy (lam=-1)
   to logistic (lam=0.0)
   to approx Normal (lam=0.14)
   to u-shape (lam = 0.5)
   to Uniform from -1 to 1 (lam = 1)
"""
                              )

# Uniform
# loc to loc + scale

class uniform_gen(rv_continuous):
    def _rvs(self):
        return rand.uniform(0.0,1.0,self._size)
    def _pdf(self, x):
        return 1.0*(x==x)
    def _cdf(self, x):
        return x
    def _ppf(self, q):
        return q
    def _stats(self):
        return 0.5, 1.0/12, 0, -1.2
    def _entropy(self):
        return 0.0
uniform = uniform_gen(a=0.0,b=1.0, name='uniform', longname="A uniform",
                      extradoc="""

Uniform distribution

   constant between loc and loc+scale
"""
                      )

# Von-Mises

# if x is not in range or loc is not in range it assumes they are angles
#   and converts them to [-pi, pi] equivalents.

eps = scipy_base.limits.double_epsilon

class vonmises_gen(rv_continuous):
    def _rvs(self, b):
        return rv._inst._von_Mises(b,size=(self._size,))
    def _pdf(self, x, b):
        x = arr(angle(exp(1j*x)))
        Px = where(b < 100, exp(b*cos(x)) / (2*pi*special.i0(b)),
                   norm.pdf(x, 0.0, sqrt(1.0/b)))
        return Px
    def _cdf(self, x, b):
        x = arr(angle(exp(1j*x)))
        eps2 = sqrt(eps)

        c_xsimple = atleast_1d((b==0)&(x==x))
        c_xiter = atleast_1d((b<100)&(b > 0)&(x==x))
        c_xnormal = atleast_1d((b>=100)&(x==x))
        c_bad = atleast_1d((b<=0) | (x != x))
    
        indxiter = nonzero(c_xiter)
        xiter = take(x, indxiter)

        vals = ones(len(c_xsimple),Num.Float)
        putmask(vals, c_bad, nan)
        putmask(vals, c_xsimple, x / 2.0/pi)
        st = sqrt(b-0.5)
        st = where(isnan(st),0.0,st)
        putmask(vals, c_xnormal, norm.cdf(x, scale=st))
        
        biter = take(atleast_1d(b)*(x==x), indxiter)
        if len(xiter) > 0:
            fac = special.i0(biter)
            x2 = xiter
            val = x2 / 2.0/ pi
            for j in range(1,501):
                trm1 = special.iv(j,biter)/j/fac
                trm2 = sin(j*x2)/pi
                val += trm1*trm2
                if all(trm1 < eps2):
                    break
            if (j == 500):
                print "Warning: did not converge..."
            put(vals, indxiter, val)        
        return vals + 0.5
    def _stats(self, b):
        return 0, None, 0, None
vonmises = vonmises_gen(name='vonmises', longname="A Von Mises",
                        shapes="b", extradoc="""

Von Mises distribution

  if x is not in range or loc is not in range it assumes they are angles
     and converts them to [-pi, pi] equivalents.

"""
                        )


## Wald distribution (Inverse Normal with shape parameter mu=1.0)

class wald_gen(invnorm_gen):
    def _rvs(self):
        return invnorm_gen._rvs(self, 1.0)
    def _pdf(self, x):
        return invnorm.pdf(x,1.0)
    def _cdf(self, x):
        return invnorm.cdf(x,1,0)
    def _stats(self):
        return 1.0, 1.0, 3.0, 15.0
wald = wald_gen(a=0.0, name="wald", longname="A Wald",
                extradoc="""

Wald distribution
"""
                )
    
## Weibull
## See Frechet

# Wrapped Cauchy

class wrapcauchy_gen(rv_continuous):
    def _argcheck(self, c):
        return (c > 0) & (c < 1)
    def _pdf(self, x, c):
        return (1.0-c*c)/(2*pi*(1+c*c-2*c*cos(x)))
    def _cdf(self, x, c):
        output = 0.0*x
        val = (1.0+c)/(1.0-c)
        c1 = x<pi
        c2 = 1-c1
        xp = extract( c1,x)
        valp = extract(c1,val)
        xn = extract( c2,x)
        valn = extract(c2,val)
        if (any(xn)):
            xn = 2*pi - xn
            yn = tan(xn/2.0)
            on = 1.0-1.0/pi*arctan(valn*yn)
            insert(output, c2, on)
        if (any(xp)):
            yp = tan(xp/2.0)
            op = 1.0/pi*arctan(valp*yp)
            insert(output, c1, op)
        return output
    def _ppf(self, q, c):
        val = (1.0-c)/(1.0+c)
        rcq = 2*arctan(val*tan(pi*q))
        rcmq = 2*pi-2*arctan(val*tan(pi*(1-q)))
        return where(q < 1.0/2, rcq, rcmq)
    def _entropy(self, c):
        return log(2*pi*(1-c*c))
wrapcauchy = wrapcauchy_gen(a=0.0,b=2*pi, name='wrapcauchy',
                            longname="A wrapped Cauchy",
                            shapes="c", extradoc="""

Wrapped Cauchy distribution
"""
                            )

### DISCRETE DISTRIBUTIONS
###

def entropy(pk,qk=None):
    """S = entropy(pk,qk=None)

    calculate the entropy of a distribution given the p_k values
    S = -sum(pk * log(pk))

    If qk is not None, then compute a relative entropy
    S = -sum(pk * log(pk / qk))

    Routine will normalize pk and qk if they don't sum to 1
    """
    pk = arr(pk)
    pk = 1.0* pk / sum(pk)
    if qk is None:
        vec = where(pk == 0, 0.0, pk*log(pk))
    else:
        qk = arr(qk)
        if len(qk) != len(pk):
            raise ValueError, "qk and pk must have same length."
        qk = 1.0*qk / sum(qk)
        # If qk is zero anywhere, then unless pk is zero at those places
        #   too, the relative entropy is infinite.
        if any(take(pk,nonzero(qk==0.0))!=0.0):
            return inf
        vec = where (pk == 0, 0.0, pk*log(pk / qk))
    return -sum(vec)
    

## Handlers for generic case where xk and pk are given

def _drv_pdf(self, xk, *args):
    try:
        return self.P[xk]
    except KeyError:
        return 0.0

def _drv_cdf(self, xk, *args):
    indx = argmax((self.xk>xk))-1
    return self.F[self.xk[indx]]

def _drv_ppf(self, q, *args):
    indx = argmax((self.qvals>=q)) 
    return self.Finv[self.qvals[indx]]

def _drv_nonzero(self, k, *args):
    return 1

def _drv_moment(self, n, *args):
    n = arr(n)
    return sum(self.xk**n[NewAxis,...] * self.pk, axis=0)

def _drv_moment_gen(self, t, *args):
    t = arr(t)
    return sum(exp(self.xk * t[NewAxis,...]) * self.pk, axis=0)

def _drv2_moment(self, n, *args):
    tot = 0.0
    diff = 1e100
    pos = self.a
    count = 0
    while (pos <= self.b) and ((pos >= (self.b + self.a)/2.0) and \
                               (diff > self.moment_tol)):
        diff = pos**n * self._pdf(pos,*args)
        tot += diff
        pos += self.inc

def _drv2_ppfsingle(self, q, *args):  # Use basic bisection algorithm
    b = self.invcdf_b
    a = self.invcdf_a
    if isinf(b):            # Be sure ending point is > q
        b = max(100*q,10)
        while 1:
            if b >= self.b: qb = 1.0; break
            qb = self._cdf(b,*args)
            if (qb < q): b += 10
            else: break
    else:
        qb = 1.0
    if isinf(a):    # be sure starting point < q
        a = min(-100*q,-10)
        while 1:
            if a <= self.a: qb = 0.0; break
            qa = self._cdf(a,*args)
            if (qa > q): a -= 10
            else: break
    else:
        qa = self._cdf(a, *args)
        
    while 1:
        if (qa == q):
            return a
        if (qb == q):
            return b
        if b == a+1:
            return b
        c = int((a+b)/2.0)
        qc = self._cdf(c, *args)            
        if (qc < q):
            a = c
            qa = qc
        elif (qc > q):
            b = c
            qb = qc
        else:
            return c                

def reverse_dict(dict):
    newdict = {}
    for key in dict.keys():
        newdict[dict[key]] = key
    return newdict

def make_dict(keys, values):
    d = {}
    for key, value in zip(keys, values):
        d[key] = value
    return d

# Must over-ride one of _pmf or _cdf or pass in
#  x_k, p(x_k) lists in initialization

class rv_discrete:
    """A generic discrete random variable.

    Discrete random variables are defined from a standard form.
    The standard form may require some other parameters to complete
    its specification.  The distribution methods also take an optional location
    parameter using loc= keyword.  The default is loc=0.  The calling form
    of the methods follow:

    generic.rvs(<shape(s)>,loc=0)
        - random variates 

    generic.pmf(x,<shape(s)>,loc=0)
        - probability mass function

    generic.cdf(x,<shape(s)>,loc=0)
        - cumulative density function

    generic.sf(x,<shape(s)>,loc=0)
        - survival function (1-cdf --- sometimes more accurate)

    generic.ppf(q,<shape(s)>,loc=0)
        - percent point function (inverse of cdf --- percentiles)

    generic.isf(q,<shape(s)>,loc=0)
        - inverse survival function (inverse of sf)

    generic.stats(<shape(s)>,loc=0,moments='mv')
        - mean('m'), variance('v'), skew('s'), and/or kurtosis('k')

    generic.entropy(<shape(s)>,loc=0)
        - entropy of the RV

    Alternatively, the object may be called (as a function) to fix
       the shape and location parameters returning a
       "frozen" discrete RV object:

    myrv = generic(<shape(s)>,loc=0)
        - frozen RV object with the same methods but holding the
            given shape and location fixed.

    You can construct an aribtrary discrete rv where P{X=xk} = pk
    by passing to the rv_discrete initialization method (through the values=
    keyword) a tuple of sequences (xk,pk) which describes only those values of
    X (xk) that occur with nonzero probability (pk).  
    """
    def __init__(self, a=0, b=inf, name=None, badvalue=None,
                 moment_tol=1e-8,values=None,inc=1,longname=None,
                 shapes=None, extradoc=None):
        if badvalue is None:
            badvalue = nan
        self.badvalue = badvalue
        self.a = a
        self.b = b
        self.invcdf_a = a
        self.invcdf_b = b
        self.name = name
        self.moment_tol = moment_tol
        self.inc = inc
        self._cdfvec = sgf(self._cdfsingle,otypes='d')
        self.return_integers = 1
        self.vecentropy = vectorize(self._entropy)

        if values is not None:
            self.xk, self.pk = values
            self.return_integers = 0
            indx = argsort(ravel(self.xk))
            self.xk = take(ravel(self.xk),indx)
            self.pk = take(ravel(self.pk),indx)
            self.a = self.xk[0]
            self.b = self.xk[-1]
            self.P = make_dict(self.xk, self.pk)
            self.qvals = scipy_base.cumsum(self.pk)
            self.F = make_dict(self.xk, self.qvals)
            self.Finv = reverse_dict(self.F)
            self._ppf = new.instancemethod(sgf(_drv_ppf,otypes='d'),
                                           self, rv_discrete)
            self._pmf = new.instancemethod(sgf(_drv_pmf,otypes='d'),
                                           self, rv_discrete)
            self._cdf = new.instancemethod(sgf(_drv_cdf,otypes='d'),
                                           self, rv_discrete)
            self._nonzero = new.instancemethod(_drv_nonzero, self, rv_discrete)
            self.generic_moment = new.instancemethod(_drv_moment,
                                                     self, rv_discrete)
            self.moment_gen = new.instancemethod(_drv_moment_gen,
                                                 self, rv_discrete)
            self.numargs=0
        else:
            self._vecppf = new.instancemethod(sgf(_drv2_ppfsingle,otypes='d'),
                                              self, rv_discrete)
            self.generic_moment = new.instancemethod(sgf(_drv2_moment,
                                                         otypes='d'),
                                                     self, rv_discrete)
            cdf_signature = inspect.getargspec(self._cdf.im_func)
            numargs1 = len(cdf_signature[0]) - 2
            pmf_signature = inspect.getargspec(self._pmf.im_func)
            numargs2 = len(pmf_signature[0]) - 2
            self.numargs = max(numargs1, numargs2)

        if longname is None:
            if name[0] in ['aeiouAEIOU']: hstr = "An "
            else: hstr = "A "
            longname = hstr + name
        if self.__doc__ is None:
            self.__doc__ = rv_discrete.__doc__
        if self.__doc__ is not None:
            self.__doc__ = self.__doc__.replace("A Generic",longname)
            if name is not None:
                self.__doc__ = self.__doc__.replace("generic",name)
            if shapes is None:
                self.__doc__ = self.__doc__.replace("<shape(s)>,","")
            else:
                self.__doc__ = self.__doc__.replace("<shape(s)>",shapes)
            ind = self.__doc__.find("You can construct an arbitrary")
            self.__doc__ = self.__doc__[:ind].strip()
            if extradoc is not None:
                self.__doc__ = self.__doc__ + extradoc

    def _rvs(self, *args):
        return self._ppf(rand.sample(self._size),*args)

    def __fix_loc(self, args, loc):
        N = len(args)
        if N > self.numargs:
            if N == self.numargs + 1 and loc is None:  # loc is given without keyword
                loc = args[-1]
            args = args[:self.numargs]
        if loc is None:
            loc = 0
        return args, loc

    def _nonzero(self, k, *args):
        return floor(k)==k
    
    def _argcheck(self, *args):
        cond = 1
        for arg in args:
            cond &= (arg > 0)
        return cond

    def _pmf(self, k, *args):
        return self.cdf(k,*args) - self.cdf(k-1,*args)

    def _cdfsingle(self, k, *args):
        m = arange(int(self.a),k+1)
        return sum(self._pmf(m,*args))

    def _cdf(self, x, *args):
        k = floor(x)
        return self._cdfvec(k,*args)
    
    def _sf(self, x, *args):
        return 1.0-self._cdf(x,*args)
        
    def _ppf(self, q, *args):
        return self._vecppf(q, *args)

    def _isf(self, q, *args):
        return self._ppf(1-q,*args)

    def _stats(self, *args):
        return None, None, None, None

    def _munp(self, n, *args):
        return self.generic_moment(n)


    def rvs(self, *args, **kwds):
        loc,size=map(kwds.get,['loc','size'])
        args, loc = self.__fix_loc(args, loc)
        cond = self._argcheck(*args)
        if not all(cond):
            raise ValueError, "Domain error in arguments."

        if size is None:
            size = 1
        else:
            self._size = product(size)
        if scipy.isscalar(size):
            self._size = size
            size = (size,)
            
        vals = reshape(self._rvs(*args),size)
        if self.return_integers:
            vals = arr(vals)
            if vals.typecode() not in scipy.typecodes['AllInteger']:
                vals = vals.astype(Num.Int)
        return vals + loc

    def pmf(self, k,*args, **kwds):
        """Probability mass function at k of the given RV.

        *args
        =====
        The shape parameter(s) for the distribution (see docstring of the
           instance object for more information)
        
        **kwds
        ======
        loc   - location parameter (default=0)
        """
        loc = kwds.get('loc')
        args, loc  = self.__fix_loc(args, loc)
        k,loc = map(arr,(k,loc))
        args = tuple(map(arr,args))
        k = arr((k-loc))
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k <= self.b) & self._nonzero(k,*args)
        cond = cond0 & cond1
        output = zeros(shape(cond),'d')
        insert(output,(1-cond0)*(cond1==cond1),self.badvalue)
        goodargs = argsreduce(cond, *((k,)+args))
        insert(output,cond,self._pmf(*goodargs))
        return output
        
    def cdf(self, k, *args, **kwds):
        """Cumulative distribution function at k of the given RV

        *args
        =====
        The shape parameter(s) for the distribution (see docstring of the
           instance object for more information)
        
        **kwds
        ======
        loc   - location parameter (default=0)
        """        
        loc = kwds.get('loc')
        args, loc = self.__fix_loc(args, loc)
        k,loc = map(arr,(k,loc))
        args = tuple(map(arr,args))
        k = arr((k-loc))
        cond0 = self._argcheck(*args)
        cond1 = (k >= self.a) & (k < self.b)
        cond2 = (k >= self.b)
        cond = cond0 & cond1
        output = zeros(shape(cond),'d')
        insert(output,(1-cond0)*(cond1==cond1),self.badvalue)
        insert(output,cond2*(cond0==cond0), 1.0)
        goodargs = argsreduce(cond, *((k,)+args))
        insert(output,cond,self._cdf(*goodargs))
        return output

    def sf(self,k,*args,**kwds):
        """Survival function (1-cdf) at k of the given RV

        *args
        =====
        The shape parameter(s) for the distribution (see docstring of the
           instance object for more information)
        
        **kwds
        ======
        loc   - location parameter (default=0)
        """        
        loc= kwds.get('loc')
        args, loc = self.__fix_loc(args, loc)
        k,loc = map(arr,(k,loc))
        args = tuple(map(arr,args))
        k = arr(k-loc)
        cond0 = self._argcheck(*args) 
        cond1 = (k >= self.a) & (k <= self.b)
        cond2 = (k < self.a) & cond0
        cond = cond0 & cond1
        output = zeros(shape(cond),'d')
        insert(output,(1-cond0)*(cond1==cond1),self.badvalue)
        insert(output,cond2,1.0)
        goodargs = argsreduce(cond, *((k,)+args))
        insert(output,cond,self._sf(*goodargs))
        return output

    def ppf(self,q,*args,**kwds):
        """Percent point function (inverse of cdf) at q of the given RV

        *args
        =====
        The shape parameter(s) for the distribution (see docstring of the
           instance object for more information)
        
        **kwds
        ======
        loc   - location parameter (default=0)
        """                
        loc = kwds.get('loc')
        args, loc = self.__fix_loc(args, loc)
        q,loc  = map(arr,(q,loc))
        args = tuple(map(arr,args))
        cond0 = self._argcheck(*args) & (loc == loc)
        cond1 = (q > 0) & (q < 1)
        cond2 = (q==1) & cond0
        cond = cond0 & cond1
        output = valarray(shape(cond),value=self.a-1)
        insert(output,(1-cond0)*(cond1==cond1), self.badvalue)
        insert(output,cond2,self.b)
        goodargs = argsreduce(cond, *((q,)+args+(loc,)))
        loc, goodargs = goodargs[-1], goodargs[:-1]
        insert(output,cond,self._ppf(*goodargs) + loc)
        return output
        
    def isf(self,q,*args,**kwds):
        """Inverse survival function (1-sf) at q of the given RV

        *args
        =====
        The shape parameter(s) for the distribution (see docstring of the
           instance object for more information)
        
        **kwds
        ======
        loc   - location parameter (default=0)
        """        

        loc = kwds.get('loc')
        args, loc = self.__fix_loc(args, loc)
        q,loc  = map(arr,(q,loc))
        args = tuple(map(arr,args))
        cond0 = self._argcheck(*args) & (loc == loc)
        cond1 = (q > 0) & (q < 1)
        cond2 = (q==1) & cond0
        cond = cond0 & cond1
        output = valarray(shape(cond),value=self.b)
        insert(output,(1-cond0)*(cond1==cond1), self.badvalue)
        insert(output,cond2,self.a-1)
        goodargs = argsreduce(cond, *((q,)+args+(loc,)))
        loc, goodargs = goodargs[-1], goodargs[:-1]
        insert(output,cond,self._ppf(*goodargs) + loc)
        return output

    def stats(self, *args, **kwds):
        """Some statistics of the given discrete RV

        *args
        =====
        The shape parameter(s) for the distribution (see docstring of the
           instance object for more information)

        **kwds
        ======
        loc     - location parameter (default=0)
        moments - a string composed of letters ['mvsk'] specifying
                   which moments to compute (default='mv')
                   'm' = mean,
                   'v' = variance,
                   's' = (Fisher's) skew,
                   'k' = (Fisher's) kurtosis.
        """        
        loc,moments=map(kwds.get,['loc','moments'])
        N = len(args)
        if N > self.numargs:
            if N == self.numargs + 1 and loc is None:  # loc is given without keyword
                loc = args[-1]
            if N == self.numargs + 2 and moments is None: # loc, scale, and moments
                loc, moments = args[-2:]
            args = args[:self.numargs]
        if loc is None: loc = 0.0
        if moments is None: moments = 'mv'
                        
        loc = arr(loc)
        args = tuple(map(arr,args))
        cond = self._argcheck(*args) & (loc==loc)

        signature = inspect.getargspec(self._stats.im_func)
        if (signature[2] is not None) or ('moments' in signature[0]):
            mu, mu2, g1, g2 = self._stats(*args,**{'moments':moments})
        else:
            mu, mu2, g1, g2 = self._stats(*args)
        if g1 is None:
            mu3 = None
        else:
            mu3 = g1*(mu2**1.5)
        default = valarray(shape(cond), self.badvalue)
        output = []

        # Use only entries that are valid in calculation
        goodargs = argsreduce(cond, *(args+(loc,)))
        loc, goodargs = goodargs[-1], goodargs[:-1]

        if 'm' in moments:
            if mu is None:
                mu = self._munp(1.0,*goodargs)
            out0 = default.copy()
            insert(out0,cond,mu+loc)
            output.append(out0)
            
        if 'v' in moments:
            if mu2 is None:
                mu2p = self._munp(2.0,*goodargs)
                if mu is None:
                    mu = self._munp(1.0,*goodargs)
                mu2 = mu2p - mu*mu
            out0 = default.copy()
            insert(out0,cond,mu2)
            output.append(out0)
            
        if 's' in moments:
            if g1 is None:
                mu3p = self._munp(3.0,*goodargs)
                if mu is None:
                    mu = self._munp(1.0,*goodargs)                    
                if mu2 is None:
                    mu2p = self._munp(2.0,*goodargs)
                    mu2 = mu2p - mu*mu
                mu3 = mu3p - 3*mu*mu2 - mu**3
                g1 = mu3 / mu2**1.5
            out0 = default.copy()
            insert(out0,cond,g1)
            output.append(out0)
                
        if 'k' in moments:
            if g2 is None:
                mu4p = self._munp(4.0,*goodargs)
                if mu is None:
                    mu = self._munp(1.0,*goodargs)                    
                if mu2 is None:
                    mu2p = self._munp(2.0,*goodargs)
                    mu2 = mu2p - mu*mu
                if mu3 is None:
                    mu3p = self._munp(3.0,*goodargs)
                    mu3 = mu3p - 3*mu*mu2 - mu**3 
                mu4 = mu4p - 4*mu*mu3 - 6*mu*mu*mu2 - mu**4
                g2 = mu4 / mu2**2.0 - 3.0
            out0 = default.copy()
            insert(out0,cond,g2)
            output.append(out0)

        if len(output) == 1:
            return output[0]
        else:
            return tuple(output)

    def moment(self, n, *args, **kwds):   # Non-central moments in standard form.
        if (floor(n) != n):
            raise ValueError, "Moment must be an integer."
        if (n < 0): raise ValueError, "Moment must be positive."
        if (n == 0): return 1.0
        if (n > 0) and (n < 5):
            signature = inspect.getargspec(self._stats.im_func)
            if (signature[2] is not None) or ('moments' in signature[0]):
                dict = {'moments':{1:'m',2:'v',3:'vs',4:'vk'}[n]}
            else:
                dict = {}
            mu, mu2, g1, g2 = self._stats(*args,**dict)
            if (n==1):
                if mu is None: return self._munp(1,*args)
                else: return mu
            elif (n==2):
                if mu2 is None: return self._munp(2,*args)
                else: return mu
            elif (n==3):
                if g1 is None or mu2 is None: return self._munp(3,*args)
                else: return g1*(mu2**1.5)
            else: # (n==4)
                if g2 is None or mu2 is None: return self._munp(4,*args)
                else: return (g2+3.0)*(mu2**2.0)
        else:
            return self._munp(n,*args)

    def freeze(self, *args, **kwds):
        return rv_frozen(self, *args, **kwds)

    def _entropy(self, *args):
        if hasattr(self,'pk'):
            return entropy(self.pk)
        else:
            mu = int(self.stats(*args, **{'moments':'m'}))
            val = self.pmf(mu,*args)
            if (val==0.0): ent = 0.0
            else: ent = -val*log(val)
            k = 1
            term = 1.0
            while (abs(term) > eps):
                val = self.pmf(mu+k,*args)
                if val == 0.0: term = 0.0
                else: term = -val * log(val)
                val = self.pmf(mu-k,*args)
                if val != 0.0: term -= val*log(val)
                k += 1
                ent += term
            return ent

    def entropy(self, *args, **kwds):
        loc= kwds.get('loc')
        args, loc = self.__fix_loc(args, loc)
        loc = arr(loc)
        args = map(arr,args)
        cond0 = self._argcheck(*args) & (loc==loc)
        output = zeros(shape(cond0),'d')
        insert(output,(1-cond0),self.badvalue)
        goodargs = argsreduce(cond0, *args)
        insert(output,cond0,self.vecentropy(*goodargs))
        return output

    def __call__(self, *args, **kwds):
        return self.freeze(*args,**kwds)
    
# Binomial

class binom_gen(rv_discrete):
    def _rvs(self, n, pr):
        return rand.binomial(n,pr,self._size)
    def _argcheck(self, n, pr):
        self.b = n
        return (n>=0) & (pr >= 0) & (pr <= 1)
    def _cdf(self, x, n, pr):
        k = floor(x)
        vals = special.bdtr(k,n,pr)
        return vals
    def _sf(self, x, n, pr):
        k = floor(x)
        return special.bdtrc(k,n,pr)
    def _ppf(self, q, n, pr):
        vals = ceil(special.bdtrik(q,n,pr))
        vals1 = vals-1
        temp = special.bdtr(vals1,n,pr)
        return where(temp >= q, vals1, vals)
    def _stats(self, n, pr):
        q = 1.0-pr
        mu = n * pr
        var = n * pr * q
        g1 = (q-pr) / sqrt(n*pr*q)
        g2 = (1.0-6*pr*q)/(n*pr*q)
        return mu, var, g1, g2
    def _entropy(self, n, pr):
        k = r_[0:n+1]
        vals = self._pmf(k,n,pr)
        lvals = where(vals==0,0.0,log(vals))
        return -sum(vals*lvals)
binom = binom_gen(name='binom',shapes="n,pr",extradoc="""

Binomial distribution

   Counts the number of successes in *n* independent
   trials when the probability of success each time is *pr*.
""")

# Bernoulli distribution

class bernoulli_gen(binom_gen):
    def _rvs(self, pr):
        return binom_gen._rvs(self, 1, pr)
    def _argcheck(self, pr):
        return (pr >=0 ) & (pr <= 1)
    def _cdf(self, x, pr):
        return binom_gen._cdf(self, x, 1, pr)
    def _sf(self, x, pr):
        return binom_gen._sf(self, x, 1, pr)
    def _ppf(self, q, pr):
        return binom_gen._ppf(self, q, 1, pr)
    def _stats(self, pr):
        return binom_gen._stats(self, 1, pr)
    def _entropy(self, pr):
        return -pr*log(pr)-(1-pr)*log(1-pr)
bernoulli = bernoulli_gen(b=1,name='bernoulli',shapes="pr",extradoc="""

Bernoulli distribution

   1 if binary experiment succeeds, 0 otherwise.  Experiment
   succeeds with probabilty *pr*.
"""
)

# Negative binomial
class nbinom_gen(rv_discrete):
    def _rvs(self, n, pr):
        return rand.negative_binomial(n, pr, self._size)
    def _argcheck(self, n, pr):
        return (n >= 0) & (pr >= 0) & (pr <= 1)
    def _cdf(self, x, n, pr):
        k = floor(x)
        return special.nbdtr(k,n,pr)
    def _sf(self, x, n, pr):
        k = floor(x)
        return special.nbdtrc(k,n,pr)
    def _ppf(self, q, n, pr):
        vals = ceil(special.nbdtrik(q,n,pr))
        vals1 = vals-1
        temp = special.nbdtr(vals1,n,pr)
        return where(temp >= q, vals1, vals)
    def _stats(self, n, pr):
        Q = 1.0 / pr
        P = Q - 1.0
        mu = n*P
        var = n*P*Q
        g1 = (Q+P)/sqrt(n*P*Q)
        g2 = (1.0 + 6*P*Q) / (n*P*Q)
        return mu, var, g1, g2
nbinom = nbinom_gen(name='nbinom', longname="A negative binomial",
                    shapes="n,pr", extradoc="""

Negative binomial distribution
"""
                    )

## Geometric distribution

class geom_gen(rv_discrete):
    def _rvs(self, pr):
        return rv._inst._geom(pr,size=(self._size,))
    def _argcheck(self, pr):
        return (pr<=1) & (pr >= 0)
    def _pmf(self, k, pr):
        return (1-pr)**k * pr
    def _cdf(self, x, pr):
        k = floor(x)
        return (1.0-(1.0-pr)**k)
    def _sf(self, x, pr):
        k = floor(x)
        return (1.0-pr)**k
    def _ppf(self, q, pr):
        vals = ceil(log(1.0-q)/log(1-pr))
        temp = 1.0-(1.0-pr)**(vals-1)
        return where((temp >= q) & (vals > 0), vals-1, vals)
    def _stats(self, pr):        
        mu = 1.0/pr
        qr = 1.0-pr
        var = qr / pr / pr
        g1 = (2.0-pr) / sqrt(qr)
        g2 = scipy.polyval([1,-6,6],pr)/(1.0-pr)
        return mu, var, g1, g2
geom = geom_gen(a=1,name='geom', longname="A geometric",
                shapes="pr", extradoc="""

Geometric distribution
"""
                )

## Hypergeometric distribution

class hypergeom_gen(rv_discrete):
    def _rvs(self, M, n, N):
        return rv._inst._hypergeom(M,n,N,size=(self._size,))
    def _argcheck(self, M, n, N):
        cond = rv_discrete._argcheck(self,M,n,N)
        cond &= (n <= M) & (N <= M)
        self.a = N-(M-n)
        self.b = min(n,N)
        return cond
    def _pmf(self, k, M, n, N):
        tot, good = M, n
        comb = scipy.comb
        bad = tot - good
        return comb(good,k) * comb(bad,N-k) / comb(tot,N)
    def _stats(self, M, n, N):
        tot, good = M, n
        n = good*1.0
        m = (tot-good)*1.0
        N = N*1.0
        tot = m+n
        p = n/tot
        mu = N*p
        var = m*n*N*(tot-N)*1.0/(tot*tot*(tot-1))
        g1 = (m - n)*(tot-2*N) / (tot-2.0)*sqrt((tot-1.0)/(m*n*N*(tot-N)))
        m2, m3, m4, m5 = m**2, m**3, m**4, m**5
        n2, n3, n4, n5 = n**2, n**2, n**4, n**5
        g2 = m3 - m5 + n*(3*m2-6*m3+m4) + 3*m*n2 - 12*m2*n2 + 8*m3*n2 + n3 \
             - 6*m*n3 + 8*m2*n3 + m*n4 - n5 - 6*m3*N + 6*m4*N + 18*m2*n*N \
             - 6*m3*n*N + 18*m*n2*N - 24*m2*n2*N - 6*n3*N - 6*m*n3*N \
             + 6*n4*N + N*N*(6*m2 - 6*m3 - 24*m*n + 12*m2*n + 6*n2 + \
                             12*m*n2 - 6*n3)
        return mu, var, g1, g2
    def _entropy(self, M, n, N):
        k = r_[N-(M-n):min(n,N)+1]
        vals = self.pmf(k,M,n,N)
        lvals = where(vals==0.0,0.0,log(vals))
        return -sum(vals*lvals)
hypergeom = hypergeom_gen(name='hypergeom',longname="A hypergeometric",
                          shapes="M,n,N", extradoc="""

Hypergeometric distribution

   Models drawing objects from a bin. 
   M is total number of objects, n is total number of Type I objects.
   RV counts number of Type I objects in N drawn without replacement from
   population.
"""
                          )

## Logarithmic (Log-Series), (Series) distribution

class logser_gen(rv_discrete):
    def _rvs(self, pr):
        return rv._inst._logser(pr,size=(self._size,))
    def _argcheck(self, pr):
        return (pr > 0) & (pr < 1)
    def _pmf(self, k, pr):
        return -pr**k * 1.0 / k / log(1-pr)
    def _stats(self, pr):
        r = log(1-pr)
        mu = pr / (pr - 1.0) / r
        mu2p = -pr / r / (pr-1.0)**2
        var = mu2p - mu*mu
        mu3p = -pr / r * (1.0+pr) / (1.0-pr)**3
        mu3 = mu3p - 3*mu*mu2p + 2*mu**3
        g1 = mu3 / var**1.5

        mu4p = -pr / r * (1.0/(pr-1)**2 - 6*pr/(pr-1)**3 + \
                          6*pr*pr / (pr-1)**4)
        mu4 = mu4p - 4*mu3p*mu + 6*mu2p*mu*mu - 3*mu**4
        g2 = mu4 / var**2 - 3.0
        return mu, var, g1, g2
logser = logser_gen(a=1,name='logser', longname='A logarithmic',
                    shapes='pr', extradoc="""

Logarithmic (Log-Series, Series) distribution
"""
                    )

## Poisson distribution

class poisson_gen(rv_discrete):
    def _rvs(self, mu):
        return rand.poisson(mu, self._size)
    def _pmf(self, k, mu):
        Pk = mu**k * exp(-mu) / arr(special.gamma(k+1))
        return Pk
    def _cdf(self, x, mu):
        k = floor(x)
        return special.pdtr(k,mu)
    def _sf(self, x, mu):
        k = floor(x)
        return special.pdtrc(k,mu)
    def _ppf(self, q, mu):
        vals = ceil(special.pdtrik(q,mu))
        temp = special.pdtr(vals-1,mu)
        return where((temp >= q), vals1, vals)
    def _stats(self, mu):
        var = mu
        g1 = 1.0/arr(sqrt(mu))
        g2 = 1.0 / arr(mu)
        return mu, var, g1, g2
poisson = poisson_gen(name="poisson", longname='A Poisson',
                      shapes="mu", extradoc="""

Poisson distribution
"""
                      )

## (Planck) Discrete Exponential 

class planck_gen(rv_discrete):
    def _argcheck(self, lambda_):
        if (lambda_ > 0):
            self.a = 0
            self.b = inf
            return 1
        elif (lambda_ < 0):
            self.a = -inf
            self.b = 0
            return 1
        return 0  # lambda_ = 0
    def _pmf(self, k, lambda_):
        fact = (1-exp(-lamba))
        return fact*exp(-lambda_(k))
    def _cdf(self, x, lambda_):
        k = floor(x)
        return 1-exp(-lambda_*(k+1))
    def _ppf(self, q, lambda_):
        val = ceil(-1.0/lambda_ * log(1-q)-1)
        return val
    def _stats(self, lambda_):
        mu = 1/(exp(lambda_)-1)
        var = exp(-lambda_)/(1-exp(-lambda_))**2
        g1 = 2*cosh(lambda_/2.0)
        g2 = 4+2*cosh(lambda_)
        return mu, var, g1, g2
    def _entropy(self, lambda_):
        l = lambda_
        C = (1-exp(-l))
        return l*exp(-l)/C - log(C)
planck = planck_gen(name='planck',longname='A discrete exponential ',
                    shapes="lambda_",
                    extradoc="""

Planck (Discrete Exponential)

    pmf is p(k; b) = (1-exp(-b))*exp(-b*k) for kb >= 0
"""
                      )

class boltzmann_gen(rv_discrete):
    def _pmf(self, k, lambda_, N):
        fact = (1-exp(-lambda_))/(1-exp(-lambda_*N))
        return fact*exp(-lambda_(k))
    def _cdf(self, x, lambda_, N):
        k = floor(x)
        return (1-exp(-lambda_*(k+1)))/(1-exp(-lambda_*N))
    def _ppf(self, q, lambda_, N):
        qnew = q*(1-exp(-lambda_*N))
        val = ceil(-1.0/lambda_ * log(1-qnew)-1)
        return val
    def _stats(self, lambda_, N):
        z = exp(-lambda_)
        zN = exp(-lambda_*N)
        mu = z/(1.0-z)-N*zN/(1-zN)
        var = z/(1.0-z)**2 - N*N*zN/(1-zN)**2
        trm = (1-zN)/(1-z)
        trm2 = (z*trm**2 - N*N*zN)
        g1 = z*(1+z)*trm**3 - N**3*zN*(1+zN)
        g1 = g1 / trm2**(1.5)
        g2 = z*(1+4*z+z*z)*trm**4 - N**4 * zN*(1+4*zN+zN*zN)
        g2 = g2 / trm2 / trm2
        return mu, var, g1, g2

boltzmann = boltzmann_gen(name='boltzmann',longname='A truncated discrete exponential ',
                    shapes="lambda_,N",
                    extradoc="""

Boltzmann (Truncated Discrete Exponential)

    pmf is p(k; b, N) = (1-exp(-b))*exp(-b*k)/(1-exp(-b*N)) for k=0,..,N-1
"""
                      )




## Discrete Uniform

class randint_gen(rv_discrete):
    def _argcheck(self, min, max):
        self.a = min
        self.b = max-1
        return (max > min)
    def _pmf(self, k, min, max):
        fact = 1.0 / (max - min)
        return fact
    def _cdf(self, x, min, max):
        k = floor(x)
        return (k-min+1)*1.0/(max-min)
    def _ppf(self, q, min, max):
        val = ceil(q*(max-min)+min)
        return val
    def _stats(self, min, max):
        m2, m1 = arr(max), arr(min)
        mu = (m2 + m1 - 1.0) / 2
        d = m2 - m1
        var = (d-1)*(d+1.0)/12.0
        g1 = 0.0
        g2 = -6.0/5.0*(d*d+1.0)/(d-1.0)*(d+1.0)
        return mu, var, g1, g2
    def rvs(self, min, max=None, size=None):
        """An array of *size* random integers >= min and < max.
    
        If max is None, then range is >=0  and < min
        """
        if max is None:
            max = min
            min = 0
        U = random(size=size)
        val = floor((max-min)*U + min)
        return arr(val).astype(Num.Int)
    def _entropy(self, min, max):
        return log(max-min)
randint = randint_gen(name='randint',longname='A discrete uniform '\
                      '(random integer)', shapes="min,max",
                      extradoc="""

Discrete Uniform

    Random integers >=min and <max.
"""
                      )


# Zipf distribution

class zipf_gen(rv_discrete):
    def _rvs(self, a):
        return rv._inst._Zipf(a, size=(self._size,))
    def _argcheck(self, a):
        return a > 1
    def _pmf(self, k, a):
        Pk = 1.0 / arr(special.zeta(a,1) * k**a)
        return Pk
    def _munp(self, n, a):
        return special.zeta(a-n,1) / special.zeta(a,1)
    def _stats(self, a):
        sv = errp(0)
        fac = arr(special.zeta(a,1))
        mu = special.zeta(a-1.0,1)/fac
        mu2p = special.zeta(a-2.0,1)/fac
        var = mu2p - mu*mu
        mu3p = special.zeta(a-3.0,1)/fac
        mu3 = mu3p - 3*mu*mu2p + 2*mu**3
        g1 = mu3 / arr(var**1.5)
        
        mu4p = special.zeta(a-4.0,1)/fac
        sv = errp(sv)
        mu4 = mu4p - 4*mu3p*mu + 6*mu2p*mu*mu - 3*mu**4
        g2 = mu4 / arr(var**2) - 3.0
        return mu, var, g1, g2
zipf = zipf_gen(a=1,name='zipf', longname='A Zipf',
                shapes="a", extradoc="""

Zipf distribution
"""
                )
                

# Discrete Laplacian

class dlaplace_gen(rv_discrete):
    def _pmf(self, k, a):
        return tanh(a/2.0)*exp(-a*abs(k))
    def _cdf(self, x, a):
        k = floor(x)
        ind = (k >= 0)
        const = exp(a)+1
        return where(ind, 1.0-exp(-a*k)/const, exp(a*(k+1))/const)
    def _ppf(self, q, a):
        const = 1.0/(1+exp(-a))
        cons2 = 1+exp(a)
        ind = q < const
        return ceil(1.0/a*where(ind, log(q*cons2)-1, -log((1-q)*cons2)))
    def _stats(self, a):
        ea = exp(-a)
        e2a = exp(-2*a)
        e3a = exp(-3*a)
        e4a = exp(-4*a)
        mu2 = 2* (e2a + ea) / (1-ea)**3.0
        mu4 = 2* (e4a + 11*e3a + 11*e2a + ea) / (1-ea)**5.0
        return 0.0, mu2, 0.0, mu4 / mu2**2.0 - 3
    def _entropy(self, a):
        return a / sinh(a) - log(tanh(a/2.0))
dlaplace = dlaplace_gen(a=-inf,
                        name='dlaplace', longname='A discrete Laplacian',
                        shapes="a", extradoc="""
                        
Discrete Laplacian distribution. 
"""
                        )

