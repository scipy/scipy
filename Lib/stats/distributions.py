# Functions to implement several important functions for 
#   for various Continous and Discrete Probability Distributions
#
# Author:  Travis Oliphant  2002
# 

from __future__ import nested_scopes
import scipy
import scipy.special as special
import Numeric
from Numeric import alltrue, where, arange, put, putmask, nonzero, \
     ravel, compress, angle, take, ones, sum
from scipy_base.fastumath import *

from scipy_base import atleast_1d, polyval
errp = special.errprint
select = scipy.select
arr = Numeric.asarray
import types
import stats as st

all = alltrue
## Special defines some of these distributions
##  using backwards argument order.  This
##  is because the underlying C-library uses this order.

## We'd also like to control what happens on domain errors in this
##   wrapper

_quantstr = "Quantile must be in [0,1]."
#_posstr = "Parameter must be > 0."
#_nonnegstr = "Parameter must be >= 0."

# A function to return nan on condition failure.
def _wc(cond, val):
    return where (cond, val, scipy.nan)

# A moment computation function -- assumes
#   parameters are valid.
def generic_moment(m, pdfunc, a, b, *args):
    def _integrand1(x):
        return x**m * apply(pdfunc, (x,)+args)    
    return scipy.integrate.quad(_integrand1, a, b)[0]

def moment_ppf(m, ppffunc, *args):
    def _integrand2(q):
        return apply(ppffunc, (q,)+args)**m
    return scipy.integrate.quad(_integrand2, 0, 1)[0]

## Internal class to compute a ppf given a distribution.
##  (needs cdf and stats function) and uses fsolve from scipy.optimize
##  to compute ppf from cdf.
class general_cont_ppf:
    def __init__(self, dist):
        self.dist = dist
        self.cdf = eval('%scdf'%dist)
        try:    
            self.stats = eval('%sstats'%dist)
        except NameError:
            self.stats = None
        self.x0 = None
        self.vecfunc = scipy.special.general_function(self._single_call)
    def _tosolve(self, x, q, *args):
        return apply(self.cdf, (x, )+args) - q
    def _single_call(self, q, *args):
        if self.x0 is None:
            if self.stats is None:
                self.stats = eval('%sstats'%self.dist)
            self.x0 = self.stats(*args)[0]
            if not isfinite(self.x0):
                self.x0 = 1.0
        return scipy.optimize.fsolve(self._tosolve, self.x0, args=(q,)+args)
    def __call__(self, q, *args):
        return self.vecfunc(q, *args)
            
### Each distribution has up to 10 functions defined plus one function
##    to return random variates following the distribution ---
##    these functions are in (this is in rv2.py or rv.py).

##  if <dist> is the name of the function to return random variates, then
##  <dist>pdf --- PDF -- Probability density function
##  <dist>cdf --- CDF -- Cumulative distribution Function
##  <dist>sf --- Survival Function 1-CDF
##  <dist>ppf --- Percent Point Function (Inverse of CDF, quantiles)
##  <dist>isf --- Inverse Survival Function (inverse of SF)
##  <dist>stats --- Return mean, variance and optionally (Fisher's) skew
##                           and kurtosis of the distribution.

##   Other things to think about supporting
##  <dist>psf --- Probability sparsity function (reciprocal of the pdf) in
##                units of percent-point-function (as a function of q).
##                Also, the derivative of the percent-point function.
##  <dist>fit --- Model-fitting (with flag for method, maximum-likelihood,
##                                  MVUB, least-squares, etc.)
##  <dist>like --- negative log-likelihood function for use with ML estimation.
##  <dist>hf -- Hazard function (PDF / SF)
##  <dist>chf -- Cumulative hazard function (-log(1-CDF)) 

##  NANs are returned for unsupported parameters.
##    location and scale parameters can be defined for each distribution.
##    The shape parameters are generally required
##
##    The order is shape parameters (if needed), loc=0.0, scale=1.0
##    These are related to the common symbols in the docs.

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


_EULER = 0.5772156649015328606   # -special.psi(1)
_ZETA3 = special.zeta(3,1)

## Kolmogorov-Smirnov one-sided and two-sided test statistics

def ksonesf(x,n):
    return special.smirnov(n,x)

def ksoneisf(p,n):
    return special.smirnovi(n,p)

def ksonecdf(x,n):
    return 1-special.smirnov(n,x)

def ksoneppf(q,n):
    assert all((0<=q) & (q<=1)), _quantstr
    return special.smirnovi(n,1-q)

def kstwosf_largen(y):
    return special.kolmogorov(y)

def kstwop_largen(p):
    assert all((0<=p)&(p<=1)), _quantstr
    return special.kolmogi(p)

def kstwocdf_largen(y):
    assert(y>=0)
    return 1-special.kolmogorov(y)

def kstwoq_largen(q):
    assert(all((0<=q) & (q<=1)))
    return special.kolmogi(1-q)

## Normal distributions

def stnormpdf(x):
    return 1.0/sqrt(2*pi)*exp(-x**2/2.0)

def stnormcdf(x):
    return special.ndtr(x)

def stnormsf(x):
    return 1-special.ndtr(x)

def stnormppf(q):
    q = arr(q)
    sv = errp(0)
    vals = where((0<=q) & (q<=1),special.ndtri(q),scipy.nan)
    sv = errp(sv)
    return vals

def stnormisf(p):
    p = arr(p)
    sv = errp(0)
    vals = where((0<=p)&(p<=1),special.ndtri(1-p),scipy.nan)
    sv = errp(sv)
    return vals

def stnormstats(full=0):
    if not full:
        return 0, 1
    return 0, 1, 0, 0

# loc = mu, scale = std

def normpdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr, (x, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Px = exp(-x*x/2.0)/sqrt(2*pi)
    return select([scale <= 0],[scipy.nan], Px/scale)

def normcdf(x, loc=0.0, scale=1.0):
    sv = errp(0)
    vals = special.ndtr((x-loc)*1.0/scale)
    sv = errp(sv)
    return vals

def normppf(q, loc=0.0, scale=1.0):
    q = arr(q)
    sv = errp(0)
    vals = where((0<=q) & (q<=1),special.ndtri(q)*scale+loc,scipy.nan)
    sv = errp(sv)    
    return vals

def normsf(x, loc=0.0, scale=1.0):
    return 1-special.ndtr((x-loc)*1.0/scale)

def normisf(p, loc=0.0, scale=1.0):
    return normppf(1-p,loc,scale)

# See Papoulis pg. 247--253
def normfit(z, loc=None, scale=None, alpha=0.05, conf=1):
    if loc is not None and scale is not None:
        raise ValueError, "Must fit at least one of location and/or scale."

    z = ravel(z)
    z = compress(isfinite(z),z)

    if loc is None and scale is None:
        muhat = st.mean(z)
        stdhat = st.std(z)
        rettup = (muhat, stdhat)
    elif loc is None:
        muhat = st.mean(z)
        rettup = (muhat,)
    else: # scale is None
        stdhat = sqrt(st.mean((z-loc)**2))
        rettup = (stdhat,)

    if conf:
        N = len(z)*1.0
        adiv2 = alpha/2.0
        if (loc is None) and (scale is None):
            critmu = tppf([adiv2, 1-adiv2], N-1)
            muconf = muhat + critmu*stdhat/sqrt(N)
            critstd = chi2ppf([1-adiv2,adiv2], N-1)
            stdconf = stdhat*sqrt((N-1)/critstd)
            rettup += (muconf, stdconf)

        elif (loc is None):   # known variance
            muconf = normppf([adiv2, 1-adiv2],loc=muhat, scale=scale/sqrt(N))
            rettup += (muconf,)
            
        else:  # scale is None,  known mean
            stdconf = stdhat*sqrt(N/chi2ppf([1-adiv2,adiv2],N))
            rettup += (stdconf,)                       
    return rettup
        
    
    

# full also returns skewness and kurtosis
def normstats(loc=0.0, scale=1.0, full=0):
    if not full:
        return loc, scale**2
    else:
        return loc, scale**2, 0, 0


## Alpha distribution
##

def alphapdf(x, a, loc=0.0, scale=1.0):
    x, a, loc, scale = map(arr,(x,a,loc,scale))
    x = arr((x-loc)/scale)
    sv = errp(0)
    Px = 1.0/arr(x**2)/special.ndtr(a)*stnormpdf(a-1.0/x)
    sv = errp(sv)
    return select([(scale <=0)|(a<=0),x>0], [scipy.nan, Px/scale]) 

def alphacdf(x, a, loc=0.0, scale=1.0):
    x, a, loc, scale = map(arr,(x,a,loc,scale))
    x = arr((x-loc)/scale)
    sv = errp(0)
    Cx = special.ndtr(a-1.0/x) / special.ndtr(a)
    sv = errp(sv)
    return select([(scale <=0)|(a<=0),x>0], [scipy.nan, Cx]) 

def alphappf(q, a, loc=0.0, scale=1.0):
    q, a, loc, scale = map(arr,(q,a,loc,scale))
    sv = errp(0)
    vals = 1.0/arr(a-special.ndtri(q*special.ndtr(a)))
    sv = errp(sv)
    cond = (a>0) & (scale > 0) & (q>=0) & (q<=1)
    return select([1-cond,q==0,q==1],[scipy.nan, 0, scipy.inf], vals*scale + loc)

def alphasf(x, a, loc=0.0, scale=1.0):    
    return 1.0 - alphacdf(x, a, loc, scale)

def alphaisf(q, a, loc=0.0, scale=1.0):
    return alphappf(1-arr(q), a, loc, scale)

def alphastats(a, loc=0.0, scale=1.0,full=0):
    # No moments
    ret = [scipy.inf]*2
    if full:
        ret += [scipy.nan]*2
    return ret

## Anglit distribution
##

def anglitpdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc)/scale)
    Px = sin(2*x + pi/2)
    return select([scale <=0,(x >= -pi/4) & (x <= pi/4)], [scipy.nan, Px/scale])

def anglitcdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc)/scale)
    Cx = sin(x+pi/4)**2
    return select([scale<=0, x>=pi/4,x > -pi/4],[scipy.nan, 1, Cx])

def anglitppf(q, loc=0.0, scale=1.0):
    q, loc, scale = map(arr, (q, loc, scale))
    return (arcsin(sqrt(q))-pi/4) * scale + loc

def anglitsf(x, loc=0.0, scale=1.0):
    return 1.0 - anglitcdf(x, loc, scale)

def anglitisf(q, loc=0.0, scale=1.0):
    return anglitppf(1.0-q, loc, scale)

def anglitstats(loc=0.0, scale=1.0, full=0):
    scale = arr(scale)
    cond = (scale>0)
    mu2 = pi*pi/16-0.5
    mn = _wc(cond, loc)
    var = _wc(cond, mu2*scale*scale)
    if not full:
        return mn, var
    g1 = 0
    g2 = -2*(pi**4 - 96) / (pi*pi-8)**2
    return mn, var, _wc(cond, g1), _wc(cond, g2)


## Arcsine distribution
##
    
def arcsinepdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc)/scale)
    Px = 1.0/pi/sqrt(x*(1-x))
    return select([scale <=0,(x>0) & (x<1)], [scipy.nan, Px/scale])

def arcsinecdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc)/scale)
    Cx = 2.0/pi*arcsin(sqrt(x))
    return select([scale <=0, x>=1,x>0],[scipy.nan, 1, Cx])

def arcsineppf(q, loc=0.0, scale=1.0):
    q, loc, scale = map(arr, (q, loc, scale))
    return sin(pi/2.0*q)**2 * scale + loc

def arcsinesf(x, loc=0.0, scale=1.0):
    return 1.0 - arcsinecdf(x, loc, scale)

def arcsineisf(q, loc=0.0, scale=1.0):
    return arcsineppf(1.0-q, loc, scale)

def arcsinestats(loc=0.0, scale=1.0, full=0):
    scale = arr(scale)
    cond = (scale>0)
    mu = 0.5
    mu2 = 3.0/8 - mu*mu
    mn = _wc(cond, mu*scale + loc)
    var = _wc(cond, mu2*scale*scale)
    if not full:
        return mn, var
    mu3 = 15.0/48 - 3*mu*mu2-mu**3
    g1 = mu3 / mu2**1.5
    mu4 = 35.0/128 - 4*mu*mu3 - 6*mu*mu*mu2 - mu**4
    g2 = mu4 / mu2**2 - 3
    return mn, var, _wc(cond, g1), _wc(cond, g2)


## Beta distribution
## 

def betapdf(x, a, b, loc=0.0, scale=1.0):
    x, loc, scale, a, b = map(arr,(x,loc,scale, a, b))
    x = arr((x-loc)/scale)
    sv = errp(0)
    Px = (1.0-x)**(b-1.0) * x**(a-1.0)
    Px /= special.beta(a,b)
    sv = errp(sv)
    vals = select([(a<=0)|(b<=0)|(scale<=0),(0<x)&(x<1)],[scipy.nan,Px/scale],0)
    return vals
    
def betacdf(x, a, b, loc=0.0, scale=1.0):
    x, loc, scale, a, b = map(arr,(x,loc,scale, a, b))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.btdtr(a,b,x)
    sv = errp(sv)
    return select([(a<=0)|(b<=0)|(scale<=0),x>1,x>0],[scipy.nan,1.0,Cx],0.)

def betappf(q, a, b, loc=0.0, scale=1.0):
    q, loc, scale, a, b = map(arr,(q,loc,scale, a, b))    
    sv = errp(0)
    vals = special.btdtri(a,b,q)
    sv = errp(sv)    
    return vals*scale + loc

def betasf(x, a, b, loc=0.0, scale=1.0):
    return 1-betacdf(x,a,b,loc,scale)

def betaisf(q, a, b, loc=0.0, scale=1.0):
    return betappf(1-arr(q),a,b,loc,scale)

def betastats(a, b, loc=0.0, scale=1.0, full=0):
    a, b, loc, scale = map(arr, (a, b, loc, scale))
    cond = (arr(a)>0) & (arr(b) > 0) & (arr(scale)>0)
    mn = _wc(cond,loc + scale*a*1.0 / (a+b))
    var = _wc(cond, scale*scale*(a*b)*1.0 / ((a+b)**2 * (a+b+1)))
    if not full:
        return mn, var
    g1 = 2.0*(b-a)*sqrt((1.0+a+b) / (a*b)) / (2+a+b)
    g2 = 6.0*(a**3 + a**2*(1-2*b) + b**2*(1+b) - 2*a*b*(2+b))
    g2 /= a*b*(a+b+2)*(a+b+3)
    return mn, var, _wc(cond, g1), _wc(cond, g2)

def betafit(z, alpha=0.05, conf=1):
    z = ravel(z)
    N = len(z)*1.0
    adiv2 = alpha/2.0
    z = compress(isfinite(z),z)
    mini = Numeric.minimum.reduce(z)
    maxi = Numeric.maximum.reduce(z)
    if (mini <=0) or (maxi >=1):
        raise ValueError, "All values must be in [0,1] range."
    rettup = (None,)
    return rettup
    

## Beta Prime

def betaprimepdf(x, a, b, loc=0.0, scale=1.0):
    x, a, b, loc, scale = map(arr, (x,a,b,loc,scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Px = 1.0/special.beta(a,b)*x**(a-1.0)/(1+x)**(a+b)
    sv = errp(sv)
    return select([(a<=0)|(b<=0)|(scale<=0),x>0],[scipy.nan, Px/scale])

def betaprimecdf(x, a, b, loc=0.0, scale=1.0):
    x, a, b, loc, scale = map(arr, (x,a,b,loc,scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    x = where(x==1.0,1.0-1e-6,x)  # hyp2f1 doesn't work at -1.0
    Cx = pow(x,a)*special.hyp2f1(a+b,a,1+a,-x)/a/special.beta(a,b)
    sv = errp(sv)
    return select([(a<=0)|(b<=0)|(scale<=0),x>0],[scipy.nan, Cx])

_betap_ppf = general_cont_ppf('betaprime')
def betaprimeppf(q, a, b, loc=0.0, scale=1.0):
    return _betap_ppf(q, a, b, loc, scale)

def betaprimesf(x, a, b, loc=0.0, scale=1.0):
    return 1.0-betaprimecdf(x,a,b,loc,scale)

def betaprimeisf(q, a, b, loc=0.0, scale=1.0):
    return betaprimeppf(1-arr(q), a, b, loc, scale)

def betaprimestats(a, b, loc=0.0, scale=1.0, full=0):
    cond = (scale <= 0) | (a <= 0)
    cond1 = (b > 1)
    cond2 = (b > 2)
    cond3 = (b > 3)
    cond4 = (b > 4)
    mu = a/(b-1.0)
    mn = select([cond,cond1], [scipy.nan, mu*scale + loc], scipy.inf)
    mu2p = mu*(a+1.0)/(b-2.0)
    mu2 = mu2p - mu*mu
    var = select([cond,cond2], [scipy.nan, mu2*scale*scale], scipy.inf)
    if not full:
        return mn, var

    mu3p = mu2p*(a+2.0)/(b-3.0)
    mu3 = mu3p - 3*mu*mu2 - mu**3
    g1 = select([cond,cond3],[scipy.nan, mu3 / mu2**1.5], scipy.inf)
    mu4 = mu3p*(a+3.0)/(b-4) -  4*mu*mu3 - 6*mu*mu*mu2 - mu**4
    g2 = select([cond,cond4],[scipy.nan, mu4 / mu2**2.0 - 3.0], scipy.inf)
    return mn, var, g1, g2
     
## Bradford
##

def bradfordpdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x,c,loc,scale))
    x = arr((x-loc*1.0)/scale)
    Px = c / (c*x + 1.0) / log(1.0+c)
    return select([(c<=0)|(scale<=0),(0<x)&(x<1)],
                  [scipy.nan,Px/scale],0)

def bradfordcdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x,c,loc,scale))
    x = arr((x-loc*1.0)/scale)
    Cx = log(1.0+c*x) / log(c+1.0)
    return select([(c<=0)|(scale<=0),x>1,x>0],
                  [scipy.nan,1,Cx],0)

def bradfordppf(al, c, loc=0.0, scale=1.0):
    al, c, loc, scale = map(arr, (al,c,loc,scale))
    val = loc + scale*((1.0+c)**al-1)/c
    return select([(c<=0)|(scale<=0)|(al>1)|(al<0)],[scipy.nan],val)

def bradfordsf(x, c, loc=0.0, scale=1.0):
    return 1.0-bradfordcdf(x,c,loc,scale)

def bradfordisf(al, c, loc=0.0, scale=1.0):
    return bradfordppf(1-al,c,loc,scale)

def bradfordstats(c, loc=0.0, scale=1.0, full=0):
    cond = (arr(c)>0)&(arr(scale)>0)
    k = log(1.0+c)
    mu = where(cond, loc + scale*(c-k)/c/k, scipy.nan)
    var = where(cond, scale**2*((c+2.0)*k-2.0*c)/(2.0*c*k**2), scipy.nan)
    if not full:
        return mu, var
    g1 = sqrt(2)*(12*c*c-9*c*k*(c+2)+2*k*k*(c*(c+3)+3))
    g1 /= sqrt(c*(c*(k-2)+2*k))*(3*c*(k-2)+6*k)
    g2 = c**3*(k-3)*(k*(3*k-16)+24)+12*k*c*c*(k-4)*(k-3) \
         + 6*c*k*k*(3*k-14) + 12*k**3
    g2 /= 3*c*(c*(k-2)+2*k)**2
    return mu, var, _wc(cond, g1),\
           _wc(cond, g2)

## Burr

# burr with d=1 is called the fisk distribution
def burrpdf(x,c,d,loc=0.0,scale=1.0):
    x, c, d, loc, scale = map(arr, (x, c, d, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Px = c*d*(x**(-c-1.0))*((1+x**(-c*1.0))**(-d-1.0))
    return select([(c<=0)|(d<=0)|(scale<=0),(x>0)],[scipy.nan,Px/scale])

def burrcdf(x,c,d,loc=0.0, scale=1.0):
    x, c, d, loc, scale = map(arr, (x, c, d, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Cx = (1+x**(-c*1.0))**(-d**1.0)
    return select([(c<=0)|(d<=0)|(scale<=0),x>0],[scipy.nan,Cx])

def burrppf(al,c,d,loc=0.0, scale=1.0):    
    al, c, d, loc, scale = map(arr, (al, c, d, loc, scale))
    cond = (al>=0) & (al<=1) & (c>0) & (d>0) & (scale>0)
    vals = arr((al**(-1.0/d)-1))**(-1.0/c)
    return where(cond,vals,scipy.nan)

def burrsf(x,c,d,loc=0.0, scale=1.0):
    return 1.0-burrcdf(x,c,d,loc,scale)

def burrisf(al,c,d,loc=0.0,scale=1.0):
    return burrppf(1.0-al,c,d,loc,scale)

def burrstats(c,d,loc=0.0,scale=1.0,full=0):
    gam = special.gamma
    cond = (arr(c)>0) & (arr(d)>0) & (arr(scale)>0)
    g2c, g2cd = gam(1-2.0/c), gam(2.0/c+d)
    g1c, g1cd = gam(1-1.0/c), gam(1.0/c+d)
    gd = gam(d)
    k = gd*g2c*g2cd - g1c**2 * g1cd**2
    mu = where(cond, g1c*g1cd / gd*scale+loc, scipy.nan)
    var = where(cond, k/gd**2*scale*scale, scipy.nan)
    if not full:
        return mu, var
    g3c, g3cd = gam(1-3.0/c), gam(3.0/c+d)
    g4c, g4cd = gam(1-4.0/c), gam(4.0/c+d)
    g1 = 2*g1c**3 * g1cd**3 + gd*gd*g3c*g3cd - 3*gd*g2c*g1c*g1cd*g2cd
    g1 /= sqrt(k**3)
    g2 = 6*gd*g2c*g2cd * g1c**2 * g1cd**2 + gd**3 * g4c*g4cd
    g2 -= 3*g1c**4 * g1cd**4 -4*gd**2*g3c*g1c*g1cd*g3cd
    return mu, var, _wc(cond, g1), _wc(cond, g2)

# Fisk distribution
# burr is a generalization

def fiskpdf(x, c, loc=0.0, scale=1.0):
    return burrpdf(x, c, 1.0, loc, scale)

def fiskcdf(x, c, loc=0.0, scale=1.0):
    return burrcdf(x, c, 1.0, loc, scale)

def fisksf(x, c, loc=0.0, scale=1.0):
    return burrsf(x, c, 1.0, loc, scale)

def fiskppf(x, c, loc=0.0, scale=1.0):
    return burrppf(x, c, 1.0, loc, scale)

def fiskisf(x, c, loc=0.0, scale=1.0):
    return burrisf(x, c, 1.0, loc, scale)

def fiskstats(x, c, loc=0.0, scale=1.0):
    return burrstats(x, c, 1.0, loc, scale)


## Cauchy

# median = loc
def cauchypdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc*1.0)/scale)
    Px = 1.0/pi/(1.0+x*x)
    return select([scale>0], [Px/scale], scipy.nan)

def cauchycdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc*1.0)/scale)
    return select([scale>0], [0.5 + 1.0/pi*arctan(x)],scipy.nan)

def cauchysf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc*1.0)/scale)
    return select([scale>0], [0.5 - 1.0/pi*arctan(x)],scipy.nan)

def cauchyppf(al, loc=0.0, scale=1.0):
    al, loc, scale = map(arr,(al,loc,scale))
    cond = ((0<al) & (al<1)) & (scale > 0)
    return select([cond,al==1,al==0], [scale*tan(pi*(al-0.5))+loc,
                                       scipy.inf,-scipy.inf], scipy.nan)

def cauchyisf(al, loc=0.0, scale=1.0):
    return cauchyppf(1-al, loc, scale)
    
def cauchystats(loc=0.0, scale=1.0, full=0):
    if not full:
        return scipy.nan, scipy.nan
    else:
        return scipy.nan, scipy.nan, scipy.nan, scipy.nan

## Chi
##   (positive square-root of chi-square)
##   chi(1, loc, scale) = halfnormal
##   chi(2, 0, scale) = Rayleigh
##   chi(3, 0, scale) = MaxWell

def chipdf(x,df,loc=0.0,scale=1.0):
    x, df, loc, scale = map(arr, (x, df, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Px = x**(df-1.)*exp(-x*x*0.5)/(2.0)**(df*0.5-1)/special.gamma(df*0.5)
    sv = errp(sv)
    return select([(df<=0)|(scale<=0),x>0],[scipy.nan,Px/scale])

def chicdf(x,df,loc=0.0,scale=1.0):
    x, df, loc, scale = map(arr, (x, df, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.gammainc(df*0.5,0.5*x*x)
    sv = errp(sv)
    return select([(df<=0)|(scale<=0),x>0],[scipy.nan,Cx])

def chippf(al,df,loc=0.0,scale=1.0):
    al, df, loc, scale = map(arr, (al, df, loc, scale))
    sv = errp(0)
    vals = sqrt(2*special.gammaincinv(df*0.5,al))*scale + loc
    sv = errp(sv)
    cond = (al>=0)&(al<=1)&(df>0)&(scale>0)
    return where(cond, vals, scipy.nan)

def chisf(x,df,loc=0.0,scale=1.0):
    return 1.0-chicdf(x,df,loc,scale)

def chiisf(al,df,loc=0.0,scale=1.0):
    return chippf(1-al,df,loc,scale)

def chistats(df,loc=0.0,scale=1.0,full=0):
    cond = (arr(df) > 0) & (arr(scale) > 0)
    sv = errp(0)
    df2 = special.gamma(df*0.5)
    mu = sqrt(2)*special.gamma((df+1.0)/2)/df2
    mu = where(cond, mu, scipy.nan)
    mu2 = (df - mu*mu)
    mn = mu*scale + loc
    var = mu2*scale*scale
    sv = errp(sv)
    if not full:
        return mn, var
    g1 = 2*mu**3 + mu*(1.0-2*df)
    g1 /= arr(mu2**1.5)
    g2 = 2*df*(1-df)-6*mu**4+4*(2*df-1)*mu*mu
    g2 /= arr(mu2**2.0)
    return mn, var, _wc(cond, g1), _wc(cond, g2)
    
## Chi-squared (gamma-distributed with loc=0 and scale=2 and shape=df/2)

def chi2pdf(x, df, loc=0.0, scale=1.0):
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc)/scale)
    df = arr(df)
    sv = errp(0)
    Px = x**(df/2.0-1)*exp(-x/2.0)
    Px /= special.gamma(df/2.0)* 2**(df/2.0)
    sv = errp(sv)
    return select([(df<=0)&(scale<=0),x>0],[scipy.nan,Px/scale])

def chi2cdf(x, df, loc=0.0, scale=1.0):
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc)/scale)
    df = arr(df)
    sv = errp(0)
    vals = select([(df<=0)&(scale<=0),x>0],[scipy.nan,special.chdtr(df, x)])
    sv = errp(sv)
    return vals

def chi2sf(x, df, loc=0.0, scale=1.0):
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc)/scale)
    df = arr(df)
    sv = errp(0)
    vals = select([(df<=0)&(scale<=0), x>0],[scipy.nan,special.chdtrc(df,x)])
    sv = errp(sv)
    return vals

def chi2isf(p, df, loc=0.0, scale=1.0):
    p, loc, scale = map(arr,(p,loc,scale))
    sv = errp(0)
    vals = where((0<=p)&(p<=1)&(df>0)&(scale>0),
                 special.chdtri(df, p),scipy.nan)
    sv = errp(sv)
    return vals*scale + loc

def chi2ppf(q, df, loc=0.0, scale=1.0):
    return chi2isf(1-arr(q), df, loc, scale)

def chi2stats(df, loc=0.0, scale=1.0, full=0):
    cond = arr(df) > 0
    mn = where(cond,df*scale+loc,scipy.nan)
    var = where(cond,2*df*scale**2,scipy.nan)
    if not full:
        return mn, var
    g1 = where(cond,2*sqrt(2.0/df),scipy.nan)
    g2 = where(cond, 12.0/df, scipy.nan)
    return mn, var, g1, g2

## Cosine
def cosinepdf(x, loc=0.0, scale=1.0):
    loc, B, x = arr(loc), arr(scale), arr(x)
    x = arr((x-loc*1.0)/B)
    Px = 1.0/2.0/pi*(1+cos(x))
    return select([B<=0,(x>=-pi)&(x<=pi)],[scipy.nan, Px/B])

def cosinecdf(x, loc=0.0, scale=1.0):
    A, B, x = arr(loc), arr(scale), arr(x)
    x = arr((x-A*1.0)/B)
    Cx = 1.0/2/pi*(pi + x + sin(x))
    return select([B<=0, x>pi, x>=-pi],[scipy.nan, 1, Cx])

def cosinesf(x, loc=0.0, scale=1.0):
    return 1.0-cosinecdf(x, loc, scale)

_cosppf = general_cont_ppf('cosine')
def cosineppf(q, loc=0.0, scale=1.0):
    return _cosppf(q, loc, scale)

def cosineisf(q, loc=0.0, scale=1.0):
    return cosineppf(1-arr(q), loc, scale)

def cosinestats(loc=0.0, scale=1.0, full=0):
    B = arr(scale)
    mu = where(B>0,loc,scipy.nan)
    var = where(B>0,(pi*pi/3.0-2)*scale*scale, scipy.nan)
    if not full:
        return mu, var
    g1 = where(B>0,0.0,scipy.nan)
    g2 = where(B>0,-6*(pi**4-90)/(5.0*(pi*pi-6)**2),scipy.nan)    
    return mu, var, g1, g2

## Double Gamma distribution

def dgammapdf(x, a, loc=0.0, scale=1.0):
    A, B, C, x = arr(loc), arr(scale), arr(a), arr(x)
    sv = errp(0)
    y = arr((x-A*1.0)/B)
    Px = 1.0/(2*B*special.gamma(C))*abs(y)**(C-1) * exp(-abs(y))
    sv = errp(sv)
    return select([(B>0) & (C>0)], [Px], scipy.nan)

def dgammacdf(x, a, loc=0.0, scale=1.0):
    A, B, C, x = arr(loc), arr(scale), arr(a), arr(x)
    z = arr((x-A*1.0)/B)
    sv = errp(0)
    fac1 = special.gammainc(C,abs(z))
    sv = errp(sv)
    return select([(B<=0)|(C<=0),x<=A], [scipy.nan, 0.5-0.5*fac1],
                  0.5+0.5*fac1)

def dgammasf(x, a, loc=0.0, scale=1.0):
    A, B, C, x = arr(loc), arr(scale), arr(a), arr(x)
    z = arr((x-A*1.0)/B)
    sv = errp(0)
    fac1 = special.gammainc(C,abs(z))
    sv = errp(sv)
    return select([(B<=0)|(C<=0),x<=A], [scipy.nan, 0.5+0.5*fac1],
                  0.5-0.5*fac1)

def dgammappf(q, a, loc=0.0, scale=1.0):
    A, B, C, q = arr(loc), arr(scale), arr(a), arr(q)
    sv = errp(0)
    fac = special.gammainccinv(C,1-abs(2*q-1))  # note the complementary inv.
    sv = errp(sv)
    return select([(B<=0) | (C<=0), q<=0.5], [scipy.nan, A-B*fac], A+B*fac)

def dgammaisf(q, a, loc=0.0, scale=1.0):
    return dgammappf(1.0-q, a, loc, scale)

def dgammastats(shape, loc=0.0, scale=1.0, full=0):
    A, B, C = arr(loc), arr(scale), arr(shape)
    cond = (B>0) & (C>0)
    mu = select([cond], [A], scipy.nan)
    var = select([cond],[C*(C+1)*B*B], scipy.nan)
    if not full:
        return mu, var
    g1 = where(cond,0.0,scipy.nan)
    g2 = (C+2)*(C+3.0)/C/(C+1.0) - 3.0
    return mu, var, g1, _wc(cond, g2)


## Double Weibull distribution
##

def dweibullpdf(x, a, loc=0.0, scale=1.0):
    A, B, C, x = arr(loc), arr(scale), arr(a), arr(x)
    x = arr((x-A*1.0)/B)
    Px = C/2.0*abs(x)**(C-1.0)*exp(-abs(x)**C)
    return select([(B>0) & (C>0)], [Px/B], scipy.nan)

def dweibullcdf(x, a, loc=0.0, scale=1.0):
    A, B, C, x = arr(loc), arr(scale), arr(a), arr(x)
    x = arr((x-A*1.0)/B)
    Cx1 = 0.5*exp(-abs(x)**C)
    return select([(B<=0)|(C<=0),x<=0], [scipy.nan, Cx1],1-Cx1)

def dweibullsf(x, a, loc=0.0, scale=1.0):
    return 1.0 - dweibullcdf(x, a, loc, scale)
    
def dweibullppf(q, a, loc=0.0, scale=1.0):
    A, B, C, q = arr(loc), arr(scale), arr(a), arr(q)
    fac = where(q<=0.5,2*q,2*q-1)
    fac = pow(arr(log(1.0/fac)),1.0/C)
    return select([(B<=0) | (C<=0), q<=0.5], [scipy.nan, A-B*fac], A+B*fac)

def dweibullisf(q, a, loc=0.0, scale=1.0):
    return dweibullppf(1.0-q, a, loc, scale)

def dweibullstats(a, loc=0.0, scale=1.0, full=0):
    A, B, C = arr(loc), arr(scale), arr(a)
    cond = (B>0) & (C>0)
    mu = select([cond], [A], scipy.nan)
    var = select([cond],[special.gamma(1+2.0/C)*B*B], scipy.nan)
    if not full:
        return mu, var
    g1 = where(cond,0.0,scipy.nan)
    gam = special.gamma
    ic = 1.0/C
    u = gam(1+ic)
    u2p = gam(1+2*ic)
    u3 = 2*u**3 - 3*u*u2p + gam(1+3*ic)
    u3p = u3 + 3*u2p*u-2*u**3
    u4 = -3*u**4 + 6*u*u*u2p - 4*u*gam(1+3*ic) + gam(1+4*ic)
    u4p = u4 + 4*u3p*u - 6*u2p*u*u + 3*u**4
    g2 = u4p / u2p**2 - 3.0
    return mu, var, _wc(cond,g1), _wc(cond,g2)


## ERLANG
##
## Special case of the Gamma distribution with shape parameter an integer.
##

def erlangpdf(x, n, loc=0.0, scale=1.0):
    x, loc, n, scale = arr(x), arr(loc), arr(n), arr(scale)
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Px = (x)**(n-1.0)*exp(-x)/special.gamma(n)
    sv = errp(sv)
    return select([(scale<=0) | (n<=0) | (floor(n)!=n), x>0],[scipy.nan, Px/scale])

def erlangcdf(x, n, loc=0.0, scale=1.0):
    x, loc, n, scale = arr(x), arr(loc), arr(n), arr(scale)
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.gdtr(1.0,n,x)
    sv = errp(sv)
    return select([(scale<=0) | (n<=0) | (floor(n)!=n), x>0],[scipy.nan,Cx])

def erlangsf(x, n, loc=0.0, scale=1.0):
    x, loc, n, scale = arr(x), arr(loc), arr(n), arr(scale)
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.gdtrc(1.0,n,x)
    sv = errp(sv)
    return select([(scale<=0) | (n<=0) | (floor(n)!=n), x>0],[scipy.nan,Cx])
    
def erlangppf(al, n, loc=0.0, scale=1.0):
    al, loc, n, scale = map(arr, (al, loc, n, scale))
    sv = errp(0)
    vals = special.gdtrix(1.0,n,al)*scale + loc
    sv = errp(sv)
    return select([(scale<=0) | (n<=0) | (floor(n)!=n) | (al<0) | (al>1),al!=1],[scipy.nan,vals], scipy.inf)

def erlangisf(al, n, loc=0.0, scale=1.0):
    return gammappf(1-al, n, loc, scale)

def erlangstats(n, loc=0.0, scale=1.0, full=0):
    a, b = arr(n)*1.0, arr(scale)
    cond = (arr(scale)>0) & (arr(a) > 0) & (floor(a) == a)
    mn = where(cond, a*b + loc, scipy.nan)
    var = where(cond, a*b*b, scipy.nan)
    if not full:
        return mn, var
    g1 = 2.0/sqrt(a)
    g2 = 6.0/a
    return mn, var, where(cond, g1,scipy.nan), where(cond,g2,scipy.nan)

       
## Exponential (gamma distributed with a=1.0, loc=loc and scale=scale)
## scale == 1.0 / lambda

def exponpdf(x, loc=0.0, scale=1.0): 
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc*1.0)/scale)
    return select([scale<=0,x>=0], [scipy.nan, exp(-x)/scale])

def exponcdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc*1.0)/scale)
    return select([scale<=0,x>=0], [scipy.nan, 1.0-exp(-x)])

def exponsf(x, loc=0.0, scale=1.0):
    return 1.0-exponcdf(x,loc, scale)

def exponppf(q, loc=0.0, scale=1.0):
    q, loc, scale = map(arr, (q, loc, scale))
    assert all((0<=q) & (q<=1)), _quantstr
    vals = -log(1.0-q)
    return select([(scale<=0)|(q>1)|(q<0)],[scipy.nan],vals*scale+loc)

def exponisf(q, loc=0.0, scale=1.0):
    return exponppf(1-arr(q), loc, scale)

def exponstats(loc=0.0, scale=1.0, full=0):
    cond = (arr(scale) > 0)
    mn = _wc(cond, scale+loc)
    var = _wc(cond, scale*scale)
    if not full:
        return mn, var
    return mn, var, _wc(cond, 2), _wc(cond, 6)


## Extreme Value Type III  (Weibull-Type)
#
#  Type III as defined by JKB

def extreme3pdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x,c,loc,scale))
    x = arr((x-loc*1.0)/scale)
    Px = c*(-x)**(c-1)*exp(-(-x)**c)
    return select([(c<=0) | (scale<=0), x<=0], [scipy.nan, Px/scale])

def extreme3cdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x,c,loc,scale))
    x = arr((x-loc*1.0)/scale)
    Cx = exp(-(-x)**c)
    return select([(c<=0) | (scale<=0), x<=0], [scipy.nan, Cx], 1)

def extreme3sf(x, c, loc=0.0, scale=1.0):
    return 1.0-extreme3cdf(x, c, loc, scale)

def extreme3ppf(q, c, loc=0.0, scale=1.0):
    q, c, loc, scale = map(arr, (q,c,loc,scale))
    vals = -(arr(-log(q)))**(1.0/c)
    cond = (q<=1) & (q>=0) & (scale>0) & (c>0)
    return _wc(cond, vals*scale + loc)

def extreme3isf(q, c, loc=0.0, scale=1.0):
    return extreme3ppf(1-arr(q), c, loc, scale)

def extreme3stats(c, loc=0.0, scale=1.0, full=0):
    c, loc, scale = map(arr, (c, loc, scale))
    cond = (c>0) & (scale>0)
    sv = errp(0)
    gm1 = special.gamma(1+1.0/c)
    gm2 = special.gamma(1+2.0/c)
    mu = -gm1
    mn = _wc(cond, mu*scale + loc)
    mu2 = gm2-mu*mu
    var = _wc(cond, mu2*scale*scale)
    if not full:
        sv = errp(sv)
        return mn, var
    gm3 = special.gamma(1+3.0/c)
    gm4 = special.gamma(1+4.0/c)
    mu3 = -gm3 - 3*mu*mu2 - mu**3
    g1 = mu3 / mu2**1.5
    mu4 = gm4 - 4*mu*mu3 - 6*mu*mu*mu2 - mu**4
    g2 = mu4 / mu2**2 - 3.0
    sv = errp(sv)
    return mn, var, _wc(cond, g1), _wc(cond, g2)

## Exponentiated Weibull

def exponweibpdf(x, a, c, loc=0.0, scale=1.0):
    x, a, c, loc, scale = map(arr, (x, a, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    exc = exp(-x**c)
    Px = a*c*(1-exc)**arr(a-1) * exc * x**arr(c-1)
    return select([(a<=0) | (c<=0) | (scale<=0), x>0],
                  [scipy.nan, Px/scale])

def exponweibcdf(x, a, c, loc=0.0, scale=1.0):
    x, a, c, loc, scale = map(arr, (x, a, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    exc = exp(-x**c)
    Cx = arr(1-exc)**a
    return select([(a<=0) | (c<=0) | (scale<=0), x>0],[scipy.nan, Cx])

def exponweibppf(q, a, c, loc=0.0, scale=1.0):
    q, a, c, loc, scale = map(arr, (q, a, c, loc, scale))
    vals = (-log(1-q**(1.0/a)))**arr(1.0/c)
    cond = (q >= 0) & (q <=1) & (scale > 0) & (a > 0) & (c > 0)
    return _wc(cond, vals*scale + loc)

def exponweibsf(x, a, c, loc=0.0, scale=1.0):
    return 1.0 - exponweibcdf(x, a, c, loc, scale)

def exponweibisf(q, a, c, loc=0.0, scale=1.0):
    return exponweibppf(1.0-q, a, c, loc, scale)

def exponweibstats(a, c, loc=0.0, scale=1.0, full=0):
    a, c, loc, scale = map(arr, (a, c, loc, scale))
    cond = (a>0) & (c>0) & (scale > 0)
    _vecfunc = special.general_function(generic_moment)
    mu = _vecfunc(1, exponweibpdf, 0, scipy.inf, a, c)
    mu2 = _vecfunc(2, exponweibpdf, 0, scipy.inf, a, c)
    mn = _wc(cond, mu*scale + loc)
    var = _wc(cond, mu2*scale*scale)
    if not full:
        return mn, var
    mu3 = _vecfunc(3, exponweibpdf, 0, scipy.inf, a, c)
    mu4 = _vecfunc(4, exponweibpdf, 0, scipy.inf, a, c)
    g1 = mu3 / mu2**1.5
    g2 = mu4 / mu2**2.0 - 3.0
    return mn, var, _wc(cond, g1), _wc(cond, g2)
                   
## Exponential Power

def exponpowpdf(x, b, loc=0.0, scale=1.0):
    x, b, loc, scale = map(arr, (x, b, loc, scale))
    x = arr((x-loc*1.0)/scale)
    xbm1 = arr(x**(b-1.0))
    xb = xbm1 * x
    Px = exp(1)*b*xbm1 * exp(xb - exp(xb))
    return select([(b<=0) | (scale<=0), x>=0], [scipy.nan, Px/scale])

def exponpowcdf(x, b, loc=0.0, scale=1.0):
    x, b, loc, scale = map(arr, (x, b, loc, scale))
    x = arr((x-loc*1.0)/scale)
    xb = arr(x**b)
    Cx = 1.0-exp(1-exp(xb))
    return select([(b<=0) | (scale<=0), x>=0], [scipy.nan, Cx])

def exponpowppf(q, b, loc=0.0, scale=1.0):
    q, b, loc, scale = map(arr, (q, b, loc, scale))
    cond = (q>=0) & (q<=1) & (b>0) & (scale > 0)
    vals = pow(log(1.0-log(1.0-q)), 1.0/b)
    return _wc(cond, vals*scale + loc)

def exponpowsf(x, b, loc=0.0, scale=1.0):
    return 1.0 - exponpowcdf(x, b, loc, scale)

def exponpowisf(q, b, loc=0.0, scale=1.0):
    return exponpowppf(1.0-q, b, loc, scale)

def exponpowstats(b, loc=0.0, scale=1.0, full=0):
    b, loc, scale = map(arr, (b, loc, scale))
    cond = (b>0) & (scale > 0)
    _vecfunc = special.general_function(generic_moment)
    mu = _vecfunc(1, exponpowpdf, 0, scipy.inf, b)
    mu2 = _vecfunc(2, exponpowpdf, 0, scipy.inf, b)
    mn = _wc(cond, mu*scale + loc)
    var = _wc(cond, mu2*scale*scale)
    if not full:
        return mn, var
    mu3 = _vecfunc(3, exponpowpdf, 0, scipy.inf, b)
    mu4 = _vecfunc(4, exponpowpdf, 0, scipy.inf, b)
    g1 = mu3 / mu2**1.5
    g2 = mu4 / mu2**2.0 - 3.0
    return mn, var, _wc(cond, g1), _wc(cond, g2)
        

## Faigue-Life (Birnbaum-Sanders)

def fatiguelifepdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Px = (x+1)/arr(2*c*sqrt(2*pi*x**3))*exp(-(x-1)**2/arr((2.0*x*c**2)))
    return select([(c<=0)|(scale<=0), x>0],[scipy.nan, Px/scale])

def fatiguelifecdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.ndtr(1.0/c*(sqrt(x)-1.0/arr(sqrt(x))))
    sv = errp(sv)
    return select([(c<=0) | (scale <=0), x>0], [scipy.nan, Cx])

def fatiguelifeppf(q, c, loc=0.0, scale=1.0):
    q, c, loc, scale = map(arr, (q, c, loc, scale))
    cond = (q >=0) & (q<=1) & (c > 0) & (scale > 0)
    sv = errp(0)
    tmp = c*special.ndtri(q)
    sv = errp(sv)
    vals = 0.25*(tmp + sqrt(tmp**2 + 4))**2
    return _wc(cond, vals*scale + loc)

def fatiguelifesf(x, c, loc=0.0, scale=1.0):
    return 1.0 - fatiguelifecdf(x, c, loc, scale)

def fatiguelifeisf(q, c, loc=0.0, scale=1.0):
    return fatiguelifeppf(1.0-q, c, loc, scale)

def fatiguelifestats(c, loc=0.0, scale=1.0, full=0):
    c, loc, scale = map(arr, (c, loc, scale))
    cond = (c > 0) & (scale > 0)
    mu = c*c/2 + 1
    mn = _wc(cond, mu*scale + loc)
    c2 = c*c
    mu2 = c2*(5.0/4*c2+1)
    var = _wc(cond, mu2*scale*scale)
    if not full:
        return mn, var
    g1 = 4*c*sqrt(11*c2+6) / (5*c2 + 4)**1.5
    g2 = 6*c2*(93*c2+41.0) / (5*c2 + 4)**2.0
    return mn, var, _wc(cond, g1), _wc(cond, g2)


## Folded Cauchy

def foldcauchypdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Px = 1.0/pi*(1.0/(1+(x-c)**2) + 1.0/(1+(x+c)**2))
    return select([(c<0) |(scale<=0), x>=0], [scipy.nan, Px/scale])

def foldcauchycdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Cx = 1.0/pi*(arctan(x-c) + arctan(x+c))
    return select([(c<0) |(scale<=0), x>=0], [scipy.nan, Cx])

_foldcauchyppf = general_cont_ppf('foldcauchy')
def foldcauchyppf(q, c, loc=0.0, scale=1.0):
    return _foldcauchyppf(q, c, loc, scale)

def foldcauchystats(c, loc=0.0, scale=1.0, full=0):
    # no moments
    if not full:
        return (scipy.inf,)*2
    else:
        return (scipy.inf,)*4

        
## F

def fpdf(x, dfn, dfd, loc=0.0, scale=1.0):
    x, dfn, dfd, loc, scale = map(arr, (x, dfn, dfd, loc, scale))
    x = arr((x-loc*1.0)/scale)
    n = arr(1.0*dfn)
    m = arr(1.0*dfd)
    sv = errp(0)
    Px = m**(m/2) * n**(n/2) * x**(n/2-1)
    Px /= (m+n*x)**((n+m)/2)*special.beta(n/2,m/2)
    sv = errp(sv)
    return select([(dfn<=0)|(dfd<=0)|(scale<=0),x>0],[scipy.nan, Px/scale])

def fcdf(x, dfn, dfd, loc=0.0, scale=1.0):
    x, dfn, dfd, loc, scale = map(arr, (x, dfn, dfd, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.fdtr(dfn, dfd, x)
    sv = errp(sv)
    return select([(dfn<=0)|(dfd<=0)|(scale<=0),x>0],[scipy.nan, Cx])

def fsf(x, dfn, dfd, loc=0.0, scale=1.0):
    x, dfn, dfd, loc, scale = map(arr, (x, dfn, dfd, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.fdtrc(dfn, dfd, x)
    sv = errp(sv)
    return select([(dfn<=0)|(dfd<=0)|(scale<=0),x>0],[scipy.nan, Cx])

def fisf(q, dfn, dfd, loc=0.0, scale=1.0):
    q, dfn, dfd, loc, scale = map(arr, (q, dfn, dfd, loc, scale))    
    cond = (0<=q) & (q<=1) & (scale > 0) & (dfn > 0) & (dfd > 0)
    sv = errp(0)
    vals = special.fdtri(dfn, dfd, p)
    sv = errp(sv)
    return _wc(cond, vals*scale + loc)
        
def fppf(q, dfn, dfd, loc=0.0, scale=1.0):
    return fisf(1-arr(q), dfn, dfd, loc, scale)

def fstats(dfn, dfd, loc=0.0, scale=1.0, full=0):
    dfn, dfd, loc, scale = map(arr, (dfn, dfd, loc, scale))
    cond1 = (dfd > 2)
    cond2 = (dfd > 4)
    v2 = arr(dfd*1.0)
    v1 = arr(dfn*1.0)
    mn = _wc(cond1, v2 / (v2-2)*scale + loc)
    var = 2*v2*v2*(v2+v1-2)/(v1*(v2-2)**2 * (v2-4))
    var = _wc(cond2, var*scale*scale)
    if not full:
        return mn, var
    cond3 = (dfd > 6)
    cond4 = (dfd > 8)
    g1 = 2*(v2+2*v1-2)/(v2-6)*sqrt((2*v2-4)/(v1*(v2+v1-2)))
    g2 = 3/(2*v2-16)*(8+g1*g1*(v2-6))
#    g2 = 12*(-16+20*m-8*m*m + m**3 + 44*n \
#             -32*m*n + 5*m*m*n-22*n*n + 5*m*n*n)
#    g2 /= n*(m-6)*(m-8)*(n+m-2)
    return mn, var, _wc(cond3,g1), _wc(cond4,g2)


## Folded Normal  
##   abs(Z) where (Z is normal with mu=L and std=S so that c=abs(L)/S)
##
##  note: regress docs have scale parameter correct, but first parameter
##    he gives is a shape parameter A = c * scale

##  Half-normal is folded normal with shape-parameter c=0.

def foldnormpdf(x, c=0.0, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Px = sqrt(2.0/pi)*cosh(c*x)*exp(-(x*x+c*c)/2.0)
    return select([(c<0) | (scale <=0), x>=0],[scipy.nan, Px/scale])

def foldnormcdf(x, c=0.0, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Cx = special.ndtr(x-c) + special.ndtr(x+c) - 1.0
    return select([(c<0) | (scale <=0), x>0],[scipy.nan, Cx])

_foldnormppf = general_cont_ppf('foldnorm')
def foldnormppf(al, c=0.0, loc=0.0, scale=1.0):
    return _foldnormppf(al, c, loc, scale)

def foldnormsf(x, c=0.0, loc=0.0, scale=1.0):
    return 1.0-foldnormcdf(x, c, loc, scale)

def foldnormisf(al, c=0.0, loc=0.0, scale=1.0):
    return foldnormppf(1.0-al, c, loc, scale)

def foldnormstats(c=0.0, loc=0.0, scale=1.0, full=0):
    cond = (arr(c)>=0) & (arr(scale) > 0)
    fac = special.erf(c/sqrt(2))
    mu = sqrt(2.0/pi)*exp(-0.5*c*c)+c*fac
    mu2 = c*c + 1 - mu*mu
    mn = _wc(cond, mu*scale + loc)
    var =_wc(cond, mu2*scale*scale)
    if not full:
        return mn, var

    c2 = c*c
    g1 = sqrt(2/pi)*exp(-1.5*c2)*(4-pi*exp(c2)*(2*c2+1.0))
    g1 += 2*c*fac*(6*exp(-c2) + 3*sqrt(2*pi)*c*exp(-c2/2.0)*fac + \
                   pi*c*(fac*fac-1))
    g1 /= pi*mu2**1.5

    g2 = c2*c2+6*c2+3+6*(c2+1)*mu*mu - 3*mu**4
    g2 -= 4*exp(-c2/2.0)*mu*(sqrt(2.0/pi)*(c2+2)+c*(c2+3)*exp(c2/2.0)*fac)
    g2 /= mu2**2.0
    return mn, var, _wc(cond, g1), _wc(cond, g2-3.0)
                   

## Extreme Value Type II or Frechet
## (defined in Regress+ documentation as Extreme LB) as
##   a limiting value distribution.
##

def frechetpdf(x, c, left=0, loc=0.0, scale=1.0):
    x, c, loc, scale, left = map(arr, (x, c, loc,scale, left))
    x = arr((x-loc*1.0)/scale)
    x = arr(where(left, -x, x))
    Px = c*pow(x,-(c+1))*exp(-x**(-c))
    return select([(scale <=0)|(c<=0),((x>0)&(left==0))|((x<0)&(left!=0))],[scipy.nan, Px/scale])

def frechetcdf(x, c, left=0, loc=0.0, scale=1.0):
    x, c, loc, scale, left = map(arr, (x, c, loc,scale, left))
    x = arr((x-loc*1.0)/scale)
    x = arr(where(left, -x, x))
    Cx = exp(-x**(-c))
    return select([(scale <=0)|(c<=0),(x>0)&(left==0), (x<0)&(left!=0), left!=0],[scipy.nan, Cx, 1-Cx, 1])

def frechetppf(q, c, left=0, loc=0.0, scale=1.0):
    q, c, loc, scale, left = map(arr, (q, c, loc,scale, left))
    q = arr(where(left,1-q,q))
    vals = pow(-log(q),-1.0/c)
    cond = (q>=0)&(q<=1)&(scale >0)&(c >0)
    return select([1-cond, left==0], [scipy.nan, vals*scale + loc], -vals*scale + loc)
    
def frechetsf(x, c, left=0, loc=0.0, scale=1.0):
    return 1.0-frechetcdf(x, c, left, loc, scale)

def frechetisf(q, c, left=0, loc=0.0, scale=1.0):
    return frechetppf(1-arr(q), c, left, loc, scale)

def frechetstats(c, left=0, loc=0.0, scale=1.0, full=0):
    c, loc, scale, left = map(arr, (c, loc, scale, left))
    cond = (c > 0) & (scale > 0)
    ic = 1.0/c
    mu = special.gamma(1-ic)
    g1c = mu
    g2c = special.gamma(1-2*ic)
    mu2 = g2c - mu*mu
    mn = select([1-cond,left==0], [scipy.nan, mu*scale + loc], -mu*scale+loc)
    var = _wc(cond, mu2*scale*scale)
    if not full:
        return mn, var

    g3c = special.gamma(1-3*ic)
    g4c = special.gamma(1-4*ic)
    g1 = (g3c - 3*mu*mu2 - mu**3) / mu2**1.5
    g2 = (g4c - 4*g3c*g1c + 3*g2c*g2c) / mu2**2 - 6
    return mn, var, _wc(cond, g1), _wc(cond, g2)
    

## Generalized Logistic
##
def genlogisticpdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Px = c*exp(-x)/(1+exp(-x))**(c+1.0)
    return select([(scale<=0)|(c<=0)],[scipy.nan], Px/scale)

def genlogisticcdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Cx = (1+exp(-x))**(-c)
    return select([(scale<=0)|(c<=0)],[scipy.nan], Cx)
    
def genlogisticsf(x, c, loc=0.0, scale=1.0):
    return 1.0-genlogisticcdf(x, c, loc, scale)

def genlogisticppf(q, c, loc=0.0, scale=1.0):
    q, c, loc, scale = map(arr, (q, c, loc, scale))
    vals = -log(pow(q,-1.0/c)-1)
    return where((c>0) & (scale>0) & (q>=0) & (q<=1), vals, scipy.nan)

def genlogisticisf(q, c, loc=0.0, scale=1.0):
    return genlogsisticppf(1-arr(q), c, loc=0.0, scale=1.0)

def genlogisticstats(c, loc=0.0, scale=1.0, full=0):
    cond = (arr(c) > 0) & (arr(scale) > 0)
    zeta = special.zeta
    mu = _EULER + special.psi(c)
    mu2 = pi*pi/6.0 + zeta(2,c)
    mn = _wc(cond, mu*scale+loc)
    var = _wc(cond, mu2*scale*scale)
    if not full:
        return mn, var
    g1 = -2*zeta(3,c) + 2*_ZETA3
    g1 /= mu2**1.5
    g2 = pi**4/15.0 + 6*zeta(4,c)
    g2 /= mu2**2.0
    return mn, var, _wc(cond, g1), _wc(cond, g2)

## Generalized Pareto

def genparetopdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Px = pow(1+c*x,arr(-1.0-1.0/c))
    upper = 1.0 / arr(abs(c))
    return select([(c==0)|(scale<=0), (x>=0) & ( ((c<0)&(x<upper)) | (c>0))],
                  [scipy.nan, Px/scale])

def genparetocdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Cx = 1.0 - pow(1+c*x,arr(-1.0/c))
    upper = 1.0 / arr(abs(c))
    return select([(c==0)|(scale<=0),
                   (x>=0)&( (c>0) | ((c<0) & (x<upper))),
                   (c<0) & (x >= upper)],
                  [scipy.nan, Cx, 1])

def genparetoppf(q, c, loc=0.0, scale=1.0):
    q, c, loc, scale = map(arr, (q, c, loc, scale))
    vals = 1.0/c * (pow(1-q, -c)-1)
    cond = (q>=0) & (q<=1) & (c!=0) & (scale > 0)
    return _wc(cond, vals*scale + loc)

def genparetosf(x, c, loc=0.0, scale=1.0):
    return 1.0-genparetocdf(x, c, loc, scale)

def genparetoisf(q, c, loc=0.0, scale=1.0):
    return genparetoppf(1-arr(q), c, loc, scale)

def genparetostats(c, loc=0.0, scale=1.0, full=0):
    mu = 1.0 / arr(1-c)
    mu2 = 2.0/ arr((1-2*c)*(1-c)) - mu*mu
    mn = select([(c==0) | (scale <=0), c<1],
                [scipy.nan, mu*scale+loc], scipy.inf)
    var = select([(c==0) | (scale <=0), c<0.5],
                [scipy.nan, mu2*scale*scale], scipy.inf)
    if not full:
        return mn, var
    mu3p = 6.0 / arr((1-c)*(1-2*c)*(1-3*c))
    mu4p = 24.0 / arr((1-c)*(1-2*c)*(1-3*c)*(1-4*c))
    mu3 = (mu3p - 3*mu*mu2 - mu**3)
    g1 = mu3 / mu2**1.5
    g2 = (mu4p - 4*mu*mu3 - 6*mu*mu*mu2 - mu**4) / mu2**2.0 - 3.0
    g1 = select([(c==0) | (scale <=0), c<1.0/3],
                [scipy.nan, g1], scipy.inf)
    g2 = select([(c==0) | (scale <=0), c<1.0/4],
                [scipy.nan, g2], scipy.inf)    
    return mn, var, g1, g2


## Generalized Exponential

def genexponpdf(x, a, b, c, loc=0.0, scale=1.0):
    x, a, b, c, loc, scale = map(arr, (x, a, b, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Px = 0
    pass



## Generalized Extreme Value
##  c=0 is just gumbel distribution. 

def genextremepdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    tmp0 = arr(1-c*x)
    tmp = pow(tmp0,1.0/c-1)
    tmp2 = tmp*tmp0
    Px = exp(-tmp2)*tmp
    Px2 = exp(-exp(-x)-x)
    limit = 1.0/c
    return select([(scale <=0), c==0,
                   ((c>0)&(x<limit)) | ((c<0) & (x > limit))],
                  [scipy.nan, Px2/scale, Px / scale])

def genextremecdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    tmp0 = arr(1-c*x)
    limit = 1.0/c
    tmp = pow(tmp0,limit)
    Cx = exp(-tmp)
    Cx2 = exp(-exp(-x))
    return select([(scale <=0), c==0,
                   ((c>0)&(x<limit)) | ((c<0) & (x>limit)),
                   (c>0)],
                  [scipy.nan, Cx2, Cx, 1])

def genextremeppf(q, c, loc=0.0, scale=1.0):
    q, c, loc, scale = map(arr, (q, c, loc, scale))
    vals = 1.0/c*(1-(-log(q))**c)
    vals2 = -log(-log(q))
    return select([(scale <=0) | (q<0) | (q>1), c!=0],
                  [scipy.nan, vals*scale+loc], vals2*scale+loc)

def genextremesf(x, c, loc=0.0, scale=1.0):
    return 1.0 - genextremecdf(x, c, loc, scale)

def genextremeisf(q, c, loc=0.0, scale=1.0):
    return genextremeppf(1-arr(q), c, loc, scale)

def genextremestats(c, loc=0.0, scale=1.0, full=0):
    c, loc, scale = map(arr, (c, loc, scale))
    sv = errp(0)
    g1c = special.gamma(1+c)
    g2c = special.gamma(1+2*c)                        
    mu = 1.0/c * (1-g1c)
    mu2 = 1.0/arr(c*c) * (1-2*g1c + g2c)
    mn = select([scale <=0, c<=-1, c!=0],
                [scipy.nan, scipy.inf, mu*scale+loc],
                _EULER*scale + loc)
    var = select([scale <=0, c<=-0.5, c!=0],
                [scipy.nan, scipy.inf, mu2*scale*scale],
                pi*pi/6*scale*scale)
    if not full:
        sv = errp(sv)
        return mn, var

    g3c = special.gamma(1+3*c)
    g4c = special.gamma(1+4*c)
    mu3p = 1.0/arr(c**3)*(1-3*g1c+3*g2c-g3c)
    mu4p = 1.0/arr(c**4)*(1-4*g1c+6*g2c-4*g3c+g4c)
    g1 = select([scale <=0, c<=-1.0/3.0, c!=0],
                [scipy.nan, scipy.inf, mu3p / mu2**1.5],
                12*sqrt(6)/pi**3 * _ZETA3)
    g2 = select([scale<=0, c<=-0.25, c!=0],
                [scipy.nan, scipy.inf, mu4p / mu2**2],
                12.0/5.0)
    sv = errp(sv)
    return mn, var, g1, g2
    

## Gamma (Use MATLAB and MATHEMATICA (b=theta=scale, a=alpha=shape) definition)

## gamma(a, loc, scale)  with a an integer is the Erlang distribution
## gamma(1, loc, scale)  is the Exponential distribution
## gamma(df/2, 0, 2) is the chi2 distribution with df degrees of freedom.

def gammapdf(x, a, loc=0.0, scale=1.0):
    x, loc, a, scale = arr(x), arr(loc), arr(a), arr(scale)
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Px = (x)**(a-1.0)*exp(-x)/special.gamma(a)
    sv = errp(sv)
    return select([(scale<=0) | (a<=0), x>0],[scipy.nan, Px/scale])

def gammacdf(x, a, loc=0.0, scale=1.0):
    x, loc, a, scale = arr(x), arr(loc), arr(a), arr(scale)
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.gdtr(1.0,a,x)
    sv = errp(sv)
    return select([(scale<=0) | (a<=0), x>0],[scipy.nan,Cx])

def gammasf(x, a, loc=0.0, scale=1.0):
    x, loc, a, scale = arr(x), arr(loc), arr(a), arr(scale)
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.gdtrc(1.0,a,x)
    sv = errp(sv)
    return select([(scale<=0) | (a<=0), x>0],[scipy.nan,Cx])
    
def gammappf(al, a, loc=0.0, scale=1.0):
    al, loc, a, scale = map(arr, (al, loc, a, scale))
    sv = errp(0)
    vals = special.gdtrix(1.0,a,al)*scale + loc
    sv = errp(sv)
    return select([(scale<=0) | (a<=0) | (al<0) | (al>1),al!=1],[scipy.nan,vals], scipy.inf)

def gammaisf(al, a, loc=0.0, scale=1.0):
    return gammappf(1-al, a, loc, scale)

def gammastats(a, loc=0.0, scale=1.0, full=0):
    a, b = arr(a)*1.0, arr(scale)
    cond = (arr(scale)>0) & (arr(a) > 0)
    mn = where(cond, a*b + loc, scipy.nan)
    var = where(cond, a*b*b, scipy.nan)
    if not full:
        return mn, var
    g1 = 2.0/sqrt(a)
    g2 = 6.0/a
    return mn, var, where(cond, g1,scipy.nan), where(cond,g2,scipy.nan)

# Generalized Gamma

def gengammapdf(x, a, c, loc=0.0, scale=1.0):
    x, a, c, loc, scale = map(arr, (x, a, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Px = abs(c)* x**(c*a-1) / special.gamma(a) * exp(-x**c)
    sv = errp(sv)
    return select([(a<=0) | (c==0) | (scale <=0), x>0],[scipy.nan, Px/scale])

def gengammacdf(x, a, c, loc=0.0, scale=1.0):
    x, a, c, loc, scale = map(arr, (x, a, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.gammainc(a, x**c) / special.gamma(a)
    sv = errp(sv)
    return select([(a<=0) | (c==0) | (scale <=0), x>0],[scipy.nan, Cx])

def gengammappf(q, a, c, loc=0.0, scale=1.0):
    q, a, c, loc, scale = map(arr, (q, a, c, loc, scale))
    sv = errp(0)
    vals = pow(special.gammaincinv(a, special.gamma(a)*q),1.0/c)
    sv = errp(sv)
    cond = (c!=0) & (a>0) & (scale > 0) & (q>=0) & (q<=1)
    return _wc(cond, vals * scale + loc)

def gengammasf(x, a, c, loc=0.0, scale=1.0):
    return 1.0 - gengammacdf(x, a, c, loc, scale)

def gengammaisf(q, a, c, loc=0.0, scale=1.0):
    return gengammacdf(1.0-q, a, c, loc, scale)

def gengammastats(a, c, loc=0.0, scale=1.0, full=0):
    cond = (c!=0) & (a>0) & (scale > 0)
    ga = special.gamma(a)
    g1c = special.gamma(a+1.0/c)
    g2c = special.gamma(a+2.0/c)
    mu = g1c / ga
    mu2 = g2c/ga - mu*mu
    mn = _wc(cond, mu*scale + loc)
    var = _wc(cond, mu2*scale*scale)
    if not full:
        return mn, var
    g3c = special.gamma(a+3.0/c)
    g4c = special.gamma(a+4.0/c)    
    mu3 = g3c/ga - 3*mu*mu2-mu**3
    mu4 = g4c/ga - 4*mu*mu3-6*mu*mu*mu2 - mu**4
    g1 = mu3/mu2**1.5
    g2 = mu4/mu2**2.0 - 3
    return mn, var, _wc(cond, g1), _wc(cond, g2)

##  Generalized Half-Logistic
##

def genhalflogisticpdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    limit = 1.0/c
    tmp = arr(1-c*x)
    tmp0 = tmp**(limit-1)
    tmp2 = tmp0*tmp
    Px = 2*tmp0 / (1+tmp2)**2
    return select([(c<=0) | (scale <=0), (x>=0)&(x<=limit)],
                  [scipy.nan, Px/scale])

def genhalflogisticcdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    limit = 1.0/c
    tmp = arr(1-c*x)
    tmp2 = tmp**(limit)
    Cx = (1.0-tmp2) / (1+tmp2)
    return select([(c<=0) | (scale <=0), x>=limit, x>=0],
                  [scipy.nan, 1, Cx])

def genhalflogisticppf(q, c, loc=0.0, scale=1.0):
    q, c, loc, scale = map(arr, (q, c, loc, scale))
    vals = 1.0/c*(1-((1.0-q)/(1.0+q))**c)
    cond = (q>=0) & (q<=1) & (c>0) & (scale > 0)
    return _wc(cond,vals*scale + loc)

def genhalflogisticsf(x, c, loc=0.0, scale=1.0):
    return 1.0-genhalflogisticcdf(x, c, loc, scale)

def genhalflogisticisf(q, c, loc=0.0, scale=1.0):
    return genhalflogisticppf(1.0-q, c, loc, scale)

def _genhalflogstats(c, loc=0.0, scale=1.0, full=0):
    if (c<=0) or (scale <=0):
        if full:
            return (scipy.nan,)*2
        else:
            return (scipy.nan,)*4
    mu = generic_moment(1, genhalflogisticpdf, 0, 1.0/c, c)
    mu2p = generic_moment(2, genhalflogisticpdf, 0, 1.0/c, c)
    mu2 = mu2p - mu*mu
    mn = mu*scale + loc
    var = mu2*scale*scale
    if not full:
        return mn, var
    mu3p = generic_moment(3, genhalflogisticpdf, 0, 1.0/c, c)
    mu4p = generic_moment(4, genhalflogisticpdf, 0, 1.0/c, c)
    mu3 = mu3p - 3*mu*mu2 - mu**3
    g1 = mu3 / mu2**1.5
    mu4 = mu4p -  4*mu*mu3 - 6*mu*mu*mu2 - mu**4
    g2 = mu4 / mu2**2.0 - 3.0
    return mn, var, g1, g2

_vgenhalflogst = special.general_function(_genhalflogstats)    

def genhalflogisticstats(c, loc=0.0, scale=1.0, full=0):
    if not isinstance(full, types.IntType):
        raise ValueError, "Full must be an integer."
    return _vgenhalflogst(c, loc, scale, full)
    


## Gompertz (Truncated Gumbel)
##  Defined for x>=0

def gompertzpdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    ex = exp(x)
    Px = c*ex*exp(-c*(ex-1))
    return select([(c<=0) | (scale <=0), (x>=0)&(x<100)],
                  [scipy.nan, Px/scale])

def gompertzcdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    ex = exp(x)
    Cx = 1.0-exp(-c*(ex-1))
    return select([(c<=0) | (scale <=0), (x>=0)],
                  [scipy.nan, Cx])

def gompertzppf(q, c, loc=0.0, scale=1.0):
    q, c, loc, scale = map(arr, (q, c, loc, scale))
    cond = (q>=0) & (q<=1) & (c > 0) & (scale > 0)
    vals = log(1.0-1.0/c*log(1.0-q))
    return _wc(cond, vals*scale + loc)

def gompertzsf(x, c, loc=0.0, scale=1.0):
    return 1.0-gompertzcdf(x, c, loc, scale)

def gompertzisf(q, c, loc=0.0, scale=1.0):
    return gompertzppf(1.0-q, c, loc, scale)

def _gompstats(c, loc, scale, full):
    if (c<=0) or (scale<=0):
        return (scipy.nan,)*(2*(full!=0)+2)
    mu = generic_moment(1, gompertzpdf, 0, scipy.inf, c)
    mu2p = generic_moment(2, gompertzpdf, 0, scipy.inf, c)
    mu2 = mu2p - mu*mu
    mn = mu*scale + loc
    var = mu2*scale*scale
    if not full:
        return mn, var
    mu3p = generic_moment(3, gompertzpdf, 0, scipy.inf, c)
    mu4p = generic_moment(4, gompertzpdf, 0, scipy.inf, c)
    mu3 = mu3p - 3*mu*mu2 - mu**3
    g1 = mu3 / mu2**1.5
    mu4 = mu4p -  4*mu*mu3 - 6*mu*mu*mu2 - mu**4
    g2 = mu4 / mu2**2.0 - 3.0
    return mn, var, g1, g2

_vgompstats = special.general_function(_gompstats)
    
def gompertzstats(c, loc=0.0, scale=1.0, full=0):
    if not isinstance(full, types.IntType):
        raise ValueError, "Full must be an integer."
    return _vgompstats(c, loc, scale, full)
    
    

## Gumbel, Log-Weibull, Fisher-Tippett, Gompertz
## if left is non-zero, reconstruct the left-skewed gumbel distribution. 

def gumbelpdf(x,left=0,loc=0.0,scale=1.0):
    x, a, b = map(arr, (x, loc, scale))
    x = (x-a*1.0)/b
    x = arr(where(left, -x,x))
    fac = x+exp(-x)
    return select([scale>0],[exp(-fac)/b],scipy.nan)

def gumbelcdf(x,left=0,loc=0.0,scale=1.0):
    x, a, b = map(arr, (x, loc, scale))
    x = (x-a*1.0)/b
    return select([scale<=0, left==0],[scipy.nan, exp(-exp(-x))],1-exp(-exp(x)))
    
def gumbelsf(x,left=0,loc=0.0,scale=1.0):
    return 1.0-gumbelcdf(x,left,loc,scale)

def gumbelppf(q,left=0,loc=0.0,scale=1.0):
    q = arr(q)
    cond = (arr(scale)>0) & (q >= 0) & (q <=1)
    q = arr(where(left, 1-q, q))
    vals = -log(-log(q))
    return select([1-cond,left==0],[scipy.nan, loc+scale*vals], loc-scale*vals)

def gumbelisf(p,left=0,loc=0.0,scale=1.0):
    return gumbelppf(1-p,left,loc,scale)

def gumbelstats(left=0, loc=0.0,scale=1.0,full=0):
    a, b, left = map(arr, (loc, scale, left))
    cond = (b > 0)
    mn = select([1-cond,left==0],[scipy.nan, a + b*_EULER],a-b*_EULER)
    var = _wc(cond, pi*pi/6*b*b)
    if not full:
        return mn, var
    g1 = 12*sqrt(6)/pi**3 * _ZETA3
    g2 = 12.0/5.0
    return mn, var, _wc(cond, g1), _wc(cond, g2)

# Half-Cauchy

def halfcauchypdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc*1.0)/scale)
    Px = 2.0/pi/(1.0+x*x)
    return select([scale<=0, x>=0], [scipy.nan, Px/scale])

def halfcauchycdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr,(x,loc,scale))
    x = arr((x-loc*1.0)/scale)
    return select([scale<=0, x>=0], [scipy.nan, 2.0/pi*arctan(x)])

def halfcauchysf(x, loc=0.0, scale=1.0):
    return 1.0 - halfcauchycdf(x, loc, scale)

def halfcauchyppf(al, loc=0.0, scale=1.0):
    al, loc, scale = map(arr,(al,loc,scale))
    cond = ((0<al) & (al<1)) & (scale > 0)
    return select([cond,al==1,al==0], [scale*tan(pi/2*al)+loc,
                                       scipy.inf,-scipy.inf], scipy.nan)

def halfcauchyisf(al, loc=0.0, scale=1.0):
    return halfcauchyppf(1-al, loc, scale)
    
def halfcauchystats(loc=0.0, scale=1.0, full=0):
    if not full:
        return scipy.nan, scipy.nan
    else:
        return scipy.nan, scipy.nan, scipy.nan, scipy.nan


## Half-Logistic
##  

def halflogisticpdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr, (x, loc, scale))
    x = arr((x-loc*1.0)/scale)
    ex = exp(-x)
    Px = 2*ex / (1+ex)**2
    return select([scale <=0, x>=0], [scipy.nan, Px/scale])

def halflogisticcdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr, (x, loc, scale))
    x = arr((x-loc*1.0)/scale)
    ex = exp(-x)
    Cx = (1-ex)/(1+ex)
    return select([scale <=0, x>=0], [scipy.nan, Cx])

def halflogisticppf(q, loc=0.0, scale=1.0):
    q, loc, scale = map(arr, (q, loc, scale))
    vals = log((1.0+q)/(1-q))
    cond = (scale > 0) & (q<=1) & (q>=0)
    return _wc(cond, vals*scale+loc)

def halflogisticsf(x, loc=0.0, scale=1.0):
    return 1.0-halflogisticcdf(x, loc, scale)

def halflogisticisf(q, loc=0.0, scale=1.0):
    return halflogisticppf(1.0-q, loc, scale)

def halflogisticstats(loc=0.0, scale=1.0, full=0):
    loc, scale = arr(loc), arr(scale)
    mu = 2*log(2)
    mu2p = pi*pi/3.0
    mu2 = mu2p - mu*mu
    cond = scale > 0
    mn = select([1-cond],[scipy.nan],mu*scale + loc)
    var = select([1-cond], [scipy.nan], mu2*scale*scale)
    if not full:
        return mn, var
    mu3p = 9*(special.zetac(3.0)+1)
    mu4p = 7*pi**4/15.0
    mu3 = mu3p - 3*mu*mu2 - mu**3
    mu4 = mu4p - 4*mu*mu3 - 6*mu**2 * mu2 - mu**4
    g1 = mu3 / mu2**1.5
    g2 = mu4 / mu2**2.0 - 3.0
    return mn, var, _wc(cond, g1), _wc(cond, g2)


## Half-normal = chi(1, loc, scale)

def halfnormpdf(x, loc=0.0, scale=1.0):
    x, loc, scale = arr(x), arr(loc), arr(scale)
    x = arr((x-loc*1.0)/scale)
    Px = sqrt(2.0/pi)*exp(-x*x/2.0)
    return select([scale<=0,x>0],[scipy.nan, Px/scale])

def halfnormcdf(x, loc=0.0, scale=1.0):
    x, loc, scale = arr(x), arr(loc), arr(scale)
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.ndtr(x)*2-1.0
    sv = errp(sv)
    return select([scale<=0,x>0],[scipy.nan, Cx])

def halfnormppf(al, loc=0.0, scale=1.0):
    al, loc, scale = arr(al), arr(loc), arr(scale)
    sv = errp(0)
    vals = special.ndtri((1+al)/2.0)
    sv = errp(sv)
    cond = (al>=0)&(al<=1)&(scale>0)
    return where(cond, vals*scale + loc, scipy.nan)

def halfnormsf(x, loc=0.0, scale=1.0):
    return 1.0-halfnormcdf(x, loc, scale)

def halfnormisf(al, loc=0.0, scale=1.0):
    return halfnormppf(1.0-al,loc,scale)

def halfnormstats(loc=0.0, scale=1.0, full=0):
    cond = arr(scale) > 0
    mu = where(cond, scale*sqrt(2.0/pi)+loc, scipy.nan)
    var = where(cond, scale**2*(1-2.0/pi), scipy.nan)
    if not full:
        return mu, var
    g1 = sqrt(2)*(4-pi)/(pi-2.0)**1.5
    g2 = 8*(pi-3)/(pi-2.0)**2
    return mu, var, _wc(cond, g1), _wc(cond, g2)


## Hyperbolic Secant

def hypsecantpdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr, (x, loc, scale))
    x = arr((x-loc*1.0)/scale)
    return select([scale>0],[1.0/(pi*cosh(x))/scale],scipy.nan)

def hypsecantcdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr, (x, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Cx = 2.0/pi * arctan(exp(x))
    return select([scale>0],[Cx],scipy.nan)

def hypsecantppf(q, loc=0.0, scale=1.0):
    q = arr(q)
    cond = (q>=0) & (q<1) & (scale > 0)
    return select([cond,q==1],[log(tan(pi*q/2.0))*scale+loc,scipy.inf],scipy.nan)

def hypsecantsf(x, loc=0.0, scale=1.0):
    return 1.0-hypsecantcdf(x, loc, scale)

def hypsecantisf(q, loc=0.0, scale=1.0):
    return hypsecantppf(1-arr(q), loc, scale)

def hypsecantstats(loc=0.0, scale=1.0, full=0):
    cond = (arr(scale)>0)
    mn = _wc(cond, loc)
    var = _wc(cond, pi*pi/4*scale*scale)
    if not full:
        return mn, var
    g1, g2 = 0, 2
    return mn, var, _wc(cond, g1), _wc(cond, g2)

##  Inverted Gamma
#     special case of generalized gamma with c=-1
# 

def invgammapdf(x, a, loc=0.0, scale=1.0):
    return gengammapdf(x, a, -1, loc, scale)

def invgammacdf(x, a, loc=0.0, scale=1.0):
    return gengammacdf(x, a, -1, loc, scale)

def invgammappf(q, a, loc=0.0, scale=1.0):
    return gengammappf(q, a, -1, loc, scale)

def invgammasf(x, a, loc=0.0, scale=1.0):
    return gengammasf(x, a, -1, loc, scale)

def invgammaisf(q, a, loc=0.0, scale=1.0):
    return gengammaisf(q, a, -1, loc, scale)

def invgammastats(a, loc=0.0, scale=1.0, full=0):
    return gengammastats(a, -1, loc, scale, full=full)


## Inverse Normal Distribution
# scale is gamma from DATAPLOT and B from Regress

def invnormpdf(x, mu, loc=0.0, scale=1.0):
    x, mu, loc, scale = map(arr,(x, mu, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Px = 1.0/sqrt(2*pi*x**3.0)*exp(-1.0/(2*x)*((x-mu)/mu)**2)
    return select([(mu<=0)|(scale<=0),x>0],[scipy.nan,Px/scale])

def invnormcdf(x, mu, loc=0.0, scale=1.0):
    x, mu, loc, scale = map(arr,(x, mu, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = special.errprint(0)
    fac = sqrt(1.0/x)
    C1 = stnormcdf(fac*(x-mu)/mu)
    C1 += exp(2.0/mu)*stnormcdf(-fac*(x+mu)/mu)
    sv = special.errprint(sv)
    return select([(mu<=0)|(scale<=0),x>0],[scipy.nan, C1])

def invnormsf(x, mu, loc=0.0, scale=1.0):
    return 1-invnormcdf(x, mu, loc, scale)

def _invnormqfunc(x,q,mu,loc,scale):
    return invnormcdf(x,mu,loc,scale)-q

def _invnormppf(q,mu,loc,scale,x0):
    import scipy.optimize as optimize
    return optimize.fsolve(_invnormqfunc,x0,args=(q,mu,loc,scale))
_vec_invnormppf = special.general_function(_invnormppf,'d')

def invnormppf(q, mu, loc=0.0, scale=1.0, x0=None):
    mu, loc, scale, q = map(arr,(mu, loc, scale, q))
    cond = (mu > 0) & (scale > 0) & (0<=q) & (q<=1)
    if x0 is None:
        x0 = mu
    qvals = _vec_invnormppf(q, mu, loc, scale, x0)
    qvals = where(qvals<0,0,qvals)
    return where(cond, qvals, scipy.nan)

def invnormisf(q, mu, loc=0.0, scale=1.0):
    return invnormppf(1-arr(q), mu, loc, scale)

def invnormstats(mu, loc=0.0, scale=1.0, full=0):
    cond = (arr(mu)>0) & (arr(scale) > 0)
    mn = _wc(cond, mu*scale + loc)
    var = _wc(cond, mu**3 * scale**2)
    if not full:
        return mn, var
    g1 = 3*sqrt(mu)
    g2 = 15*mu
    return mn, var, _wc(cond,g1), _wc(cond,g2)


## Inverted Weibull

def invweibullpdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    xc1 = x**(-c-1.0)
    xc2 = xc1*x
    Px = c*xc1*xc2
    return select([(c<=0) | (scale<=0), x>0],
                  [Px/scale])

def invweibullcdf(x, c, loc=0.0,  scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    xc1 = x**(-c)
    Cx = exp(-xc1)
    return select([(c<=0) | (scale<=0), x>0],
                  [Cx])

def invweibullppf(q, c, loc=0.0, scale=1.0):
    q, c, loc, scale = map(arr, (q, c, loc, scale))
    vals = pow(-log(q),arr(-1.0/c))
    cond = (q >=0) & (q<=1) & (c > 0 ) & (scale > 0)
    return _wc(cond, vals*scale+loc)

def invweibullisf(q, c, loc=0.0, scale=1.0):
    return invweibullppf(1.0-q, c, loc, scale)

def invweibullsf(x, c, loc=0.0, scale=1.0):
    return 1.0-invweibullcdf(x, c, loc, scale)

def invweibullstats(x, c, loc=0.0, scale=1.0):
    raise NotImplementedError, "Not Implemented."
    

## Johnson SB

def johnsonsbpdf():
    pass

## Johnson SU

def johnsonsupdf():
    pass

## Laplace Distribution

def laplacepdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr, (x, loc, scale))
    x = arr((x-loc*1.0)/scale)
    return select([scale>0],[0.5*exp(-abs(x))/scale],scipy.nan)
    
def laplacecdf(x, loc=0.0, scale=1.0):
    x, loc, scale = map(arr, (x, loc, scale))
    x = arr((x-loc*1.0)/scale)
    return select([scale<=0,x<0],[scipy.nan, 0.5*exp(x)],1-0.5*exp(-x))

def laplacesf(x, loc=0.0, scale=1.0):
    return 1.0-laplacecdf(x, loc, scale)

def laplaceppf(q, loc=0.0, scale=1.0):
    q, loc, scale = map(arr, (q, loc, scale))
    cond = (q<0) | (q>1) | (scale<=0)
    vals1 = log(2*q)*scale + loc
    vals2 = -log(2-2*q)*scale + loc
    return select([cond,q<0.5],[scipy.nan,vals1],vals2)

def laplaceisf(q, loc=0.0, scale=1.0):
    return laplaceppf(1-arr(q), loc, scale)

def laplacestats(loc=0.0, scale=1.0, full=0):
    loc, scale = arr(loc), arr(scale)
    cond = (scale>0)
    mn = _wc(cond, loc)
    var = _wc(cond, 2.0*scale**2)
    if not full:
        return mn, var
    g1 = 0
    g2 = 3
    return mn, var, _wc(cond, g1), _wc(cond, g2)

## Logistic

def logisticpdf(x, loc=0.0, scale=1.0):
    iscale = 1.0/scale
    fac1 = exp((x-loc)*iscale)
    Px = iscale*fac1
    Px /= (1+fac1)**2
    return select([scale>0],[Px],scipy.nan)

def logisticcdf(x, loc=0.0, scale=1.0):
    scale = arr(scale)
    fac1 = exp((loc-x)*1.0/scale)
    return select([scale>0],[1.0/(1+fac1)],scipy.nan)

def logisticsf(x, loc=0.0, scale=1.0):
    return 1.0-logisticcdf(x,loc=loc,scale=scale)

def logisticppf(q, loc=0.0, scale=1.0):
    q, scale = arr(q), arr(scale)
    cond = (0<=q) & (q<=1) & (scale > 0)
    vals = loc - scale*log(1.0/q-1)
    return where(cond, vals, scipy.nan)

def logisticisf(q, loc=0.0, scale=1.0):
    return logisticppf(1-arr(q), loc, scale)

def logisticstats(loc=0.0, scale=1.0, full=0):
    mn = loc
    var = pi*pi/3.0*scale*scale
    if not full:
        return mn, var
    g1 = 0.0
    g2 = 6.0/5
    return mn, var, g1, g2

## Log Gamma
#

def loggammapdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Px = exp(c*x-exp(x))/ special.gamma(c)
    sv = errp(sv)
    return select([(scale <=0) | (c<=0)],[scipy.nan], Px / scale)

def loggammacdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.gammainc(c, exp(x))/ special.gamma(c)
    sv = errp(sv)
    return select([(scale <=0) | (c<=0)],[scipy.nan], Cx)

def loggammappf(q, c, loc=0.0, scale=1.0):
    q, c, loc, scale = map(arr, (q, c, loc, scale))
    sv = errp(0)
    vals = log(special.gammaincinv(c,q*special.gamma(c)))
    sv = errp(sv)
    return vals*scale + loc

def loggammaisf(q, c, loc=0.0, scale=1.0):
    return loggammappf(1.0-q, c, loc, scale)

def loggammasf(x, c, loc=0.0, scale=1.0):
    return 1.0-loggammacdf(x, c, loc, scale)

def loggammastats(c, loc=0.0, select=1.0, full=1):
    pass
        

## Log-Laplace
##

def loglaplacepdf():
    pass


## Lognormal
## std is a shape parameter and is the variance of the underlying
##    distribution.
## the mean of the underlying distribution is log(scale)

def lognormpdf(x, std=1.0, loc=0.0, scale=1.0):
    x, std, loc, scale = map(arr, (x, std, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Px = exp(-log(x)**2 / (2*std**2))
    Px /= std*x*sqrt(2*pi)
    return select([(std<=0)|(scale<=0),x>0],[scipy.nan,Px/scale])

def lognormcdf(x, std=1.0, loc=0.0, scale=1.0):
    x, std, loc, scale = map(arr, (x, std, loc, scale))
    x = arr((x-loc*1.0)/scale)
    cond = (scale > 0) & (std > 0)
    sv = errp(0)
    Cx = normcdf(log(x)/std)
    sv = errp(sv)
    return select([cond],[Cx],scipy.nan)

def lognormsf(x, std=1.0, loc=0.0, scale=1.0):
    return 1-lognormcdf(x, std, loc, scale)

def lognormppf(q, std=1.0, loc=0.0, scale=1.0):
    sv = errp(0)
    vals = exp(normppf(q)*std)*scale + loc
    sv = errp(sv)
    return vals

def lognormisf(p, std=1.0, loc=0.0, scale=1.0):
    sv = errp(0)
    vals = exp(normisf(p)*std)*scale + loc
    sv = errp(sv)
    return vals

def lognormstats(std=1.0, loc=0.0, scale=1.0, full=0):
    cond = (arr(std)>0) & (arr(scale)>0)
    s2 = std*std
    mn = _wc(cond, exp(s2/2.0)*scale + loc)
    fac1 = exp(s2)
    fac2 = fac1-1
    var = _wc(cond, fac1*fac2*scale*scale)
    if not full:
        return mn, var
    g1 = sqrt(fac2)*(2+fac1)
    g2 = polyval([1,2,3,0,-6],fac1)
    return mn, var, _wc(cond, g1), _wc(cond, g2)
    
# Gibrat's distribution is just default lognormal

def gilbratpdf(x):
    return lognormpdf(x)

def gilbratcdf(x):
    return lognormcdf(x)

def gilbratsf(x):
    return lognormsf(x)

def gilbratppf(x):
    return lognormppf(x)

def gilbratisf(x):
    return lognormisf(x)

def gilbratstats(full=0):
    return lognormstats(full=full)

# LOMAX (Pareto of the second kind.)
#  Special case of Pareto of the first kind.

def lomaxpdf(x, c, scale=1.0):
    return paretopdf(x, c, -1, scale)

def lomaxcdf(x, c, scale=1.0):
    return paretocdf(x, c, -1, scale)

def lomaxppf(q, c, scale=1.0):
    return paretoppf(q, c, -1, scale)


def lomaxsf(x, c, scale=1.0):
    return paretosf(x, c, -1, scale)


def lomaxisf(x, c, scale=1.0):
    return paretoisf(x, c, -1, scale)


def lomaxstats(c, scale=1.0, full=0):
    return paretostats(c, -1, scale, full)



# MAXWELL
#  a special case of chi with df = 3, loc=0.0, and given scale = 1.0/sqrt(a)
#    where a is the parameter used in mathworld description

def maxwellpdf(x, scale=1.0):
    x, scale = arr(x), arr(scale)
    x = arr(x*1.0/scale)
    Px = sqrt(2.0/pi)*x*x*exp(-x*x/2.0)
    return select([scale<=0,x>0],[scipy.nan, Px/scale])

def maxwellcdf(x, scale=1.0):
    x, scale = arr(x), arr(scale)
    x = arr(x*1.0/scale)
    sv = errp(0)
    Cx = special.gammainc(1.5,x*x/2.0)
    sv = errp(sv)
    return select([scale<=0,x>0],[scipy.nan, Cx])

def maxwellppf(al, scale=1.0):
    al, scale = arr(al), arr(scale)
    sv = errp(0)
    vals = sqrt(2*special.gammaincinv(1.5,al))
    sv = errp(sv)
    cond = (al>=0)&(al<=1)&(scale>0)
    return where(cond, vals, scipy.nan)

def maxwellsf(x, scale=1.0):
    return 1.0-maxwellcdf(x, scale)

def maxwellisf(al, scale=1.0):
    return maxwellppf(1.0-al,scale)

def maxwellstats(scale=1.0, full=0):
    cond = arr(scale) > 0
    mu = where(cond, scale*2*sqrt(2.0/pi), scipy.nan)
    var = where(cond, scale**2*(3-8.0/pi), scipy.nan)
    if not full:
        return mu, var
    g1 = sqrt(2)*(32-10*pi)/(3*pi-8.0)**1.5
    g2 = (-12*pi*pi+160*pi-384)/(3*pi-8.0)**2
    return mu, var, _wc(cond, g1), _wc(cond, g2)

# Mielke's Beta-Kappa

def mielkepdf():
    pass

# Nakagami (cf Chi)

def nakagamipdf(x, df, loc=0.0, scale=1.0):
    x, loc, scale = map(arr, (x, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Px = 2*(df**df)*x**(2*df-1.0)/special.gamma(df)*exp(-df*x*x)
    sv = errp(sv)
    return select([(df<=0)|(scale<=0),x>0],[scipy.nan,Px/scale])

def nakagamicdf(x, df, loc=0.0, scale=1.0): 
    x, loc, scale = map(arr, (x, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.gammainc(df, df*x*x)
    sv = errp(sv)
    return select([(df<=0)|(scale<=0),x>0],[scipy.nan, Cx])

def nakagamippf(q, df, loc=0.0, scale=1.0):
    sv = errp(0)
    vals = sqrt(special.gammaincinv(df,q)/arr(df))
    sv = errp(sv)
    return vals*scale + loc

def nakagamisf(x, df, loc=0.0, scale=1.0):
    return 1.0-nakagamicdf(x, df, loc, scale)

def nakagamiisf(q, df, loc=0.0, scale=1.0):
    return nakagamippf(1-arr(q), df, loc, scale)

def nakagamistats(df, loc=0.0, scale=1.0, full=0):
    df, scale = arr(df), arr(scale)
    cond = (df > 0) & (scale > 0)
    mu = special.gamma(df+0.5)/sqrt(df)/special.gamma(df)
    mu2 = 1-mu*mu
    mn = _wc(cond, mu*scale + loc)
    var = _wc(cond, mu2*scale*scale)
    if not full:
        return mn, var

    g1 = mu*(1-4*df*mu2)/(2*df*mu2**1.5)
    g2 = (-6*mu**4 * df + (8*df-2)*mu*mu - 2*df + 1)/(df*mu2*mu2)
    return mn, var, _wc(cond, g1), _wc(cond, g2)
    

# Non-central chi-squared
# nc is lambda of definition, df is nu
def ncx2pdf(x, df, nc, loc=0.0, scale=1.0):
    x, df, nc, loc, scale = map(arr, (x, df, nc, loc, scale))
    x = arr((x-loc*1.0)/scale)
    df = where(df==1.0,0.9999999999999,df)
    sv = errp(0)
    a = arr(df/2.0)
    Px = exp(-nc/2.0)*special.hyp0f1(a,nc*x/4.0)
    Px *= exp(-x/2.0)*x**(a-1) / arr(2**a * special.gamma(a))
    sv = errp(sv)
    return select([(nc<0)|(df<=0)|(scale<=0),x>0],[scipy.nan,Px/scale])

def ncx2cdf(x, df, nc, loc=0.0, scale=1.0):
    x, df, nc, loc, scale = map(arr, (x, df, nc, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.chndtr(x,df,nc)
    sv = errp(sv)
    return select([(nc<0)|(df<=0)|(scale<=0),x>0],[scipy.nan,Cx])

def ncx2sf(x,df,nc, loc=0.0, scale=1.0):
    return 1-ncx2cdf(x,df,nc, loc, scale)

def ncx2ppf(q,df,nc, loc=0.0, scale=1.0):
    q, df, nc, loc, scale = map(arr, (q, df, nc, loc, scale))
    cond = (0<=q) & (q<=1) & (scale > 0) & (df > 0) & (nc > 0)
    sv = errp(0)
    vals = special.chndtrix(q,df,nc)
    sv = errp(sv)
    return _wc(cond, vals * scale + loc)
    
def ncx2isf(p,df,nc, loc=0.0, scale=1.0):
    return ncx2ppf(1-p,df,nc,loc,scale)

def ncx2stats(df,nc,loc=0.0, scale=1.0, full=0):
    df, nc, loc, scale = map(arr, (df, nc, loc, scale))
    cond = (df > 0) & (scale > 0) & (nc > 0)
    mn = _wc(cond, (nc+df)*scale + loc)
    var = _wc(cond, 2*(2*nc+df)*scale*scale)
    if not full:
        return mn, var
    g1 = 2*sqrt(2)*(3*nc+df)/(2*nc+df)**1.5
    g2 = 12.0*(4*nc+df)/(2*nc+df)**2.0
    return mn, var, _wc(cond, g1), _wc(cond, g2)

# Non-central F

def ncfpdf(x,dfn,dfd,nc, loc=0.0, scale=1.0):
    x,n1,n2,nc,loc,scale = map(arr, (x, dfn, dfd, nc, loc, scale))
    x = arr((x-loc*1.0)/scale)
    y = where(x<0,1.0,x)
    sv = errp(0)
    Px = exp(-nc/2+nc*n1*y/(2*(n2+n1*y)))
    Px *= n1**(n1/2) * n2**(n2/2) * y**(n1/2-1)
    Px *= (n2+n1*y)**(-(n1+n2)/2)
    Px *= special.gamma(n1/2)*special.gamma(1+n2/2)
    Px *= special.assoc_laguerre(-nc*n1*y/(2.0*(n2+n1*y)),n2/2,n1/2-1)
    Px /= special.beta(n1/2,n2/2)*special.gamma((n1+n2)/2.0)
    sv = errp(sv)
    cond = (n1<=0) | (n2<=0) | (nc<0) | (scale<=0)
    return select([cond, x>=0],[scipy.nan, Px])

def ncfcdf(x,dfn,dfd,nc, loc=0.0, scale=1.0):
    x,dfn,dfd,nc,loc,scale = map(arr, (x, dfn, dfd, nc, loc, scale))
    x = arr((x-loc*1.0)/scale)
    badcond = (dfn<=0) | (dfd<=0) | (nc<0) | (scale<=0)
    sv = errp(0)
    vals = special.ncfdtr(dfn,dfd,nc,x)
    sv = errp(sv)
    return select([badcond, x>=0],[scipy.nan, vals])
    
def ncfsf(x,dfn,dfd,nc, loc=0.0, scale=1.0):
    return 1-ncfcdf(x,dfn,dfd,nc, loc=loc, scale=scale)

def ncfppf(q, dfn, dfd, nc, loc=0.0, scale=1.0):
    q,dfn,dfd,nc,loc,scale = map(arr, (q, dfn, dfd, nc, loc, scale))    
    cond = ((0<=q) & (q<=1) & (dfn>0) & (dfd>0) & (nc>0) & (scale > 0))
    vals = special.ncfdtri(dfn, dfd, nc, q)
    return _wc(cond, vals*scale + loc)

def ncfisf(p,dfn,dfd,nc, loc=0.0, scale=1.0):
    return ncfppf(1-p,dfn,dfd,nc, loc, scale)

def ncfstats(dfn,dfd,nc,loc=0.0,scale=1.0,full=0):
    dfn = arr(dfn)*1.0
    dfd = arr(dfd)*1.0
    nc = arr(nc)*1.0
    cond = (dfd > 0) & (nc > 0) & (scale > 0)
    mn = where((1-cond) & (dfd<=2),scipy.nan,dfd/(dfd-2)*(1+nc/dfn)*scale + loc)
    var1 = 2*(dfd/dfn)**2 * ((dfn+nc/2)**2+(dfn+nc)*(dfd-2))
    var1 /= (dfd-2)**2 * (dfd-4)
    var = where((1-cond) & (dfd<=4),scipy.nan,var1*scale*scale)       
    if not full:
        return mn, var
    _vecfunc = special.general_function(generic_moment)
    mu3 = _vecfunc(3,  ncfpdf, 0, scipy.inf, dfn, dfd, nc, loc, scale)
    mu4 = _vecfunc(4,  ncfpdf, 0, scipy.inf, dfn, dfd, nc, loc, scale)
    g1 = mu3 / var1**1.5
    g2 = mu4 / var1**2 - 3.0
    return mn, var, _wc(cond & (dfd>6), g1), _wc(cond & (dfd>8), g2)


## Student t distribution

def tpdf(x, df, loc=0.0, scale=1.0):
    x, df, loc, scale = map(arr,(x, df, loc, scale))
    x = arr((x-loc*1.0)/scale)
    r = arr(df*1.0)
    sv = errp(0)
    Px = exp(special.gammaln((r+1)/2)-special.gammaln(r/2))
    Px /= sqrt(r*pi)*(1+(x**2)/r)**((r+1)/2)
    sv = errp(sv)
    return select([scale <=0],[scipy.nan], Px/scale)

def tcdf(x, df, loc=0.0, scale=1.0):
    x, df, loc, scale = map(arr,(x, df, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.stdtr(df, x)
    sv = errp(sv)
    return select([scale <=0],[scipy.nan], Cx)

def tsf(x, df, loc=0.0, scale=1.0):
    return 1-tcdf(x,df, loc, scale)

def tppf(q, df, loc=0.0, scale=1.0):
    q, df, loc, scale = map(arr, (q, df, loc, scale))    
    cond = (q<=1) & (q>=0) & (scale > 0)
    return _wc(cond, special.stdtrit(df, q)*scale + loc)

def tisf(q, df, loc=0.0, scale=1.0):
    return tppf(1-arr(q),df,loc,scale)

def tstats(df, loc=0.0, scale=1.0, full=0):
    df, loc, scale = map(arr, (df, loc, scale))
    r = arr(df*1.0)
    mn = loc
    var = _wc(r>2,r/(r-2)*scale*scale)
    if not full:
        return mn, var
    return mn, var, 0.0, _wc(r>4,6.0/(r-4.0))

## Non-central T distribution

def nctpdf(x, df, nc, loc=0.0, scale=1.0):
    x, df, nc, loc, scale = map(arr, (x, df, nc, loc, scale))
    x = arr((x-loc*1.0)/scale)
    n = df*1.0
    nc = nc*1.0
    x2 = x*x
    ncx2 = nc*nc*x2
    fac1 = n + x2
    sv = errp(0)
    Px = n**(n/2) * special.gamma(n+1)
    Px /= arr(2.0**n*exp(nc*nc/2)*fac1**(n/2)*special.gamma(n/2))
    valF = ncx2 / (2*fac1)
    trm1 = sqrt(2)*nc*x*special.hyp1f1(n/2+1,1.5,valF)
    trm1 /= arr(fac1*special.gamma((n+1)/2))
    trm2 = special.hyp1f1((n+1)/2,0.5,valF)
    trm2 /= arr(sqrt(fac1)*special.gamma(n/2+1))
    Px *= trm1+trm2
    sv = errp(sv)
    return select([(df<=0) | (nc<=0) | (scale<=0)],[scipy.nan],Px/scale)
    
def nctcdf(x,df,nc, loc=0.0, scale=1.0):
    x, df, nc, loc, scale = map(arr, (x, df, nc, loc, scale))
    x = arr((x-loc*1.0)/scale)
    sv = errp(0)
    Cx = special.nctdtr(df, nc, x)
    sv = errp(sv)
    return select([(df<=0) | (nc<=0) | (scale<=0)],[scipy.nan],Cx)

def nctsf(x,df,nc, loc=0.0, scale=1.0):
    return 1-nctcdf(x,df,nc, loc, scale)

def nctppf(q,df,nc, loc=0.0, scale=1.0):
    q, df, nc, loc, scale = map(arr, (q, df, nc, loc, scale))
    cond = (0<=q) & (q<=1) & (df > 0) & (nc > 0)
    sv = errp(0)
    vals = special.nctdtrit(df, nc, q)
    sv = errp(sv)
    return _wc(cond, vals*scale + loc)

def nctisf(p,df,nc, loc=0.0, scale=1.0):
    return nctppf(1-p,df,nc, loc, scale)

def nctstats(df,nc,loc=0.0, scale=1.0, full=0):
    df, nc, loc, scale = map(arr, (df, nc, loc, scale))
    cond = (df > 0) & (nc > 0) & (scale > 0)
    nc = nc*1.0
    df = df*1.0
    gam = special.gamma
    val1 = gam((df-1.0)/2.0)
    val2 = gam(df/2.0)
    mn = _wc(cond, nc*sqrt(df/2.0)*val1/val2*scale + loc)
    var = (nc*nc+1.0)*df/(df-2.0)
    var -= nc*nc*df* val1**2 / 2.0 / val2**2
    var = _wc(cond, var*scale*scale)
    if not full:
        return mn, var
    g1n = 2*nc*sqrt(df)*val1*((nc*nc*(2*df-7)-3)*val2**2 \
                             -nc*nc*(df-2)*(df-3)*val1**2)
    g1d = (df-3)*sqrt(2*df*(nc*nc+1)/(df-2) - nc*nc*df*(val1/val2)**2) \
          * val2 * (nc*nc*(df-2)*val1**2 - 2*(nc*nc+1)*val2**2)
    g2n = 2*(-3*nc**4*(df-2)**2 *(df-3) *(df-4)*val1**4 + \
             2**(6-2*df) * nc*nc*(df-2)*(df-4)*(nc*nc*(2*df-7)-3)*pi* \
             gam(df+1)**2 - 4*(nc**4*(df-5)-6*nc*nc-3)*(df-3)*val2**4)
    g2d = (df-3)*(df-4)*(nc*nc*(df-2)*val1**2 - 2*(nc*nc+1)*val2)**2
    return mn, var, _wc(cond, g1n/g1d), _wc(cond, g2n/g2d)



# Pareto

def paretopdf(x, b, loc=0.0, scale=1.0):
    x, b, loc, scale = map(arr, (x,b,loc,scale))
    x = arr((x-loc*1.0)/scale)
    Px = b * x**(-b-1)
    return select([(b<=0)&(scale<=0.0)],[scipy.nan],Px/scale)

def paretocdf(x, b, loc=0.0, scale=1.0):
    x, b, loc, scale = map(arr, (x,b,loc,scale))
    x = arr((x-loc*1.0)/scale)
    Cx = 1 -  x**(-b)
    return select([(b<=0)&(scale<=0.0)],[scipy.nan],Cx)

def paretosf(x, b, loc=0.0, scale=1.0):
    return 1-paretocdf(x, b, loc, scale)

def paretoppf(q, b, loc=0.0, scale=1.0):
    q, b, loc, scale = map(arr, (q,b,loc,scale))
    cond = (q>=0)&(q<=1)&(scale>0)&(b>0)
    vals = (1-q, -1.0/b)
    return _wc(cond, vals*scale + loc)

def paretoisf(q, b, loc=0.0, scale=1.0):
    q, b, loc, scale = map(arr, (q,b,loc,scale))
    cond = (q>=0)&(q<=1)&(scale>0)&(b>0)
    vals = (q, -1.0/b)
    return _wc(cond, vals*scale + loc)

def paretostats(b, loc=0.0, scale=1.0, full=0):
    b, loc, scale = map(arr, (b,loc,scale))
    cond (scale > 0) & (b>0)
    mu = b / (b-1.0)
    mn = select([cond, b<=1],[scipy.nan, scipy.inf], mu*scale + loc)
    mu2 = b / (b-1.0)**2 / (b-2.0)
    var = select([cond, b<=2],[scipy.nan, scipy.inf], mu2*scale*scale)
    if not full:
        return mn, var
    g1 = sqrt((b-2.0)/b)*2.0*(b+1.0)/(b-3.0)
    g2 = 6*scipy.polyval([1,1,-6,-2],b)/(b*(b-3.0)*(b-4))
    g1 = select([cond, b<=3],[scipy.nan, scipy.inf], g1)
    g2 = select([cond, b<=4],[scipy.nan, scipy.inf], g2)
    return mn, var, g1, g2

## Power distribution
##   Special case of beta dist. with d =1.0

def powerpdf(x, a, loc=0.0, scale=1.0):
    x, a, loc, scale = map(arr, (x, a, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Px = a*x**(a-1.0)
    return select([(a<=0)|(scale<=0),(x>=0)&(x<=1)],[scipy.nan, Px/scale])

def powercdf(x, a, loc=0.0, scale=1.0):
    x, a, loc, scale = map(arr, (x, a, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Cx = x**a
    return select([(a<=0)|(scale<=0),(x>1),(x>=0)],[scipy.nan, 1, Cx])

def powersf(x, a, loc=0.0, scale=1.0):
    return 1.0-powercdf(x, a, loc, scale)

def powerppf(q, a, loc=0.0, scale=1.0):
    q, a, loc, scale = map(arr, (q, a, loc, scale))
    vals = pow(q, 1.0/a)
    cond = (q>=0)&(q<=1)&(a>0)&(scale>0)
    return _wc(cond, vals)

def powerisf(q, a, loc=0.0, scale=1.0):
    return powerppf(1-arr(q), a, loc, scale)

def powerstats(a, loc=0.0, scale=1.0, full=0):
    a, scale = arr(a), arr(scale)
    cond = (a>0) & (scale>0)
    mn = _wc(cond, a/(a+1.0)*scale + loc)
    var = _wc(cond, a*(a+1.0)/(a+1.0)**2 * scale**2)
    if not full:
        return mn, var
    g1 = 2*(1.0-a)*sqrt((a+2.0)/a/(a+3.0))
    g2 = 6*polyval([1,-1,-6,2.],a)/(a*(a+3)*(a+4.0))
    return mn, var, _wc(cond, g1), _wc(cond, g2)

# Power log normal

def powerlognormpdf():
    pass

# Power Normal

def powernormpdf():
    pass

# Reciprocal Distribution

def reciprocalpdf(x, a, b, loc=0.0, scale=1.0):
    x, a, b, loc, scale = map(arr, (x, a, b, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Px = 1.0/(x*log(b/a))
    return select([(a<=0)|(b<=a)|(scale<=0),(x>=a)&(x<=b)],
                  [scipy.nan, Px/scale])

def reciprocalcdf(x, a, b, loc=0.0, scale=1.0):
    x, a, b, loc, scale = map(arr, (x, a, b, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Cx = log(a/x)/log(a/b)
    return select([(a<=0)|(b<=a)|(scale<=0),x>b,x>=a],
                  [scipy.nan, 1.0, Cx])

def reciprocalppf(q, a, b, loc=0.0, scale=1.0):
    q, a, b, loc, scale = map(arr, (q, a, b, loc, scale))
    vals = a*(arr(b*1.0/a)**q)
    cond = (a>0)&(b>a)&(scale>0)&(q>=0)&(q<=1.0)
    return where(cond, vals*scale + loc, scipy.nan)

def reciprocalsf(x, a, b, loc=0.0, scale=1.0):
    return 1.0-reciprocalcdf(x, a, b, loc, scale)

def reciprocalisf(q, a, b, loc=0.0, scale=1.0):
    return reciprocalppf(1-arr(q), a, b, loc, scale)

def reciprocalstats(a, b, loc=0.0, scale=1.0, full=0):
    a, b, loc, scale = map(arr, (a, b, loc, scale))
    cond = (a>0) & (b>a) & (scale > 0)
    b = arr(b*1.0)
    d = log(a/b)
    mu = (a-b)/d
    mu2 = (a+b)*mu/2 - mu*mu
    mn = _wc(cond, mu*scale + loc)
    var = _wc(cond, mu2*scale*scale)
    if not full:
        return mn, var
    fac = 2*d*d*mu2/(a-b)
    g1 = sqrt(2)*(12*d*(a-b)**2 + d*d*(a*a*(2*d-9) + 2*a*b*d + b*b*(2*d+9)))
    g1 /= 3*d*sqrt((a-b)*fac)*fac
    g2 = -36*(a-b)**3 + 36*d*(a-b)**2*(a+b) - 16*d*d*(a**3-b**3) + \
         3*d**3*(a*a+b*b)*(a+b)
    g2 /= 3*(a-b)*fac*fac
    g2 -= 3
    return mn, var, _wc(cond, g1), _wc(cond, g2)

# Reciprocal Inverse Gaussian

def recipinvgauss():
    pass

# Rayleigh distribution (this is chi with df=2 and loc=0.0)
# scale is the mode.
def rayleighpdf(r, scale=1.0):
    r, scale = map(arr, (r,scale))
    r = arr(r *1.0 / scale)
    Px = r*exp(-r*r/2.0)
    return select([(scale<=0),r>=0],[scipy.nan,Px / scale])

def rayleighcdf(r, scale=1.0):
    r, scale = map(arr, (r,scale))
    r = arr(r *1.0 / scale)
    return select([scale<=0,r>=0],[scipy.nan, 1-exp(-r*r/2.0)])

def rayleighsf(r, scale=1.0):
    r, scale = map(arr, (r,scale))
    r = arr(r *1.0 / scale)
    return select([scale<=0,r>=0],[scipy.nan,exp(-r*r/2.0)],1)

def rayleighppf(q, scale=1.0):
    q, scale = map(arr, (q, scale))
    cond = (0<=q) & (q<=1) & (scale>0)
    vals = scale*sqrt(2*log(1.0/arr(1.0-q)))
    return select([1-cond],[scipy.nan],vals)

def rayleighisf(q, scale=1.0):
    q, scale = map(arr, (q, scale))
    cond = (0<=q) & (q<=1) & (scale>0)
    vals = scale*sqrt(2*log(1.0/arr(1.0-q)))
    return select([1-cond],[scipy.nan],vals)


def rayleighstats(scale=1.0, full=0):
    s = scale
    mn = s*sqrt(pi/2)
    var = (4-pi)/2*s*s
    if not full:
        return mn, var
    g1 = 2*(pi-3)*sqrt(pi)/(4-pi)**1.5
    g2 = (24*pi-6*pi*pi-16)/(pi-4)**2
    return mn, var, g1, g2

# Semicircular

def semicircularpdf(x, loc=0.0, scale=1.0):
    x, loc, shape = map(arr, (x, loc, scale))
    x = arr((x-1.0*loc)/scale)
    Px = 2.0/pi*sqrt(1-x*x)
    return select([scale <=0,(x>=-1)&(x<=1)],[scipy.nan,Px/scale])

def semicircularcdf(x, loc=0.0, scale=1.0):
    x, loc, shape = map(arr, (x, loc, scale))
    x = arr((x-1.0*loc)/scale)
    Cx = 0.5 + 1.0/pi*(x*sqrt(1-x*x) + arcsin(x))
    return select([scale <=0,x>1,x>0],[scipy.nan,1,Cx])

_semicircppf = general_cont_ppf('semicircular')
def semicircularppf(q, loc=0.0, scale=1.0):
    return _semicircppf(q, loc, scale)

def semicircularsf(x, loc=0.0, scale=1.0):
    return 1.0-semicircularcdf(x, loc, scale)

def semicircularisf(q, loc=0.0, scale=1.0):
    return semicircularppf(1-arr(q), loc, scale)

def semicircularstats(loc=0.0, scale=1.0, full=0):
    loc, scale = arr(loc), arr(scale)
    cond = (scale > 0)
    mu2 = 0.25
    mn = _wc(cond, loc)
    var = _wc(cond, mu2*scale*scale)
    if not full:
        return mn, var
    return mn, var, _wc(cond, 0.0), _wc(cond, -1.0)

# Triangular
# up-sloping line from loc to (loc + c) and then downsloping line from
#    loc + c to loc + shape

# _trstr = "Left must be <= mode which must be <= right with left < right"

def triangpdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    P1x = 2.0*x / c
    P2x = 2.0*(1-x)/arr(1-c)
    return select([(c<0)|(c>1)|(scale <=0), x>1, x>=c, x>=0],
                  [scipy.nan, 0, P2x/scale,P1x/scale])
    
def triangcdf(x, c, loc=0.0, scale=1.0):
    x, c, loc, scale = map(arr, (x, c, loc, scale))
    x = arr((x-loc*1.0)/scale)
    C1x = x*x / c
    C2x = (x*x-2*x+c) / arr(c-1)
    return select([(c<0)|(c>1)|(scale <=0), x>1, x>=c, x>=0],
                  [scipy.nan, 1, C2x, C1x])

def triangsf(x, c, loc=0.0, scale=1.0):
    return 1.0-triangcdf(x, c, loc, scale)

def triangppf(q, c, loc=0.0, scale=1.0):
    q, c, loc, scale = map(arr, (q, c, loc, scale))
    cond = (0<=q) & (q<=1) & (scale > 0) & (c>=0) & (c<=1)
    x1 = sqrt(c*q)
    x2 = 1-sqrt((1-c)*(1-q))
    select([1-cond, q<c],[scipy.nan, x1*scale + loc], x2*scale + loc)

def triangisf(p, c, loc=0.0, scale=1.0):
    return triangppf(1.0-p, c, loc, scale)

def triangstats(c, loc=0.0, scale=1.0, full=0):
    cond = (c >=0) & (c <= 1) & (scale > 0)
    mn = _wc(cond, (c+1.0) / 3 * scale + loc)
    fac = (1.0-c+c*c)
    var = _wc(cond, fac / 18.0*scale*scale)
    if not full:
        return mn, var
    g1 = sqrt(2)*(2*c-1)*(c+1)*(c-2) / 5.0 / fac**1.5
    g2 = -3.0/5
    return mn, var, _wc(cond, g1), _wc(cond, g2)

# Truncated Exponential

def truncexponpdf():
    pass

def truncnormpdf():
    pass

# Tukey-Lambda
# A flexible distribution ranging from Cauchy (lam=-1)
#   to logistic (lam=0.0)
#   to approx Normal (lam=0.14)
#   to u-shape (lam = 0.5)
#   to Uniform from -1 to 1 (lam = 1)

def tukeylambdapdf(x, lam, loc=0.0, scale=1.0):
    x, lam, loc, scale = map(arr, (x, lam, loc, scale))
    cond = (scale > 0)&(x==x)
    x = where(cond, arr((x-loc*1.0)/scale), 0.0)
    Fx = arr(special.tklmbda(x,lam))
    Px = Fx**(lam-1.0) + (arr(1-Fx))**(lam-1.0)
    Px = 1.0/arr(Px)
    maglim = arr(1.0/lam)
    return select([scale<=0, (lam>0)&((x<-maglim)|(x>maglim))],
                  [scipy.nan, 0.0], Px/scale)

def tukeylambdacdf(x, lam, loc=0.0, scale=1.0):
    x, lam, loc, scale = map(arr, (x, lam, loc, scale))
    cond = (scale > 0)&(x==x)
    x = where(cond, arr((x-loc*1.0)/scale), 0.0)
    return _wc(cond, special.tklmbda(x, lam))
                
def tukeylambdappf(q, lam, loc=0.0, scale=1.0):
    q, lam, loc, scale = map(arr, (q, lam, loc, scale))
    q = arr(q*1.0)
    vals1 = (q**lam - (1-q)**lam)/lam
    vals2 = log(q/(1-q))
    return select([scale<=0, lam==0], [scipy.nan, vals2*scale+loc],
                  vals1*scale + loc)

def tukeylambdasf(x, lam, loc=0.0, scale=1.0):
    return 1.0-tukeylambdacdf(x, lam, loc, scale)

def tukeylambdaisf(q, lam, loc=0.0, scale=1.0):
    return tukeylambdappf(1.0-q, lam, loc, scale)

def tukeylambdastats(lam, loc=0.0, scale=1.0, full=0):
    lam, loc, scale = arr(lam), arr(loc), arr(scale)
    mn = _wc(scale>0,loc)
    _vecfunc = special.general_function(moment_ppf)
    mu2 = _vecfunc(2, tukeylambdappf, lam)
    var = _wc(scale>0, mu2*scale*scale)  
    if not full:
        return mn, var
    g1 = 0
    g2 = _vecfunc(4, tukeylambdappf, lam) / mu2**2 - 3.0
    return mn, var, _wc(scale>0, g1), _wc(scale>0, g2)


# Uniform

# loc to loc + shape
def uniformpdf(x, loc=0.0, shape=1.0):
    x, loc, shape = map(arr, (x, loc, shape))
    x = arr((x-loc*1.0)/shape)
    return select([shape <= 0],[scipy.nan],1.0*shape*(x>0)*(x<1))

def uniformcdf(x, loc=0.0, shape=1.0):
    x, loc, shape = map(arr, (x, loc, shape))
    x = arr((x-loc*1.0)/shape)
    return select([shape<=0, x>=1, x>0],[scipy.nan, 1.0, x])
    
def uniformsf(x, loc=0.0, shape=1.0):
    return 1-uniformcdf(x,loc,shape)

def uniformppf(q, loc=0.0, shape=1.0):
    q, loc, shape = map(arr, (q, loc, shape))
    cond (shape > 0) & (q<=1) & (q>=0)
    return _wc(cond, q*shape + loc)

def uniformisf(q, loc=0.0, shape=1.0):
    return uniformppf(1-arr(q),loc,shape)

def uniformstats(loc=0.0, shape=1.0, full=0):
    loc, shape = arr(loc), arr(shape)
    cond = (shape > 0)
    mn = _wc(cond, 0.5*shape+loc)
    var = _wc(cond, shape**2.0 / 12.0)
    if not full:
        return mn, var
    return mn, var, _wc(cond, 0.0), _wc(cond, -6.0/5.0)

# Von-Mises

# if x is not in range or loc is not in range it assumes they are angles
#   and converts them to [-pi, pi] equivalents. 
def von_misespdf(x,b,loc=0.0):
    x, b, loc = map(arr, (x, b, loc))
    loc = arr(angle(exp(1j*loc)))
    x = arr(angle(exp(1j*(x-loc))))
    Px1 = exp(b*cos(x)) / (2.0*pi*special.i0(b))
    Px2 = normpdf(x, 0.0, sqrt(1.0/b))
    return select([b < 0, b < 100],[scipy.nan, Px1], Px2)

def von_misescdf(x, b, loc=0.0):
    x, b, loc = map(arr, (x, b, loc))
    if len(x.shape) > 1 or len(b.shape) > 1 or len(loc.shape) > 1:
        raise ValueError, "Only works for 1-d arrays."
    loc = arr(angle(exp(1j*loc)))
    x = atleast_1d(angle(exp(1j*(x-loc))))*(b==b)
    from scipy_base.limits import double_epsilon as eps
    eps2 = sqrt(eps)

    c_xsimple = atleast_1d((b==0)&(x==x))
    c_xiter = atleast_1d((b<100)&(b > 0)&(x==x))
    c_xnormal = atleast_1d((b>=100)&(x==x))
    c_bad = atleast_1d((b<=0) | (x != x))

    indxiter = nonzero(c_xiter)
    xiter = take(x, indxiter)

    vals = ones(len(c_xsimple),Numeric.Float)
    putmask(vals, c_bad, scipy.nan)
    putmask(vals, c_xsimple, x / 2.0/pi)
    st = sqrt(b-0.5)
    st = where(isnan(st),0.0,st)
    putmask(vals, c_xnormal, normcdf(x, scale=st))
        
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

def _vmqfunc(x,q,b,loc):
    return von_misescdf(x,b,loc)-q

def _vmppf(q,b,loc,x0):
    import scipy.optimize as optimize
    return optimize.fsolve(_vmqfunc,x0,args=(q,b, loc))
_vec_vmppf = special.general_function(_vmppf,'d')

def von_misesppf(q, b, loc=0.0, x0=None):
    if x0 is None:
        x0 = loc
    assert all((0<=q) & (q<=1)), _quantstr
    return _vec_vmppf(q, b, loc,x0)

def von_misesisf(p, b, loc=0.0, x0=None):
    return von_misesppf(1-p, b, loc, x0)

def von_misesstats(b, loc=0.0, full=0):
    b, loc = arr(b), arr(loc)
    mn = _wc(b>=0,loc)
    _vecfunc = special.general_function(generic_moment)
    mu2 = _vecfunc(2, von_misespdf, -pi, pi, b, loc)
    var = _wc(b>=0, mu2)    
    if not full:
        return mn, var
    g1 = 0
    g2 = _vecfunc(4, von_misespdf, -pi, pi, b, loc) / mu2**2 - 3.0
    return mn, var, _wc(b>=0, g1), _wc(b>=0, g2)

## Wald distribution (Inverse Normal with shape parameter mu=1.0)

def waldpdf(x, loc=0.0, scale=1.0):
    return invnormpdf(x, 1.0, loc, scale)

def waldcdf(x, loc=0.0, scale=1.0):
    return invnormcdf(x, 1.0, loc, scale)

def waldsf(x, loc=0.0, scale=1.0):
    return invnormsf(x, 1.0, loc, scale)

def waldppf(q, loc=0.0, scale=1.0, x0=None):
    return invnormppf(q, 1.0, loc, scale)

def waldisf(q, loc=0.0, scale=1.0):
    return invnormisf(q, 1.0, loc, scale)

def waldstats(loc=0.0, scale=1.0, full=0):
    return invnormstats(1.0, loc, scale, full=full)

## Weibull

def weibullpdf(x, shape, left=0, loc=0.0, scale=1.0):
    c, b, A, x, left = map(arr,(shape, scale, loc, x, left))
    x = arr((x-A*1.0)/b)
    x = arr(where(left,-x,x))
    Px = c * x**(c-1.0) * exp(-x**c)
    return select([(c<=0)|(b<=0),((x>0)&(left==0))|((x<0)&(left!=0))],[scipy.nan,Px/b])

def weibullcdf(x, shape, left=0, loc=0.0, scale=1.0):
    c, b, A, x, left = map(arr,(shape, scale, loc, x, left))
    x = arr((x-A*1.0)/b)
    sv =errp(0)
    Cx = -special.expm1(-x**c)
    Cx2 = exp(-(-x)**c)
    sv = errp(sv)
    return select([(c<=0)|(b<=0),(x>0)&(left==0), (x<0)&(left!=0), left!=0],[scipy.nan,Cx, Cx2,1])

def weibullsf(x, shape, left=0, loc=0.0, scale=1.0):
    c, b, A, x, left = map(arr,(shape, scale, loc, x, left))
    x = arr((x-A*1.0)/b)
    sv = errp(0)
    Cx = exp(-x**c)
    Cx2 = -special.expm1(-(-x)**c)
    sv = errp(sv)
    return select([(c<=0)|(b<=0),(x>0)&(left==0),(x<0)&(left!=0), left==0],[scipy.nan,Cx,Cx2,1])

def weibullppf(q, shape, left=0, loc=0.0, scale=1.0):
    a, b, loc, q, left = map(arr,(shape, scale, loc, q, left))
    cond1 = (a>0) & (b>0) & (0<=q) & (q<=1)
    q = arr(where(left, 1-q, q))
    vals = pow(arr(log(1.0/arr(1-q))),1.0/a)
    return select([1-cond1,left==0], [scipy.nan, b*vals+loc], -b*vals+loc)

def weibullisf(q, shape, left=0, loc=0.0, scale=1.0):
    a, b, loc, q, left = map(arr,(shape, scale, loc, q, left))
    cond1 = (a>0) & (b>0) & (0<=q) & (q<=1)
    q = arr(where(left, 1-q, q))
    vals = pow(arr(log(1.0/q)),1.0/a)
    return select([1-cond1,left==0], [scipy.nan, b*vals+loc], -b*vals+loc)

def weibullstats(shape, left=0, loc=0.0, scale=1.0, full=0):
    a, loc, b, left = map(arr,(shape, loc, scale, left))
    cond = (a>0) & (b>0)
    gam = special.gamma
    ia = 1.0/a
    mn = select([1-cond, left==0], [scipy.nan, b*gam(1+ia)+loc], -b*gam(1+ia)+loc)
    var = _wc(cond, b*b*(gam(1+2*ia)-gam(1+ia)**2))
    if not full:
        return mn, var
    den = (gam(1+2*ia)-gam(1+ia)**2)
    g1 = 2*gam(1+ia)**3 - 3*gam(1+ia)*gam(1+2*ia) + gam(1+3*ia)
    g1 /= den**1.5
    g2 = 12*gam(1+ia)**2 * gam(1+2*ia) - 6*gam(1+ia)**4
    g2 += gam(1+4*ia) - 4*gam(1+ia)*gam(1+3*ia) - 3*gam(1+2*ia)**2
    g2 /= den**2.0
    return mn, var, _wc(cond, g1), _wc(cond, g2)


# Wrapped Cauchy

def wrapcauchypdf():
    pass

        
### DISCRETE DISTRIBUTIONS
###


# Binomial 
    
def binompdf(k, n, pr=0.5):
    k = arr(k)
    assert (0<pr<1)
    cond = arr((k > 0) & (k == floor(k)))
    return _wc(cond, scipy.comb(n,k)* pr**k * (1-pr)**(n-k))

def binomcdf(k, n, pr=0.5):
    return special.bdtr(k,n,pr)

def binomsf(k, n, pr=0.5):
    return special.bdtrc(k,n,pr)

def binomppf(q, n, pr=0.5):
    return special.bdtrik(q,n,pr)

def binomisf(p, n, pr=0.5):
    return special.bdtrik(1-p,n,pr)

def binomstats(n, pr=0.5, full=0):
    q = 1.0-pr
    mu = n * pr
    var = n * pr * q
    if not full:
        return mu, var
    g1 = (q-pr) / sqrt(n*pr*q)
    g2 = (1.0-6*pr*q)/(n*pr*q)
    return mu, var, g1, g2


# Bernoulli distribution

def bernoullipdf(k, pr=0.5):
    return binompdf(k, 1, pr)

def bernoullicdf(k, pr=0.5):
    return binomcdf(k, 1, pr)

def bernoullisf(k, pr=0.5):
    return binomsf(k, 1, pr)

def bernoullippf(k, pr=0.5):
    return binomppf(k, 1, pr)

def bernoulliisf(k, pr=0.5):
    return binomp(k, 1, pr)

def bernoullistats(pr=0.5, full=0):
    return binomstats(1, pr, full)

# Negative binomial

def nbinompdf(k, n, pr=0.5):
    return scipy.comb(n+k-1,k)* pr**n * (1-pr)**k

def nbinomcdf(k, n, pr=0.5):
    return special.nbdtr(k,n,pr)

def nbinomsf(k, n, pr=0.5):
    return special.nbdtrc(k,n,pr)

def nbinomppf(q, n, pr=0.5):
    return special.nbdtrik(q,n,pr)

def nbinomisf(p, n, pr=0.5):
    return special.nbdtrik(1-p,n,pr)

def nbinomstats(n, pr=0.5, full=0):
    Q = 1.0 / pr
    P = Q - 1.0
    mu = n*P
    var = n*P*Q
    if not full:
        return mu, var
    g1 = (Q+P)/sqrt(n*P*Q)
    g2 = (1.0 + 6*P*Q) / (n*P*Q)
    return mu, var, g1, g2 


## Geometric distribution

def geompdf(k, pr=0.5):
    return (1-pr)**k * pr

def geomcdf(k, pr=0.5):
    return 1.0-(1.0-pr)**(k+1)

def geomsf(k, pr=0.5):
    return (1.0-pr)**(k+1)

def geomppf(q, pr=0.5):
    return log(1.0-q)/log(1-pr)-1

def geomisf(p, pr=0.5):
    return log(p)/log(1.0-pr)-1

def geomstats(pr=0.5,full=0):
    mu = 1.0/pr
    qr = 1.0-pr
    var = qr / pr / pr
    if not full:
        return mu, var
    g1 = (2.0-pr) / sqrt(qr)
    g2 = scipy.polyval([1,-6,6],pr)/(1.0-pr)
    return mu, var, g1, g2


## Hypergeometric distribution

def hypergeompdf(k, tot=35, good=25, N=10):
    bad = tot - good
    return scipy.comb(good,k) * scipy.comb(bad,N-k) / scipy.comb(tot,N)

def _hypergeomcdf(k, tot, good, N):
    j = arange(0,k+1)
    return sum(hypergeompdf(j, tot, good, N),axis=-1)
_vhypergeomcdf = special.general_function(_hypergeomcdf,'d')

def hypergeomcdf(k, tot=35, good=25, N=10):
    return _vhypergeomcdf(k, tot, good, N)

def hypergeomsf(k, tot=35, good=25, N=10):
    return 1.0 - hypergeomcdf(k, tot, good, N)

def hypergeomppf(q, tot=35, good=25, N=10):
    pass

def hypergeomisf(p, tot=35, good=25, N=10):
    pass

def hypergeomstats(tot=35, good=25, N=10, full=0):
    n = good*1.0
    m = (tot-good)*1.0
    N = N*1.0
    tot = m+n
    p = n/tot
    mu = N*p
    var = m*n*N*(tot-N)*1.0/(tot*tot*(tot-1))
    if not full:
        return mu, var
    g1 = (m - n)*(tot-2*N) / (tot-2.0)*sqrt((tot-1.0)/(m*n*N*(tot-N)))
    m2, m3, m4, m5 = m**2, m**3, m**4, m**5
    n2, n3, n4, n5 = n**2, n**2, n**4, n**5
    g2 = m3 - m5 + n*(3*m2-6*m3+m4) + 3*m*n2 - 12*m2*n2 + 8*m3*n2 + n3 \
         - 6*m*n3 + 8*m2*n3 + m*n4 - n5 - 6*m3*N + 6*m4*N + 18*m2*n*N \
         - 6*m3*n*N + 18*m*n2*N - 24*m2*n2*N - 6*n3*N
    return mu, var, g1, g2


## Logarithmic (Log-Series), (Series) distribution

def logserpdf(k,pr=0.5):
    k = arr(k)
    assert (0<pr<1)
    cond = arr((k > 0) & (k == floor(k)))
    return where(cond,-pr**k / k / log(1-pr),0)

def logsercdf(k, pr=0.5):
    j = arange(0,k+1)
    pr = arr(pr)
    return sum(logserpdf(j, pr),axis=-1)

def logsersf(k, pr=0.5):
    return 1.0-logsercdf(k, pr=pr)

def logserppf(q, pr=0.5):
    pass

def logserisf(p, pr=0.5):
    pass

def logserstats(pr=0.5, full=0):
    r = log(1-pr)
    mu = pr / (pr - 1.0) / r
    mu2p = -pr / r / (pr-1.0)**2
    var = mu2p - mu*mu
    if not full:
        return mu, var
    mu3p = -pr / r * (1.0+pr) / (1.0-pr)**3
    mu3 = mu3p - 3*mu*mu2p + 2*mu**3
    g1 = mu3 / var**1.5

    mu4p = -pr / r * (1.0/(pr-1)**2 - 6*pr/(pr-1)**3 + 6*pr*pr / (pr-1)**4)
    mu4 = mu4p - 4*mu3p*mu + 6*mu2p*mu*mu - 3*mu**4
    g2 = mu4 / var**2 - 3.0
    return mu, var, g1, g2


## Poisson distribution

def poissonpdf(k,mu):
    k, mu = arr(k)*1.0, arr(mu)*1.0
    sv = errp(0)
    Pk = mu**k * exp(-mu) / arr(special.gamma(k+1))
    sv = errp(sv)
    return select([mu<0,(k>=0) & (floor(k)==k)],[scipy.nan,Pk])

def poissoncdf(k,mu):
    k, mu = arr(k), arr(mu)
    sv = errp(0)
    vals = special.pdtr(k,mu)
    sv = errp(sv)
    return select([mu<0,k>=0],[scipy.nan, vals])

def poissonsf(k,mu):
    k, mu = arr(k), arr(mu)
    sv = errp(0)
    vals = special.pdtrc(k,mu)
    sv = errp(sv)
    return select([mu<0,k>=0],[scipy.nan, vals])

def poissonppf(q,mu):
    q, mu = arr(q), arr(mu)
    sv = errp(0)
    vals = special.pdtrik(q,mu)
    sv = errp(sv)
    return where((mu<0) | (q<0) | (q>1),scipy.nan,vals)

def poissonisf(p,mu):
    return poissonppf(1-p,mu)

def poissonstats(mu,full=0):
    cond = (arr(mu) < 0)
    where(cond,scipy.nan,mu)
    var = mu
    if not full:
        return mu, var
    g1 = 1.0/arr(sqrt(mu))
    g2 = 1.0 / arr(mu)
    return mu, var, g1, g2


## Discrete Uniform

def randintpdf(k, min, max=None):
    if max is None:
        max = min
        min = 0
    cond = (arr(max) <= arr(min))
    fact = 1.0 / arr(max - min)
    return select([cond,(min<=k) & (k < max) & (k==floor(k))],[scipy.nan, fact])

def randintcdf(k, min, max=None):
    if max is None:
        max = min
        min = 0
    cond = (arr(max) <= arr(min))
    k = arr(k)
    Ck = (k-min)*1.0/(max-min)
    return select([cond,k>max,k>min],[scipy.nan,1,Ck])
    
def randintsf(k, min, max=None):
    return 1.0-randintcdf(k, min, max)

def randintppf(q, min, max=None):
    if max is None:
        max = min
        min = 0
    q = arr(q)
    cond = (arr(min) < arr(max)) & (0<=q) & (q<=1)
    return where(cond, (max-min)*q + min, scipy.nan)

def randintisf(p, min, max=None):
    return randintppf(1-p, min, max)

def randintstats(min, max=None, full=0):
    if max is None:
        max = min
        min = 0
    m2, m1 = arr(max), arr(min)
    cond = (m1 < m2)
    mu = (m2 + m1 - 1.0) / 2
    mu2p = (m2**3-m1**3)/(3.0*(m2-m1)) - (m2+m1)/2.0 + 1.0/6.0
    var = mu2p - mu*mu
    if not full:
        return mu, var
    mu3p = (m2**4-m1**4)/(4.0*(m2-m1)) - (m2**3-m1**3)/(2.0*(m2-m1)) \
           +(m2+m1)/4.0
    mu3 = mu3p - 3*mu*mu2p + 2*mu**3
    g1 = mu3 / var**1.5

    mu4p = (m2**5-m1**5)/(5.0*(m2-m1)) - (m2**4-m1**4)/(2.0*(m2-m1)) \
           + (m2**3-m1**3)/(3.0*(m2-m1)) - 1 / 30.0
    mu4 = mu4p - 4*mu3p*mu + 6*mu2p*mu*mu - 3*mu**4
    g2 = mu4 / var**2 - 3.0
    return mu, var, g1, g2

# Zipf distribution

def zipfpdf(k, a=4.0):
    k, a = arr(k)*1.0, arr(a)
    sv = errp(0)
    Pk = 1.0 / arr(special.zeta(a,1) * k**a)
    sv = errp(sv)
    return select([a<=1.0,(k>0) & (floor(k)==k)],[scipy.nan,Pk])

def _zipfcdf(k, a):
    if a<=1.0:
        return scipy.nan
    if k<1:
        return 0.0
    j = arange(1.0,k+1)
    return sum(1.0/j**a)/special.zeta(a,1)
_vzipfcdf = special.general_function(_zipfcdf,'d')

def zipfcdf(k, a=4.0):
    return _vzipfcdf(k, a)

def zipfsf(k, a=4.0):
    return 1.0-zipfcdf(k, a)

def zipfppf(q, a=4.0):
    pass

def zipfisf(p, a=4.0):
    pass

def zipfstats(a=4.0, full=0):
    sv = errp(0)
    fac = arr(special.zeta(a,1))
    mu = special.zeta(a-1.0,1)/fac
    mu2p = special.zeta(a-2.0,1)/fac
    var = mu2p - mu*mu
    if not full:
        return mu, var
    mu3p = special.zeta(a-3.0,1)/fac
    mu3 = mu3p - 3*mu*mu2p + 2*mu**3
    g1 = mu3 / arr(var**1.5)

    mu4p = special.zeta(a-4.0,1)/fac
    mu4 = mu4p - 4*mu3p*mu + 6*mu2p*mu*mu - 3*mu**4
    g2 = mu4 / arr(var**2) - 3.0
    return mu, var, g1, g2


