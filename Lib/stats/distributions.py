# Functions to implement several important functions for 
#   for various Continous and Discrete Probability Distributions
#
# Author:  Travis Oliphant  2002
# 

import scipy
import scipy.special as special
import Numeric
from Numeric import alltrue, where, arange
from fastumath import *
errp = special.errprint
select = scipy.select
arr = Numeric.asarray

all = alltrue
## Special defines some of these distributions
##  using backwards argument order.  This
##  is because the underlying C-library uses this order.

## We'd also like to control what happens on domain errors in this
##   wrapper

_quantstr = "Quantile must be in [0,1]."
_posstr = "Parameter must be > 0."
_nonnegstr = "Parameter must be >= 0."

def _wc(cond, val):
    return where (cond, val, scipy.nan)

### Each distribution has up to 8 functions defined plus one function
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
##  <dist>hf -- Hazard function (PDF / SF)
##  <dist>chf -- Cumulative hazard function (-log(1-CDF)) 
##  <dist>mle -- maximum-likelihood estimation of parameters of the
##                distribution.

##  NANs are returned for unsupported parameters.
##    location, scale, and shape parameters are used as much as possible,
##    except in certain well known cases, where the shape is calle something
##    else (e.g. degrees of freedom).
##    The order is shape parameters (if needed), loc=0.0, scale=1.0
##    These are related to the common symbols in the docs.


##  skew is third central moment / variance**(1.5)
##  kurtosis is fourth central moment / variance**2 - 3


## References::

##  Documentation for ranlib, rv2, cdflib and
## 
##  Eric Wesstein's world of Mathematics http://mathworld.wolfram.com/
##      http://mathworld.wolfram.com/topics/StatisticalDistributions.html
##
##  Documentation to Regress+ by Michael McLaughlin
##
##  Engineering and Statistics Handbook (NIST)
##      http://www.itl.nist.gov/div898/handbook/index.htm
##
##  Johnson book on Univariate Distributions


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
    assert all((0<=p)&(p<=1)), _quanstr
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
    return 1.0/sqrt(2*pi)/scale * exp(-(x-loc)**2 / (2.0*scale*scale))

def normcdf(x, loc=0.0, scale=1.0):
    return special.ndtr((x-loc)*1.0/scale)

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

# full also returns skewness and kurtosis
def normstats(loc=0.0, scale=1.0, full=0):
    if not full:
        return loc, scale**2
    else:
        return loc, scale**2, 0, 0

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
    return betappf(1-q,a,b,loc,scale)

def betastats(a, b, loc=0.0, scale=1.0, full=0):
    cond = (arr(a)>0) & (arr(b) > 0) & (arr(scale)>0)
    mn = _wc(cond,loc + scale*a*1.0 / (a+b))
    var = _wc(cond, scale*scale*(a*b)*1.0 / ((a+b)**2 * (a+b+1)))
    if not full:
        return mn, var
    g1 = 2.0*(b-a)*sqrt((1.0+a+b) / (a*b*(2+a+b)))
    g2 = 6.0*(a**3 + a**2*(1-2*b) + b**2*(1+b) - 2*a*b*(2+b))
    g2 /= a*b*(a+b+2)*(a+b+3)
    return mn, var, _wc(cond, g1), _wc(cond, g2)

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
    mu = where(cond, g1c*g1cd / gd, scipy.nan)
    var = where(cond, k/gd**2, scipy.nan)
    if not full:
        return mu, var
    g3c, g3cd = gam(1-3.0/c), gam(3.0/c+d)
    g4c, g4cd = gam(1-4.0/c), gam(4.0/c+d)
    g1 = 2*g1c**3 * g1cd**3 + gd*gd*g3c*g3cd - 3*gd*g2c*g1c*g1cd*g2cd
    g1 /= sqrt(k**3)
    g2 = 6*gd*g2c*g2cd * g1c**2 * g1cd**2 + gd**3 * g4c*g4cd
    g2 -= 3*g1c**4 * g1cd**4 -4*gd**2*g13c*g1c*g1cd*g3cd
    return mu, var, _wc(cond, g1), _wc(cond, g2)

# Fisk distribution
# burr is a generalization

def fiskpdf(x, c, loc=0.0, scale=1.0):
    return burrpdf(x, c, d, loc, scale)

def fiskcdf(x, c, loc=0.0, scale=1.0):
    return burrcdf(x, c, d, loc, scale)

def fisksf(x, c, loc=0.0, scale=1.0):
    return burrsf(x, c, d, loc, scale)

def fiskppf(x, c, loc=0.0, scale=1.0):
    return burrppf(x, c, d, loc, scale)

def fiskisf(x, c, loc=0.0, scale=1.0):
    return burrisf(x, c, d, loc, scale)

def fiskstats(x, c, loc=0.0, scale=1.0):
    return burrstats(x, c, d, loc, scale)


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
    q, loc, scale = map(arr,(p,loc,scale))    
    sv = errp(0)
    vals = where((0<=q)&(q<=1)&(df>0)&(scale>0),
                 special.chdtri(df, 1-q), scipy.nan)
    sv = errp(sv)
    return vals*scale + loc

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
def cosinepdf(x, mean=0.0, scale=1.0):
    A, B, x = arr(mean), arr(scale), arr(x)
    Px = 1.0/2.0/pi/B*(1+cos((x-A)*1.0/B))
    return select([B<=0,(y>=(A-pi*B))&(y<=(A+pi*B))],[scipy.nan, Px])

def cosinecdf(x, mean=0.0, scale=1.0):
    A, B, x = arr(mean), arr(scale), arr(x)
    y = (x-A)*1.0/B
    Cx = 1.0/2/pi*(pi + y + sin(y))
    return select([B<=0, y>(A+pi*B), y>=(A-pi*B)],[scipy.nan, 1, Cx])

def cosinesf(x, mean=0.0, scale=1.0):
    return 1.0-cosinecdf(x, mean, scale)

def cosineppf(q, mean=0.0, scale=1.0):
    pass

def cosineisf(p, mean=0.0, scale=1.0):
    pass

def cosinestats(mean=0.0, scale=1.0, full=0):
    B = arr(scale)
    mu = where(B>0,mean,scipy.nan)
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
    Px = (x)**(n-1.0)*exp(-x)/special.gamma(a)
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
    vals = -1.0/lam*log(1-q)
    return select([(scale<=0)|(q>1)|(q<0)],[scipy.nan],vals*scale+loc)

def exponisf(q, loc=0.0, scale=1.0):
    return exponppf(1-q, loc, scale)

def exponstats(loc=0.0, scale=1.0, full=0):
    cond = (arr(scale) > 0)
    mn = _wc(cond, scale+loc)
    var = _wc(cond, scale*scale)
    if not full:
        return mn, var
    return mn, var, _wc(cond, 2), _wc(cond, 6)

    
## F

def fpdf(x, dfn, dfd):
    x = arr(x)
    m = 1.0*dfd
    n = 1.0*dfn
    Px = m**(m/2) * n**(n/2) * x**(n/2-1)
    Px /= (m+n*x)**((n+m)/2)*special.beta(n/2,m/2)
    return where(x>0,Px,0)

def fcdf(x, dfn, dfd):
    x = arr(x)
    x = where(x<0, 0, x)
    return special.fdtr(dfn, dfd, x)

def fsf(x, dfn, dfd):
    x = arr(x)
    x = where(x<0, 0, x)
    return special.fdtrc(dfn, dfd, x)

def fisf(p, dfn, dfd):
    assert all((0<=p)&(p<=1)), _quantstr
    return special.fdtri(dfn, dfd, p)
    
def fppf(q, dfn, dfd):
    assert all((0<=q) & (q<=1)), _quantstr
    return special.fdtri(dfn, dfd, 1-p)

def fstats(dfn, dfd, full=0):
    m = dfd*1.0
    n = dfn*1.0
    mn = m / (m-2)
    var = 2*m*m*(m+n-2)/(n*(m-2)**2 * (m-4))
    if not full:
        return mn, var
    g1 = 2*(m+2*n-2)/(m-6)*sqrt((2*m-4)/(n*(m+n-2)))
    g2 = 12*(-16+20*m-8*m*m + m**3 + 44*n \
             -32*m*n + 5*m*m*n-22*n*n + 5*m*n*n)
    g2 /= n*(m-6)*(m-8)*(n+m-2)
    return mn, var, g1, g2


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

def foldnormppf(al, c=0.0, loc=0.0, scale=1.0):
    pass

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
    return 1.0-genlogistcdf(x, c, loc, scale)

def genlogisticppf(q, c, loc=0.0, scale=1.0):
    q, c, loc, scale = map(arr, (q, c, loc, scale))
    vals = -log(pow(q,-1.0/c)-1)
    return where((c>0) & (scale>0) & (q>=0) & (q<=1), vals, scipy.nan)

def genlogisticisf(q, c, loc=0.0, scale=1.0):
    return genlogsistppf(1-q, c, loc=0.0, scale=1.0)

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


## Gumbel, Log-Weibull, Fisher-Tippett, Gompertz

def gumbelpdf(x,loc=0.0,scale=1.0):
    x, a, b = map(arr, (x, loc, scale))
    x = (x-a*1.0)/b
    fac = x+exp(-x)
    return select([scale>0],[exp(-fac)/b],scipy.nan)

def gumbelcdf(x,loc=0.0,scale=1.0):
    x, a, b = map(arr, (x, loc, scale))
    x = (x-a*1.0)/b
    return select([scale>0],[exp(-exp(-x))],scipy.nan)
    
def gumbelsf(x,loc=0.0,scale=1.0):
    return 1.0-gumbelcdf(x,a,b)

def gumbelppf(q,loc=0.0,scale=1.0):
    q = arr(q)
    cond = (arr(scale)>0) & (q >= 0) & (q <=1)
    vals = -log(-log(q))
    return _wc(cond, loc+scale*vals)

def gumbelisf(p,loc=0.0,scale=1.0):
    return gumbelppf(1-p,a,b)

def gumbelstats(loc=0.0,scale=1.0,full=0):
    a, b = map(arr, (loc, scale))
    cond = (b > 0)
    mn = _wc(cond, a + b*_EULER)
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
    return select([cond,q==1],[log(tan(pi*q/2.0)),scipy.inf],scipy.nan)

def hypsecantsf(x, loc=0.0, scale=1.0):
    return 1.0-hypsecantcdf(x, loc, scale)

def hypsecantisf(q, loc=0.0, scale=1.0):
    return hypsecantppf(1-q, loc, scale)

def hypsecantstats(loc=0.0, scale=1.0, full=0):
    cond = (arr(scale)>0)
    mn = _wc(cond, loc)
    var = _wc(cond, pi*pi/4*scale*scale)
    if not full:
        return mn, var
    g1, g2 = 0, 2
    return mn, var, _wc(cond, g1), _wc(cond, g2)


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
    return select([(mu<=0)|(scale<=0),x>0],[scipy.nan, Cx])

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
    return invnormppf(1-q, mu, loc, scale)

def invnormstats(mu, loc=0.0, scale=1.0, full=0):
    cond = (arr(mu)>0) & (arr(scale) > 0)
    mn = _wc(cond, mu*scale + loc)
    var = _wc(cond, mu**3 * scale**2)
    if not full:
        return mn, var
    g1 = 3*sqrt(mu)
    g2 = 15*mu
    return mn, var, _wc(cond,g1), _wc(cond,g2)

    

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
    return laplaceppf(1-q, loc, scale)

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
    return logisticppf(1-q, loc, scale)

def logisticstats(loc=0.0, scale=1.0, full=0):
    mn = loc
    var = pi*pi/3.0*scale*scale
    if not full:
        return mn, var
    g1 = 0.0
    g2 = 6.0/5
    return mn, var, g1, g2

## Lognorm
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



# MAXWELL
#  a special case of chi with df = 3, loc=0.0, and given scale = 1.0/sqrt(a)
#    where a is the parameter used in Mathworld description

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


# Nakagami (cf Chi)

def nakagamipdf(x, df, loc=0.0, scale=1.0):
    x, loc, scale = map(arr, (x, loc, scale))
    x = arr((x-loc*1.0)/scale)
    Px = 2*df**df*x**(2*df-1.0)/special.gamma(df)*x**(2*df-1.0)*exp(-df*v*x)
    return select([(df<=0)|(scale<=0),x>0],[scipy.nan,Px/scale])
    

# Non-central chi-squared

def ncx2pdf(x,df,nc):
    assert all(nc>=0), _nonnegstr
    assert all(df>0), _posstr
    x = arr(x)
    y = where(x<=0,1,x)
    a = df/2.0
    Px = 0.5**a*exp(-(nc+y)/2.0)*(y*1.0)**(a-1)
    z = nc*y/4.0
    Px *= special.hyp0f1(a,z)/special.gamma(a)
    return where(x<=0,0,Px)

##def _ncx2cdf(x,df,nc):
##    from scipy.limits import double_epsilon as eps    
##    nc = nc/2.0
##    jmid = floor(nc)
##    val = 0
##    for j in range(jmid,-1,-1):
##        term = poissonpdf(j,nc)*chi2cdf(x,df+2*j)
##        val += term
##        if (term < val*eps):
##            print "Here", j
##            break
##    for j in range(jmid+1,2*jmid+2):
##        term = poissonpdf(j,nc)*chi2cdf(x,df+2*j)
##        val += term
##        if (term < val*eps):
##            break
##    if (j == 2*jmid+1):
##        for j in range(2*(jmid+1),2*(jmid+1)+1000):
##            term = poissonpdf(j,nc)*chi2cdf(x,df+2*j)
##            val += term
##            if (term < val*eps):
##                print "final break at ", j
##                break
##        if (j==(2*jmid+1002)):
##            print "Warning:  Did not converge."
##    return val
        
#def _ncx2cdf(x,df,nc):
#    import scipy.integrate as integrate
#    return integrate.quad(ncx2pdf,0,x,args=(df,nc))[0]
#_vec_ncx2cdf = special.general_function(_ncx2cdf,'d')

#def ncx2cdf(x,df,nc):
#    assert all(nc>=0), _nonnegstr
#    assert all(df>0), _posstr
#    x = where(x<0,0,x)
#    return _vec_ncx2cdf(x,df,nc)

def ncx2cdf(x,df,nc):
    assert all(nc>=0), _nonnegstr
    assert all(df>0), _posstr
    x = arr(x)
    x = where(x<0,0,x)
    return special.chndtr(x,df,nc)

def ncx2sf(x,df,nc):
    return 1-ncx2cdf(x,df,nc)
    
##def _ncx2qfunc(x,q,df,nc):
##    return _ncx2cdf(x,df,nc)-q

##def _ncx2ppf(q,df,nc):
##    import scipy.optimize as optimize
##    return optimize.fsolve(_ncx2qfunc,nc+df,args=(q,df,nc))
##_vec_ncx2q = special.general_function(_ncx2q,'d')

##def ncx2ppf(q,df,nc):
##    assert all((0<=q) & (q<=1)), _quantstr
##    return _vec_ncx2ppf(q, df, nc)

def ncx2ppf(q,df,nc):
    assert all((0<=q) & (q<=1)), _quantstr
    return special.chndtrix(q,df,nc)
    
def ncx2isf(p,df,nc):
    assert all((0<=p)&(p<=1)), _quantstr
    return ncx2ppf(1-p,df,nc)

def ncx2stats(df,nc,full=0):
    mn = nc+df
    var = 2*(2*nc+df)
    if not full:
        return mn, var
    g1 = 2*sqrt(2)*(3*nc+df)/(2*nc+df)**1.5
    g2 = 12.0*(4*nc+df)/(2*nc+df)**2.0
    return mn, var, g1, g2

# Non-central F

def ncfpdf(x,n1,n2,nc):
    assert all((n1>0) & (n2>0)), _posstr
    assert all(nc>=0), _nonnegstr
    x = arr(x)
    n1 = n1*1.0
    n2 = n2*1.0
    nc = nc*1.0
    y = where(x<0,1.0,x)
    Px = exp(-nc/2+nc*n1*y/(2*(n2+n1*y)))
    Px *= n1**(n1/2) * n2**(n2/2) * y**(n1/2-1)
    Px *= (n2+n1*y)**(-(n1+n2)/2)
    Px *= special.gamma(n1/2)*special.gamma(1+n2/2)
    Px *= special.assoc_laguerre(-nc*n1*y/(2.0*(n2+n1*y)),n2/2,n1/2-1)
    Px /= special.beta(n1/2,n2/2)*special.gamma((n1+n2)/2.0)
    return where(x<0,0.0,Px)

##def _ncfcdf(x,n1,n2,nc):
##    from scipy.limits import double_epsilon as eps
##    eps2 = sqrt(eps)
##    eps4 = sqrt(eps2)
##    n1, n2, nc = n1/2.0, n2/2.0, nc/2.0
##    val = 0
##    j = 0
##    valm1 = 100
##    term = 80
##    gam = special.gamma
##    bI = special.betainc
##    bval =  n1*x/(n2+n1*x)
##    jmid = floor(nc)
##    for j in range(jmid,-1,-1):
##        term = poissonpdf(j,nc)*betacdf(bval,j+n1,n2)
##        val += term
##        if all(ravel(term / (val+eps4)) < eps2):
##            break
##    for j in range(jmid+1,jmid+2000):
##        term = poissonpdf(j,nc)*betacdf(bval,j+n1,n2)
##        val += term
##        if all(ravel(term / (val+eps4)) < eps2):
##            break
##    if (j == jmid+2000-1):
##        print "Warning: Series failed to converge."
##    return val
##_vec_ncfcdf = special.general_function(_ncfcdf,'d')

##def ncfcdf(x,dfn,dfd,nc):
##    assert all(nc>=0), _nonnegstr
##    assert all((dfn>0) & (dfd>0)), _posstr
##    x = where(x<0,0,x)
##    return _vec_ncfcdf(x,dfn,dfd,nc)

def ncfcdf(x,dfn,dfd,nc):
    assert all(nc>=0), _nonnegstr    
    assert all((dfn>0) & (dfd>0)), _posstr
    x = arr(x)
    x = where(x<0,0,x)
    return special.ncfdtr(dfn,dfd,nc,x)
    
def ncfsf(x,dfn,dfd,nc):
    return 1-ncfcdf(x,dfn,dfd,nc)

##def _ncfqfunc(x,q,dfn,dfd,nc):
##    return _ncfcdf(x,dfn,dfd,nc)-q

##def _ncfppf(q,dfn,dfd,nc,x0):
##    import scipy.optimize as optimize
##    return optimize.fsolve(_ncfqfunc,x0,args=(q,dfn,dfd,nc))
##_vec_ncfq = special.general_function(_ncfq,'d')

##def ncfppf(q,dfn,dfd,nc,x0=None):
##    assert all((0<=q) & (q<=1)), _quanstr
##    if x0 is None:
##        x0 = dfd * (dfn+nc)/(dfn*(dfd-2))
##    return _vec_ncfppf(q, dfn, dfd, nc, x0)

def ncfppf(q, dfn, dfd, nc):
    assert ((0<=q) & (q<=1))
    return special.ncfdtri(dfn, dfd, nc, q)

def ncfisf(p,dfn,dfd,nc):
    return ncfppf(1-p,dfn,dfd,nc)

def ncfstats(dfn,dfd,nc, full=0):
    dfn = arr(dfn)*1.0
    dfd = arr(dfd)*1.0
    nc = arr(nc)*1.0
    mn = where(dfd<=2,scipy.nan,dfd/(dfd-2)*(1+nc/dfn))
    var1 = 2*(dfd/dfn)**2 * ((dfn+nc/2)**2+(dfn+nc)*(dfd-2))
    var1 /= (dfd-2)**2 * (dfd-4)
    var = where(dfd<=4,scipy.nan,var1)               
    if full:
        print "Skew and kurtosis unavailable for this distribution."
    return mn, var


## Student t distribution

def tpdf(x, df):
    assert all(df > 0), _posstr
    r = df*1.0
    Px = exp(special.gammaln((r+1)/2)-special.gammaln(r/2))
    Px /= sqrt(r*pi)*(1+(x**2)/r)**((r+1)/2)
    return Px

def tcdf(x, df):
    assert all(df > 0), _posstr
    return special.stdtr(df, x)

def tsf(x, df):
    return 1-tcdf(x,df)

def tppf(q, df):
    assert all((0<=q) & (q<=1)), _quantstr
    return special.stdtrit(df, q)

def tisf(p, df):
    assert all((0<=p)&(p<=1)), _quantstr
    return special.stdtrit(df, 1-p)

def tstats(df,full=0):
    r = df*1.0
    if r <= 1:
        mn = scipy.nan
    else:
        mn = 0.0
    if r <= 2:
        var = scipy.nan
    else:
        var = r/(r-2)
    if not full:
        return mn, var

    if r <= 1:
        g1 = scipy.nan
    else:
        g1 = 0.0
    if r <= 4:
        g2 = scipy.nan
    else:
        g2 = 6 / (r-4)
        
    return mn, var, g1, g2

## Non-central T distribution

def nctpdf(x, df, nc):
    assert all((df > 0) & (nc > 0)), _posstr
    n = df*1.0
    nc = nc*1.0
    x2 = x*x
    ncx2 = nc*nc*x2
    fac1 = n + x2
    Px = n**(n/2) * special.gamma(n+1)
    Px /= 2.0**n*exp(nc*nc/2)*fac1**(n/2)*special.gamma(n/2)
    valF = ncx2 / (2*fac1)
    trm1 = sqrt(2)*nc*x*special.hyp1f1(n/2+1,1.5,valF)
    trm1 /= fac1*special.gamma((n+1)/2)
    trm2 = special.hyp1f1((n+1)/2,0.5,valF)
    trm2 /= sqrt(fac1)*special.gamma(n/2+1)
    Px *= trm1+trm2

def nctcdf(x,df,nc):
    assert all((df > 0) & (nc > 0)), _posstr
    return special.nctdtr(df, nc, x)

def nctsf(x,df,nc):
    return 1-nctcdf(x,df,nc)

##def _nctqfunc(x,q,df,nc):
##    return _nctcdf(x,dfn,dfd,nc)-q

##def _nctppf(q,df,nc,x0):
##    import scipy.optimize as optimize
##    return optimize.fsolve(_nctqfunc,x0,args=(q,dfn,dfd,nc))
##_vec_nctq = special.general_function(_nctq,'d')

##def nctppf(q,df,nc,x0=None):
##    assert all((0<=q) & (q<=1)), _quanstr
##    if x0 is None:
##        val1 = gam((df-1.0)/2.0)
##        val2 = gam(df/2.0)
##        x0 = nc*sqrt(df/2.0)*val1/val2
##    return _vec_ncfppf(q, dfn, dfd, nc, x0)

def nctppf(q,df,nc):
    assert all((0<=q) & (q<=1)), _quanstr
    assert all((df > 0) & (nc > 0)), _posstr
    return special.nctdtrit(df, nc, q)

def nctisf(p,df,nc):
    return nctppf(1-p,df,nc)

def nctstats(df,nc,full=0):
    assert all((df > 0) & (nc > 0)), _posstr
    nc = nc*1.0
    df = df*1.0
    gam = special.gamma
    val1 = gam((df-1.0)/2.0)
    val2 = gam(df/2.0)
    mn = nc*sqrt(df/2.0)*val1/val2
    var = (nc*nc+1.0)*df/(df-2.0)
    var -= nc*nc*df* val1**2 / 2.0 / val2**2
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
    return mn, var, g1n/g1d, g2n/g2d



# Pareto

def paretopdf(x, mode=1.0, shape=4.0):
    a, b = shape, mode
    assert all((a > 0) & (b > 0)), _posstr
    x = arr(x)
    return where(x<b,0,1.0*a*b**a / x**(a+1.0))

def paretocdf(x, mode=1.0, shape=4.0):
    a, b = shape, mode
    assert all((a > 0) & (b > 0)), _posstr
    x = arr(x)
    x = where(x<b,b,x)
    return 1-(b*1.0/x)**a

def paretosf(x, mode=1.0, shape=4.0):
    assert all((mode > 0) & (shape > 0)), _posstr
    return 1-paretocdf(x,mode,shape)

def paretoppf(q, mode=1.0, shape=4.0):
    a, b = shape, mode
    assert all((a > 0) & (b > 0)), _posstr
    assert all(0<=q<1), _quantstr
    return b*pow(1-q,-1.0/a)

def paretoisf(p, mode=1.0, shape=4.0):
    a, b = shape, mode
    assert all(0<=q<1), _quantstr
    assert all((a > 0) & (b > 0)), _posstr
    return b*pow(p,-1.0/a)

def paretostats(mode=1.0, shape=4.0, full=0):
    a, b = shape, mode
    assert all((a > 0) & (b > 0)), _posstr
    mn = a*b / (a-1.0)
    var = a*b*b / (a-1.0)**2 / (a-2.0)
    if not full:
        return mn, var

    g1 = sqrt((a-2.0)/a)*2.0*(a+1.0)/(a-3.0)
    g2 = 6*scipy.polyval([1,1,-6,-2],a)/(a*(a-3.0)*(a-4))
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
    return select([(a<=0)|(scale<=0),(x>1),(x>=0)],[scipy.nan, 1, Px/scale])

def powersf(x, a, loc=0.0, scale=1.0):
    return 1.0-powercdf(x, a, loc, scale)

def powerppf(q, a, loc=0.0, scale=1.0):
    q, a, loc, scale = map(arr, (q, a, loc, scale))
    vals = pow(q, 1.0/a)
    cond = (q>=0)&(q<=1)&(a>0)&(scale>0)
    return _wc(cond, vals)

def powerisf(q, a, loc=0.0, scale=1.0):
    return powerppf(1-q, a, loc, scale)

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

# Uniform

def uniformpdf(x, a=0.0, b=1.0):
    assert all(b > a)
    x = arr(x)
    return where((x>a) & (x<b), 1.0/(maxi-mini), 0.0)

def uniformcdf(x, a=0.0, b=1.0):
    assert all(b > a)
    x = arr(x)
    x = where(x<a,a,x)
    x = where(x>b,b,x)*1.0
    return (x-a)/(b-a)

def uniformsf(x, a=0.0, b=1.0):
    return 1-uniformcdf(x,a,b)

def uniformppf(q, a=0.0, b=1.0):
    assert all(b > a)
    assert all((0<=q) & (q<=1)), _quantstr
    return q*(b-a) + a

def uniformisf(p, a=0.0, b=1.0):
    assert all(b > a)
    assert all((0<=p)&(p<=1)), _quantstr
    return b - p*(b-a)

def uniformstats(a=0.0, b=1.0, full=0):
    assert all(b > a)
    mn = 0.5*(a+b)
    var = (b-a)**2.0 / 12.0
    if not full:
        return mn, var
    g1 = 0.0
    g2 = -6.0/5
    return mn, var, g1, g2


# Triangular

_trstr = "Left must be <= mode which must be <= right with left < right"

def triangpdf(x, left=0.0, mode=0.5, right=1.0):
    x = arr(x)
    a, b, c = left, right, mode
    assert all((a <= c <= b) & (a<b)), _trstr
    x = where(x<a,a,x)
    x = where(x>b,b,x)
    P1 = 2.0*(x-a) / (b-a) / (c-a)
    P2 = 2.0*(b-x) / (b-a) / (b-c)
    return where(x<c,P1,P2)
    
def triangcdf(x,left=0.0, mode=0.5, right=1.0):
    x = arr(x)
    a, b, c = left, right, mode
    assert all((a <= c <= b) & (a<b)), _trstr
    x = where(x<a,a,x)
    x = where(x>b,b,x)
    C1 = (x-a)**2.0 / (b-a) / (c-a)
    C2 = 1.0-(b-x)**2.0 / (b-a) / (b-c)
    return where(x<c,C1,C2)

def triangsf(x, left=0.0, mode=0.5, right=1.0):
    return 1.0-triangcdf(x, left, mode, right)

def triangppf(q, left=0.0, mode=0.5, right=1.0):
    a, b, c = left, right, mode
    assert all((a <= c <= b) & (a<b)), _trstr
    assert all((0<=q) & (q<=1)), _quantstr
    q = arr(q)
    Dc = (c-a)/(b-a)
    x1 = b - sqrt((1-q)*(b-a)*(b-c))
    x2 = a + sqrt(q*(b-a)*(c-a))
    return where(q > Dc, x1, x2)

def triangisf(p, left=0.0, mode=0.5, right=1.0):
    return triangppf(1.0-p, left, mode, right)

def triangstats(left=0.0, mode=0.5, right=1.0, full=0):
    a, b, c = left, right, mode
    assert all((a <= c <= b) & (a<b)), _trstr
    mu = (a+b+c)/3.0
    trm = (a*a + b*b + c*c - a*b - a*c - b*c)
    var =  trm / 18.0
    if not full:
        return mu, var
    mu3 = (a+b-2*c)*(a+c-2*b)*(b+c-2*a)/270.0
    mu4 = trm / 135.0
    g1 = mu3 / var**1.5
    g2 = mu4 / var**2.0 - 3
    return mu, var, g1, g2

# Von-Mises

def von_misespdf(x,mode=0.0, shape=1.0):
    x = arr(x)
    a, b = mode, shape
    assert (-pi<=a<=pi)
    assert (b>0), _posstr
    box = where(x<=-pi,0,1)
    box = where(x>=pi,0,box)
    Px = exp(b*cos(x+pi-a)) / (2.0*pi*special.i0(b))
    return Px*box

def von_misescdf(x, mode=0.0, shape=1.0):
    x = arr(x)
    a, b = mode, shape
    assert (-pi<=a<=pi)
    assert (b>0), _posstr
    from scipy.limits import double_epsilon as eps
    eps2 = sqrt(eps)
    x = where(x<-pi,-pi,x)
    x = where(x>pi,pi,x)
    fac = special.i0(b)
    x2 = x + pi
    val = x2 / 2.0/ pi
    for j in range(1,501):
        trm1 = special.iv(j,b)/j/fac
        trm2 = sin(j*(x2-a))/pi
        val += trm1*trm2
        if all(trm1 < eps2):
            break
    if (j == 500):
        print "Warning: did not converge..."
    return val

def _vmqfunc(x,q,mode,shape):
    return von_misescdf(x,mode,shape)-q

def _vmppf(q,mode,shape,x0):
    import scipy.optimize as optimize
    return optimize.fsolve(_vmqfunc,x0,args=(q,mode,shape))
_vec_vmppf = special.general_function(_vmppf,'d')

def von_misesppf(q,mode=0.0, shape=1.0, x0=None):
    assert all((0<=q) & (q<=1)), _quanstr
    if x0 is None:
        x0 = mode
    return _vec_vmppf(q, mode, shape,x0)

def von_misesisf(p, mode=0.0, shape=1.0, x0=None):
    return von_misesppf(1-p, mode, shape, x0)

def von_misesstats(mode=0.0, shape=1.0, full=0):
    mn = mode
    var = 1.0 - special.iv(mode, shape) / special.i0(shape)
    if not full:
        return mn, var
    print "Skew and Kertosis not available---returned as NaN"
    return mn, var, scipy.nan, scipy.nan


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

def weibullpdf(x, shape, loc=0.0, scale=1.0):
    c, b, A, x = map(arr,(shape, scale, loc, x))
    x = arr((x-A*1.0)/b)
    Px = c * x**(c-1.0) * exp(-x**c)
    return select([(c<=0)|(b<=0),x>0],[scipy.nan,Px/b])

def weibullcdf(x, shape, loc=0.0, scale=1.0):
    c, b, A, x = map(arr,(shape, scale, loc, x))
    x = arr((x-A*1.0)/b)
    Cx = -special.expm1(-x**c)
    return select([(c<=0)|(b<=0),x>0],[scipy.nan,Cx])

def weibullsf(x, shape, loc=0.0, scale=1.0):
    c, b, A, x = map(arr,(shape, scale, loc, x))
    x = arr((x-A*1.0)/b)
    Cx = exp(-x**c)
    return select([(c<=0)|(b<=0),x>0],[scipy.nan,Cx],1)

def weibullppf(q, shape, loc=0.0, scale=1.0):
    a, b, loc, q = map(arr,(shape, scale, loc, q))
    cond1 = (a>0) & (b>0) & (0<=q) & (q<=1)
    vals = b*pow(arr(log(1.0/(1-q))),1.0/a) + loc
    return where(cond1, vals, scipy.nan)

def weibullisf(q, shape, loc=0.0, scale=1.0):
    a, b, loc, q = map(arr,(shape, scale, loc, q))
    cond1 = (a>0) & (b>0) & (0<=q) & (q<=1)
    vals = b*pow(arr(log(1.0/q)),1.0/a) + loc
    return where(cond1, vals, scipy.nan)

def weibullstats(shape, loc=0.0, scale=1.0, full=0):
    a, loc, b = map(arr,(shape, loc, scale))
    cond = (a>0) & (b>0)
    gam = special.gamma
    ia = 1.0/a
    mn = _wc(cond, b*gam(1+ia)+loc)
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
        
### DISCRETE DISTRIBUTIONS
###


# Binomial 
    
def binompdf(k, n, pr=0.5):
    k = arr(k)
    assert (0<pr<1)
    cond = arr((k > 0) & (k == floor(k)))
    return scipy.comb(n,k)* pr**k * (1-pr)**(n-k)

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

def bernoullip(k, pr=0.5):
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
    return 1.0-logseriescdf(k, pr=pr)

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
    N = arr(m2 - m1)
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


