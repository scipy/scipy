# Functions to implement CDF, PDF, Quanitles, CDFC, and STATS
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

### Each distribution has up to 6 functions defined plus one function
##    to return random variates following the distribution ---
##    these functions are in (this is in rv2.py or rv.py).

##  if <dist> is the name of the function to return random variates, then
##  <dist>pdf --- PDF
##  <dist>cdf --- CDF
##  <dist>cdfc --- Complementary CDF
##  <dist>q --- Inverse of CDF (quantiles)
##  <dist>p --- Inverse of CDFC
##  <dist>stats --- Return mean, variance and optionally (Fisher's) skew and kurtosis
##                      of the distribution.

##  skew is third central moment / variance**(1.5)
##  kurtosis is fourth central moment / variance**2 - 3


## References::

##  Documentation for ranlib, rv2, cdflib and
## 
##  Eric Wesstein's world of Mathematics http://mathworld.wolfram.com/
##      http://mathworld.wolfram.com/topics/StatisticalDistributions.html
##
##  Documentation to Regress+ by Michael McLaughlin 


## Kolmogorov-Smirnov one-sided and two-sided test statistics

def ksonecdfc(x,n):
    return special.smirnov(n,x)

def ksonep(p,n):
    return special.smirnovi(n,p)

def ksonecdf(x,n):
    return 1-special.smirnov(n,x)

def ksoneq(q,n):
    assert all((0<=q) & (q<=1)), _quantstr
    return special.smirnovi(n,1-q)

def kstwocdfc_largen(y):
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

def stnormcdfc(x):
    return 1-special.ndtr(x)

def stnormq(q):
    q = arr(q)
    sv = errp(0)
    vals = where((0<=q) & (q<=1),special.ndtri(q),scipy.nan)
    sv = errp(sv)
    return vals

def stnormp(p):
    p = arr(p)
    sv = errp(0)
    vals = where((0<=p)&(p<=1),special.ndtri(1-p),scipy.nan)
    sv = errp(sv)
    return vals

def stnormstats(full=0):
    if not full:
        return 0, 1
    return 0, 1, 0, 0

def normpdf(x, mu=0.0, std=1.0):
    return 1.0/sqrt(2*pi)/std * exp(-(x-mu)**2 / (2.0*std*std))

def normcdf(x, mu=0.0, std=1.0):
    return special.ndtr((x-mu)*1.0/std)

def normq(q, mu=0.0, std=1.0):
    q = arr(q)
    sv = errp(0)
    vals = where((0<=q) & (q<=1),special.ndtri(q)*std+mu,scipy.nan)
    sv = errp(sv)    
    return vals

def normcdfc(x, mu=0.0, std=1.0):
    return 1-special.ndtr((x-mu)*1.0/std)

def normp(p, mu=0.0, std=1.0):
    return normq(1-p,mu,std)

# full also returns skewness and kurtosis
def normstats(mu=0.0, std=1.0, full=0):
    if not full:
        return mu, std**2
    else:
        return mu, std**2, 0, 0

## Beta distribution

def betapdf(x, a, b):
    x = arr(x)
    sv = errp(0)
    Px = (1.0-x)**(b-1.0) * x**(a-1.0)
    Px /= special.beta(a,b)
    sv = errp(sv)
    vals = where(0<=x<=1,Px,0)
    return vals
    
def betacdf(x, a, b):
    x = arr(x)
    x = where(x<0,0,x)
    x = where(x>1,1,x)
    return special.btdtr(a,b,x)

def betaq(q, a, b):
    sv = errp(0)
    vals = special.btdtri(a,b,q)
    sv = errp(sv)    
    return vals

def betacdfc(x, a, b):
    return 1-betacdf(x,a,b)

def betap(p, a, b):
    return betaq(1-p,a,b)

def betastats(a, b, full=0):
    cond = (arr(a)>0 & arr(b) > 0)
    mn = where(cond,a*1.0 / (a+b),scipy.nan)
    var = where(cond, (a*b)*1.0 / ((a+b)**2 * (a+b+1), scipy.nan))
    if not full:
        return mn, var
    g1 = 2.0*(b-a)*sqrt(1+a+b) / (sqrt(a*b)*(2+a+b))
    g2 = 6.0*(a**3 + a**2*(1-2*b) + b**2*(1+b) - 2*a*b*(2+b))
    g2 /= a*b*(a+b+2)*(a+b+3)
    return mn, var, where(cond, g1, scipy.nan), where(cond, g2, scipy.nan)

## Cauchy

def cauchypdf(x, median=0.0, scale=1.0):
    cond = scale>0
    Px = 1.0/pi*scale / ((x-median)**2 + scale**2)
    return where(cond, Px, scipy.nan)

def cauchycdf(x, median=0.0, scale=1.0):
    cond = scale > 0
    return where(cond, 0.5 + 1.0/pi*arctan((x-median)*1.0/scale), scipy.nan)

def cauchycdfc(x, median=0.0, scale=1.0):
    cond = scale > 0
    return where(cond, 0.5 - 1.0/pi*arctan((x-median)*1.0/scale), scipy.nan)

def cauchyq(q, median=0.0, scale=1.0):
    cond = ((0<=q) & (q<=1)) & (scale > 0)
    return where(cond, scale*tan(pi*(q-0.5))+median, scipy.nan)

def cauchyp(p, median=0.0, scale=1.0):
    return cauchyq(1-p, median, scale)
    
def cauchystats(median=0.0, scale=1.0, full=0):
    if not full:
        return scipy.nan, scipy.nan
    else:
        return scipy.nan, scipy.nan, scipy.nan, scipy.nan
        
## Chi-squared

def chi2pdf(x, df):
    x = arr(x)
    df = arr(df)
    sv = errp(0)
    Px = x**(df/2.0-1)*exp(-x/2.0)
    Px /= special.gamma(df/2.0)* 2**(df/2.0)
    sv = errp(sv)
    return select([df<=0,x>0],[scipy.nan,Px])

def chi2cdf(x, df):
    x = arr(x)
    df = arr(df)
    sv = errp(0)
    vals = select([df<=0,x>0],[scipy.nan,special.chdtr(df, x)])
    sv = errp(sv)
    return vals

def chi2cdfc(x, df):
    x = arr(x)
    df = arr(df)
    sv = errp(0)
    vals = select([df<=0, x>0],[scipy.nan,special.chdtrc(df,x)])
    sv = errp(sv)
    return vals

def chi2p(p, df):
    sv = errp(0)
    vals = where(((0<=p)&(p<=1)) & (df>0),special.chdtri(df, p),scipy.nan)
    sv = errp(sv)
    return vals

def chi2q(q, df):
    sv = errp(0)
    vals = where(((0<=q) & (q<=1)) & (df>0), special.chdtri(df, 1-q), scipy.nan)
    sv = errp(sv)
    return vals

def chi2stats(df, full=0):
    cond = arr(df) > 0
    mn = where(cond,df,scipy.nan)
    var = where(cond,2*df,scipy.nan)
    if not full:
        return mn, var
    g1 = where(cond,2*sqrt(2.0/df),scipy.nan)
    g2 = where(cond, 12.0/df, scipy.nan)
    return mn, var, g1, g2


## Exponential

def exponpdf(x, lam):
    x = arr(x)
    return where(x<0, 0, lam*exp(-lam*x))

def exponcdf(x, lam):
    x = arr(x)
    x = where(x<0, 0, x)
    return 1.0-exp(-lam*x)

def exponcdfc(x, lam):
    return 1.0-exponcdf(x,lam)

def exponq(q, lam):
    assert all((0<=q) & (q<=1)), _quantstr
    return -1.0/lam*log(1-q)

def exponp(p, lam):
    return exponq(1-p, lam)

def exponstats(lam, full=0):
    th = 1.0/lam
    mn = th
    var = th*th
    if not full:
        return mn, var
    return mn, var, 2, 6

## (Fisher-Tippett)

def fisher_tippettpdf(x,a,b):
    fac1 = (a-x)*1.0/b
    fac2 = exp(fac1)
    return exp(fac1-fac2)/b

def fisher_tippettcdf(x,a,b):
    fac1 = (a-x)*1.0/b
    fac2 = exp(fac1)
    return exp(-fac2)

def fisher_tippettcdfc(x,a,b):
    return 1.0-fisher_tippetcdf(x,a,b)

def fisher_tippettq(q,a,b):
    assert all((0<=q) & (q<=1)), _quantstr
    return a-b*log(-log(q))

def fisher_tippettp(p,a,b):
    return fisher_tippettq(1-p,a,b)

def fisher_tippettstats(a,b,full=0):
    euler = 0.5772156649015328606
    zeta3 = special.zeta(3,1)
    mn = a + b*euler
    var = pi*pi/6*b*b
    if not full:
        return mn, var
    g1 = 12*sqrt(6)/pi**3 * zeta3
    g2 = 12.0/5.0
    return mn, var, g1, g2
    
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

def fcdfc(x, dfn, dfd):
    x = arr(x)
    x = where(x<0, 0, x)
    return special.fdtrc(dfn, dfd, x)

def fp(p, dfn, dfd):
    assert all((0<=p)&(p<=1)), _quantstr
    return special.fdtri(dfn, dfd, p)
    
def fq(q, dfn, dfd):
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

## Gumbel's

def gumbelpdf(x,mode=0.0,scale=1.0):
    return fisher_tippettpdf(x,mode,scale)

def gumbelcdf(x,mode=0.0,scale=1.0):
    return fisher_tippettcdf(x,mode,scale)

def gumbelcdfc(x,mode=0.0,scale=1.0):
    return fisher_tippettcdfc(x,mode,scale)

def gumbelq(q,mode=0.0,scale=1.0):
    return fisher_tippettq(q,mode,scale)

def gumbelp(p,mode=0.0,scale=1.0):
    return fisher_tippettp(q,mode,scale)

def gumbelstats(mode=0.0,scale=1.0,full=0):
    return fisher_tippettstats(mode,scale,full=full)

## Gamma (Use MATLAB and MATHEMATICA (b=theta, a=alpha) definition)

def gammapdf(x, a, b):
    x = arr(x)
    Px = x**(a-1.0)*exp(-x*1.0/b)
    Px /= special.gamma(a) * b**a
    return where(x<0,0,Px)

def gammacdf(x, a, b):
    x = arr(x)
    x = where(x<0,0,x)
    return special.gdtr(1.0/b,a,x)

def gammacdfc(x, a, b):
    x = arr(x)
    x = where(x<0,0,x)
    return special.gdtrc(1.0/b,a,x)
    
def gammaq(q, a, b):
    assert all((0<=q) & (q<=1)), _quantstr
    return special.gdtri(1.0/b,a,q)

def gammap(p, a, b):
    assert all((0<=p)&(p<=1)), _quantstr
    return special.gdtri(1.0/b,a,1-p)

def gammastats(a, b, full=0):
    mn = a*b
    var = a*b*b
    if not full:
        return mn, var
    g1 = 2.0/sqrt(a)
    g2 = 6.0/a
    return mn, var, g1, g2

## Laplace Distribution

def laplacepdf(x, mu=0.0, scale=1.0):
    iscale = 1.0/scale
    return iscale/2.0*exp(-abs(x-mu)*iscale)
    
def laplacecdf(x, mu=0.0, scale=1.0):
    return 0.5*(1+sign(x-mu)*(1-exp(-abs(x-mu)*1.0/scale)))

def laplacecdfc(x, mu=0.0, scale=1.0):
    return 0.5*(1-sign(x-mu)*(1-exp(-abs(x-mu)*1.0/scale)))

def laplaceq(q, mu=0.0, scale=1.0):
    assert all((0<=q) & (q<=1)), _quantstr
    fac = 0.5+sign(q-0.5)*(0.5-q)
    return mu+sign(0.5-q)*scale*log(2*fac)

def laplacep(p, mu=0.0, scale=1.0):
    return laplaceq(1-p, mu, scale)

def laplacestats(mu=0.0, scale=1.0, full=0):
    mn = mu*1.0
    var = 2.0*scale**2
    if not full:
        return mn, var
    g1 = 0
    g2 = 3
    return mn, var, g1, g2

## Logistic

def logisticpdf(x, mu=0.0, scale=1.0):
    iscale = 1.0/scale
    fac1 = exp((x-mu)*iscale)
    Px = abs(iscale)*fac1
    Px /= (1+fac1)**2
    return Px

def logisticcdf(x, mu=0.0, scale=1.0):
    fac1 = exp((mu-x)*1.0/abs(scale))
    return 1.0/(1+fac1)

def logisticcdfc(x, mu=0.0, scale=1.0):
    return 1.0-logisticcdf(x,mu=mu,scale=scale)

def logisticq(q, mu=0.0, scale=1.0):
    q = arr(q)
    assert all((0<=q) & (q<=1)), _quantstr    
    return mu - abs(scale)*log((1.0-q)/q)

def logisticp(p, mu=0.0, scale=1.0):
    return logisticq(1-p, mu, scale)

def logisticstats(mu=0.0, scale=1.0, full=0):
    mn = mu
    var = pi*pi/3.0*scale*scale
    if not full:
        return mn, var
    g1 = 0.0
    g2 = 6.0/5
    return mn, var, g1, g2

## Lognorm

def lognormpdf(x, mu=0.0, std=1.0):
    x = arr(x)
    Px = exp(-(log(x)-mu)**2 / (2*std**2))
    Px /= std*x*sqrt(2*pi)
    return where(x<=0,0,Px)

def lognormcdf(x, mu=0.0, std=1.0):
    x = arr(x)
    x = where(x<0,0.0,x)
    return 0.5*(1+special.erf((log(x)-mu)/(sqrt(2)*std)))

def lognormcdfc(x, mu=0.0, std=1.0):
    return 1-lognormcdf(x, mu, std)

def lognormq(q, mu=0.0, std=1.0):
    return exp(normalq(q,mu,std))

def lognormp(p, mu=0.0, std=1.0):
    return exp(normalp(p,mu,std))

def lognormstats(mu=0.0, std=1.0, full=0):
    s2 = std*std
    mn = exp(mu+s2/2.0)
    fac1 = exp(s2)-1
    var = exp(s2+2*mu)*fac1
    if not full:
        return mn, var
    g1 = sqrt(fac1)*(2+exp(s2))
    g2 = exp(4*s2) + 2*exp(3*s2)+3*exp(2*s2)-6.0
    return mn, var, g1, g2
    
# Gibrat's distribution is just default lognormal (think of it
#   as a distribution with one argument).

gilbratpdf = lognormpdf
gilbratcdf = lognormcdf
gilbratcdfc = lognormcdfc
gilbratq = lognormq
gilbratp = lognormp
gilbratstats = lognormstats

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

def ncx2cdfc(x,df,nc):
    return 1-ncx2cdf(x,df,nc)
    
##def _ncx2qfunc(x,q,df,nc):
##    return _ncx2cdf(x,df,nc)-q

##def _ncx2q(q,df,nc):
##    import scipy.optimize as optimize
##    return optimize.fsolve(_ncx2qfunc,nc+df,args=(q,df,nc))
##_vec_ncx2q = special.general_function(_ncx2q,'d')

##def ncx2q(q,df,nc):
##    assert all((0<=q) & (q<=1)), _quantstr
##    return _vec_ncx2q(q, df, nc)

def ncx2q(q,df,nc):
    assert all((0<=q) & (q<=1)), _quantstr
    return special.chndtrix(q,df,nc)
    
def ncx2p(p,df,nc):
    assert all((0<=p)&(p<=1)), _quantstr
    return ncx2q(1-p,df,nc)

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
    
def ncfcdfc(x,dfn,dfd,nc):
    return 1-ncfcdf(x,dfn,dfd,nc)

##def _ncfqfunc(x,q,dfn,dfd,nc):
##    return _ncfcdf(x,dfn,dfd,nc)-q

##def _ncfq(q,dfn,dfd,nc,x0):
##    import scipy.optimize as optimize
##    return optimize.fsolve(_ncfqfunc,x0,args=(q,dfn,dfd,nc))
##_vec_ncfq = special.general_function(_ncfq,'d')

##def ncfq(q,dfn,dfd,nc,x0=None):
##    assert all((0<=q) & (q<=1)), _quanstr
##    if x0 is None:
##        x0 = dfd * (dfn+nc)/(dfn*(dfd-2))
##    return _vec_ncfq(q, dfn, dfd, nc, x0)

def ncfq(q, dfn, dfd, nc):
    assert ((0<=q) & (q<=1))
    return special.ncfdtri(dfn, dfd, nc, q)

def ncfp(p,dfn,dfd,nc):
    return ncfq(1-p,dfn,dfd,nc)

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

def tcdfc(x, df):
    return 1-tcdf(x,df)

def tq(q, df):
    assert all((0<=q) & (q<=1)), _quantstr
    return special.stdtrit(df, q)

def tp(p, df):
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

def nctcdfc(x,df,nc):
    return 1-nctcdf(x,df,nc)

##def _nctqfunc(x,q,df,nc):
##    return _nctcdf(x,dfn,dfd,nc)-q

##def _nctq(q,df,nc,x0):
##    import scipy.optimize as optimize
##    return optimize.fsolve(_nctqfunc,x0,args=(q,dfn,dfd,nc))
##_vec_nctq = special.general_function(_nctq,'d')

##def nctq(q,df,nc,x0=None):
##    assert all((0<=q) & (q<=1)), _quanstr
##    if x0 is None:
##        val1 = gam((df-1.0)/2.0)
##        val2 = gam(df/2.0)
##        x0 = nc*sqrt(df/2.0)*val1/val2
##    return _vec_ncfq(q, dfn, dfd, nc, x0)

def nctq(q,df,nc):
    assert all((0<=q) & (q<=1)), _quanstr
    assert all((df > 0) & (nc > 0)), _posstr
    return special.nctdtrit(df, nc, q)

def nctp(p,df,nc):
    return nctq(1-p,df,nc)

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

def paretocdfc(x, mode=1.0, shape=4.0):
    assert all((mode > 0) & (shape > 0)), _posstr
    return 1-paretocdf(x,mode,shape)

def paretoq(q, mode=1.0, shape=4.0):
    a, b = shape, mode
    assert all((a > 0) & (b > 0)), _posstr
    assert all(0<=q<1), _quantstr
    return b*pow(1-q,-1.0/a)

def paretop(p, mode=1.0, shape=4.0):
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


# Rayleigh distribution

def rayleighpdf(r, mode=1.0):
    assert all(mode>0.0), _posstr
    r = arr(r)    
    return where(r<0,0,r*exp(-r*r/(2.0*mode**2))/mode**2)

def rayleighcdf(r, mode=1.0):
    assert all(mode>0.0), _posstr
    r = arr(r)
    r = where(r<0,0,r)
    return 1-exp(-r*r/(2.0*mode*mode))

def rayleighcdfc(r, mode=1.0):
    assert all(mode>0.0), _posstr
    r = arr(r)    
    r = where(r<0,0,r)
    return exp(-r*r/(2.0*mode*mode))

def rayleighq(q, mode=1.0):
    assert all(mode>0.0), _posstr
    assert all((0<=q) & (q<=1)), _quantstr
    return mode*sqrt(2*log(1.0/(1.0-q)))

def rayleighp(p, mode=1.0):
    assert all(mode>0.0), _posstr
    assert all((0<=p)&(p<=1)), _quantstr
    return mode*sqrt(2.0*log(1.0/p))

def rayleighstats(mode=1.0, full=0):
    s = mode
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

def uniformcdfc(x, a=0.0, b=1.0):
    return 1-uniformcdf(x,a,b)

def uniformq(q, a=0.0, b=1.0):
    assert all(b > a)
    assert all((0<=q) & (q<=1)), _quantstr
    return q*(b-a) + a

def uniformp(p, a=0.0, b=1.0):
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

def triangcdfc(x, left=0.0, mode=0.5, right=1.0):
    return 1.0-triangcdf(x, left, mode, right)

def triangq(q, left=0.0, mode=0.5, right=1.0):
    a, b, c = left, right, mode
    assert all((a <= c <= b) & (a<b)), _trstr
    assert all((0<=q) & (q<=1)), _quantstr
    q = arr(q)
    Dc = (c-a)/(b-a)
    x1 = b - sqrt((1-q)*(b-a)*(b-c))
    x2 = a + sqrt(q*(b-a)*(c-a))
    return where(q > Dc, x1, x2)

def triangp(p, left=0.0, mode=0.5, right=1.0):
    return triangq(1.0-p, left, mode, right)

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

def _vmq(q,mode,shape,x0):
    import scipy.optimize as optimize
    return optimize.fsolve(_vmqfunc,x0,args=(q,mode,shape))
_vec_vmq = special.general_function(_vmq,'d')

def von_misesq(q,mode=0.0, shape=1.0, x0=None):
    assert all((0<=q) & (q<=1)), _quanstr
    if x0 is None:
        x0 = mode
    return _vec_vmq(q, mode, shape,x0)

def von_misesp(p, mode=0.0, shape=1.0, x0=None):
    return von_misesq(1-p, mode, shape, x0)

def von_misesstats(mode=0.0, shape=1.0, full=0):
    mn = mode
    var = 1.0 - special.iv(mode, shape) / special.i0(shape)
    if not full:
        return mn, var
    print "Skew and Kertosis not available---returned as NaN"
    return mn, var, scipy.nan, scipy.nan


## Wald distribution (Inverse Normal)

def waldpdf(x, mean=1.0, scale=1.0):
    A, B = map(arr,(mean, scale))
    assert all((A > 0) & (B > 0)), _posstr
    x = arr(x)*1.0
    Px = sqrt(B/(2*pi*x**3.0))*exp(-B/(2*x)*((x-A)/A)**2)
    Px = where(x<=0,0,Px)
    return Px

def waldcdf(x, mean=1.0, scale=1.0):
    A, B = map(arr,(mean, scale))    
    assert all((A > 0) & (B > 0)), _posstr
    x = arr(x)*1.0
    x = where(x<0,0,x)
    sv = special.errprint(0)
    fac = sqrt(B/x)
    C1 = stnormcdf(fac*(x-A)/A)
    C1 += exp(2.0*B/A)*stnormcdf(-fac*(x+A)/A)
    sv = special.errprint(sv)
    return where(x==0,0.0,C1)

def waldcdfc(x, mean=1.0, scale=1.0):
    return 1-waldcdf(x, mean, scale)

def _waldqfunc(x,q,mean,scale):
    return waldcdf(x,mean,scale)-q

def _waldq(q,mean,scale,x0):
    import scipy.optimize as optimize
    return optimize.fsolve(_waldqfunc,x0,args=(q,mean,scale))
_vec_waldq = special.general_function(_waldq,'d')

def waldq(q, mean=1.0, scale=1.0, x0=None):
    A, B, q = map(arr,(mean, scale, q))
    assert all((A > 0) & (B > 0)), _posstr
    assert all((0<=q) & (q<=1)), _quanstr
    if x0 is None:
        x0 = mean
    qvals = _vec_waldq(q, mean, scale, x0)
    return where(qvals<0,0,qvals)

def waldp(p, mean=1.0, scale=1.0):
    return waldq(1-p, mean, scale)

def waldstats(mean=1.0, scale=1.0, full=0):
    A, B = map(arr,(mean, scale))    
    assert all((A > 0) & (B > 0)), _posstr
    iB = 1.0/B
    mn = A
    var = A**3 * iB
    if not full:
        return mn, var
    g1 = 3*sqrt(A*iB)
    g2 = 15*A*iB
    return mn, var, g1, g2

## Weibull

def weibullpdf(x, scale=1.0, shape=0.5):
    a, b, x = map(arr,(shape, scale, x))
    assert all((a>0) & (b>0)), _posstr
    ib = 1.0/ b
    x = arr(x * ib)
    Px = a * ib * x**(a-1.0) * exp(-x**a)
    return where(x<=0,0,Px)

def weibullcdf(x, scale=1.0, shape=0.5):
    a, b, x = map(arr,(shape, scale, x))
    assert all((a>0) & (b>0)), _posstr
    x = where(x<0,0,x)
    return -special.expm1(-(x*1.0/b)**a)

def weibullcdfc(x, scale=1.0, shape=0.5):
    a, b, x = map(arr,(shape, scale, x))
    assert all((a>0) & (b>0)), _posstr
    x = where(x<0,0,x)
    return exp(-(x*1.0/b)**a)

def weibullq(q, scale=1.0, shape=0.5):
    a, b, q = map(arr,(shape, scale, q))
    assert all((a>0) & (b>0)), _posstr
    assert all((0<=q) & (q<=1)), _quantstr
    return b*pow(log(1.0/(1-q)),1.0/a)

def weibullp(p, scale=1.0, shape=0.5):
    a, b, p = map(arr,(shape, scale, p))
    assert all((a>0) & (b>0)), _posstr
    assert all((0<=p)&(p<=1)), _quantstr
    return b*pow(log(1.0/p),1.0/a)

def weibullstats(scale=1.0, shape=0.5, full=0):
    a, b = map(arr,(shape, scale))
    assert all((a>0) & (b>0)), _posstr
    gam = special.gamma
    ia = 1.0/a
    mn = b*gam(1+ia)
    var = b*b*(gam(1+2*ia)-gam(1+ia)**2)
    if not full:
        return mn, var
    den = (gam(1+2*ia)-gam(1+ia)**2)
    g1 = 2*gam(1+ia)**3 - 3*gam(1+ia)*gam(1+2*ia) + gam(1+3*ia)
    g1 /= den**1.5
    g2 = 12*gam(1+ia)**2 * gam(1+2*ia) - 6*gam(1+ia)**4
    g2 += gam(1+4*ia) - 4*gam(1+ia)*gam(1+3*ia) - 3*gam(1+2*ia)**2
    g2 /= den**2.0
    return mn, var, g1, g2
        

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

def binomcdfc(k, n, pr=0.5):
    return special.bdtrc(k,n,pr)

def binomq(q, n, pr=0.5):
    return special.bdtrik(q,n,pr)

def binomp(p, n, pr=0.5):
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

def bernoullicdfc(k, pr=0.5):
    return binomcdfc(k, 1, pr)

def bernoulliq(k, pr=0.5):
    return binomq(k, 1, pr)

def bernoullip(k, pr=0.5):
    return binomp(k, 1, pr)

def bernoullistats(pr=0.5, full=0):
    return binomstats(1, pr, full)

# Negative binomial

def nbinompdf(k, n, pr=0.5):
    return scipy.comb(n+k-1,k)* pr**n * (1-pr)**k

def nbinomcdf(k, n, pr=0.5):
    return special.nbdtr(k,n,pr)

def nbinomcdfc(k, n, pr=0.5):
    return special.nbdtrc(k,n,pr)

def nbinomq(q, n, pr=0.5):
    return special.nbdtrik(q,n,pr)

def nbinomp(p, n, pr=0.5):
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

def geomcdfc(k, pr=0.5):
    return (1.0-pr)**(k+1)

def geomq(q, pr=0.5):
    return log(1.0-q)/log(1-pr)-1

def geomp(p, pr=0.5):
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

def hypergeomcdfc(k, tot=35, good=25, N=10):
    return 1.0 - hypergeomcdf(k, tot, good, N)

def hypergeomq(q, tot=35, good=25, N=10):
    pass

def hypergeomp(p, tot=35, good=25, N=10):
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

def logsercdfc(k, pr=0.5):
    return 1.0-logseriescdf(k, pr=pr)

def logserq(q, pr=0.5):
    pass

def logserp(p, pr=0.5):
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

def poissoncdfc(k,mu):
    k, mu = arr(k), arr(mu)
    sv = errp(0)
    vals = special.pdtrc(k,mu)
    sv = errp(sv)
    return select([mu<0,k>=0],[scipy.nan, vals])

def poissonq(q,mu):
    q, mu = arr(q), arr(mu)
    sv = errp(0)
    vals = special.pdtrik(q,mu)
    sv = errp(sv)
    return where((mu<0) | (q<0) | (q>1),scipy.nan,vals)

def poissonp(p,mu):
    return poissonq(1-p,mu)

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
    
def randintcdfc(k, min, max=None):
    return 1.0-randintcdf(k, min, max)

def randintq(q, min, max=None):
    if max is None:
        max = min
        min = 0
    q = arr(q)
    cond = (arr(min) < arr(max)) & (0<=q) & (q<=1)
    return where(cond, (max-min)*q + min, scipy.nan)

def randintp(p, min, max=None):
    return randintq(1-p, min, max)

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

def zipfcdfc(k, a=4.0):
    return 1.0-zipfcdf(k, a)

def zipfq(q, a=4.0):
    pass

def zipfp(p, a=4.0):
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


