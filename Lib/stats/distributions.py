# Functions to implement CDF, PDF, Quanitles, CDFC, *P-iles, and STATS
#   for various Continous and Discrete Probability Distributions
#
# Author:  Travis Oliphant  2002
# 

import scipy
import scipy.special as special
from fastumath import *

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
##  <dist>stats --- Return mean, variance and optionally skew and kurtosis
##                      of the distribution.


## Kolmogorov-Smirnov one-sided and two-sided test statistics

def ksonecdfc(x,n):
    return special.smirnov(n,x)

def ksonep(p,n):
    return special.smirnovi(n,p)

def ksonecdf(x,n):
    return 1-special.smirnov(n,x)

def ksoneq(q,n):
    assert all(0<=q<=1), _quantstr
    return special.smirnovi(n,1-q)

def kstwocdfc_largen(y):
    return special.kolmogorov(y)

def kstwop_largen(p):
    assert all(0<=p<=1), _quanstr
    return special.kolmogi(p)

def kstwocdf_largen(y):
    assert(y>=0)
    return 1-special.kolmogorov(y)

def kstwoq_largen(q):
    assert(all(0<=q<=1))
    return special.kolmogi(1-q)

## Normal distributions

def stnormpdf(x):
    return 1.0/sqrt(2*pi)*exp(-x**2/2.0)

def stnormcdfc(x):
    return 1-special.ndtr(x)

def stnormq(q):
    assert all(0<=q<=1), _quantstr
    return special.ndtri(q)

def stnormp(p):
    assert all(0<=p<=1), _quantstr
    return special.ndtri(1-p)

def stnormcdf(x):
    return special.ndtr(x)

def stnromalq(x):
    assert(all(0<=q<=1))
    return special.ndtri(q)

def stnormstats(full=0):
    if not full:
        return 0, 1
    return 0, 1, 0, 0

def normpdf(x, mu=0.0, std=1.0):
    return 1.0/sqrt(2*pi)/std * exp(-(x-mu)**2 / (2.0*std*std))

def normcdf(x, mu=0.0, std=1.0):
    return special.ndtr((x-mu)*1.0/std)

def normq(q, mu=0.0, std=1.0):
    return special.ndtri(q)*std + mu

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
    Px = (1.0-x)**(b-1.0) * x**(a-1.0)
    Px /= special.beta(a,b)
    choice1 = where(x>1,0,Px)
    return where(x<0,0,choice1)
    
def betacdf(x, a, b):
    x = where(x<0,0,x)
    x = where(x>1,1,x)
    return special.btdtr(a,b,x)

def betaq(q, a, b):
    assert all(0<=q<=1), _quantstr
    return special.btdtri(a,b,q)

def betacdfc(x, a, b):
    return 1-betacdf(x,a,b)

def betap(p, a, b):
    return betaq(1-p,a,b)

def betastats(a, b, full=0):
    assert all((a>0) & (b>0)), _posstr
    mn = a*1.0 / (a+b)
    var = (a*b)*1.0 / ((a+b)**2 * (a+b+1))
    if not full:
        return mn, var
    g1 = 2.0*(b-a)*sqrt(1+a+b) / (sqrt(a*b)*(2+a+b))
    g2 = 6.0*(a**3 + a**2*(1-2*b) + b**2*(1+b) - 2*a*b*(2+b))
    g2 /= a*b*(a+b+2)*(a+b+3)
    return mn, var, g1, g2

## Cauchy

def cauchypdf(x, median=0.0, scale=1.0):
    assert all(scale>0), _posstr
    Px = 1.0/pi*scale / ((x-median)**2 + scale**2)
    return Px

def cauchycdf(x, median=0.0, scale=1.0):
    assert all(scale>0), _posstr
    return 0.5 + 1.0/pi*arctan((x-median)*1.0/scale)

def cauchycdfc(x, median=0.0, scale=1.0):
    assert all(scale>0), _posstr
    return 0.5 - 1.0/pi*arctan((x-median)*1.0/scale)

def cauchyq(q, median=0.0, scale=1.0):
    assert all(0<=q<=1), _quantstr
    assert all(scale>0), _posstr
    return scale*tan(pi*(q-0.5))+median

def cauchyp(p, median=0.0, scale=1.0):
    return cauchyq(1-p, median, scale)
    
def cauchystats(median=0.0, scale=1.0):
    raise ValueError, "No moments exist."

## Chi-squared

def chi2pdf(x, df):
    Px = x**(df/2.0-1)*exp(-x/2.0)
    Px /= special.gamma(df/2.0)* 2**(df/2.0)
    return where(x>=0,Px,0)

def chi2cdf(x, df):
    x = where(x<0,0,x)
    return special.chdtr(df, x)

def chi2cdfc(x, df):
    x = where(x<0,0,x)
    return special.chdtrc(df, x)

def chi2p(p, df):
    assert all(0<=p<=1), _quantstr
    return special.chdtri(df, p)

def chi2q(q, df):
    assert all(0<=q<=1), _quantstr
    return special.chdtri(df, 1-q)

def chistats(df, full=0):
    mn = df
    var = 2*df
    if not full:
        return mn, var
    g1 = 2*sqrt(2.0/r)
    g2 = 12.0/r
    return mn, var, g1, g2


## Exponential

def exponpdf(x, lam):
    return where(x<0, 0, lam*exp(-lam*x))

def exponcdf(x, lam):
    x = where(x<0, 0, x)
    return 1.0-exp(-lam*x)

def exponcdfc(x, lam):
    return 1.0-exponcdf(x,lam)

def exponq(q, lam):
    assert all(0<=q<=1), _quantstr
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
    assert all(0<=q<=1), _quantstr
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
    m = 1.0*dfd
    n = 1.0*dfn
    Px = m**(m/2) * n**(n/2) * x**(n/2-1)
    Px /= (m+n*x)**((n+m)/2)*special.beta(n/2,m/2)
    return where(x>0,Px,0)

def fcdf(x, dfn, dfd):
    x = where(x<0, 0, x)
    return special.fdtr(dfn, dfd, x)

def fcdfc(x, dfn, dfd):
    x = where(x<0, 0, x)
    return special.fdtrc(dfn, dfd, x)

def fp(p, dfn, dfd):
    assert all(0<=p<=1), _quantstr
    return special.fdtri(dfn, dfd, p)
    
def fq(q, dfn, dfd):
    assert all(0<=q<=1), _quantstr
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
    x = asarray(x)
    Px = x**(a-1.0)*exp(-x*1.0/b)
    Px /= special.gamma(a) * b**a
    return where(x<0,0,Px)

def gammacdf(x, a, b):
    x = where(x<0,0,x)
    return special.gdtr(1.0/b,a,x)

def gammacdfc(x, a, b):
    x = where(x<0,0,x)
    return special.gdtrc(1.0/b,a,x)
    
def gammaq(q, a, b):
    assert all(0<=q<=1), _quantstr
    return special.gdtri(1.0/b,a,q)

def gammap(p, a, b):
    assert all(0<=p<=1), _quantstr
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
    assert all(0<=q<=1), _quantstr
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
    q = asarray(q)
    assert all(0<=q<=1), _quantstr    
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
    x = asarray(x)
    Px = exp(-(log(x)-mu)**2 / (2*std**2))
    Px /= std*x*sqrt(2*pi)
    return where(x<=0,0,Px)

def lognormcdf(x, mu=0.0, std=1.0):
    x = asarray(x)
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
    y = where(x<=0,1,x)
    a = df/2.0
    Px = 0.5**a*exp(-(nc+y)/2.0)*(y*1.0)**(a-1)
    z = nc*y/4.0
    Px *= special.hyp0f1(a,z)/special.gamma(a)
    return where(x<=0,0,Px)

def _ncx2cdf(x,df,nc):
    from scipy.limits import double_epsilon as eps    
    nc = nc/2.0
    jmid = floor(nc)
    val = 0
    for j in range(jmid,-1,-1):
        term = poissonpdf(j,nc)*chi2cdf(x,df+2*j)
        val += term
        if (term < val*eps):
            print "Here", j
            break
    for j in range(jmid+1,2*jmid+2):
        term = poissonpdf(j,nc)*chi2cdf(x,df+2*j)
        val += term
        if (term < val*eps):
            break
    if (j == 2*jmid+1):
        for j in range(2*(jmid+1),2*(jmid+1)+1000):
            term = poissonpdf(j,nc)*chi2cdf(x,df+2*j)
            val += term
            if (term < val*eps):
                print "final break at ", j
                break
        if (j==(2*jmid+1002)):
            print "Warning:  Did not converge."
    return val
        
#def _ncx2cdf(x,df,nc):
#    import scipy.integrate as integrate
#    return integrate.quad(ncx2pdf,0,x,args=(df,nc))[0]
_vec_ncx2cdf = special.general_function(_ncx2cdf,'d')

def ncx2cdf(x,df,nc):
    assert all(nc>=0), _nonnegstr
    assert all(df>0), _posstr
    x = where(x<0,0,x)
    return _vec_ncx2cdf(x,df,nc)

def ncx2cdfc(x,df,nc):
    return 1-ncx2cdf(x,df,nc)
    
def _ncx2qfunc(x,q,df,nc):
    return _ncx2cdf(x,df,nc)-q

def _ncx2q(q,df,nc):
    import scipy.optimize as optimize
    return optimize.fsolve(_ncx2qfunc,nc+df,args=(q,df,nc))
_vec_ncx2q = special.general_function(_ncx2q,'d')

def ncx2q(q,df,nc):
    assert all(0<=q<=1), _quantstr
    return _vec_ncx2q(q, df, nc)
    
def ncx2p(p,df,nc):
    assert all(0<=p<=1), _quantstr
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

def _ncfcdf(x,n1,n2,nc):
    from scipy.limits import double_epsilon as eps
    eps2 = sqrt(eps)
    eps4 = sqrt(eps2)
    n1, n2, nc = n1/2.0, n2/2.0, nc/2.0
    val = 0
    j = 0
    valm1 = 100
    term = 80
    gam = special.gamma
    bI = special.betainc
    bval =  n1*x/(n2+n1*x)
    jmid = floor(nc)
    for j in range(jmid,-1,-1):
        term = poissonpdf(j,nc)*betacdf(bval,j+n1,n2)
        val += term
        if all(ravel(term / (val+eps4)) < eps2):
            break
    for j in range(jmid+1,jmid+2000):
        term = poissonpdf(j,nc)*betacdf(bval,j+n1,n2)
        val += term
        if all(ravel(term / (val+eps4)) < eps2):
            break
    if (j == jmid+2000-1):
        print "Warning: Series failed to converge."
    return val
_vec_ncfcdf = special.general_function(_ncfcdf,'d')

def ncfcdf(x,dfn,dfd,nc):
    assert all(nc>=0), _nonnegstr
    assert all((dfn>0) & (dfd>0)), _posstr
    x = where(x<0,0,x)
    return _vec_ncfcdf(x,dfn,dfd,nc)

def ncfcdfc(x,dfn,dfd,nc):
    return 1-ncfcdf(x,dfn,dfd,nc)

def _ncfqfunc(x,q,dfn,dfd,nc):
    return _ncfcdf(x,dfn,dfd,nc)-q

def _ncfq(q,dfn,dfd,nc,x0):
    import scipy.optimize as optimize
    return optimize.fsolve(_ncfqfunc,x0,args=(q,dfn,dfd,nc))
_vec_ncfq = special.general_function(_ncfq,'d')

def ncfq(q,dfn,dfd,nc,x0=None):
    assert all(0<=q<=1), _quanstr
    if x0 is None:
        x0 = dfd * (dfn+nc)/(dfn*(dfd-2))
    return _vec_ncfq(q, dfn, dfd, nc, x0)

def ncfp(p,dfn,dfd,nc):
    return ncfq(1-p,dfn,dfd,nc)

def ncfstats(dfn,dfd,nc, full=0):
    dfn = dfn*1.0
    dfd = dfd*1.0
    nc = nc*1.0
    mn = where(dfd<=2,scipy.nan,dfd/(dfd-2)*(1+nc/dfn))
    var1 = 2*(dfd/dfn)**2 * ((dfn+nc/2)**2+(dfn+nc)*(dfd-2))
    var1 /= (dfd-2)**2 * (dfd-4)
    var = where(dfd<=4,scipy.nan,var1)               
    if full:
        print "Skew and kurtosis unavailable for this distribution."
    return mn, var

## Non-central T distribution

def nctpdf(x, df, nc):
    pass

def nctcdf(x,df,nc):
    pass

def nctcdfc(x,df,nc):
    return 1-nctcdf(x,df,nc)

def _nctqfunc(x,q,df,nc):
    return _nctcdf(x,dfn,dfd,nc)-q

def _ncft(q,df,nc,x0):
    import scipy.optimize as optimize
    return optimize.fsolve(_nctqfunc,x0,args=(q,dfn,dfd,nc))
_vec_nctq = special.general_function(_nctq,'d')

def nctq(q,df,nc,x0=None):
    assert all(0<=q<=1), _quanstr
    if x0 is None:
        val1 = gam((df-1.0)/2.0)
        val2 = gam(df/2.0)
        x0 = nc*sqrt(df/2.0)*val1/val2
    return _vec_ncfq(q, dfn, dfd, nc, x0)

def nctp(p,df,nc):
    return nctq(1-p,df,nc)

def nctstats(df,nc,full=0):
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
    

def poissonpdf(k,mu):
    return mu**k * exp(-mu) / special.gamma(k+1)

