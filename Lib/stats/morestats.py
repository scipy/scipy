# Author:  Travis Oliphant, 2002
#

from __future__ import nested_scopes

import math
import statlib
import stats
import distributions
import inspect
from scipy_base import isscalar, r_, log, sum, around, unique, asarray
from scipy_base import zeros, arange, sort, amin, amax, any, where, \
     array, atleast_1d, sqrt, ceil, floor, array, poly1d, compress, not_equal, \
     pi, exp, ravel
import scipy
import types
import scipy.optimize as optimize
import futil
import scipy_base as sb


def find_repeats(arr):
    """Find repeats in arr and return (repeats, repeat_count)
    """    
    v1,v2, n = futil.dfreps(arr)
    return v1[:n],v2[:n]


##########################################################
###  Bayesian confidence intervals for mean, variance, std
##########################################################

##  Assumes all is known is that mean, and std (variance) exist
##   and are the same for all the data.  Uses Jeffrey's prior
##
##  Returns alpha confidence interval for the mean, variance,
##      and std.

def bayes_mvs(data,alpha=0.90):
    """Return Bayesian confidence intervals for the mean, var, and std.

    Assumes 1-d data all has same mean and variance and uses Jeffrey's prior
    for variance and std.

    alpha gives the probability that the returned interval contains
    the true parameter.

    Uses peak of conditional pdf as starting center.

    Returns (peak, (a, b)) for each of mean, variance and standard deviation.
    """
    x = ravel(data)
    n = len(x)
    assert(n > 1)
    n = float(n)
    xbar = sb.add.reduce(x)/n
    C = sb.add.reduce(x*x)/n - xbar*xbar
    #
    fac = sqrt(C/(n-1))
    tval = distributions.t.ppf((1+alpha)/2.0,n-1)
    delta = fac*tval
    ma = xbar - delta
    mb = xbar + delta
    mp = xbar
    #
    fac = n*C/2.0
    peak = 2/(n+1.)
    a = (n-1)/2.0    
    F_peak = distributions.invgamma.cdf(peak,a)
    q1 = F_peak - alpha/2.0
    q2 = F_peak + alpha/2.0
    if (q1 < 0):  # non-symmetric area
        q2 = alpha
        va = 0.0
    else:
        va = fac*distributions.invgamma.ppf(q1,a)
    vb = fac*distributions.invgamma.ppf(q2,a)
    vp = peak*fac
    #
    fac = sqrt(fac)
    peak = sqrt(2./n)
    F_peak = distributions.gengamma.cdf(peak,a,-2)
    q1 = F_peak - alpha/2.0
    q2 = F_peak + alpha/2.0
    if (q1 < 0):
        q2 = alpha
        sta = 0.0
    else:
        sta = fac*distributions.gengamma.ppf(q1,a,-2)
    stb = fac*distributions.gengamma.ppf(q2,a,-2)
    stp = peak*fac
        
    return (mp,(ma,mb)),(vp,(va,vb)),(stp,(sta,stb))
    


################################
##  K-Statistic
################################

###
## The n-th k-statistic is the unique symmetric unbiased estimator of
##   the n-th cumulant kappa_n
## The cumulants are related to central moments but are specifically defined
##   using a power series expansion of the logarithm of the characteristic
##   function (which is the Fourier transform of the PDF).
##   In particular let phi(t) be the characteristic function, then
##                     _
##         ln phi(t) = > kappa_n (it)^n / n!    (sum from n=0 to inf)
##                     -
##
##   The first few cumulants (kappa_n)  in terms of central moments (mu_n) are
##    kappa_1 = mu_1
##    kappa_2 = mu_2
##    kappa_3 = mu_3
##    kappa_4 = mu_4 - 3*mu_2**2
##    kappa_5 = mu_5 - 10*mu_2 * mu_3
##
## Source:  http://mathworld.wolfram.com/Cumulant.html
##          http://mathworld.wolfram.com/k-Statistic.html
##

def kstat(data,n=2):
    """Return the nth k-statistic (1<=n<=4 so far).

    The nth k-statistic is the unique symmetric unbiased estimator of the nth
    cumulant kappa_n
    """
    if n>4 or n<1:
        raise ValueError, "k-statistics only supported for 1<=n<=4"
    n = int(n)
    S = zeros(n+1,'d')
    data = ravel(data)
    N = len(data)
    for k in range(1,n+1):
        S[k] = sum(data**k)
    if n==1:
        return S[1]*1.0/N
    elif n==2:
        return (N*S[2]-S[1]**2.0)/(N*(N-1.0))
    elif n==3:
        return (2*S[1]**3 - 3*N*S[1]*S[2]+N*N*S[3]) / (N*(N-1.0)*(N-2.0))
    elif n==4:
        return (-6*S[1]**4 + 12*N*S[1]**2 * S[2] - 3*N*(N-1.0)*S[2]**2 - \
                4*N*(N+1)*S[1]*S[3] + N*N*(N+1)*S[4]) / \
                (N*(N-1.0)*(N-2.0)*(N-3.0))
    else:
        raise ValueError, "Should not be here."

def kstatvar(data,n=2):
    """Returns an unbiased estimator of the variance of the k-statistic:  n=1 or 2
    """
    data = ravel(data)
    N = len(data)
    if n==1:
        return kstat(data,n=2)*1.0/N
    elif n==2:
        k2 = kstat(data,n=2)
        k4 = kstat(data,n=4)
        return (2*k2*k2*N + (N-1)*k4)/(N*(N+1))
    else:
        raise ValueError, "Only n=1 or n=2 supported."


#__all__ = ['probplot','ppcc_max','ppcc_plot','boxcox','boxcox_llf',
#           'boxcox_normplot','boxcox_normmax','shapiro']

def probplot(x, sparams=(), dist='norm', fit=1, plot=None):
    """Return (osm, osr){,(scale,loc,r)} where (osm, osr) are order statistic
    medians and ordered response data respectively so that plot(osm, osr)
    is a probability plot.  If fit==1, then do a regression fit and compute the
    slope (scale), intercept (loc), and correlation coefficient (r), of the
    best straight line through the points.  If fit==0, only (osm, osr) is
    returned.
    
    sparams is a tuple of shape parameter arguments for the distribution.
    """
    N = len(x)
    Ui = zeros(N)*1.0
    Ui[-1] = 0.5**(1.0/N)
    Ui[0] = 1-Ui[-1]
    i = arange(2,N)
    Ui[1:-1] = (i-0.3175)/(N+0.365)
    try:
        ppf_func = eval('distributions.%s.ppf'%dist)
    except AttributError:
        raise dist, "is not a valid distribution with a ppf."
    if sparams is None:
        sparams = ()
    if isscalar(sparams):
        sparams = (sparams,)
    if not isinstance(sparams,types.TupleType):
        sparams = tuple(sparams)
    res = inspect.getargspec(ppf_func)
    if not ('loc' == res[0][-2] and 'scale' == res[0][-1] and \
            0.0==res[-1][-2] and 1.0==res[-1][-1]):
        raise ValueError, "Function has does not have default location", \
              "and scale parameters\n  that are 0.0 and 1.0 respectively."
    if (len(sparams) < len(res[0])-len(res[-1])-1) or \
       (len(sparams) > len(res[0])-3):
        raise ValueError, "Incorrect number of shape parameters."
    osm = ppf_func(Ui,*sparams)
    osr = sort(x)
    if fit or (plot is not None):
        # perform a linear fit.
        slope, intercept, r, prob, sterrest = stats.linregress(osm,osr)
    if plot is not None:
        try:
            import scipy.xplt as xplt
            xplt.limits()
        except: pass
        plot.plot(osm, osr, 'o', osm, slope*osm + intercept)
        plot.title('Probability Plot')
        plot.xlabel('Order Statistic Medians')
        plot.ylabel('Ordered Values')
        try: plot.expand_limits(5)
        except: pass
        xmin,xmax= amin(osm),amax(osm)
        ymin,ymax= amin(x),amax(x)
        pos = xmin+0.70*(xmax-xmin), ymin+0.01*(ymax-ymin)
        try: plot.addtext("r^2^=%1.4f" % r, xy=pos,tosys=1)
        except: pass
    if fit:
        return (osm, osr), (slope, intercept, r)
    else:
        return osm, osr

def ppcc_max(x, brack=(0.0,1.0), dist='tukeylambda'):
    """Returns the shape parameter that maximizes the probability plot
    correlation coefficient for the given data to a one-parameter
    family of distributions.

    See also ppcc_plot
    """
    try:
        ppf_func = eval('distributions.%s.ppf'%dist)
    except AttributError:
        raise dist, "is not a valid distribution with a ppf."
    res = inspect.getargspec(ppf_func)
    if not ('loc' == res[0][-2] and 'scale' == res[0][-1] and \
            0.0==res[-1][-2] and 1.0==res[-1][-1]):
        raise ValueError, "Function has does not have default location", \
              "and scale parameters\n  that are 0.0 and 1.0 respectively."
    if (1 < len(res[0])-len(res[-1])-1) or \
       (1 > len(res[0])-3):
        raise ValueError, "Must be a one-parameter family."
    N = len(x)
    # compute uniform median statistics
    Ui = zeros(N)*1.0
    Ui[-1] = 0.5**(1.0/N)
    Ui[0] = 1-Ui[-1]
    i = arange(2,N)
    Ui[1:-1] = (i-0.3175)/(N+0.365)
    osr = sort(x)
    # this function computes the x-axis values of the probability plot
    #  and computes a linear regression (including the correlation)
    #  and returns 1-r so that a minimization function maximizes the
    #  correlation
    def tempfunc(shape, mi, yvals, func):
        xvals = func(mi, shape)
        r, prob = stats.pearsonr(xvals, yvals)
        return 1-r
    return optimize.brent(tempfunc, brack=brack, args=(Ui, osr, ppf_func))

def ppcc_plot(x,a,b,dist='tukeylambda', plot=None, N=80):
    """Returns (shape, ppcc), and optionally plots shape vs. ppcc
    (probability plot correlation coefficient) as a function of shape
    parameter for a one-parameter family of distributions from shape
    value a to b.

    See also ppcc_max
    """
    svals = r_[a:b:complex(N)]
    ppcc = svals*0.0
    k=0
    for sval in svals:
        r1,r2 = probplot(x,sval,dist=dist,fit=1)
        ppcc[k] = r2[-1]
        k += 1
    if plot is not None:
        try:
            import scipy.xplt as xplt
            xplt.limits()
        except: pass
        plot.plot(svals, ppcc, 'x')
        plot.title('(%s) PPCC Plot' % dist)
        plot.xlabel('Prob Plot Corr. Coef.',deltay=-0.01)
        plot.ylabel('Shape Values',deltax=-0.01)
        try: plot.expand_limits(5)
        except: pass
    return svals, ppcc

def boxcox_llf(lmb, data):
    """The boxcox log-likelihood function.
    """
    N = len(data)
    y = boxcox(data,lmb)
    my = stats.mean(y)
    f = (lmb-1)*sum(log(data))
    f -= N/2.0*log(sum((y-my)**2.0/N))
    return f

def _boxcox_conf_interval(x, lmax, alpha):
    # Need to find the lambda for which
    #  f(x,lmbda) >= f(x,lmax) - 0.5*chi^2_alpha;1
    fac = 0.5*distributions.chi2.ppf(1-alpha,1)
    target = boxcox_llf(lmax,x)-fac
    def rootfunc(lmbda,data,target):
        return boxcox_llf(lmbda,data) - target
    # Find positive endpont
    newlm = lmax+0.5
    N = 0
    while (rootfunc(newlm,x,target) > 0.0) and (N < 500):
        newlm += 0.1
        N +=1
    if (N==500):
        raise RuntimeError, "Could not find endpoint."
    lmplus = optimize.brentq(rootfunc,lmax,newlm,args=(x,target))
    newlm = lmax-0.5
    N = 0
    while (rootfunc(newlm,x,target) > 0.0) and (N < 500):
        newlm += 0.1
        N +=1
    if (N==500):
        raise RuntimeError, "Could not find endpoint."
    lmminus = optimize.brentq(rootfunc,newlm,lmax,args=(x,target))
    return lmminus,lmplus
       
def boxcox(x,lmbda=None,alpha=None):
    """Return a positive dataset tranformed by a Box-Cox power transformation.

    If lmbda is not None, do the transformation for that value.

    If lmbda is None, find the lambda that maximizes the log-likelihood
    function and return it as the second output argument.

    If alpha is not None, return the 100(1-alpha)% confidence interval for
    lambda as the third output argument. 
    """
    if any(x < 0):
        raise ValueError, "Data must be positive."
    if lmbda is not None:  # single transformation
        lmbda = lmbda*(x==x)
        y = where(lmbda == 0, log(x), (x**lmbda - 1)/lmbda)
        return y
    # Otherwise find the lmbda that maximizes the log-likelihood function.
    def tempfunc(lmb, data):  # function to minimize
        return -boxcox_llf(lmb,data)
    lmax = optimize.brent(tempfunc, brack=(-2.0,2.0),args=(x,))
    y, lmax = boxcox(x, lmax)
    if alpha is None:
        return y, lmax
    # Otherwise find confidence interval
    interval = _boxcox_conf_interval(x, lmax, alpha)
    return y, lmax, interval

    
def boxcox_normmax(x,brack=(-1.0,1.0)):
    N = len(x)
    # compute uniform median statistics
    Ui = zeros(N)*1.0
    Ui[-1] = 0.5**(1.0/N)
    Ui[0] = 1-Ui[-1]
    i = arange(2,N)
    Ui[1:-1] = (i-0.3175)/(N+0.365)
    # this function computes the x-axis values of the probability plot
    #  and computes a linear regression (including the correlation)
    #  and returns 1-r so that a minimization function maximizes the
    #  correlation
    xvals = distributions.norm.ppf(Ui)
    def tempfunc(lmbda, xvals, samps):
        y = boxcox(samps,lmbda)
        yvals = sort(y)
        r, prob  = stats.pearsonr(xvals, yvals)
        return 1-r
    return optimize.brent(tempfunc, brack=brack, args=(xvals, x))
        

def boxcox_normplot(x,la,lb,plot=None,N=80):
    svals = r_[la:lb:complex(N)]
    ppcc = svals*0.0
    k = 0
    for sval in svals:
        r1,r2 = probplot(x,dist='norm',fit=1)
        ppcc[k] = r2[-1]
        k +=1
    if plot is not None:
        try:
            import scipy.xplt as xplt
            xplt.limits()
        except: pass
        plot.plot(svals, ppcc, 'x')
        plot.title('Box-Cox Normality Plot')
        plot.xlabel('Prob Plot Corr. Coef.',deltay=-0.01)
        plot.ylabel('Transformation parameter',deltax=-0.01)
        try: plot.expand_limits(5)
        except: pass
    return svals, ppcc

def shapiro(x,a=None,reta=0):
    """Shapiro and Wilk test for normality.

    Given random variates x, compute the W statistic and its p-value
    for a normality test.

    If p-value is high, one cannot reject the null hypothesis of normality
    with this test.  P-value is probability that the W statistic is
    as low as it is if the samples are actually from a normal distribution.

    Output:  W statistic and its p-value

              if reta is nonzero then also return the computed "a" values
                 as the third output.  If these are known for a given size
                 they can be given as input instead of computed internally. 
    
    """
    N = len(x)
    if N < 3:
        raise ValueError, "Data must be at least length 3."
    if a is None:
        a = zeros(N,'f')
        init = 0
    else:
        assert(len(a) == N/2), "a must be == len(x)/2"
        init = 1
    y = sort(x)
    a,w,pw,ifault = statlib.swilk(y,a[:N/2],init)
    if not ifault in [0,2]:
        print ifault
    if N > 5000:
        print "p-value may not be accurate for N > 5000."
    if reta:
        return w, pw, a
    else:        
        return w, pw

# Values from Stephens, M A, "EDF Statistics for Goodness of Fit and
#             Some Comparisons", Journal of he American Statistical
#             Association, Vol. 69, Issue 347, Sept. 1974, pp 730-737
_Avals_norm = array([0.576, 0.656, 0.787, 0.918, 1.092])
_Avals_expon  = array([0.922, 1.078, 1.341, 1.606, 1.957])
# From Stephens, M A, "Goodness of Fit for the Extreme Value Distribution",
#             Biometrika, Vol. 64, Issue 3, Dec. 1977, pp 583-588.
_Avals_gumbel = array([0.474, 0.637, 0.757, 0.877, 1.038])
# From Stephens, M A, "Tests of Fit for the Logistic Distribution Based
#             on the Empirical Distribution Function.", Biometrika, 
#             Vol. 66, Issue 3, Dec. 1979, pp 591-595.
_Avals_logistic = array([0.426, 0.563, 0.660, 0.769, 0.906, 1.010])
def anderson(x,dist='norm'):
    """Anderson and Darling test for normal, exponential, or Gumbel
    (Extreme Value Type I) distribution.

    Given samples x, return A2, the Anderson-Darling statistic,
    the significance levels in percentages, and the corresponding
    critical values.

    Critical values provided are for the following significance levels
    norm/expon:   15%, 10%, 5%, 2.5%, 1%
    Gumbel:       25%, 10%, 5%, 2.5%, 1%
    logistic:     25%, 10%, 5%, 2.5%, 1%, 0.5%

    If A2 is larger than these critical values then for that significance
    level, the hypothesis that the data come from a normal (exponential)
    can be rejected.
    """
    if not dist in ['norm','expon','gumbel','extreme1','logistic']:
        raise ValueError, "Invalid distribution."
    y = sort(x)
    xbar = stats.mean(x)
    N = len(y)
    if dist == 'norm':
        s = stats.std(x)
        w = (y-xbar)/s
        z = distributions.norm.cdf(w)
        sig = array([15,10,5,2.5,1])
        critical = around(_Avals_norm / (1.0 + 4.0/N - 25.0/N/N),3)
    elif dist == 'expon':
        w = y / xbar
        z = distributions.expon.cdf(w)
        sig = array([15,10,5,2.5,1])
        critical = around(_Avals_expon / (1.0 + 0.6/N),3)
    elif dist == 'logistic':
        def rootfunc(ab,xj,N):
            a,b = ab
            tmp = (xj-a)/b
            tmp2 = exp(tmp)
            val = [sum(1.0/(1+tmp2))-0.5*N,
                   sum(tmp*(1.0-tmp2)/(1+tmp2))+N]
            return array(val)
        sol0=array([xbar,stats.std(x)])
        sol = optimize.fsolve(rootfunc,sol0,args=(x,N),xtol=1e-5)
        w = (y-sol[0])/sol[1]
        z = distributions.logistic.cdf(w)
        sig = array([25,10,5,2.5,1,0.5])
        critical = around(_Avals_logistic / (1.0+0.25/N),3)
    else:
        def fixedsolve(th,xj,N):
            val = stats.sum(xj)*1.0/N
            tmp = exp(-xj/th)
            term = sum(xj*tmp)
            term /= sum(tmp)
            return val - term
        s = optimize.fixed_point(fixedsolve, 1.0, args=(x,N),xtol=1e-5)
        xbar = -s*log(sum(exp(-x/s))*1.0/N)
        w = (y-xbar)/s
        z = distributions.gumbel_l.cdf(w)
        sig = array([25,10,5,2.5,1])
        critical = around(_Avals_gumbel / (1.0 + 0.2/sqrt(N)),3)
    i = arange(1,N+1)
    S = sum((2*i-1.0)/N*(log(z)+log(1-z[::-1])))
    A2 = -N-S
    return A2, critical, sig


def _find_repeats(arr):
    """Find repeats in the array and return (repeats, repeat_count)
    """
    arr = sort(arr)
    lastval = arr[0]
    howmany = 0
    ind = 1
    N = len(arr)
    repeat = 0
    replist = []
    repnum = []
    while ind < N:
        if arr[ind] != lastval:
            if repeat:
                repnum.append(howmany+1)
                repeat = 0
                howmany = 0
        else:
            howmany += 1
            repeat = 1
            if (howmany == 1):
                replist.append(arr[ind])
        lastval = arr[ind]
        ind += 1
    if repeat:
        repnum.append(howmany+1)
    return replist, repnum

def ansari(x,y):
    """Determine if the scale parameter for two distributions with equal
    medians is the same using the Ansari-Bradley statistic.

    Specifically, compute the AB statistic and the probability of error
    that the null hypothesis is true but rejected with the computed
    statistic as the critical value.

    One can reject the null hypothesis that the ratio of variances is 1 if
    returned probability of error is small (say < 0.05)
    """
    x,y = asarray(x),asarray(y)
    n = len(x)
    m = len(y)
    if (m < 1):
        raise ValueError, "Not enough other observations."
    if (n < 1):
        raise ValueError, "Not enough test observations."
    N = m+n
    xy = r_[x,y]  # combine
    rank = stats.rankdata(xy)
    symrank = amin(array((rank,N-rank+1)),0)
    AB = sum(symrank[:n])
    uxy = unique(xy)
    repeats = (len(uxy) != len(xy))
    exact = ((m<55) and (n<55) and not repeats)
    if repeats and ((m < 55)  or (n < 55)):
        print "Ties preclude use of exact statistic."
    if exact:
        astart, a1, ifault = statlib.gscale(n,m)
        ind = AB-astart
        total = sum(a1)
        if ind < len(a1)/2.0:
            cind = int(ceil(ind))
            if (ind == cind):
                pval = 2.0*sum(a1[:cind+1])/total
            else:
                pval = 2.0*sum(a1[:cind])/total
        else:
            find = int(floor(ind))
            if (ind == floor(ind)):
                pval = 2.0*sum(a1[find:])/total
            else:
                pval = 2.0*sum(a1[find+1:])/total
        return AB, min(1.0,pval)
    
    # otherwise compute normal approximation
    if N % 2:  # N odd
        mnAB = n*(N+1.0)**2 / 4.0 / N
        varAB = n*m*(N+1.0)*(3+N**2)/(48.0*N**2)
    else:
        mnAB = n*(N+2.0)/4.0
        varAB = m*n*(N+2)*(N-2.0)/48/(N-1.0)
    if repeats:   # adjust variance estimates
        # compute sum(tj * rj**2)
        fac = sum(symrank**2)
        if N % 2: # N odd
            varAB = m*n*(16*N*fac-(N+1)**4)/(16.0 * N**2 * (N-1))
        else:  # N even
            varAB = m*n*(16*fac-N*(N+2)**2)/(16.0 * N * (N-1))
    z = (AB - mnAB)/sqrt(varAB)
    pval = (1-distributions.norm.cdf(abs(z)))*2.0
    return AB, pval

def bartlett(*args):
    """Perform Bartlett test with the null hypothesis that all input samples
    have equal variances.

    Inputs are sample vectors:  bartlett(x,y,z,...) 

    Outputs: (T, pval)

         T    -- the Test statistic
         pval -- significance level if null is rejected with this value of T
                 (prob. that null is true but rejected with this p-value.)  

    Sensitive to departures from normality.  The Levene test is
    an alternative that is less sensitive to departures from
    normality.

    References:

      http://www.itl.nist.gov/div898/handbook/eda/section3/eda357.htm

      Snedecor, George W. and Cochran, William G. (1989), Statistical
        Methods, Eighth Edition, Iowa State University Press.     
    """
    k = len(args)
    if k < 2:
        raise ValueError, "Must enter at least two input sample vectors."
    Ni = zeros(k)
    ssq = zeros(k,'d')
    for j in range(k):
        Ni[j] = len(args[j])
        ssq[j] = stats.var(args[j])
    Ntot = sum(Ni)
    spsq = sum((Ni-1)*ssq)/(1.0*(Ntot-k))
    numer = (Ntot*1.0-k)*log(spsq) - sum((Ni-1.0)*log(ssq))
    denom = 1.0 + (1.0/(3*(k-1)))*((sum(1.0/(Ni-1.0)))-1.0/(Ntot-k))
    T = numer / denom
    pval = distributions.chi2.sf(T,k-1) # 1 - cdf
    return T, pval


def levene(*args,**kwds):
    """Perform Levene test with the null hypothesis that all input samples
    have equal variances.

    Inputs are sample vectors:  bartlett(x,y,z,...)

    One keyword input, center, can be used with values
        center = 'mean', center='median' (default), center='trimmed'

    center='median' is recommended for skewed (non-normal) distributions
    center='mean' is recommended for symmetric, moderate-tailed, dist.
    center='trimmed' is recommended for heavy-tailed distributions.

    Outputs: (W, pval)

         W    -- the Test statistic
         pval -- significance level if null is rejected with this value of W
                 (prob. that null is true but rejected with this p-value.)  

    References:

       http://www.itl.nist.gov/div898/handbook/eda/section3/eda35a.htm

       Levene, H. (1960). In Contributions to Probability and Statistics:
         Essays in Honor of Harold Hotelling, I. Olkin et al. eds.,
         Stanford University Press, pp. 278-292. 
       Brown, M. B. and Forsythe, A. B. (1974), Journal of the American
         Statistical Association, 69, 364-367         
    """
    k = len(args)
    if k < 2:
        raise ValueError, "Must enter at least two input sample vectors."
    Ni = zeros(k)
    Yci = zeros(k,'d')
    if 'center' in kwds.keys():
        center = kwds['center']
    else:
        center = 'median'
    if not center in ['mean','median','trimmed']:
        raise ValueError, "Keyword argument <center> must be 'mean', 'median'"\
              + "or 'trimmed'."
    if center == 'median':
        func = stats.median
    elif center == 'mean':
        func = stats.mean
    else:
        func = stats.trim_mean
    for j in range(k):
        Ni[j] = len(args[j])
        Yci[j] = func(args[j])    
    Ntot = sum(Ni)

    # compute Zij's
    Zij = [None]*k
    for i in range(k):
        Zij[i] = abs(asarray(args[i])-Yci[i])
    # compute Zbari
    Zbari = zeros(k,'d')
    Zbar = 0.0
    for i in range(k):
        Zbari[i] = stats.mean(Zij[i])
        Zbar += Zbari[i]*Ni[i]
    Zbar /= Ntot

    numer = (Ntot-k)*sum(Ni*(Zbari-Zbar)**2)

    # compute denom_variance
    dvar = 0.0
    for i in range(k):
        dvar += sum((Zij[i]-Zbari[i])**2)

    denom = (k-1.0)*dvar

    W = numer / denom
    pval = distributions.f.sf(W,k-1,Ntot-k) # 1 - cdf
    return W, pval

def binom_test(x,n=None,p=0.5):
    """An exact (two-sided) test of the null hypothesis that the
    probability of success in a Bernoulli experiment is p.

    Inputs:

       x -- Number of successes (or a vector of length 2 giving the
              number of successes and number of failures respectively)
       n -- Number of trials (ignored if x has length 2)
       p -- Hypothesized probability of success

    Returns pval -- Probability that null test is rejected for this set
                    of x and n even though it is true.
    """
    x = atleast_1d(x)
    if len(x) == 2:
        n = x[1]+x[0]
        x = x[0]
    elif len(x) == 1:
        x = x[0]
        if n is None or n < x:
            raise ValueError, "n must be >= x"
    else:
        raise ValueError, "Incorrect length for x."

    if (p > 1.0) or (p < 0.0):
        raise ValueError, "p must be in range [0,1]"

    d = distributions.binom.pmf(x,n,p)
    rerr = 1+1e-7
    if (x*1.0/n < p):
        i = arange(x+1,n+1)
        y = sum(distributions.binom.pmf(i,n,p) <= d*rerr)
        pval = distributions.binom.cdf(x,n,p) + distributions.binom.sf(n-y,n,p)
    else:
        i = arange(0,x)
        y = sum(distributions.binom.pmf(i,n,p) <= d*rerr)
        pval = distributions.binom.cdf(y-1,n,p) + distributions.binom.sf(x-1,n,p)

    return min(1.0,pval)

def _apply_func(x,g,func):
    # g is list of indices into x
    #  separating x into different groups
    #  func should be applied over the groups
    g = unique(r_[0,g,len(x)])
    output = []
    for k in range(len(g)-1):
        output.append(func(x[g[k]:g[k+1]]))
    return asarray(output)

def fligner(*args,**kwds):
    """Perform Levene test with the null hypothesis that all input samples
    have equal variances.

    Inputs are sample vectors:  bartlett(x,y,z,...)

    One keyword input, center, can be used with values
        center = 'mean', center='median' (default), center='trimmed'

    Outputs: (Xsq, pval)

         Xsq  -- the Test statistic
         pval -- significance level if null is rejected with this value of X
                 (prob. that null is true but rejected with this p-value.)  

    References:

       http://www.stat.psu.edu/~bgl/center/tr/TR993.ps

       Fligner, M.A. and Killeen, T.J. (1976). Distribution-free two-sample
       tests for scale. 'Journal of the American Statistical Association.'
       71(353), 210-213.
    """
    k = len(args)
    if k < 2:
        raise ValueError, "Must enter at least two input sample vectors."
    if 'center' in kwds.keys():
        center = kwds['center']
    else:
        center = 'median'
    if not center in ['mean','median','trimmed']:
        raise ValueError, "Keyword argument <center> must be 'mean', 'median'"\
              + "or 'trimmed'."
    if center == 'median':
        func = stats.median
    elif center == 'mean':
        func = stats.mean
    else:
        func = stats.trim_mean

    Ni = asarray([len(args[j]) for j in range(k)])
    Yci = asarray([func(args[j]) for j in range(k)])
    Ntot = sum(Ni)
    # compute Zij's
    Zij = [abs(asarray(args[i])-Yci[i]) for i in range(k)]
    allZij = []
    g = [0]
    for i in range(k):
        allZij.extend(list(Zij[i]))
        g.append(len(allZij))

    a = distributions.norm.ppf(stats.rankdata(allZij)/(2*(Ntot+1.0)) + 0.5)

    # compute Aibar
    Aibar = _apply_func(a,g,sum) / Ni
    anbar = stats.mean(a)
    varsq = stats.var(a)

    Xsq = sum(Ni*(asarray(Aibar)-anbar)**2.0)/varsq

    pval = distributions.chi2.sf(Xsq,k-1) # 1 - cdf
    return Xsq, pval


def mood(x,y):
    """Determine if the scale parameter for two distributions with equal
    medians is the same using a Mood test.

    Specifically, compute the z statistic and the probability of error
    that the null hypothesis is true but rejected with the computed
    statistic as the critical value.

    One can reject the null hypothesis that the ratio of scale parameters is
    1 if the returned probability of error is small (say < 0.05)
    """
    n = len(x)
    m = len(y)
    xy = r_[x,y]
    N = m+n
    if (N < 3):
        raise ValueError, "Not enough observations."
    ranks = stats.rankdata(xy)
    Ri = ranks[:n]
    M = sum((Ri - (N+1.0)/2)**2)
    # Approx stat.
    mnM = n*(N*N-1.0)/12
    varM = m*n*(N+1.0)*(N+2)*(N-2)/180
    z = (M-mnM)/sqrt(varM)
    p = distributions.norm.cdf(z)
    pval = 2*min(p,1-p)
    return z, pval


def oneway(*args,**kwds):
    """Test for equal means in two or more samples from the
    normal distribution.

    If the keyword parameter <equal_var> is true then the variances
    are assumed to be equal, otherwise they are not assumed to
    be equal (default).

    Return test statistic and the p-value giving the probability
    of error if the null hypothesis (equal means) is rejected at this value.
    """
    k = len(args)
    if k < 2:
        raise ValueError, "Must enter at least two input sample vectors."
    if 'equal_var' in kwds.keys():
        if kwds['equal_var']: evar = 1
        else: evar = 0
    else:
        evar = 0

    Ni = array([len(args[i]) for i in range(k)])
    Mi = array([stats.mean(args[i]) for i in range(k)])
    Vi = array([stats.var(args[i]) for i in range(k)])
    Wi = Ni / Vi
    swi = sum(Wi)
    N = sum(Ni)
    my = sum(Mi*Ni)*1.0/N
    tmp = sum((1-Wi/swi)**2 / (Ni-1.0))/(k*k-1.0)
    if evar:
        F = ((sum(Ni*(Mi-my)**2) / (k-1.0)) / (sum((Ni-1.0)*Vi) / (N-k)))
        pval = distributions.f.sf(F,k-1,N-k)  # 1-cdf
    else:
        m = sum(Wi*Mi)*1.0/swi
        F = sum(Wi*(Mi-m)**2) / ((k-1.0)*(1+2*(k-2)*tmp))
        pval = distributions.f.sf(F,k-1.0,1.0/(3*tmp))

    return F, pval



def wilcoxon(x,y=None):
    """
Calculates the Wilcoxon signed-rank test for the null hypothesis that two samples come from the same distribution. A non-parametric T-test. (need N > 20)

Returns: t-statistic, two-tailed p-value
"""
    if y is None:
        d = x
    else:
        x, y = map(asarray, (x, y))
        if len(x) <> len(y):
            raise ValueError, 'Unequal N in wilcoxon.  Aborting.'
        d = x-y
    d = compress(not_equal(d,0),d) # Keep all non-zero differences
    count = len(d)
    if (count < 10):
        print "Warning: sample size too small for normal approximation."
    r = stats.rankdata(abs(d))
    r_plus = sum((d > 0)*r)
    r_minus = sum((d < 0)*r)
    T = min(r_plus, r_minus)
    mn = count*(count+1.0)*0.25
    se = math.sqrt(count*(count+1)*(2*count+1.0)/24)
    if (len(r) != len(unique(r))):  # handle ties in data
        replist, repnum = find_repeats(r)
        corr = 0.0
        for i in range(len(replist)):
            si = repnum[i]
            corr += 0.5*si*(si*si-1.0)
        V = se*se - corr
        se = sqrt((count*V - T*T)/(count-1.0))
    z = (T - mn)/se
    prob = 2*(1.0 -stats.zprob(abs(z)))
    return T, prob

def _hermnorm(N):
    # return the negatively normalized hermite polynomials up to order N-1
    #  (inclusive)
    #  using the recursive relationship
    #  p_n+1 = p_n(x)' - x*p_n(x)
    #   and p_0(x) = 1
    plist = [None]*N
    plist[0] = poly1d(1)
    for n in range(1,N):
        plist[n] = plist[n-1].deriv() - poly1d([1,0])*plist[n-1]
    return plist

def pdf_moments(cnt):
    """Return the Gaussian expanded pdf function given the list of central
    moments (first one is mean).
    """
    N = len(cnt)
    if N < 2:
        raise ValueError, "At least two moments must be given to" + \
              "approximate the pdf."
    totp = poly1d(1)
    sig = sqrt(cnt[1])
    mu = cnt[0]
    if N > 2:
        Dvals = _hermnorm(N+1)
    for k in range(3,N+1):
        # Find Ck
        Ck = 0.0
        for n in range((k-3)/2):
            m = k-2*n
            if m % 2: # m is odd
                momdiff = cnt[m-1]
            else:
                momdiff = cnt[m-1] - sig*sig*scipy.factorial2(m-1)
            Ck += Dvals[k][m] / sig**m * momdiff
        # Add to totp 
        totp = totp +  Ck*Dvals[k]        
        
    def thisfunc(x):
        xn = (x-mu)/sig
        return totp(xn)*exp(-xn*xn/2.0)/sqrt(2*pi)/sig
    return thisfunc
           

def pdf_fromgamma(g1,g2,g3=0.0,g4=None):
    if g4 is None:
        g4 = 3*g2*g2
    sigsq = 1.0/g2
    sig = sqrt(sigsq)
    mu = g1*sig**3.0
    p12 = _hermnorm(13)
    for k in range(13):
        p12[k] = p12[k]/sig**k

    # Add all of the terms to polynomial
    totp = p12[0] - (g1/6.0*p12[3]) + \
           (g2/24.0*p12[4] +g1*g1/72.0*p12[6]) - \
           (g3/120.0*p12[5] + g1*g2/144.0*p12[7] + g1**3.0/1296.0*p12[9]) + \
           (g4/720*p12[6] + (g2*g2/1152.0+g1*g3/720)*p12[8] +
            g1*g1*g2/1728.0*p12[10] + g1**4.0/31104.0*p12[12])
    # Final normalization
    totp = totp / sqrt(2*pi)/sig
    def thefunc(x):
        xn = (x-mu)/sig
        return totp(xn)*exp(-xn*xn/2.0)
    return thefunc

def pdfapprox(samples):
    """Return a function that approximates the pdf of a set of samples
    using a Gaussian expansion computed from the mean, variance, skewness
    and Fisher's kurtosis.
    """
    # Estimate mean, variance, skewness and kurtosis
    mu,sig,sk,kur = stats.describe(samples)[2:]
    # Get central moments
    cnt = [None]*4
    cnt[0] = mu
    cnt[1] = sig*sig
    cnt[2] = sk * sig**1.5
    cnt[3] = (kur+3.0) * sig**2.0
    return pdf_moments(cnt)
    #g2 = (1.0/sig)**2.0
    #g1 = mu / sig**3.0
    #g3 = sk / sig**3.5
    #g4 = (kur+3.0) / sig**4.0
    #return pdf_fromgamma(g1, g2, g3, g4)

def circmean(samples, high=2*pi, low=0):
    """Compute the circular mean for samples assumed to be in the range [low to high]
    """
    ang = (samples - low)*2*pi / (high-low)
    res = angle(stats.mean(exp(1j*ang)))
    if (res < 0):
        res = res + 2*pi
    return res*(high-low)/2.0/pi + low

def circvar(samples, high=2*pi, low=0):
    """Compute the circular variance for samples assumed to be in the range [low to high]
    """
    ang = (samples - low)*2*pi / (high-low)
    res = stats.mean(exp(1j*ang))
    V = 1-abs(res)
    return ((high-low)/2.0/pi)**2 * V

def circstd(samples, high=2*pi, low=0):
    """Compute the circular standard deviation for samples assumed to be in the range [low to high]
    """
    ang = (samples - low)*2*pi / (high-low)
    res = stats.mean(exp(1j*ang))
    V = 1-abs(res)
    return ((high-low)/2.0/pi) * sqrt(V)
    
    
        
#Tests to include (from R) -- some of these already in stats.
########
#X Ansari-Bradley
#X Bartlett (and Levene)
#X Binomial
#Y Pearson's Chi-squared (stats.chisquare)
#Y Association Between Paired samples (stats.pearsonr, stats.spearmanr)
#                       stats.kendalltau) -- these need work though
# Fisher's exact test
#X Fligner-Killeen Test
#Y Friedman Rank Sum (stats.friedmanchisquare?)
#Y Kruskal-Wallis
#Y Kolmogorov-Smirnov
# Cochran-Mantel-Haenszel Chi-Squared for Count
# McNemar's Chi-squared for Count
#X Mood Two-Sample
#X Test For Equal Means in One-Way Layout (see stats.ttest also)
# Pairwise Comparisons of proportions
# Pairwise t tests
# Tabulate p values for pairwise comparisons
# Pairwise Wilcoxon rank sum tests
# Power calculations two sample test of prop.
# Power calculations for one and two sample t tests
# Equal or Given Proportions
# Trend in Proportions
# Quade Test
#Y Student's T Test
#Y F Test to compare two variances
#XY Wilcoxon Rank Sum and Signed Rank Tests
