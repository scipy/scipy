# Author:  Travis Oliphant, 2002
#

import statlib
import stats
import distributions
import inspect
from scipy_base import isscalar, r_, log, sum, around
from scipy_base import zeros, arange, sort, amin, amax, any, where, array
import types
import scipy.optimize as optimize

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
        ppf_func = eval('distributions.%sppf'%dist)
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
        try: xplt.limits()
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
        ppf_func = eval('distributions.%sppf'%dist)
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
        try: xplt.limits()
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
    fac = 0.5*distributions.chi2ppf(1-alpha,1)
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
    y = boxcox(x, lmax)
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
    xvals = distributions.normppf(Ui)
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
        try: xplt.limits()
        except: pass
        plot.plot(svals, ppcc, 'x')
        plot.title('(%s) Box-Cox Normality Plot' % dist)
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
        z = distributions.normcdf(w)
        sig = array([15,10,5,2.5,1])
        critical = around(_Avals_norm / (1.0 + 4.0/N - 25.0/N/N),3)
    elif dist == 'expon':
        w = y / xbar
        z = distributions.exponcdf(w)
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
        z = distributions.logisticcdf(w)
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
        z = distributions.gumbelcdf(w)
        sig = array([25,10,5,2.5,1])
        critical = around(_Avals_gumbel / (1.0 + 0.2/sqrt(N)),3)
    i = arange(1,N+1)
    S = sum((2*i-1.0)/N*(log(z)+log(1-z[::-1])))
    A2 = -N-S
    return A2, critical, sig


def _find_repeats(arr):
    """Find repeats in the sorted array and return a list of the
    repeats, position of first one in the array and how many there were.
    """
    lastval = arr[0]
    howmany = 0
    ind = 1
    N = len(arr)
    repeat = 0
    reppos = []
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
                reppos.append(ind-1)
                replist.append(arr[ind])
        lastval = arr[ind]
        ind += 1
    if repeat:
        repnum.append(howmany+1)
    return replist, reppos, repnum
            
def ansari(x,y):
    """Determine if the scale parameter for two distributions with equal
    medians is the same.

    Specifically, compute the AB statistic and critical values for rejection
    of the null hypothesis that the ratio of variances is 1.
    """
    x,y = asarray(x),asarray(y)
    m = len(x)
    n = len(y)
    N = m+n
    xy = r_[x,y]  # combine
    sxy = sort(xy)
    replist,reppos,repnum = _find_repeats(sxy)
    if N % 2: # N odd
        ranks = r_[1:(N+1)/2+1,(N-1)/2:0:-1]*1.0
    else:
        ranks = r_[1:N/2+1,N/2:0:-1]*1.0
    repeats = 0
    if replist != []:  # repeats
        repeats = 1
        for k in range(len(replist)):
            pos = reppos[k]
            num = repnum[k]
            ranks[pos:pos+num] = sum(ranks[pos:pos+num])*1.0/num
    # Now compute AB
    AB = 0.0
    for k in n:
        # find the position in the xy array for this value
        pos = nonzero(xy==y[k])[0]
        AB += ranks[pos]
    print AB
    if N % 2:  # N odd
        mnAB = m*(N+1.0)**2 / 4.0 / N
        varAB = m*n*(N+1.0)*(3+N**2)/(48.0*N**2)
    else:
        mnAB = m*(N+2.0)/4.0
        varAB = m*n*(N+2)*(N-2.0)/48/(N-1.0)
    if repeats:
        pass
    return AB


#Tests to include (from R) -- some of these already in stats.
########
# Ansari-Bradley
# Bartlett
# Binomial
# Pearson's Chi-squared
# Association Between Paired samples
# Fisher's exact test
# Fligner-Killeen Test
# Friedman Rank Sum
# Kruskal-Wallis
# Kolmogorov-Smirnov
# Cochran-Mantel-Haenszel Chi-Squared for Count
# McNemar's Chi-squared for Count
# Mood Two-Sample
# Test For Equal Means in One-Way Layout
# Pairwise Comparisons of proportions
# Pairwise t tests
# Tabulate p values for pairwise comparisons
# Pairwise Wilcoxon rank sum tests
# Power calculations two sample test of prop.
# Power calculations for one and two sample t tests
# Equal or Given Proportions
# Trend in Proportions
# Quade Test
# Student's T Test
# F Test to compare two variances
# Wilcoxon Rank Sum and Signed Rank Tests



    


################## test functions #########################

def test(level=1):
    from scipy_base.testing import module_test
    module_test(__name__,__file__,level=level)

def test_suite(level=1):
    from scipy_base.testing import module_test_suite
    return module_test_suite(__name__,__file__,level=level)


