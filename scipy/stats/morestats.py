# Author:  Travis Oliphant, 2002
#
# Further updates and enhancements by many SciPy developers.
#

import math
import statlib
import stats
from stats import find_repeats
import distributions
from numpy import isscalar, r_, log, sum, around, unique, asarray
from numpy import zeros, arange, sort, amin, amax, any, where, \
     atleast_1d, sqrt, ceil, floor, array, poly1d, compress, not_equal, \
     pi, exp, ravel, angle
import scipy
import numpy as np
import types
import scipy.optimize as optimize
import scipy.special as special
import futil
from numpy.testing.decorators import setastest
import warnings

__all__ = ['mvsdist',
           'bayes_mvs', 'kstat', 'kstatvar', 'probplot', 'ppcc_max', 'ppcc_plot',
           'boxcox_llf', 'boxcox', 'boxcox_normmax', 'boxcox_normplot',
           'shapiro', 'anderson', 'ansari', 'bartlett', 'levene', 'binom_test',
           'fligner', 'mood', 'oneway', 'wilcoxon',
           'pdf_fromgamma', 'circmean', 'circvar', 'circstd',
          ]


def bayes_mvs(data, alpha=0.90):
    """Bayesian confidence intervals for the mean, var, and std.

    Parameters
    ----------
    data : array_like
       Input data, if multi-dimensional it is flattened to 1-D by `bayes_mvs`.
       Requires 2 or more data points.
    alpha : float, optional
       Probability that the returned confidence interval contains
       the true parameter.

    Returns
    -------
    Returns a 3 output arguments for each of mean, variance, and standard deviation.
       Each of the outputs is a pair:
          (center, (lower, upper))
       with center the mean of the conditional pdf of the value given the data
       and (lower, upper) is a confidence interval centered on the median,
       containing the estimate to a probability alpha.

    mctr, (ma, mb) :
       Estimates for mean
    vctr, (va, vb) :
       Estimates for variance
    sctr, (sa, sb) :
       Estimates for standard deviation

    Notes
    -----
    Converts data to 1-D and assumes all data has the same mean and variance.
    Uses Jeffrey's prior for variance and std.

    Equivalent to tuple((x.mean(), x.interval(alpha)) for x in mvsdist(dat))

    References
    ----------
    T.E. Oliphant, "A Bayesian perspective on estimating mean, variance, and
    standard-deviation from data", http://hdl.handle.net/1877/438, 2006.

    """
    res = mvsdist(data)
    if alpha >= 1 or alpha <= 0:
        raise ValueError("0 < alpha < 1 is required, but alpha=%s was given." % alpha)
    return tuple((x.mean(), x.interval(alpha)) for x in res)

def mvsdist(data):
    """
    'Frozen' distributions for mean, variance, and standard deviation of data.

    Parameters
    ----------
    data : array_like
        Input array. Converted to 1-D using ravel.
        Requires 2 or more data-points.

    Returns
    -------
    mdist : "frozen" distribution object
        Distribution object representing the mean of the data
    vdist : "frozen" distribution object
        Distribution object representing the variance of the data
    sdist : "frozen" distribution object
        Distribution object representing the standard deviation of the data

    Notes
    -----
    The return values from bayes_mvs(data) is equivalent to
    ``tuple((x.mean(), x.interval(0.90)) for x in mvsdist(data))``.

    In other words, calling ``<dist>.mean()`` and ``<dist>.interval(0.90)``
    on the three distribution objects returned from this function will give
    the same results that are returned from `bayes_mvs`.

    Examples
    --------
    >>> from scipy.stats import mvsdist
    >>> data = [6, 9, 12, 7, 8, 8, 13]
    >>> mean, var, std = mvsdist(data)

    We now have frozen distribution objects "mean", "var" and "std" that we can
    examine:

    >>> mean.mean()
    9.0
    >>> mean.interval(0.95)
    (6.6120585482655692, 11.387941451734431)
    >>> mean.std()
    1.1952286093343936

    """
    x = ravel(data)
    n = len(x)
    if (n < 2):
        raise ValueError("Need at least 2 data-points.")
    xbar = x.mean()
    C = x.var()
    if (n > 1000): # gaussian approximations for large n
        mdist = distributions.norm(loc=xbar, scale=math.sqrt(C/n))
        sdist = distributions.norm(loc=math.sqrt(C), scale=math.sqrt(C/(2.*n)))
        vdist = distributions.norm(loc=C, scale=math.sqrt(2.0/n)*C)
    else:
        nm1 = n-1
        fac = n*C/2.
        val = nm1/2.
        mdist = distributions.t(nm1,loc=xbar,scale=math.sqrt(C/nm1))
        sdist = distributions.gengamma(val,-2,scale=math.sqrt(fac))
        vdist = distributions.invgamma(val,scale=fac)
    return mdist, vdist, sdist


def kstat(data,n=2):
    """
    Return the nth k-statistic (1<=n<=4 so far).

    The nth k-statistic is the unique symmetric unbiased estimator of the nth
    cumulant kappa_n.

    Parameters
    ----------
    data : array_like
        Input array.
    n : int, {1, 2, 3, 4}, optional
        Default is equal to 2.

    Returns
    -------
    kstat : float
        The nth k-statistic.

    See Also
    --------
    kstatvar: Returns an unbiased estimator of the variance of the k-statistic.

    Notes
    -----
    The cumulants are related to central moments but are specifically defined
    using a power series expansion of the logarithm of the characteristic
    function (which is the Fourier transform of the PDF).
    In particular let phi(t) be the characteristic function, then::

        ln phi(t) = > kappa_n (it)^n / n!    (sum from n=0 to inf)

    The first few cumulants (kappa_n)  in terms of central moments (mu_n) are::

        kappa_1 = mu_1
        kappa_2 = mu_2
        kappa_3 = mu_3
        kappa_4 = mu_4 - 3*mu_2**2
        kappa_5 = mu_5 - 10*mu_2 * mu_3

    References
    ----------
    http://mathworld.wolfram.com/k-Statistic.html

    http://mathworld.wolfram.com/Cumulant.html

    """
    if n > 4 or n < 1:
        raise ValueError("k-statistics only supported for 1<=n<=4")
    n = int(n)
    S = zeros(n+1,'d')
    data = ravel(data)
    N = len(data)
    for k in range(1,n+1):
        S[k] = sum(data**k,axis=0)
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
        raise ValueError("Should not be here.")

def kstatvar(data,n=2):
    """
    Returns an unbiased estimator of the variance of the k-statistic.

    See `kstat` for more details of the k-statistic.

    Parameters
    ----------
    data : array_like
        Input array.
    n : int, {1, 2}, optional
        Default is equal to 2.

    Returns
    -------
    kstatvar : float
        The nth k-statistic variance.

    See Also
    --------
    kstat

    """
    data = ravel(data)
    N = len(data)
    if n == 1:
        return kstat(data,n=2)*1.0/N
    elif n == 2:
        k2 = kstat(data,n=2)
        k4 = kstat(data,n=4)
        return (2*k2*k2*N + (N-1)*k4)/(N*(N+1))
    else:
        raise ValueError("Only n=1 or n=2 supported.")


def probplot(x, sparams=(), dist='norm', fit=True, plot=None):
    """
    Calculate quantiles for a probability plot of sample data against a
    specified theoretical distribution.

    `probplot` optionally calculates a best-fit line for the data and plots the
    results using Matplotlib or a given plot function.

    Parameters
    ----------
    x : array_like
        Sample/response data from which `probplot` creates the plot.
    sparams : tuple, optional
        Distribution-specific shape parameters (location(s) and scale(s)).
    dist : str, optional
        Distribution function name. The default is 'norm' for a normal
        probability plot.
    fit : bool, optional
        Fit a least-squares regression (best-fit) line to the sample data if
        True (default).
    plot : object, optional
        If given, plots the quantiles and least squares fit.
        `plot` is an object with methods "plot", "title", "xlabel", "ylabel"
        and "text". The matplotlib.pyplot module or a Matplotlib axes object can
        be used, or a custom object with the same methods.
        By default, no plot is created.

    Notes
    -----
    Even if `plot` is given, the figure is not shown or saved by `probplot`;
    ``plot.show()`` or ``plot.savefig('figname.png')`` should be used after
    calling `probplot`.

    Returns
    -------
    (osm, osr) : tuple of ndarrays
        Tuple of theoretical quantiles (osm, or order statistic medians) and
        ordered responses (osr).
    (slope, intercept, r) : tuple of floats, optional
        Tuple  containing the result of the least-squares fit, if that is
        performed by `probplot`. `r` is the square root of the coefficient of
        determination.  If ``fit=False`` and ``plot=None``, this tuple is not
        returned.

    Examples
    --------
    >>> import scipy.stats as stats
    >>> nsample = 100
    >>> np.random.seed(7654321)

    A t distribution with small degrees of freedom:

    >>> ax1 = plt.subplot(221)
    >>> x = stats.t.rvs(3, size=nsample)
    >>> res = stats.probplot(x, plot=plt)

    A t distribution with larger degrees of freedom:

    >>> ax2 = plt.subplot(222)
    >>> x = stats.t.rvs(25, size=nsample)
    >>> res = stats.probplot(x, plot=plt)

    A mixture of 2 normal distributions with broadcasting:

    >>> ax3 = plt.subplot(223)
    >>> x = stats.norm.rvs(loc=[0,5], scale=[1,1.5], size=(nsample/2.,2)).ravel()
    >>> res = stats.probplot(x, plot=plt)

    A standard normal distribution:

    >>> ax4 = plt.subplot(224)
    >>> x = stats.norm.rvs(loc=0, scale=1, size=nsample)
    >>> res = stats.probplot(x, plot=plt)

    """
    N = len(x)
    Ui = zeros(N) * 1.0
    Ui[-1] = 0.5**(1.0 /N)
    Ui[0] = 1 - Ui[-1]
    i = arange(2, N)
    Ui[1:-1] = (i - 0.3175) / (N + 0.365)
    try:
        ppf_func = eval('distributions.%s.ppf' % dist)
    except AttributeError:
        raise ValueError("%s is not a valid distribution with a ppf." % dist)
    if sparams is None:
        sparams = ()
    if isscalar(sparams):
        sparams = (sparams,)
    if not isinstance(sparams, types.TupleType):
        sparams = tuple(sparams)
    """
    res = inspect.getargspec(ppf_func)
    if not ('loc' == res[0][-2] and 'scale' == res[0][-1] and \
            0.0==res[-1][-2] and 1.0==res[-1][-1]):
        raise ValueError("Function has does not have default location "
              "and scale parameters\n  that are 0.0 and 1.0 respectively.")
    if (len(sparams) < len(res[0])-len(res[-1])-1) or \
       (len(sparams) > len(res[0])-3):
        raise ValueError("Incorrect number of shape parameters.")
    """
    osm = ppf_func(Ui, *sparams)
    osr = sort(x)
    if fit or (plot is not None):
        # perform a linear fit.
        slope, intercept, r, prob, sterrest = stats.linregress(osm, osr)
    if plot is not None:
        plot.plot(osm, osr, 'o', osm, slope*osm + intercept)
        plot.title('Probability Plot')
        plot.xlabel('Quantiles')
        plot.ylabel('Ordered Values')

        xmin = amin(osm)
        xmax = amax(osm)
        ymin = amin(x)
        ymax = amax(x)
        posx = xmin + 0.70 * (xmax - xmin)
        posy = ymin + 0.01 * (ymax - ymin)
        plot.text(posx, posy, "r^2=%1.4f" % r)
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
    except AttributeError:
        raise ValueError("%s is not a valid distribution with a ppf." % dist)
    """
    res = inspect.getargspec(ppf_func)
    if not ('loc' == res[0][-2] and 'scale' == res[0][-1] and \
            0.0==res[-1][-2] and 1.0==res[-1][-1]):
        raise ValueError("Function has does not have default location "
              "and scale parameters\n  that are 0.0 and 1.0 respectively.")
    if (1 < len(res[0])-len(res[-1])-1) or \
       (1 > len(res[0])-3):
        raise ValueError("Must be a one-parameter family.")
    """
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
        plot.plot(svals, ppcc, 'x')
        plot.title('(%s) PPCC Plot' % dist)
        plot.xlabel('Prob Plot Corr. Coef.')#,deltay=-0.01)
        plot.ylabel('Shape Values')#,deltax=-0.01)
    return svals, ppcc

def boxcox_llf(lmb, data):
    """The boxcox log-likelihood function.
    """
    N = len(data)
    y = boxcox(data,lmb)
    my = np.mean(y, axis=0)
    f = (lmb-1)*sum(log(data),axis=0)
    f -= N/2.0*log(sum((y-my)**2.0/N,axis=0))
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
    if N == 500:
        raise RuntimeError("Could not find endpoint.")
    lmplus = optimize.brentq(rootfunc,lmax,newlm,args=(x,target))
    newlm = lmax-0.5
    N = 0
    while (rootfunc(newlm,x,target) > 0.0) and (N < 500):
        newlm += 0.1
        N +=1
    if N == 500:
        raise RuntimeError("Could not find endpoint.")
    lmminus = optimize.brentq(rootfunc, newlm, lmax, args=(x,target))
    return lmminus, lmplus

def boxcox(x,lmbda=None,alpha=None):
    """Return a positive dataset tranformed by a Box-Cox power transformation.

    If lmbda is not None, do the transformation for that value.

    If lmbda is None, find the lambda that maximizes the log-likelihood
    function and return it as the second output argument.

    If alpha is not None, return the 100(1-alpha)% confidence interval for
    lambda as the third output argument.
    """
    if any(x < 0):
        raise ValueError("Data must be positive.")
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
        #JP: this doesn't use sval, creates constant ppcc, and horizontal line
        z = boxcox(x,sval)  #JP: this was missing
        r1,r2 = probplot(z,dist='norm',fit=1)
        ppcc[k] = r2[-1]
        k +=1
    if plot is not None:
        plot.plot(svals, ppcc, 'x')
        plot.title('Box-Cox Normality Plot')
        plot.xlabel('Prob Plot Corr. Coef.')
        plot.ylabel('Transformation parameter')
    return svals, ppcc

def shapiro(x,a=None,reta=False):
    """
    Perform the Shapiro-Wilk test for normality.

    The Shapiro-Wilk test tests the null hypothesis that the
    data was drawn from a normal distribution.

    Parameters
    ----------
    x : array_like
        Array of sample data.
    a : array_like, optional
        Array of internal parameters used in the calculation.  If these
        are not given, they will be computed internally.  If x has length
        n, then a must have length n/2.
    reta : bool, optional
        Whether or not to return the internally computed a values.  The
        default is False.

    Returns
    -------
    W : float
        The test statistic.
    p-value : float
        The p-value for the hypothesis test.
    a : array_like, optional
        If `reta` is True, then these are the internally computed "a"
        values that may be passed into this function on future calls.

    See Also
    --------
    anderson : The Anderson-Darling test for normality

    References
    ----------
    .. [1] http://www.itl.nist.gov/div898/handbook/prc/section2/prc213.htm

    """
    N = len(x)
    if N < 3:
        raise ValueError("Data must be at least length 3.")
    if a is None:
        a = zeros(N,'f')
        init = 0
    else:
        if len(a) != N//2:
            raise ValueError("len(a) must equal len(x)/2")
        init = 1
    y = sort(x)
    a, w, pw, ifault = statlib.swilk(y, a[:N//2], init)
    if not ifault in [0,2]:
        warnings.warn(str(ifault))
    if N > 5000:
        warnings.warn("p-value may not be accurate for N > 5000.")
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
    """
    Anderson-Darling test for data coming from a particular distribution

    The Anderson-Darling test is a modification of the Kolmogorov-
    Smirnov test kstest_ for the null hypothesis that a sample is
    drawn from a population that follows a particular distribution.
    For the Anderson-Darling test, the critical values depend on
    which distribution is being tested against.  This function works
    for normal, exponential, logistic, or Gumbel (Extreme Value
    Type I) distributions.

    Parameters
    ----------
    x : array_like
        array of sample data
    dist : {'norm','expon','logistic','gumbel','extreme1'}, optional
        the type of distribution to test against.  The default is 'norm'
        and 'extreme1' is a synonym for 'gumbel'

    Returns
    -------
    A2 : float
        The Anderson-Darling test statistic
    critical : list
        The critical values for this distribution
    sig : list
        The significance levels for the corresponding critical values
        in percents.  The function returns critical values for a
        differing set of significance levels depending on the
        distribution that is being tested against.

    Notes
    -----
    Critical values provided are for the following significance levels:

    normal/exponenential
        15%, 10%, 5%, 2.5%, 1%
    logistic
        25%, 10%, 5%, 2.5%, 1%, 0.5%
    Gumbel
        25%, 10%, 5%, 2.5%, 1%

    If A2 is larger than these critical values then for the corresponding
    significance level, the null hypothesis that the data come from the
    chosen distribution can be rejected.

    References
    ----------
    .. [1] http://www.itl.nist.gov/div898/handbook/prc/section2/prc213.htm
    .. [2] Stephens, M. A. (1974). EDF Statistics for Goodness of Fit and
           Some Comparisons, Journal of the American Statistical Association,
           Vol. 69, pp. 730-737.
    .. [3] Stephens, M. A. (1976). Asymptotic Results for Goodness-of-Fit
           Statistics with Unknown Parameters, Annals of Statistics, Vol. 4,
           pp. 357-369.
    .. [4] Stephens, M. A. (1977). Goodness of Fit for the Extreme Value
           Distribution, Biometrika, Vol. 64, pp. 583-588.
    .. [5] Stephens, M. A. (1977). Goodness of Fit with Special Reference
           to Tests for Exponentiality , Technical Report No. 262,
           Department of Statistics, Stanford University, Stanford, CA.
    .. [6] Stephens, M. A. (1979). Tests of Fit for the Logistic Distribution
           Based on the Empirical Distribution Function, Biometrika, Vol. 66,
           pp. 591-595.

    """
    if not dist in ['norm','expon','gumbel','extreme1','logistic']:
        raise ValueError("Invalid distribution; dist must be 'norm', "
                            "'expon', 'gumbel', 'extreme1' or 'logistic'.")
    y = sort(x)
    xbar = np.mean(x, axis=0)
    N = len(y)
    if dist == 'norm':
        s = np.std(x, ddof=1, axis=0)
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
            val = [sum(1.0/(1+tmp2),axis=0)-0.5*N,
                   sum(tmp*(1.0-tmp2)/(1+tmp2),axis=0)+N]
            return array(val)
        sol0=array([xbar,np.std(x, ddof=1, axis=0)])
        sol = optimize.fsolve(rootfunc,sol0,args=(x,N),xtol=1e-5)
        w = (y-sol[0])/sol[1]
        z = distributions.logistic.cdf(w)
        sig = array([25,10,5,2.5,1,0.5])
        critical = around(_Avals_logistic / (1.0+0.25/N),3)
    else:  # (dist == 'gumbel') or (dist == 'extreme1'):
        #the following is incorrect, see ticket:1097
##        def fixedsolve(th,xj,N):
##            val = stats.sum(xj)*1.0/N
##            tmp = exp(-xj/th)
##            term = sum(xj*tmp,axis=0)
##            term /= sum(tmp,axis=0)
##            return val - term
##        s = optimize.fixed_point(fixedsolve, 1.0, args=(x,N),xtol=1e-5)
##        xbar = -s*log(sum(exp(-x/s),axis=0)*1.0/N)
        xbar, s = distributions.gumbel_l.fit(x)
        w = (y-xbar)/s
        z = distributions.gumbel_l.cdf(w)
        sig = array([25,10,5,2.5,1])
        critical = around(_Avals_gumbel / (1.0 + 0.2/sqrt(N)),3)

    i = arange(1,N+1)
    S = sum((2*i-1.0)/N*(log(z)+log(1-z[::-1])),axis=0)
    A2 = -N-S
    return A2, critical, sig


def ansari(x,y):
    """
    Perform the Ansari-Bradley test for equal scale parameters

    The Ansari-Bradley test is a non-parametric test for the equality
    of the scale parameter of the distributions from which two
    samples were drawn.

    Parameters
    ----------
    x, y : array_like
        arrays of sample data

    Returns
    -------
    p-value : float
        The p-value of the hypothesis test

    See Also
    --------
    fligner : A non-parametric test for the equality of k variances
    mood : A non-parametric test for the equality of two scale parameters

    Notes
    -----
    The p-value given is exact when the sample sizes are both less than
    55 and there are no ties, otherwise a normal approximation for the
    p-value is used.

    References
    ----------
    .. [1] Sprent, Peter and N.C. Smeeton.  Applied nonparametric statistical
           methods.  3rd ed. Chapman and Hall/CRC. 2001.  Section 5.8.2.

    """
    x,y = asarray(x),asarray(y)
    n = len(x)
    m = len(y)
    if m < 1:
        raise ValueError("Not enough other observations.")
    if n < 1:
        raise ValueError("Not enough test observations.")
    N = m+n
    xy = r_[x,y]  # combine
    rank = stats.rankdata(xy)
    symrank = amin(array((rank,N-rank+1)),0)
    AB = sum(symrank[:n],axis=0)
    uxy = unique(xy)
    repeats = (len(uxy) != len(xy))
    exact = ((m<55) and (n<55) and not repeats)
    if repeats and ((m < 55)  or (n < 55)):
        warnings.warn("Ties preclude use of exact statistic.")
    if exact:
        astart, a1, ifault = statlib.gscale(n,m)
        ind = AB-astart
        total = sum(a1,axis=0)
        if ind < len(a1)/2.0:
            cind = int(ceil(ind))
            if (ind == cind):
                pval = 2.0*sum(a1[:cind+1],axis=0)/total
            else:
                pval = 2.0*sum(a1[:cind],axis=0)/total
        else:
            find = int(floor(ind))
            if (ind == floor(ind)):
                pval = 2.0*sum(a1[find:],axis=0)/total
            else:
                pval = 2.0*sum(a1[find+1:],axis=0)/total
        return AB, min(1.0,pval)

    # otherwise compute normal approximation
    if N % 2:  # N odd
        mnAB = n*(N+1.0)**2 / 4.0 / N
        varAB = n*m*(N+1.0)*(3+N**2)/(48.0*N**2)
    else:
        mnAB = n*(N+2.0)/4.0
        varAB = m*n*(N+2)*(N-2.0)/48/(N-1.0)
    if repeats:   # adjust variance estimates
        # compute sum(tj * rj**2,axis=0)
        fac = sum(symrank**2,axis=0)
        if N % 2: # N odd
            varAB = m*n*(16*N*fac-(N+1)**4)/(16.0 * N**2 * (N-1))
        else:  # N even
            varAB = m*n*(16*fac-N*(N+2)**2)/(16.0 * N * (N-1))
    z = (AB - mnAB)/sqrt(varAB)
    pval = distributions.norm.sf(abs(z)) * 2.0
    return AB, pval

def bartlett(*args):
    """
    Perform Bartlett's test for equal variances

    Bartlett's test tests the null hypothesis that all input samples
    are from populations with equal variances.  For samples
    from significantly non-normal populations, Levene's test
    `levene`_ is more robust.

    Parameters
    ----------
    sample1, sample2,... : array_like
        arrays of sample data.  May be different lengths.

    Returns
    -------
    T : float
        The test statistic.
    p-value : float
        The p-value of the test.

    References
    ----------
    .. [1]  http://www.itl.nist.gov/div898/handbook/eda/section3/eda357.htm

    .. [2]  Snedecor, George W. and Cochran, William G. (1989), Statistical
              Methods, Eighth Edition, Iowa State University Press.

    """
    k = len(args)
    if k < 2:
        raise ValueError("Must enter at least two input sample vectors.")
    Ni = zeros(k)
    ssq = zeros(k,'d')
    for j in range(k):
        Ni[j] = len(args[j])
        ssq[j] = np.var(args[j], ddof=1)
    Ntot = sum(Ni,axis=0)
    spsq = sum((Ni-1)*ssq,axis=0)/(1.0*(Ntot-k))
    numer = (Ntot*1.0-k)*log(spsq) - sum((Ni-1.0)*log(ssq),axis=0)
    denom = 1.0 + (1.0/(3*(k-1)))*((sum(1.0/(Ni-1.0),axis=0))-1.0/(Ntot-k))
    T = numer / denom
    pval = distributions.chi2.sf(T,k-1) # 1 - cdf
    return T, pval


def levene(*args,**kwds):
    """
    Perform Levene test for equal variances.

    The Levene test tests the null hypothesis that all input samples
    are from populations with equal variances.  Levene's test is an
    alternative to Bartlett's test `bartlett` in the case where
    there are significant deviations from normality.

    Parameters
    ----------
    sample1, sample2, ... : array_like
        The sample data, possibly with different lengths
    center : {'mean', 'median', 'trimmed'}, optional
        Which function of the data to use in the test.  The default
        is 'median'.
    proportiontocut : float, optional
        When `center` is 'trimmed', this gives the proportion of data points
        to cut from each end. (See `scipy.stats.trim_mean`.)
        Default is 0.05.

    Returns
    -------
    W : float
        The test statistic.
    p-value : float
        The p-value for the test.

    Notes
    -----
    Three variations of Levene's test are possible.  The possibilities
    and their recommended usages are:

      * 'median' : Recommended for skewed (non-normal) distributions>
      * 'mean' : Recommended for symmetric, moderate-tailed distributions.
      * 'trimmed' : Recommended for heavy-tailed distributions.

    References
    ----------
    .. [1]  http://www.itl.nist.gov/div898/handbook/eda/section3/eda35a.htm
    .. [2]   Levene, H. (1960). In Contributions to Probability and Statistics:
               Essays in Honor of Harold Hotelling, I. Olkin et al. eds.,
               Stanford University Press, pp. 278-292.
    .. [3]  Brown, M. B. and Forsythe, A. B. (1974), Journal of the American
              Statistical Association, 69, 364-367

    """
    # Handle keyword arguments.
    center = 'median'
    proportiontocut = 0.05
    for kw, value in kwds.items():
        if kw not in ['center', 'proportiontocut']:
            raise TypeError("levene() got an unexpected keyword argument '%s'" % kw)
        if kw == 'center':
            center = value
        else:
            proportiontocut = value

    k = len(args)
    if k < 2:
        raise ValueError("Must enter at least two input sample vectors.")
    Ni = zeros(k)
    Yci = zeros(k,'d')

    if not center in ['mean','median','trimmed']:
        raise ValueError("Keyword argument <center> must be 'mean', 'median'"
              + "or 'trimmed'.")

    if center == 'median':
        func = lambda x: np.median(x, axis=0)
    elif center == 'mean':
        func = lambda x: np.mean(x, axis=0)
    else: # center == 'trimmed'
        args = tuple(stats.trimboth(arg, proportiontocut) for arg in args)
        func = lambda x: np.mean(x, axis=0)

    for j in range(k):
        Ni[j] = len(args[j])
        Yci[j] = func(args[j])
    Ntot = sum(Ni,axis=0)

    # compute Zij's
    Zij = [None]*k
    for i in range(k):
        Zij[i] = abs(asarray(args[i])-Yci[i])
    # compute Zbari
    Zbari = zeros(k,'d')
    Zbar = 0.0
    for i in range(k):
        Zbari[i] = np.mean(Zij[i], axis=0)
        Zbar += Zbari[i]*Ni[i]
    Zbar /= Ntot

    numer = (Ntot-k)*sum(Ni*(Zbari-Zbar)**2,axis=0)

    # compute denom_variance
    dvar = 0.0
    for i in range(k):
        dvar += sum((Zij[i]-Zbari[i])**2,axis=0)

    denom = (k-1.0)*dvar

    W = numer / denom
    pval = distributions.f.sf(W,k-1,Ntot-k) # 1 - cdf
    return W, pval

@setastest(False)
def binom_test(x,n=None,p=0.5):
    """
    Perform a test that the probability of success is p.

    This is an exact, two-sided test of the null hypothesis
    that the probability of success in a Bernoulli experiment
    is `p`.

    Parameters
    ----------
    x : integer or array_like
        the number of successes, or if x has length 2, it is the
        number of successes and the number of failures.
    n : integer
        the number of trials.  This is ignored if x gives both the
        number of successes and failures
    p : float, optional
        The hypothesized probability of success.  0 <= p <= 1. The
        default value is p = 0.5

    Returns
    -------
    p-value : float
        The p-value of the hypothesis test

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Binomial_test

    """
    x = atleast_1d(x).astype(np.integer)
    if len(x) == 2:
        n = x[1]+x[0]
        x = x[0]
    elif len(x) == 1:
        x = x[0]
        if n is None or n < x:
            raise ValueError("n must be >= x")
        n = np.int_(n)
    else:
        raise ValueError("Incorrect length for x.")

    if (p > 1.0) or (p < 0.0):
        raise ValueError("p must be in range [0,1]")

    d = distributions.binom.pmf(x,n,p)
    rerr = 1+1e-7
    if (x < p*n):
        i = np.arange(np.ceil(p*n),n+1)
        y = np.sum(distributions.binom.pmf(i,n,p) <= d*rerr,axis=0)
        pval = distributions.binom.cdf(x,n,p) + distributions.binom.sf(n-y,n,p)
    else:
        i = np.arange(np.floor(p*n))
        y = np.sum(distributions.binom.pmf(i,n,p) <= d*rerr,axis=0)
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
    """
    Perform Fligner's test for equal variances.

    Fligner's test tests the null hypothesis that all input samples
    are from populations with equal variances.  Fligner's test is
    non-parametric in contrast to Bartlett's test `bartlett` and
    Levene's test `levene`.

    Parameters
    ----------
    sample1, sample2, ... : array_like
        arrays of sample data.  Need not be the same length
    center : {'mean', 'median', 'trimmed'}, optional
        keyword argument controlling which function of the data
        is used in computing the test statistic.  The default
        is 'median'.
    proportiontocut : float, optional
        When `center` is 'trimmed', this gives the proportion of data points
        to cut from each end. (See `scipy.stats.trim_mean`.)
        Default is 0.05.

    Returns
    -------
    Xsq : float
        the test statistic
    p-value : float
        the p-value for the hypothesis test

    Notes
    -----
    As with Levene's test there are three variants
    of Fligner's test that differ by the measure of central
    tendency used in the test.  See `levene` for more information.

    References
    ----------
    .. [1] http://www.stat.psu.edu/~bgl/center/tr/TR993.ps

    .. [2] Fligner, M.A. and Killeen, T.J. (1976). Distribution-free two-sample
           tests for scale. 'Journal of the American Statistical Association.'
           71(353), 210-213.

    """
    # Handle keyword arguments.
    center = 'median'
    proportiontocut = 0.05
    for kw, value in kwds.items():
        if kw not in ['center', 'proportiontocut']:
            raise TypeError("fligner() got an unexpected keyword argument '%s'" % kw)
        if kw == 'center':
            center = value
        else:
            proportiontocut = value

    k = len(args)
    if k < 2:
        raise ValueError("Must enter at least two input sample vectors.")

    if not center in ['mean','median','trimmed']:
        raise ValueError("Keyword argument <center> must be 'mean', 'median'"
              + "or 'trimmed'.")

    if center == 'median':
        func = lambda x: np.median(x, axis=0)
    elif center == 'mean':
        func = lambda x: np.mean(x, axis=0)
    else: # center == 'trimmed'
        args = tuple(stats.trimboth(arg, proportiontocut) for arg in args)
        func = lambda x: np.mean(x, axis=0)

    Ni = asarray([len(args[j]) for j in range(k)])
    Yci = asarray([func(args[j]) for j in range(k)])
    Ntot = sum(Ni,axis=0)
    # compute Zij's
    Zij = [abs(asarray(args[i])-Yci[i]) for i in range(k)]
    allZij = []
    g = [0]
    for i in range(k):
        allZij.extend(list(Zij[i]))
        g.append(len(allZij))

    ranks = stats.rankdata(allZij)
    a = distributions.norm.ppf(ranks/(2*(Ntot+1.0)) + 0.5)

    # compute Aibar
    Aibar = _apply_func(a,g,sum) / Ni
    anbar = np.mean(a, axis=0)
    varsq = np.var(a,axis=0, ddof=1)
    Xsq = sum(Ni*(asarray(Aibar)-anbar)**2.0,axis=0)/varsq
    pval = distributions.chi2.sf(Xsq,k-1) # 1 - cdf
    return Xsq, pval


def mood(x,y):
    """
    Perform Mood's test for equal scale parameters.

    Mood's two-sample test for scale parameters is a non-parametric
    test for the null hypothesis that two samples are drawn from the
    same distribution with the same scale parameter.

    Parameters
    ----------
    x, y : array_like
        Arrays of sample data.

    Returns
    -------
    p-value : float
        The p-value for the hypothesis test.

    See Also
    --------
    fligner : A non-parametric test for the equality of k variances
    ansari : A non-parametric test for the equality of 2 variances
    bartlett : A parametric test for equality of k variances in normal samples
    levene : A parametric test for equality of k variances

    Notes
    -----
    The data are assumed to be drawn from probability distributions f(x) and
    f(x/s)/s respectively, for some probability density function f.  The
    null hypothesis is that s = 1.

    """
    n = len(x)
    m = len(y)
    xy = r_[x,y]
    N = m+n
    if N < 3:
        raise ValueError("Not enough observations.")
    ranks = stats.rankdata(xy)
    Ri = ranks[:n]
    M = sum((Ri - (N+1.0)/2)**2,axis=0)
    # Approx stat.
    mnM = n*(N*N-1.0)/12
    varM = m*n*(N+1.0)*(N+2)*(N-2)/180
    z = (M-mnM)/sqrt(varM)

    # Numerically better than p = norm.cdf(x); p = min(p, 1 - p)
    if z > 0:
        pval = distributions.norm.sf(z)
    else:
        pval = distributions.norm.cdf(z)

    # Account for two-sidedness
    pval *= 2.
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
        raise ValueError("Must enter at least two input sample vectors.")
    if 'equal_var' in kwds.keys():
        if kwds['equal_var']: evar = 1
        else: evar = 0
    else:
        evar = 0

    Ni = array([len(args[i]) for i in range(k)])
    Mi = array([np.mean(args[i], axis=0) for i in range(k)])
    Vi = array([np.var(args[i]) for i in range(k)])
    Wi = Ni / Vi
    swi = sum(Wi,axis=0)
    N = sum(Ni,axis=0)
    my = sum(Mi*Ni,axis=0)*1.0/N
    tmp = sum((1-Wi/swi)**2 / (Ni-1.0),axis=0)/(k*k-1.0)
    if evar:
        F = ((sum(Ni*(Mi-my)**2,axis=0) / (k-1.0)) / (sum((Ni-1.0)*Vi,axis=0) / (N-k)))
        pval = distributions.f.sf(F,k-1,N-k)  # 1-cdf
    else:
        m = sum(Wi*Mi,axis=0)*1.0/swi
        F = sum(Wi*(Mi-m)**2,axis=0) / ((k-1.0)*(1+2*(k-2)*tmp))
        pval = distributions.f.sf(F,k-1.0,1.0/(3*tmp))

    return F, pval


def wilcoxon(x,y=None):
    """
    Calculate the Wilcoxon signed-rank test.

    The Wilcoxon signed-rank test tests the null hypothesis that two
    related samples come from the same distribution. It is a a
    non-parametric version of the paired T-test.

    Parameters
    ----------
    x : array_like
        The first set of measurements.
    y : array_like, optional
        The second set of measurements.  If y is not given, then the x array
        is considered to be the differences between the two sets of
        measurements.

    Returns
    -------
    z-statistic : float
        The test statistic under the large-sample approximation that the
        signed-rank statistic is normally distributed.
    p-value : float
        The two-sided p-value for the test.

    Notes
    -----
    Because the normal approximation is used for the calculations, the
    samples used should be large.  A typical rule is to require that
    n > 20.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test

    """
    if y is None:
        d = x
    else:
        x, y = map(asarray, (x, y))
        if len(x) <> len(y):
            raise ValueError('Unequal N in wilcoxon.  Aborting.')
        d = x-y
    d = compress(not_equal(d,0),d,axis=-1) # Keep all non-zero differences
    count = len(d)
    if (count < 10):
        warnings.warn("Warning: sample size too small for normal approximation.")
    r = stats.rankdata(abs(d))
    r_plus = sum((d > 0)*r,axis=0)
    r_minus = sum((d < 0)*r,axis=0)
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
    prob = 2 * distributions.norm.sf(abs(z))
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

def circmean(samples, high=2*pi, low=0):
    """
    Compute the circular mean for samples assumed to be in the range
    [low to high].

    Parameters
    ----------
    samples : array_like
        Input array.
    low : float or int, optional
        Low boundary for circular mean range.  Default is 0.
    high : float or int, optional
        High boundary for circular mean range.  Default is 2*pi.

    Returns
    -------
    circmean : float
        Circular mean.

    """
    ang = (samples - low)*2*pi / (high-low)
    res = angle(np.mean(exp(1j*ang), axis=0))
    if (res < 0):
        res = res + 2*pi
    return res*(high-low)/2.0/pi + low

def circvar(samples, high=2*pi, low=0):
    """
    Compute the circular variance for samples assumed to be in the range
    [low to high].

    Parameters
    ----------
    samples : array_like
        Input array.
    low : float or int, optional
        Low boundary for circular variance range.  Default is 0.
    high : float or int, optional
        High boundary for circular variance range.  Default is 2*pi.

    Returns
    -------
    circvar : float
        Circular variance.

    """
    ang = (samples - low)*2*pi / (high-low)
    res = np.mean(exp(1j*ang), axis=0)
    V = 1-abs(res)
    return ((high-low)/2.0/pi)**2 * V

def circstd(samples, high=2*pi, low=0):
    """
    Compute the circular standard deviation for samples assumed to be in the
    range [low to high].

    Parameters
    ----------
    samples : array_like
        Input array.
    low : float or int, optional
        Low boundary for circular standard deviation range.  Default is 0.
    high : float or int, optional
        High boundary for circular standard deviation range.  Default is 2*pi.

    Returns
    -------
    circstd : float
        Circular standard deviation.

    """
    ang = (samples - low)*2*pi / (high-low)
    res = np.mean(exp(1j*ang), axis=0)
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
