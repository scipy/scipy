Statistics
==========

.. sectionauthor:: Travis E. Oliphant

Introduction
------------

SciPy has a tremendous number of basic statistics routines with more
easily added by the end user (if you create one please contribute it).
All of the statistics functions are located in the sub-package
:mod:`scipy.stats` and a fairly complete listing of these functions
can be had using ``info(stats)``.

Random Variables
^^^^^^^^^^^^^^^^

There are two general distribution classes that have been implemented
for encapsulating
:ref:`continuous random variables <continuous-random-variables>`
and
:ref:`discrete random variables <discrete-random-variables>`
. Over 80 continuous random variables  and 10 discrete random
variables have been implemented using these classes. The list of the
random variables available is in the docstring for the stats sub-
package.


Note: The following is work in progress

Distributions
-------------


First some imports

    >>> import numpy as np
    >>> from scipy import stats
    >>> import warnings
    >>> warnings.simplefilter('ignore', DeprecationWarning)

We can obtain the list of available distribution through introspection:

    >>> dist_continu = [d for d in dir(stats) if
    ...                 isinstance(getattr(stats,d), stats.rv_continuous)]
    >>> dist_discrete = [d for d in dir(stats) if
    ...                  isinstance(getattr(stats,d), stats.rv_discrete)]
    >>> print 'number of continuous distributions:', len(dist_continu)
    number of continuous distributions: 84
    >>> print 'number of discrete distributions:  ', len(dist_discrete)
    number of discrete distributions:   12




Distributions can be used in one of two ways, either by passing all distribution
parameters to each method call or by freezing the parameters for the instance
of the distribution. As an example, we can get the median of the distribution by using
the percent point function, ppf, which is the inverse of the cdf:

    >>> print stats.nct.ppf(0.5, 10, 2.5)
    2.56880722561
    >>> my_nct = stats.nct(10, 2.5)
    >>> print my_nct.ppf(0.5)
    2.56880722561

``help(stats.nct)`` prints the complete docstring of the distribution. Instead
we can print just some basic information::

    >>> print stats.nct.extradoc #contains the distribution specific docs
    Non-central Student T distribution

                                     df**(df/2) * gamma(df+1)
    nct.pdf(x,df,nc) = --------------------------------------------------
                       2**df*exp(nc**2/2)*(df+x**2)**(df/2) * gamma(df/2)
    for df > 0, nc > 0.


    >>> print 'number of arguments: %d, shape parameters: %s'% (stats.nct.numargs,
    ...                                                         stats.nct.shapes)
    number of arguments: 2, shape parameters: df,nc
    >>> print 'bounds of distribution lower: %s, upper: %s' % (stats.nct.a,
    ...                                                        stats.nct.b)
    bounds of distribution lower: -1.#INF, upper: 1.#INF

We can list all methods and properties of the distribution with
``dir(stats.nct)``. Some of the methods are private methods, that are
not named as such, i.e. no leading underscore, for example veccdf or
xa and xb are for internal calculation. The main methods we can see
when we list the methods of the frozen distribution:

    >>> print dir(my_nct) #reformatted
        ['__class__', '__delattr__', '__dict__', '__doc__', '__getattribute__',
        '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__',
        '__repr__', '__setattr__', '__str__', '__weakref__', 'args', 'cdf', 'dist',
        'entropy', 'isf', 'kwds', 'moment', 'pdf', 'pmf', 'ppf', 'rvs', 'sf', 'stats']


The main public methods are:

* rvs:   Random Variates
* pdf:   Probability Density Function
* cdf:   Cumulative Distribution Function
* sf:    Survival Function (1-CDF)
* ppf:   Percent Point Function (Inverse of CDF)
* isf:   Inverse Survival Function (Inverse of SF)
* stats: Return mean, variance, (Fisher's) skew, or (Fisher's) kurtosis
* moment: non-central moments of the distribution

The main additional methods of the not frozen distribution are related to the estimation
of distrition parameters:

* fit:   maximum likelihood estimation of distribution parameters, including location
         and scale
* fit_loc_scale: estimation of location and scale when shape parameters are given
* nnlf:  negative log likelihood function
* expect: Calculate the expectation of a function against the pdf or pmf

All continuous distributions take `loc` and `scale` as keyword
parameters to adjust the location and scale of the distribution,
e.g. for the standard normal distribution location is the mean and
scale is the standard deviation. The standardized distribution for a
random variable `x` is obtained through ``(x - loc) / scale``.

Discrete distribution have most of the same basic methods, however
pdf is replaced the probability mass function `pmf`, no estimation
methods, such as fit, are available, and scale is not a valid
keyword parameter. The location parameter, keyword `loc` can be used
to shift the distribution.

The basic methods, pdf, cdf, sf, ppf, and isf are vectorized with
``np.vectorize``, and the usual numpy broadcasting is applied. For
example, we can calculate the critical values for the upper tail of
the t distribution for different probabilites and degrees of freedom.

    >>> stats.t.isf([0.1, 0.05, 0.01], [[10], [11]])
    array([[ 1.37218364,  1.81246112,  2.76376946],
           [ 1.36343032,  1.79588482,  2.71807918]])

Here, the first row are the critical values for 10 degrees of freedom and the second row
is for 11 d.o.f., i.e. this is the same as

    >>> stats.t.isf([0.1, 0.05, 0.01], 10)
    array([ 1.37218364,  1.81246112,  2.76376946])
    >>> stats.t.isf([0.1, 0.05, 0.01], 11)
    array([ 1.36343032,  1.79588482,  2.71807918])

If both, probabilities and degrees of freedom have the same array shape, then element
wise matching is used. As an example, we can obtain the 10% tail for 10 d.o.f., the 5% tail
for 11 d.o.f. and the 1% tail for 12 d.o.f. by

    >>> stats.t.isf([0.1, 0.05, 0.01], [10, 11, 12])
    array([ 1.37218364,  1.79588482,  2.68099799])



Performance and Remaining Issues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The performance of the individual methods, in terms of speed, varies
widely by distribution and method. The results of a method are
obtained in one of two ways, either by explicit calculation or by a
generic algorithm that is independent of the specific distribution.
Explicit calculation, requires that the method is directly specified
for the given distribution, either through analytic formulas or
through special functions in scipy.special or numpy.random for
`rvs`. These are usually relatively fast calculations. The generic
methods are used if the distribution does not specify any explicit
calculation. To define a distribution, only one of pdf or cdf is
necessary, all other methods can be derived using numeric integration
and root finding. These indirect methods can be very slow. As an
example, ``rgh = stats.gausshyper.rvs(0.5, 2, 2, 2, size=100)`` creates
random variables in a very indirect way and takes about 19 seconds
for 100 random variables on my computer, while one million random
variables from the standard normal or from the t distribution take
just above one second.


The distributions in scipy.stats have recently been corrected and improved
and gained a considerable test suite, however a few issues remain:

* skew and kurtosis, 3rd and 4th moments and entropy are not thoroughly
  tested and some coarse testing indicates that there are still some
  incorrect results left.
* the distributions have been tested over some range of parameters,
  however in some corner ranges, a few incorrect results may remain.
* the maximum likelihood estimation in `fit` does not work with
  default starting parameters for all distributions and the user
  needs to supply good starting parameters. Also, for some
  distribution using a maximum likelihood estimator might
  inherently not be the best choice.


The next example shows how to build our own discrete distribution,
and more examples for the usage of the distributions are shown below
together with the statistical tests.



Example: discrete distribution rv_discrete
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the following we use stats.rv_discrete to generate a discrete distribution
that has the probabilites of the truncated normal for the intervalls
centered around the integers.


    >>> npoints = 20 # number of integer support points of the distribution minus 1
    >>> npointsh = npoints / 2
    >>> npointsf = float(npoints)
    >>> nbound = 4 # bounds for the truncated normal
    >>> normbound = (1+1/npointsf) * nbound # actual bounds of truncated normal
    >>> grid = np.arange(-npointsh, npointsh+2, 1) # integer grid
    >>> gridlimitsnorm = (grid-0.5) / npointsh * nbound # bin limits for the truncnorm
    >>> gridlimits = grid - 0.5
    >>> grid = grid[:-1]
    >>> probs = np.diff(stats.truncnorm.cdf(gridlimitsnorm, -normbound, normbound))
    >>> gridint = grid
    >>> normdiscrete = stats.rv_discrete(values = (gridint,
    ...              np.round(probs, decimals=7)), name='normdiscrete')

From the docstring of rv_discrete:
 "You can construct an aribtrary discrete rv where P{X=xk} = pk by
 passing to the rv_discrete initialization method (through the values=
 keyword) a tuple of sequences (xk, pk) which describes only those
 values of X (xk) that occur with nonzero probability (pk)."

There are some requirements for this distribution to work. The
keyword `name` is required. The support points of the distribution
xk have to be integers. Also, I needed  to limit the number of
decimals. If the last two requirements are not satisfied an
exception may be raised or the resulting numbers may be incorrect.

After defining the distribution, we obtain access to all methods of
discrete distributions.

    >>> print 'mean = %6.4f, variance = %6.4f, skew = %6.4f, kurtosis = %6.4f'% \
    ...       normdiscrete.stats(moments =  'mvsk')
    mean = -0.0000, variance = 6.3302, skew = 0.0000, kurtosis = -0.0076

    >>> nd_std = np.sqrt(normdiscrete.stats(moments =  'v'))

**Generate a random sample and compare observed frequencies with probabilities**

    >>> n_sample = 500
    >>> np.random.seed(87655678) # fix the seed for replicability
    >>> rvs = normdiscrete.rvs(size=n_sample)
    >>> rvsnd = rvs
    >>> f, l = np.histogram(rvs, bins=gridlimits)
    >>> sfreq = np.vstack([gridint, f, probs*n_sample]).T
    >>> print sfreq
    [[ -1.00000000e+01   0.00000000e+00   2.95019349e-02]
     [ -9.00000000e+00   0.00000000e+00   1.32294142e-01]
     [ -8.00000000e+00   0.00000000e+00   5.06497902e-01]
     [ -7.00000000e+00   2.00000000e+00   1.65568919e+00]
     [ -6.00000000e+00   1.00000000e+00   4.62125309e+00]
     [ -5.00000000e+00   9.00000000e+00   1.10137298e+01]
     [ -4.00000000e+00   2.60000000e+01   2.24137683e+01]
     [ -3.00000000e+00   3.70000000e+01   3.89503370e+01]
     [ -2.00000000e+00   5.10000000e+01   5.78004747e+01]
     [ -1.00000000e+00   7.10000000e+01   7.32455414e+01]
     [  0.00000000e+00   7.40000000e+01   7.92618251e+01]
     [  1.00000000e+00   8.90000000e+01   7.32455414e+01]
     [  2.00000000e+00   5.50000000e+01   5.78004747e+01]
     [  3.00000000e+00   5.00000000e+01   3.89503370e+01]
     [  4.00000000e+00   1.70000000e+01   2.24137683e+01]
     [  5.00000000e+00   1.10000000e+01   1.10137298e+01]
     [  6.00000000e+00   4.00000000e+00   4.62125309e+00]
     [  7.00000000e+00   3.00000000e+00   1.65568919e+00]
     [  8.00000000e+00   0.00000000e+00   5.06497902e-01]
     [  9.00000000e+00   0.00000000e+00   1.32294142e-01]
     [  1.00000000e+01   0.00000000e+00   2.95019349e-02]]


.. plot:: examples/normdiscr_plot1.py
   :align: center
   :include-source: 0


.. plot:: examples/normdiscr_plot2.py
   :align: center
   :include-source: 0


Next, we can test, whether our sample was generated by our normdiscrete
distribution. This also verifies, whether the random numbers are generated
correctly

The chisquare test requires that there are a minimum number of observations
in each bin. We combine the tail bins into larger bins so that they contain
enough observations.

    >>> f2 = np.hstack([f[:5].sum(), f[5:-5], f[-5:].sum()])
    >>> p2 = np.hstack([probs[:5].sum(), probs[5:-5], probs[-5:].sum()])
    >>> ch2, pval = stats.chisquare(f2, p2*n_sample)

    >>> print 'chisquare for normdiscrete: chi2 = %6.3f pvalue = %6.4f' % (ch2, pval)
    chisquare for normdiscrete: chi2 = 12.466 pvalue = 0.4090

The pvalue in this case is high, so we can be quite confident that
our random sample was actually generated by the distribution.



Analysing One Sample
--------------------

First, we create some random variables. We set a seed so that in each run
we get identical results to look at. As an example we take a sample from
the Student t distribution:

    >>> np.random.seed(282629734)
    >>> x = stats.t.rvs(10, size=1000)

Here, we set the required shape parameter of the t distribution, which
in statistics corresponds to the degrees of freedom, to 10. Using size=100 means
that our sample consists of 1000 independently drawn (pseudo) random numbers.
Since we did not specify the keyword arguments `loc` and `scale`, those are
set to their default values zero and one.

Descriptive Statistics
^^^^^^^^^^^^^^^^^^^^^^

`x` is a numpy array, and we have direct access to all array methods, e.g.

    >>> print x.max(), x.min()  # equivalent to np.max(x), np.min(x)
    5.26327732981 -3.78975572422
    >>> print x.mean(), x.var() # equivalent to np.mean(x), np.var(x)
    0.0140610663985 1.28899386208


How do the some sample properties compare to their theoretical counterparts?

    >>> m, v, s, k = stats.t.stats(10, moments='mvsk')
    >>> n, (smin, smax), sm, sv, ss, sk = stats.describe(x)

    >>> print 'distribution:',
    distribution:
    >>> sstr = 'mean = %6.4f, variance = %6.4f, skew = %6.4f, kurtosis = %6.4f'
    >>> print sstr %(m, v, s ,k)
    mean = 0.0000, variance = 1.2500, skew = 0.0000, kurtosis = 1.0000
    >>> print 'sample:      ',
    sample:
    >>> print sstr %(sm, sv, ss, sk)
    mean = 0.0141, variance = 1.2903, skew = 0.2165, kurtosis = 1.0556

Note: stats.describe uses the unbiased estimator for the variance, while
np.var is the biased estimator.


For our sample the sample statistics differ a by a small amount from
their theoretical counterparts.


T-test and KS-test
^^^^^^^^^^^^^^^^^^

We can use the t-test to test whether the mean of our sample differs
in a statistcally significant way from the theoretical expectation.

    >>> print 't-statistic = %6.3f pvalue = %6.4f' %  stats.ttest_1samp(x, m)
    t-statistic =  0.391 pvalue = 0.6955

The pvalue is 0.7, this means that with an alpha error of, for
example, 10%, we cannot reject the hypothesis that the sample mean
is equal to zero, the expectation of the standard t-distribution.


As an exercise, we can calculate our ttest also directly without
using the provided function, which should give us the same answer,
and so it does:

    >>> tt = (sm-m)/np.sqrt(sv/float(n))  # t-statistic for mean
    >>> pval = stats.t.sf(np.abs(tt), n-1)*2  # two-sided pvalue = Prob(abs(t)>tt)
    >>> print 't-statistic = %6.3f pvalue = %6.4f' % (tt, pval)
    t-statistic =  0.391 pvalue = 0.6955

The Kolmogorov-Smirnov test can be used to test the hypothesis that
the sample comes from the standard t-distribution

    >>> print 'KS-statistic D = %6.3f pvalue = %6.4f' % stats.kstest(x, 't', (10,))
    KS-statistic D =  0.016 pvalue = 0.9606

Again the p-value is high enough that we cannot reject the
hypothesis that the random sample really is distributed according to the
t-distribution. In real applications, we don't know what the
underlying distribution is. If we perform the Kolmogorov-Smirnov
test of our sample against the standard normal distribution, then we
also cannot reject the hypothesis that our sample was generated by the
normal distribution given that in this example the p-value is almost 40%.

    >>> print 'KS-statistic D = %6.3f pvalue = %6.4f' % stats.kstest(x,'norm')
    KS-statistic D =  0.028 pvalue = 0.3949

However, the standard normal distribution has a variance of 1, while our
sample has a variance of 1.29. If we standardize our sample and test it
against the normal distribution, then the p-value is again large enough
that we cannot reject the hypothesis that the sample came form the
normal distribution.

    >>> d, pval = stats.kstest((x-x.mean())/x.std(), 'norm')
    >>> print 'KS-statistic D = %6.3f pvalue = %6.4f' % (d, pval)
    KS-statistic D =  0.032 pvalue = 0.2402

Note: The Kolmogorov-Smirnov test assumes that we test against a
distribution with given parameters, since in the last case we
estimated mean and variance, this assumption is violated, and the
distribution of the test statistic on which the p-value is based, is
not correct.

Tails of the distribution
^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, we can check the upper tail of the distribution. We can use
the percent point function ppf, which is the inverse of the cdf
function, to obtain the critical values, or, more directly, we can use
the inverse of the survival function

    >>> crit01, crit05, crit10 = stats.t.ppf([1-0.01, 1-0.05, 1-0.10], 10)
    >>> print 'critical values from ppf at 1%%, 5%% and 10%% %8.4f %8.4f %8.4f'% (crit01, crit05, crit10)
    critical values from ppf at 1%, 5% and 10%   2.7638   1.8125   1.3722
    >>> print 'critical values from isf at 1%%, 5%% and 10%% %8.4f %8.4f %8.4f'% tuple(stats.t.isf([0.01,0.05,0.10],10))
    critical values from isf at 1%, 5% and 10%   2.7638   1.8125   1.3722

    >>> freq01 = np.sum(x>crit01) / float(n) * 100
    >>> freq05 = np.sum(x>crit05) / float(n) * 100
    >>> freq10 = np.sum(x>crit10) / float(n) * 100
    >>> print 'sample %%-frequency at 1%%, 5%% and 10%% tail %8.4f %8.4f %8.4f'% (freq01, freq05, freq10)
    sample %-frequency at 1%, 5% and 10% tail   1.4000   5.8000  10.5000

In all three cases, our sample has more weight in the top tail than the
underlying distribution.
We can briefly check a larger sample to see if we get a closer match. In this
case the empirical frequency is quite close to the theoretical probability,
but if we repeat this several times the fluctuations are still pretty large.

    >>> freq05l = np.sum(stats.t.rvs(10, size=10000) > crit05) / 10000.0 * 100
    >>> print 'larger sample %%-frequency at 5%% tail %8.4f'% freq05l
    larger sample %-frequency at 5% tail   4.8000

We can also compare it with the tail of the normal distribution, which
has less weight in the tails:

    >>> print 'tail prob. of normal at 1%%, 5%% and 10%% %8.4f %8.4f %8.4f'% \
    ...       tuple(stats.norm.sf([crit01, crit05, crit10])*100)
    tail prob. of normal at 1%, 5% and 10%   0.2857   3.4957   8.5003

The chisquare test can be used to test, whether for a finite number of bins,
the observed frequencies differ significantly from the probabilites of the
hypothesized distribution.

    >>> quantiles = [0.0, 0.01, 0.05, 0.1, 1-0.10, 1-0.05, 1-0.01, 1.0]
    >>> crit = stats.t.ppf(quantiles, 10)
    >>> print crit
    [       -Inf -2.76376946 -1.81246112 -1.37218364  1.37218364  1.81246112
      2.76376946         Inf]
    >>> n_sample = x.size
    >>> freqcount = np.histogram(x, bins=crit)[0]
    >>> tprob = np.diff(quantiles)
    >>> nprob = np.diff(stats.norm.cdf(crit))
    >>> tch, tpval = stats.chisquare(freqcount, tprob*n_sample)
    >>> nch, npval = stats.chisquare(freqcount, nprob*n_sample)
    >>> print 'chisquare for t:      chi2 = %6.3f pvalue = %6.4f' % (tch, tpval)
    chisquare for t:      chi2 =  2.300 pvalue = 0.8901
    >>> print 'chisquare for normal: chi2 = %6.3f pvalue = %6.4f' % (nch, npval)
    chisquare for normal: chi2 = 64.605 pvalue = 0.0000

We see that the standard normal distribution is clearly rejected while the
standard t-distribution cannot be rejected. Since the variance of our sample
differs from both standard distribution, we can again redo the test taking
the estimate for scale and location into account.

The fit method of the distributions can be used to estimate the parameters
of the distribution, and the test is repeated using probabilites of the
estimated distribution.

    >>> tdof, tloc, tscale = stats.t.fit(x)
    >>> nloc, nscale = stats.norm.fit(x)
    >>> tprob = np.diff(stats.t.cdf(crit, tdof, loc=tloc, scale=tscale))
    >>> nprob = np.diff(stats.norm.cdf(crit, loc=nloc, scale=nscale))
    >>> tch, tpval = stats.chisquare(freqcount, tprob*n_sample)
    >>> nch, npval = stats.chisquare(freqcount, nprob*n_sample)
    >>> print 'chisquare for t:      chi2 = %6.3f pvalue = %6.4f' % (tch, tpval)
    chisquare for t:      chi2 =  1.577 pvalue = 0.9542
    >>> print 'chisquare for normal: chi2 = %6.3f pvalue = %6.4f' % (nch, npval)
    chisquare for normal: chi2 = 11.084 pvalue = 0.0858

Taking account of the estimated parameters, we can still reject the
hypothesis that our sample came from a normal distribution (at the 5% level),
but again, with a p-value of 0.95, we cannot reject the t distribution.



Special tests for normal distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since the normal distribution is the most common distribution in statistics,
there are several additional functions available to test whether a sample
could have been drawn from a normal distribution

First we can test if skew and kurtosis of our sample differ significantly from
those of a normal distribution:

    >>> print 'normal skewtest teststat = %6.3f pvalue = %6.4f' % stats.skewtest(x)
    normal skewtest teststat =  2.785 pvalue = 0.0054
    >>> print 'normal kurtosistest teststat = %6.3f pvalue = %6.4f' % stats.kurtosistest(x)
    normal kurtosistest teststat =  4.757 pvalue = 0.0000

These two tests are combined in the normality test

    >>> print 'normaltest teststat = %6.3f pvalue = %6.4f' % stats.normaltest(x)
    normaltest teststat = 30.379 pvalue = 0.0000

In all three tests the p-values are very low and we can reject the hypothesis
that the our sample has skew and kurtosis of the normal distribution.

Since skew and kurtosis of our sample are based on central moments, we get
exactly the same results if we test the standardized sample:

    >>> print 'normaltest teststat = %6.3f pvalue = %6.4f' % \
    ...                      stats.normaltest((x-x.mean())/x.std())
    normaltest teststat = 30.379 pvalue = 0.0000

Because normality is rejected so strongly, we can check whether the
normaltest gives reasonable results for other cases:

    >>> print 'normaltest teststat = %6.3f pvalue = %6.4f' % stats.normaltest(stats.t.rvs(10, size=100))
    normaltest teststat =  4.698 pvalue = 0.0955
    >>> print 'normaltest teststat = %6.3f pvalue = %6.4f' % stats.normaltest(stats.norm.rvs(size=1000))
    normaltest teststat =  0.613 pvalue = 0.7361

When testing for normality of a small sample of t-distributed observations
and a large sample of normal distributed observation, then in neither case
can we reject the null hypothesis that the sample comes from a normal
distribution. In the first case this is because the test is not powerful
enough to distinguish a t and a normally distributed random variable in a
small sample.


Comparing two samples
---------------------

In the following, we are given two samples, which can come either from the
same or from different distribution, and we want to test whether these
samples have the same statistical properties.

Comparing means
^^^^^^^^^^^^^^^

Test with sample with identical means:

    >>> rvs1 = stats.norm.rvs(loc=5, scale=10, size=500)
    >>> rvs2 = stats.norm.rvs(loc=5, scale=10, size=500)
    >>> stats.ttest_ind(rvs1, rvs2)
    (-0.54890361750888583, 0.5831943748663857)


Test with sample with different means:

    >>> rvs3 = stats.norm.rvs(loc=8, scale=10, size=500)
    >>> stats.ttest_ind(rvs1, rvs3)
    (-4.5334142901750321, 6.507128186505895e-006)



Kolmogorov-Smirnov test for two samples ks_2samp
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the example where both samples are drawn from the same distribution,
we cannot reject the null hypothesis since the pvalue is high

    >>> stats.ks_2samp(rvs1, rvs2)
    (0.025999999999999995, 0.99541195173064878)

In the second example, with different location, i.e. means, we can
reject the null hypothesis since the pvalue is below 1%

    >>> stats.ks_2samp(rvs1, rvs3)
    (0.11399999999999999, 0.0027132103661283141)
