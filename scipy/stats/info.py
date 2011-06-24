"""
Statistical Functions
=====================

This module contains a large number of probability distributions as
well as a growing library of statistical functions.

Each included distribution is an instance of the class rv_continous.
For each given name the following methods are available.  See docstring
for rv_continuous for more information

:rvs:
   random variates with the distribution
:pdf:
   probability density function
:cdf:
   cumulative distribution function
:sf:
   survival function (1.0 - cdf)
:ppf:
   percent-point function (inverse of cdf)
:isf:
   inverse survival function
:stats:
   mean, variance, and optionally skew and kurtosis

Calling the instance as a function returns a frozen pdf whose shape,
location, and scale parameters are fixed.

Distributions
---------------

The distributions available with the above methods are:


Continuous  (Total == 81 distributions)
---------------------------------------

.. autosummary::
   :toctree: generated/

   norm              Normal (Gaussian)
   alpha             Alpha
   anglit            Anglit
   arcsine           Arcsine
   beta              Beta
   betaprime         Beta Prime
   bradford          Bradford
   burr              Burr
   cauchy            Cauchy
   chi               Chi
   chi2              Chi-squared
   cosine            Cosine
   dgamma            Double Gamma
   dweibull          Double Weibull
   erlang            Erlang
   expon             Exponential
   exponweib         Exponentiated Weibull
   exponpow          Exponential Power
   f                 F (Snecdor F)
   fatiguelife       Fatigue Life (Birnbaum-Sanders)
   fisk              Fisk
   foldcauchy        Folded Cauchy
   foldnorm          Folded Normal
   frechet_r         Frechet Right Sided, Extreme Value Type II (Extreme LB) or weibull_min
   frechet_l         Frechet Left Sided, Weibull_max
   genlogistic       Generalized Logistic
   genpareto         Generalized Pareto
   genexpon          Generalized Exponential
   genextreme        Generalized Extreme Value
   gausshyper        Gauss Hypergeometric
   gamma             Gamma
   gengamma          Generalized gamma
   genhalflogistic   Generalized Half Logistic
   gompertz          Gompertz (Truncated Gumbel)
   gumbel_r          Right Sided Gumbel, Log-Weibull, Fisher-Tippett, Extreme Value Type I
   gumbel_l          Left Sided Gumbel, etc.
   halfcauchy        Half Cauchy
   halflogistic      Half Logistic
   halfnorm          Half Normal
   hypsecant         Hyperbolic Secant
   invgamma          Inverse Gamma
   invgauss          Inverse Gaussian
   invweibull        Inverse Weibull
   johnsonsb         Johnson SB
   johnsonsu         Johnson SU
   ksone             Kolmogorov-Smirnov one-sided (no stats)
   kstwobign         Kolmogorov-Smirnov two-sided test for Large N (no stats)
   laplace           Laplace
   logistic          Logistic
   loggamma          Log-Gamma
   loglaplace        Log-Laplace (Log Double Exponential)
   lognorm           Log-Normal
   gilbrat           Gilbrat
   lomax             Lomax (Pareto of the second kind)
   maxwell           Maxwell
   mielke            Mielke's Beta-Kappa
   nakagami          Nakagami
   ncx2              Non-central chi-squared
   ncf               Non-central F
   nct               Non-central Student's T
   pareto            Pareto
   powerlaw          Power-function
   powerlognorm      Power log normal
   powernorm         Power normal
   rdist             R distribution
   reciprocal        Reciprocal
   rayleigh          Rayleigh
   rice              Rice
   recipinvgauss     Reciprocal Inverse Gaussian
   semicircular      Semicircular
   t                 Student's T
   triang            Triangular
   truncexpon        Truncated Exponential
   truncnorm         Truncated Normal
   tukeylambda       Tukey-Lambda
   uniform           Uniform
   vonmises          Von-Mises (Circular)
   wald              Wald
   weibull_min       Minimum Weibull (see Frechet)
   weibull_max       Maximum Weibull (see Frechet)
   wrapcauchy        Wrapped Cauchy


=============== ==============================================================
Discrete    (Total == 10 distributions)
==============================================================================
binom             Binomial
bernoulli         Bernoulli
nbinom            Negative Binomial
geom              Geometric
hypergeom         Hypergeometric
logser            Logarithmic (Log-Series, Series)
poisson           Poisson
planck            Planck (Discrete Exponential)
boltzmann         Boltzmann (Truncated Discrete Exponential)
randint           Discrete Uniform
zipf              Zipf
dlaplace          Discrete Laplacian
=============== ==============================================================

Statistical Functions (adapted from Gary Strangman)
-----------------------------------------------------

================= ==============================================================
gmean             Geometric mean
hmean             Harmonic mean
mean              Arithmetic mean
cmedian           Computed median
median            Median
mode              Modal value
tmean             Truncated arithmetic mean
tvar              Truncated variance
tmin              _
tmax              _
tstd              _
tsem              _
moment            Central moment
variation         Coefficient of variation
skew              Skewness
kurtosis          Fisher or Pearson kurtosis
describe          Descriptive statistics
skewtest          _
kurtosistest      _
normaltest        _
================= ==============================================================

================= ==============================================================
itemfreq          _
scoreatpercentile _
percentileofscore _
histogram2        _
histogram         _
cumfreq           _
relfreq           _
================= ==============================================================

================= ==============================================================
obrientransform   _
signaltonoise     _
bayes_mvs         _
sem               _
zmap              _
================= ==============================================================

================= ==============================================================
threshold         _
trimboth          _
trim1             _
================= ==============================================================

================= ==============================================================
f_oneway          _
paired            _
pearsonr          _
spearmanr         _
pointbiserialr    _
kendalltau        _
linregress        _
================= ==============================================================

================= ==============================================================
ttest_1samp       _
ttest_ind         _
ttest_rel         _
kstest            _
chisquare         _
ks_2samp          _
meanwhitneyu      _
tiecorrect        _
ranksums          _
wilcoxon          _
kruskal           _
friedmanchisquare _
================= ==============================================================

================= ==============================================================
ansari            _
bartlett          _
levene            _
shapiro           _
anderson          _
binom_test        _
fligner           _
mood              _
oneway            _
================= ==============================================================

================= =======================================================================
fisher_exact      Performs a Fisher exact test on a 2x2 contingency table.
chi2_contingency  Chi-square test of independence of observations in a contingency table.
expected_freq     Compute the expected frequencies from a contingency table.
margins           Return a list of the marginal sums of the array `a`.
================= =======================================================================

================= ==============================================================
glm               _
================= ==============================================================


================= ==============================================================
Plot-tests
================================================================================
probplot          _
ppcc_max          _
ppcc_plot         _
================= ==============================================================


For many more stat related functions install the software R and the
interface package rpy.

"""

postpone_import = 1
global_symbols = ['find_repeats']

depends  = ['linalg','special']
ignore = False # importing stats causes a segfault
