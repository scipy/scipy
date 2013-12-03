"""
==========================================
Statistical functions (:mod:`scipy.stats`)
==========================================

.. module:: scipy.stats

This module contains a large number of probability distributions as
well as a growing library of statistical functions.

Each included distribution is an instance of the class rv_continuous:
For each given name the following methods are available:

.. autosummary::
   :toctree: generated/

   rv_continuous
   rv_continuous.pdf
   rv_continuous.logpdf
   rv_continuous.cdf
   rv_continuous.logcdf
   rv_continuous.sf
   rv_continuous.logsf
   rv_continuous.ppf
   rv_continuous.isf
   rv_continuous.moment
   rv_continuous.stats
   rv_continuous.entropy
   rv_continuous.fit
   rv_continuous.expect

Calling the instance as a function returns a frozen pdf whose shape,
location, and scale parameters are fixed.

Similarly, each discrete distribution is an instance of the class
rv_discrete:

.. autosummary::
   :toctree: generated/

   rv_discrete
   rv_discrete.rvs
   rv_discrete.pmf
   rv_discrete.logpmf
   rv_discrete.cdf
   rv_discrete.logcdf
   rv_discrete.sf
   rv_discrete.logsf
   rv_discrete.ppf
   rv_discrete.isf
   rv_discrete.stats
   rv_discrete.moment
   rv_discrete.entropy
   rv_discrete.expect

Continuous distributions
========================

.. autosummary::
   :toctree: generated/

   alpha             -- Alpha
   anglit            -- Anglit
   arcsine           -- Arcsine
   beta              -- Beta
   betaprime         -- Beta Prime
   bradford          -- Bradford
   burr              -- Burr
   cauchy            -- Cauchy
   chi               -- Chi
   chi2              -- Chi-squared
   cosine            -- Cosine
   dgamma            -- Double Gamma
   dweibull          -- Double Weibull
   erlang            -- Erlang
   expon             -- Exponential
   exponweib         -- Exponentiated Weibull
   exponpow          -- Exponential Power
   f                 -- F (Snecdor F)
   fatiguelife       -- Fatigue Life (Birnbaum-Saunders)
   fisk              -- Fisk
   foldcauchy        -- Folded Cauchy
   foldnorm          -- Folded Normal
   frechet_r         -- Frechet Right Sided, Extreme Value Type II (Extreme LB) or weibull_min
   frechet_l         -- Frechet Left Sided, Weibull_max
   genlogistic       -- Generalized Logistic
   genpareto         -- Generalized Pareto
   genexpon          -- Generalized Exponential
   genextreme        -- Generalized Extreme Value
   gausshyper        -- Gauss Hypergeometric
   gamma             -- Gamma
   gengamma          -- Generalized gamma
   genhalflogistic   -- Generalized Half Logistic
   gilbrat           -- Gilbrat
   gompertz          -- Gompertz (Truncated Gumbel)
   gumbel_r          -- Right Sided Gumbel, Log-Weibull, Fisher-Tippett, Extreme Value Type I
   gumbel_l          -- Left Sided Gumbel, etc.
   halfcauchy        -- Half Cauchy
   halflogistic      -- Half Logistic
   halfnorm          -- Half Normal
   hypsecant         -- Hyperbolic Secant
   invgamma          -- Inverse Gamma
   invgauss          -- Inverse Gaussian
   invweibull        -- Inverse Weibull
   johnsonsb         -- Johnson SB
   johnsonsu         -- Johnson SU
   ksone             -- Kolmogorov-Smirnov one-sided (no stats)
   kstwobign         -- Kolmogorov-Smirnov two-sided test for Large N (no stats)
   laplace           -- Laplace
   logistic          -- Logistic
   loggamma          -- Log-Gamma
   loglaplace        -- Log-Laplace (Log Double Exponential)
   lognorm           -- Log-Normal
   lomax             -- Lomax (Pareto of the second kind)
   maxwell           -- Maxwell
   mielke            -- Mielke's Beta-Kappa
   nakagami          -- Nakagami
   ncx2              -- Non-central chi-squared
   ncf               -- Non-central F
   nct               -- Non-central Student's T
   norm              -- Normal (Gaussian)
   pareto            -- Pareto
   pearson3          -- Pearson type III
   powerlaw          -- Power-function
   powerlognorm      -- Power log normal
   powernorm         -- Power normal
   rdist             -- R-distribution
   reciprocal        -- Reciprocal
   rayleigh          -- Rayleigh
   rice              -- Rice
   recipinvgauss     -- Reciprocal Inverse Gaussian
   semicircular      -- Semicircular
   t                 -- Student's T
   triang            -- Triangular
   truncexpon        -- Truncated Exponential
   truncnorm         -- Truncated Normal
   tukeylambda       -- Tukey-Lambda
   uniform           -- Uniform
   vonmises          -- Von-Mises (Circular)
   wald              -- Wald
   weibull_min       -- Minimum Weibull (see Frechet)
   weibull_max       -- Maximum Weibull (see Frechet)
   wrapcauchy        -- Wrapped Cauchy

Multivariate distributions
==========================

.. autosummary::
   :toctree: generated/

   multivariate_normal   -- Multivariate normal distribution
   dirichlet             -- Dirichlet

Discrete distributions
======================

.. autosummary::
   :toctree: generated/

   bernoulli         -- Bernoulli
   binom             -- Binomial
   boltzmann         -- Boltzmann (Truncated Discrete Exponential)
   dlaplace          -- Discrete Laplacian
   geom              -- Geometric
   hypergeom         -- Hypergeometric
   logser            -- Logarithmic (Log-Series, Series)
   nbinom            -- Negative Binomial
   planck            -- Planck (Discrete Exponential)
   poisson           -- Poisson
   randint           -- Discrete Uniform
   skellam           -- Skellam
   zipf              -- Zipf

Statistical functions
=====================

Several of these functions have a similar version in scipy.stats.mstats
which work for masked arrays.

.. autosummary::
   :toctree: generated/

   describe          -- Descriptive statistics
   gmean             -- Geometric mean
   hmean             -- Harmonic mean
   kurtosis          -- Fisher or Pearson kurtosis
   kurtosistest      --
   mode              -- Modal value
   moment            -- Central moment
   normaltest        --
   skew              -- Skewness
   skewtest          --
   tmean             -- Truncated arithmetic mean
   tvar              -- Truncated variance
   tmin              --
   tmax              --
   tstd              --
   tsem              --
   nanmean           -- Mean, ignoring NaN values
   nanstd            -- Standard deviation, ignoring NaN values
   nanmedian         -- Median, ignoring NaN values
   variation         -- Coefficient of variation

.. autosummary::
   :toctree: generated/

   cumfreq           _
   histogram2        _
   histogram         _
   itemfreq          _
   percentileofscore _
   scoreatpercentile _
   relfreq           _

.. autosummary::
   :toctree: generated/

   binned_statistic     -- Compute a binned statistic for a set of data.
   binned_statistic_2d  -- Compute a 2-D binned statistic for a set of data.
   binned_statistic_dd  -- Compute a d-D binned statistic for a set of data.

.. autosummary::
   :toctree: generated/

   obrientransform
   signaltonoise
   bayes_mvs
   sem
   zmap
   zscore

.. autosummary::
   :toctree: generated/

   sigmaclip
   threshold
   trimboth
   trim1

.. autosummary::
   :toctree: generated/

   f_oneway
   pearsonr
   spearmanr
   pointbiserialr
   kendalltau
   linregress
   theilslopes

.. autosummary::
   :toctree: generated/

   ttest_1samp
   ttest_ind
   ttest_rel
   kstest
   chisquare
   power_divergence
   ks_2samp
   mannwhitneyu
   tiecorrect
   rankdata
   ranksums
   wilcoxon
   kruskal
   friedmanchisquare
   stouffers_method
   fishers_method

.. autosummary::
   :toctree: generated/

   ansari
   bartlett
   levene
   shapiro
   anderson
   anderson_ksamp
   binom_test
   fligner
   median_test
   mood

.. autosummary::
   :toctree: generated/

   boxcox
   boxcox_normmax
   boxcox_llf

   entropy

Contingency table functions
===========================

.. autosummary::
   :toctree: generated/

   chi2_contingency
   contingency.expected_freq
   contingency.margins
   fisher_exact

Plot-tests
==========

.. autosummary::
   :toctree: generated/

   ppcc_max
   ppcc_plot
   probplot
   boxcox_normplot


Masked statistics functions
===========================

.. toctree::

   stats.mstats


Univariate and multivariate kernel density estimation (:mod:`scipy.stats.kde`)
==============================================================================

.. autosummary::
   :toctree: generated/

   gaussian_kde

For many more stat related functions install the software R and the
interface package rpy.

"""
from __future__ import division, print_function, absolute_import

from .stats import *
from .distributions import *
from .rv import *
from .morestats import *
from ._binned_statistic import *
from .kde import gaussian_kde
from . import mstats
from .contingency import chi2_contingency
from ._multivariate import *

#remove vonmises_cython from __all__, I don't know why it is included
__all__ = [s for s in dir() if not (s.startswith('_') or s.endswith('cython'))]

from numpy.testing import Tester
test = Tester().test
