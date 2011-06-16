==========================================
Statistical functions (:mod:`scipy.stats`)
==========================================

.. module:: scipy.stats

This module contains a large number of probability distributions as
well as a growing library of statistical functions.

Each included continuous distribution is an instance of the class rv_continous:

.. autosummary::
   :toctree: generated/

   rv_continuous
   rv_continuous.rvs
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
   rv_continuous.median
   rv_continuous.mean
   rv_continuous.var
   rv_continuous.std
   rv_continuous.interval

Each discrete distribution is an instance of the class rv_discrete:

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
   rv_discrete.median
   rv_discrete.mean
   rv_discrete.var
   rv_discrete.std
   rv_discrete.interval

Continuous distributions
========================

.. autosummary::
   :toctree: generated/

   norm
   alpha
   anglit
   arcsine
   beta
   betaprime
   bradford
   burr
   fisk
   cauchy
   chi
   chi2
   cosine
   dgamma
   dweibull
   erlang
   expon
   exponweib
   exponpow
   fatiguelife
   foldcauchy
   f
   foldnorm
   genlogistic
   genpareto
   genexpon
   genextreme
   gausshyper
   gamma
   gengamma
   genhalflogistic
   gompertz
   gumbel_r
   gumbel_l
   halfcauchy
   halflogistic
   halfnorm
   hypsecant
   invgamma
   invgauss
   invweibull
   johnsonsb
   johnsonsu
   laplace
   logistic
   loggamma
   loglaplace
   lognorm
   gilbrat
   lomax
   maxwell
   mielke
   nakagami
   ncx2
   ncf
   t
   nct
   pareto
   powerlaw
   powerlognorm
   powernorm
   rdist
   reciprocal
   rayleigh
   rice
   recipinvgauss
   semicircular
   triang
   truncexpon
   truncnorm
   tukeylambda
   uniform
   vonmises
   wald
   weibull_min
   weibull_max
   wrapcauchy
   ksone
   kstwobign

Discrete distributions
======================

.. autosummary::
   :toctree: generated/

   binom
   bernoulli
   nbinom
   geom
   hypergeom
   logser
   poisson
   planck
   boltzmann
   randint
   zipf
   dlaplace

Statistical functions
=====================

Several of these functions have a similar version in scipy.stats.mstats
which work for masked arrays.

.. autosummary::
   :toctree: generated/

   gmean
   hmean
   cmedian
   mode
   tmean
   tvar
   tmin
   tmax
   tstd
   tsem
   moment
   variation
   skew
   kurtosis
   describe
   skewtest
   kurtosistest
   normaltest


.. autosummary::
   :toctree: generated/

   itemfreq
   scoreatpercentile
   percentileofscore
   histogram2
   histogram
   cumfreq
   relfreq

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

   threshold
   trimboth
   trim1

.. autosummary::
   :toctree: generated/

   f_oneway
   pearsonr
   spearmanr
   fisher_exact
   pointbiserialr
   kendalltau
   linregress

.. autosummary::
   :toctree: generated/

   ttest_1samp
   ttest_ind
   ttest_rel
   kstest
   chisquare
   ks_2samp
   mannwhitneyu
   tiecorrect
   ranksums
   wilcoxon
   kruskal
   friedmanchisquare

.. autosummary::
   :toctree: generated/

   ansari
   bartlett
   levene
   shapiro
   anderson
   binom_test
   fligner
   mood
   oneway

Contingency Tables
------------------

.. autosummary::
   :toctree: generated/
   
   chi2_contingency
   expected_freq
   margins
   
General Linear Model
--------------------

.. autosummary::
   :toctree: generated/

   glm

Plot-tests
==========

.. autosummary::
   :toctree: generated/

   probplot
   ppcc_max
   ppcc_plot


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
