""" Statistical Functions

This module contains a large number of probability distributions as
well as a growing library of statistical functions.

Each included distribution has associated functions to compute the
following

<XXX>      --  random variates with distribution <XXX>
<XXX>pdf   --  probability density function
<XXX>cdf   --  cummulative distribution function
<XXX>sf    --  survival function (1.0 - cdf)
<XXX>ppf   --  percent-point function (inverse of cdf)
<XXX>isf   --  inverse survival function
<XXX>stats --  mean, variance, and optionally skew and kurtosis

The distributions that can be used to replace the <XXX> are

CONTINUOUS
===========
stnorm          --  Standard normal
norm            --  Normal (Gaussian)
alpha           --  Alpha
anglit          --  Anglit
arcsine         --  Arcsine
beta            --  Beta
betaprime       --  Beta Prime
bradford        --  Bradford
burr            --  Burr
fisk            --  Fisk
cauchy          --  Cauchy
chi             --  Chi
chi2            --  Chi-squared
cosine          --  Cosine
dgamma          --  Double Gamma
dweibull        --  Double Weibull
erlang          --  Erlang
expon           --  Exponential
extreme3        --  Extreme Value Type III
exponweib       --  Exponentiated Weibull
fatiguelife     --  Fatigue Life (Birnbaum-Sanders)
foldcauchy      --  Folded Cauchy
f               --  F (Snecdor F)
foldnorm        --  Folded Normal
frechet         --  Frechet or Extreme Value Type II (Exterme LB)
genlogistic     --  Generalized Logistic
genpareto       --  Generalized Pareto
genextreme      --  Generalized Extreme Value
gamma           --  Gamma 
gengamma        --  Generalized gamma
genhalflogistic --  Generalized Half Logistic
gompertz        --  Gompertz (Truncated Gumbel)
gumbel          --  Gumbel, Log-Weibull, Fisher-Tippett
halfcauchy      --  Half Cauchy
halflogistic    --  Half Logistic
halfnorm        --  Half Normal
hypsecant       --  Hyperbolic Secant
invgamma        --  Inverse Gamma
invnorm         --  Inverse Normal
invweibull      --  Inverse Weibull
laplace         --  Laplace
logistic        --  Logistic
loggamma        --  Log-Gamma
lognorm         --  Log-Normal
gilbrat         --  Gilbrat
lomax           --  Lomax (Pareto of the second kind)
maxwell         --  Maxwell
nakagami        --  Nakagami
ncx2            --  Non-central chi-squared
ncf             --  Non-central F
t               --  Student's T
nct             --  Non-central Student's T
pareto          --  Pareto
power           --  Power
reciprocal      --  Reciprocal
rayleigh        --  Rayleigh
semicircular    --  Semicircular
triang          --  Triangular
tukeylambda     --  Tukey-Lambda
uniform         --  Uniform
von_mises       --  Von-Mises (Circular)
wald            --  Wald
weibull         --  Weibull
ksone           --  Kolmogorov-Smirnov one-sided (no stats)
kstwobign       --  Kolmogorov-Smirnov two-sided test for Large N (no stats)

DISCRETE
============
binom           --  Binomial
bernoulli       --  Bernoulli
nbinom          --  Negative Binomial
geom            --  Geometric
hypergeom       --  Hypergeometric
logser          --  Logarithmic (Log-Series, Series)
poisson         --  Poisson
randint         --  Discrete Uniform
zipf            --  Zipf


Statistical Functions (adapted from Gary Strangman)

gmean
hmean
mean
cmedian
median
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

itemfreq
scoreatpercentile
percentileofscore
histogram2
histogram
cumfreq
relfreq

obrientransform
samplevar
samplestd
signaltonoise
var
std
stderr
sem
z
zs
zmap

threshold
trimboth
trim1
cov
corrcoef

f_oneway
paired
pearsonr
spearmanr
pointbiserialr
kendalltau
linregress

ttest_1samp
ttest_ind
ttest_rel
kstest
chisquare
ks_2samp
meanwhitneyu
tiecorrect
ranksums
wilcoxont
kruskalwallish
friedmanchisquare

glm
anova

Plot-tests

probplot
ppcc_max
ppcc_plot
"""

import pstats
from stats import *
from rv2 import *
from distributions import *
from rv import *
from morestats import *


#---- testing ----#
def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level=level))
    return runner

def test_suite(level=1):
    import scipy_base.testing
    import scipy.stats
    this_mod = scipy.stats
    return scipy_base.testing.harvest_test_suites(this_mod,level=level)
