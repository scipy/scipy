""" Statistical Functions

This module contains a large number of probability distributions as
well as a growing library of statistical functions.

Each included distribution is an instance of a class.  The instance
name is given each instance has the following methods:

rvs   --  random variates with the distribution (also available by calling the instance)
pdf   --  probability density function
cdf   --  cummulative distribution function
sf    --  survival function (1.0 - cdf)
ppf   --  percent-point function (inverse of cdf)
isf   --  inverse survival function
stats --  mean, variance, and optionally skew and kurtosis

The distributions available with the above methods are:

CONTINUOUS  (Total == 81 distributions)
===========
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
exponweib       --  Exponentiated Weibull
exponpow        --  Exponential Power
fatiguelife     --  Fatigue Life (Birnbaum-Sanders)
foldcauchy      --  Folded Cauchy
f               --  F (Snecdor F)
foldnorm        --  Folded Normal
frechet_r       --  Frechet Right Sided, Extreme Value Type II (Extreme LB) or weibull_min
frechet_l       --  Frechet Left Sided, Weibull_max
genlogistic     --  Generalized Logistic
genpareto       --  Generalized Pareto
genexpon        --  Generalized Exponential
genextreme      --  Generalized Extreme Value
gausshyper      --  Gauss Hypergeometric
gamma           --  Gamma 
gengamma        --  Generalized gamma
genhalflogistic --  Generalized Half Logistic
gompertz        --  Gompertz (Truncated Gumbel)
gumbel_r        --  Right Sided Gumbel, Log-Weibull, Fisher-Tippett, Extreme Value Type I
gumbel_l        --  Left Sided Gumbel, etc.
halfcauchy      --  Half Cauchy
halflogistic    --  Half Logistic
halfnorm        --  Half Normal
hypsecant       --  Hyperbolic Secant
invgamma        --  Inverse Gamma
invnorm         --  Inverse Normal
invweibull      --  Inverse Weibull
johnsonsb       --  Johnson SB
johnsonsu       --  Johnson SU
laplace         --  Laplace
logistic        --  Logistic
loggamma        --  Log-Gamma
loglaplace      --  Log-Laplace (Log Double Exponential)
lognorm         --  Log-Normal
gilbrat         --  Gilbrat
lomax           --  Lomax (Pareto of the second kind)
maxwell         --  Maxwell
mielke          --  Mielke's Beta-Kappa
nakagami        --  Nakagami
ncx2            --  Non-central chi-squared
ncf             --  Non-central F
t               --  Student's T
nct             --  Non-central Student's T
pareto          --  Pareto
powerlaw        --  Power-function
powerlognorm    --  Power log normal
powernorm       --  Power normal
rdist           --  R distribution
reciprocal      --  Reciprocal
rayleigh        --  Rayleigh
rice            --  Rice
recipinvgauss   --  Reciprocal Inverse Gaussian
semicircular    --  Semicircular
triang          --  Triangular
truncexpon      --  Truncated Exponential
truncnorm       --  Truncated Normal
tukeylambda     --  Tukey-Lambda
uniform         --  Uniform
von_mises       --  Von-Mises (Circular)
wald            --  Wald
weibull_min     --  Minimum Weibull (see Frechet)
weibull_max     --  Maximum Weibull (see Frechet)
wrapcauchy      --  Wrapped Cauchy
ksone           --  Kolmogorov-Smirnov one-sided (no stats)
kstwobign       --  Kolmogorov-Smirnov two-sided test for Large N (no stats)


DISCRETE    (Total == 10 distributions)
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
dlaplace        --  Discrete Laplacian

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
wilcoxon
kruskal
friedmanchisquare

ansari
bartlett
levene
shapiro
anderson
binom_test
fligner
mood
oneway


glm
anova

Plot-tests

probplot
ppcc_max
ppcc_plot

For many more stat related functions install the software R and the 
interface package rpy.
"""

import pstats
from stats import *
from distributions import *
from rv import *
from morestats import *

try:  # use R functions if installed.
    import rpy
    from rfuncs import *
except ImportError:
    pass


#---- testing ----#
def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level=level))
    return runner

def test_suite(level=1):
    import scipy_test.testing
    import scipy.stats
    this_mod = scipy.stats
    return scipy_test.testing.harvest_test_suites(this_mod,level=level)
