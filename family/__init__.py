'''
This module contains the one-parameter exponential families used
for fitting GLMs and GAMs.

These families are described in

   P. McCullagh and J. A. Nelder.  "Generalized linear models."
   Monographs on Statistics and Applied Probability.
   Chapman & Hall, London, 1983.

'''

from scipy.stats.models.family.family import Gaussian, Family, \
     Poisson, Gamma, InverseGaussian, Binomial
