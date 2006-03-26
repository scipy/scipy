#!/usr/bin/env python

""" Example use of the maximum entropy module fit a model on a continuous
    sample space using simulation.
    
    We fit a model in 4 dimensions with constraints on the first two
    moments in each dimension.  The maximum entropy model subject to
    these constraints will be a 4-dimensional Gaussian distribution with
    diagonal covariance matrix.
    
    Computing the moments requires integration.  This could be done
    numerically here, but it might be quite expensive with 4 dimensions,
    and with any more dimensions numerical integration would quickly
    become infeasible.  The 'bigmodel' class allows the moments to be
    estimated with importance sampling instead, which, like other Monte
    Carlo methods, scales very well to high-dimensional spaces.
"""

__author__  =  'Ed Schofield'
__version__ =  '2.1'


import sys
from scipy import maxentropy
from numpy import rand

try:
    algorithm = sys.argv[1]
except IndexError:
    algorithm = 'CG'
else:
    assert algorithm in ['CG', 'BFGS', 'LBFGSB', 'Powell', 'Nelder-Mead']

def identity(x):
    return x

def square(x):
    return x**2

f = [identity, square]

model = maxentropy.bigmodel()

dims = 4

# Now set the desired feature expectations
K = [0.0, 1.0, 2.0, 3.0, 1.0, 3.0, 6.0, 9.0]
# mean is constrained to be at point (0, 1, 2, 3)  in R^4
# variances (diagonal covariance elements) are constrained to be (1, 3, 6, 9)

# Define a uniform auxiliary distribution on [-20, 20] for sampling
def featuresampler(batchsize=1000):
    lower = -20
    upper = 20
    while True:
        x = numpy.rand(dims, batchsize) * (upper - lower) - lower
        # log pdf value of this point x
        logprob = numpy.log(numpy.ones(dims, float) / (upper - lower))
        # Compute f(x)
        fx = numpy.array(numpy.flatten([fi(x) for fi in f]))
        yield fx, logprob

# NOT FINISHED!


