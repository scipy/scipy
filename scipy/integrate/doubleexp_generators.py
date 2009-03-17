from mpmath import *

# generators for abcissas and weights for double exponential
# transformations.  Note that w = d ab(x) / dx (except for
# abw_0_inf_expfalloff, where the weights include an exponential
# falloff term.
#
# The proper way to use these for integration is to transform the
# function to be integrated to match one of the domains below, and
# then compute the transformed integral as:
# 
# sum = 0.0
# for x in numpy.linspace(-X, X, steps):
#    ab, weight = abw_???(x)
#    sum += f(inverse_transform(ab)) * weight * 2 * X / steps
# 
# The resulting sum then must be adjusted for the transformation of
# the domain.
#
# In practice, weight is quite small and the abcissae are very close
# to the integration bounds (or very large if the bounds are infinite)
# for X = 3 or 4.
# 
# See "Discovery of the Double Exponential Transformation and its
# Developments", Mori, 2005.


# integrals from -1 to 1
def abw_m1_1(x):
    mp.dps=500
    ab = tanh((pi / 2) * sinh(x))
    w = (pi/2) * cosh(x) / (cosh((pi/2) * sinh(x)))**2
    return ab, w

# integrals from 0 to infinity
def abw_0_inf(x):
    mp.dps=500
    ab = exp((pi/2) * sinh(x))
    w = (pi/2) * exp((pi/2) * sinh(x)) * cosh(x)
    return ab, w

# integrals from -infinity to infinity
def abw_minf_inf(x):
    mp.dps=500
    ab = sinh((pi/2) * sinh(2))
    w = (pi/2) * cosh(x) * cosh((pi/2) * sinh(x))
    return ab, w

# integrals from 0 to infinity of the form exp(-x) * f(x).  
#
# NB: the exp(-x) term is included in the weights directly.
def abw_0_inf_expfalloff(x):
    mp.dps=500
    ab = exp(x - exp(-x))
    w = exp(- x - exp(-x)) * (1 + exp(x))
    return ab, w
