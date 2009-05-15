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

from mpmath import *
mp.dps=500

# integrals from -1 to 1
def abw_m1_1(x):
    ab = tanh((pi / 2) * sinh(x))
    w = (pi/2) * cosh(x) / (cosh((pi/2) * sinh(x)))**2
    return ab, w

# integrals from 0 to infinity
def abw_0_inf(x):
    ab = exp((pi/2) * sinh(x))
    w = (pi/2) * exp((pi/2) * sinh(x)) * cosh(x)
    return ab, w

# integrals from -infinity to infinity
def abw_minf_inf(x):
    ab = sinh((pi/2) * sinh(2))
    w = (pi/2) * cosh(x) * cosh((pi/2) * sinh(x))
    return ab, w

# integrals from 0 to infinity of the form exp(-x) * f(x).  
#
# NB: the exp(-x) term is included in the weights directly.
def abw_0_inf_expfalloff(x):
    ab = exp(x - exp(-x))
    w = exp(- x - exp(-x)) * (1 + exp(x))
    return ab, w

# print out the _abscissas_raw and _weights_raw for use in
# doubleexp.py
def print_abw_raw(lo=0.0, hi=3.155, initial_subdiv=4, levels=9, abw_func=abw_m1_1):
    abcissas = []
    weights = []
    step = (mpf(hi) - mpf(lo)) / initial_subdiv
    pts = [mpf(lo) + i * step for i in range(initial_subdiv + 1)]
    prev_pts = []
    for k in range(levels):
        a, w = zip(*[abw_func(p) for p in pts])
        abcissas += [a]
        weights += [w]
        prev_pts = prev_pts + pts
        prev_pts.sort()
        step = step / 2
        pts = [p + step for p in prev_pts[:-1]]

    print "_initial_step_size = ", float((mpf(hi) - mpf(lo)) / initial_subdiv)
    print "_abcissas_raw = ", [[float(v) for v in a] for a in abcissas]
    print "_weights_raw = ", [[float(v) for v in w] for w in weights]

if __name__ == '__main__':    
    print_abw_raw()
