# This script generates the Python code for the dictionary
# of polynomial coefficients that is used to implement the method
# _munp(order, a) for the skewnorm distribution in scipy.stats.
#
# This script is not part of the build system.  This script was
# run, and the small snippet of code that was generated was copied
# by hand into _continuous_distns.py in scipy/stats and then edited
# a bit for PEP 8 compliance.

from sympy import symbols, simplify, exp, sqrt, diff, pi, poly
from sympy.stats import Normal, cdf


num_moments = 20

x, t, δ = symbols('x t δ')

X = Normal('x', 0, 1)

# M is the moment generating function for the skew-normal distribution.
# See, for example, https://en.wikipedia.org/wiki/Skew_normal_distribution;
# here we set the location xi to 0 and the scale omega to 1.
M = 2*exp(t**2/2)*cdf(X)(δ*t)

moments = [1]
dM = M
for k in range(1, num_moments):
    dM = simplify(diff(dM, t))
    mom = simplify(dM.subs(t, 0))
    moments.append(mom)
    print(k, mom)


with open('skewnorm_odd_moment_coefficients.py', 'w') as f:
    f.write('from numpy.polynomial import Polynomial\n')
    f.write('\n\n')
    f.write('_skewnorm_odd_moments = {\n')
    for k in range(1, num_moments, 2):
        coeffs = poly(moments[k]*sqrt(pi)/sqrt(2)).all_coeffs()[::-1][1::2]
        f.write(f'    {k}: Polynomial({coeffs}),\n')
    f.write('}\n\n\n')
