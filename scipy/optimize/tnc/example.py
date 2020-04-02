# Python TNC example
# @(#) $Jeannot: example.py,v 1.4 2004/04/02 18:51:04 js Exp $

import tnc

# A function to minimize
# Must return a tuple with the function value and the gradient (as a list)
# or None to abort the minimization


def function(x):
    f = pow(x[0],2.0)+pow(abs(x[1]),3.0)
    g = [0,0]
    g[0] = 2.0*x[0]
    g[1] = 3.0*pow(abs(x[1]),2.0)
    if x[1] < 0:
        g[1] = -g[1]
    return f, g


# Optimizer call
rc, nf, x = tnc.minimize(function, [-7, 3], [-10, 1], [10, 10])

print("After", nf, "function evaluations, TNC returned:", tnc.RCSTRINGS[rc])
print("x =", x)
print("exact value = [0, 1]")
