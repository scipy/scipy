"""optimize.py

A collection of general-purpose optimization routines using Numeric

fmin        ---      Nelder-Mead Simplex algorithm (uses only function calls)
fminBFGS    ---      Quasi-Newton method (uses function and gradient)
fminNCG     ---      Line-search Newton Conjugate Gradient (uses function, 
                     gradient and hessian (if it's provided))

"""

from optimize import fmin, fminBFGS, fminNCG
