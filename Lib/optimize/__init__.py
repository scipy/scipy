""" Optimization Tools

 A collection of general-purpose optimization routines.

   fmin --      Nelder-Mead Simplex algorithm
                (uses only function calls)
   fminBFGS --  Quasi-Newton method (can use function and gradient)
   fminNCG --   Line-search Newton Conjugate Gradient (can use
                function, gradient and hessian).
   leastsq --   Minimize the sum of squares of M equations in
                 N unknowns given a starting estimate.
   fsolve --    Non-linear equation solver.                 
"""

from optimize import fmin, fminBFGS, fminNCG
from minpack import fsolve, leastsq

