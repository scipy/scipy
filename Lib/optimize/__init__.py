""" Optimization Tools

 A collection of general-purpose optimization routines.

   fmin --       Nelder-Mead Simplex algorithm
                 (uses only function calls)
   fmin_bfgs --  Quasi-Newton method (can use function and gradient)
   fmin_ncg --   Line-search Newton Conjugate Gradient (can use
                 function, gradient and hessian).
   leastsq --    Minimize the sum of squares of M equations in
                 N unknowns given a starting estimate.
   fminbound --  Bounded minimization of a scalar function.
   fsolve --     Non-linear multi-variable equation solver.   
   newton --     Single-variable function root finder.
   fixed_point -- Single-variable fixed-point solver.
   bisection --  Single-variable function root finder using bisection method.
"""

_moddict = {'optimize' : ['fmin', 'fmin_bfgs', 'fmin_ncg', 'fminbound',
                          'rosen','rosen_der', 'rosen_hess',
                          'rosen_hess_prod'],
            'minpack' : ['fsolve', 'leastsq', 'newton', 'fixed_point','bisection']
            }
__all__ = []

import scipy
scipy.somenames2all(__all__, _moddict, globals())
del scipy

