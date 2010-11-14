"""
Optimization Tools
==================

General-purpose Optimization Routines
-------------------------------------

.. autosummary::
   :toctree: generated/

   fmin - Nelder-Mead Simplex algorithm #
   fmin_powell - Powell's (modified) level set method #
   fmin_cg - Non-linear (Polak-Ribiere) conjugate gradient algorithm ##
   fmin_bfgs - Quasi-Newton method (Broydon-Fletcher-Goldfarb-Shanno) ##
   fmin_ncg - Line-search Newton Conjugate Gradient ###
   leastsq - Minimize the sum of squares of M equations in N unknowns

Constrained Optimizers (Multivariate)
-------------------------------------

.. autosummary::
   :toctree: generated/

   fmin_l_bfgs_b - Zhu, Byrd, and Nocedal's constrained optimizer %
   fmin_tnc - Truncated Newton code %%
   fmin_cobyla - Constrained optimization by linear approximation
   fmin_slsqp - Minimization using sequential least-squares programming
   nnls - Linear least-squares problem with non-negativity constraint

Global Optimizers
-----------------

.. autosummary::
   :toctree: generated/

   anneal - Simulated annealing
   brute - Brute force searching optimizer

Scalar Function Minimizers
--------------------------

.. autosummary::
   :toctree: generated/

   fminbound - Bounded minimization of a scalar function
   brent - 1-D function minimization using Brent method
   golden - 1-D function minimization using Golden Section method
   bracket - Bracket a minimum, given two starting points

General-purpose Root-finding Routines
-------------------------------------

.. autosummary::
   :toctree: generated/

   fsolve - Non-linear multi-variable equation solver

Scalar Function Solvers
-----------------------

.. autosummary::
   :toctree: generated/

   brentq - quadratic interpolation Brent method
   brenth - Brent method, modified by Harris with hyperbolic extrapolation
   ridder - Ridder's method
   bisect - Bisection method
   newton - Secant method or Newton's method
   fixed_point - Single-variable fixed-point solver

General-purpose Non-linear Multidimensional Solvers
---------------------------------------------------

.. autosummary::
   :toctree: generated/

   broyden1 - Broyden's first method $
   broyden2 - Broyden's second method $$
   broyden3 - Broyden's third method $$$
   broyden_generalized - Generalized Broyden's method &
   anderson - Extended Anderson method &&
   anderson2 - The Anderson method &&&

Utility Functions
-----------------

.. autosummary::
   :toctree: generated/

   line_search - Return a step that satisfies the strong Wolfe conditions
   check_grad - Check the supplied derivative using finite differences

Related Software
----------------

OpenOpt - A BSD-licensed optimization framework (see `<http://openopt.org>`_) that includes: a number of constrained and
unconstrained solvers from and beyond the scipy.optimize module; unified
text and graphical output of convergence information; and automatic
differentiation.

Notes
-----
# Uses only function calls

## Can use function and gradient

### Can use function, gradient, and Hessian

% If you use fmin_l_bfgs_b, please cite Zhu, Byrd, and Nocedal's papers; see the function's docstring for references.

%% Originally written by Stephen Nash, adapted to C by Jean-Sebastien Roy.

$ broyden1 is a quasi-Newton-Raphson method for updating an approximate
Jacobian and then inverting it.

$$ broyden2 is the same as broyden1, but updates the inverse Jacobian
directly.

$$$ broyden3 is the same as broyden2, but instead of directly computing
the inverse Jacobian, it remembers how to construct it using vectors, and
when computing inv(J)*F, it uses those vectors to compute this product,
thus avoding the expensive NxN matrix multiplication.

& broyden_generalized is the same as broyden2, but instead of
approximating the full NxN Jacobian, it constructs it at every iteration
in a way that avoids the NxN matrix multiplication.  This is not as
precise as broyden3.

&& anderson is the same as broyden_generalized, but (w_0**2)*I is added
before taking inversion to improve the stability.

&&& anderson2 is the same as anderson, but formulated differently.

"""
#
# optimize - Optimization Tools
#

from info import __doc__

from optimize import *
from minpack import *
from zeros import *
from anneal import *
from lbfgsb import fmin_l_bfgs_b
from tnc import fmin_tnc
from cobyla import fmin_cobyla
from nonlin import *
from slsqp import fmin_slsqp
from nnls import nnls

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
