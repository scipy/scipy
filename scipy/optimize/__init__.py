"""
=====================================================
Optimization and root finding (:mod:`scipy.optimize`)
=====================================================

.. currentmodule:: scipy.optimize

Optimization
============

General-purpose
---------------

.. autosummary::
   :toctree: generated/

   minimize - Unified interface for minimizers of multivariate functions
   fmin - Nelder-Mead Simplex algorithm
   fmin_powell - Powell's (modified) level set method
   fmin_cg - Non-linear (Polak-Ribiere) conjugate gradient algorithm
   fmin_bfgs - Quasi-Newton method (Broydon-Fletcher-Goldfarb-Shanno)
   fmin_ncg - Line-search Newton Conjugate Gradient
   leastsq - Minimize the sum of squares of M equations in N unknowns

Constrained (multivariate)
--------------------------

.. autosummary::
   :toctree: generated/

   fmin_l_bfgs_b - Zhu, Byrd, and Nocedal's constrained optimizer
   fmin_tnc - Truncated Newton code
   fmin_cobyla - Constrained optimization by linear approximation
   fmin_slsqp - Minimization using sequential least-squares programming
   nnls - Linear least-squares problem with non-negativity constraint

Global
------

.. autosummary::
   :toctree: generated/

   anneal - Simulated annealing
   basinhopping - Basinhopping stochastic optimizer
   brute - Brute force searching optimizer

Scalar function minimizers
--------------------------

.. autosummary::
   :toctree: generated/

   minimize_scalar - Unified interface for minimizers of univariate functions
   fminbound - Bounded minimization of a scalar function
   brent - 1-D function minimization using Brent method
   golden - 1-D function minimization using Golden Section method
   bracket - Bracket a minimum, given two starting points

Rosenbrock function
-------------------

.. autosummary::
   :toctree: generated/

   rosen - The Rosenbrock function.
   rosen_der - The derivative of the Rosenbrock function.
   rosen_hess - The Hessian matrix of the Rosenbrock function.
   rosen_hess_prod - Product of the Rosenbrock Hessian with a vector.

Fitting
=======

.. autosummary::
   :toctree: generated/

   curve_fit -- Fit curve to a set of points

Root finding
============

Scalar functions
----------------
.. autosummary::
   :toctree: generated/

   brentq - quadratic interpolation Brent method
   brenth - Brent method, modified by Harris with hyperbolic extrapolation
   ridder - Ridder's method
   bisect - Bisection method
   newton - Secant method or Newton's method

Fixed point finding:

.. autosummary::
   :toctree: generated/

   fixed_point - Single-variable fixed-point solver

Multidimensional
----------------

General nonlinear solvers:

.. autosummary::
   :toctree: generated/

   root - Unified interface for nonlinear solvers of multivariate functions
   fsolve - Non-linear multi-variable equation solver
   broyden1 - Broyden's first method
   broyden2 - Broyden's second method

Large-scale nonlinear solvers:

.. autosummary::
   :toctree: generated/

   newton_krylov
   anderson

Simple iterations:

.. autosummary::
   :toctree: generated/

   excitingmixing
   linearmixing
   diagbroyden

:mod:`Additional information on the nonlinear solvers <scipy.optimize.nonlin>`

Utility Functions
=================

.. autosummary::
   :toctree: generated/

   line_search - Return a step that satisfies the strong Wolfe conditions
   check_grad - Check the supplied derivative using finite differences

   show_options - Show specific options optimization solvers

"""

from __future__ import division, print_function, absolute_import

from .optimize import *
from ._minimize import *
from ._root import *
from .minpack import *
from .zeros import *
from .anneal import *
from .lbfgsb import fmin_l_bfgs_b
from .tnc import fmin_tnc
from .cobyla import fmin_cobyla
from .nonlin import *
from .slsqp import fmin_slsqp
from .nnls import nnls
from ._basinhopping import basinhopping

__all__ = [s for s in dir() if not s.startswith('_')]
from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
