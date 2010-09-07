=====================================================
Optimization and root finding (:mod:`scipy.optimize`)
=====================================================

.. module:: scipy.optimize

Optimization
============

General-purpose
---------------

.. autosummary::
   :toctree: generated/

   fmin
   fmin_powell
   fmin_cg
   fmin_bfgs
   fmin_ncg
   leastsq


Constrained (multivariate)
--------------------------

.. autosummary::
   :toctree: generated/

   fmin_l_bfgs_b
   fmin_tnc
   fmin_cobyla
   fmin_slsqp
   nnls

Global
------

.. autosummary::
   :toctree: generated/

   anneal
   brute

Scalar function minimizers
--------------------------

.. autosummary::
   :toctree: generated/

   fminbound
   golden
   bracket
   brent

Fitting
=======

.. autosummary::
   :toctree: generated/

   curve_fit

Root finding
============

Scalar functions
----------------

.. autosummary::
   :toctree: generated/

   brentq
   brenth
   ridder
   bisect
   newton

Fixed point finding:

.. autosummary::
   :toctree: generated/

   fixed_point

Multidimensional
----------------

.. toctree::
   :maxdepth: 1

   optimize.nonlin

General nonlinear solvers:

.. autosummary::
   :toctree: generated/

   fsolve
   broyden1
   broyden2

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


Utility Functions
=================

.. autosummary::
   :toctree: generated/

   line_search
   check_grad
