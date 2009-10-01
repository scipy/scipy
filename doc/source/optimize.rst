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

.. autosummary::
   :toctree: generated/

   fsolve

Scalar function solvers
-----------------------

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

General-purpose nonlinear (multidimensional)
--------------------------------------------

.. autosummary::
   :toctree: generated/

   broyden1
   broyden2
   broyden3
   broyden_generalized
   anderson
   anderson2

Utility Functions
=================

.. autosummary::
   :toctree: generated/

   line_search
   check_grad
