"""
=================================================
Orthogonal distance regression (:mod:`scipy.odr`)
=================================================

.. currentmodule:: scipy.odr

Package Content
===============

.. autosummary::
   :toctree: generated/

   odr -- Perform orthogonal distance regression

   ODR           -- Gathers all info & manages the main fitting routine.
   Data          -- Stores the data to fit.
   Model         -- Stores information about the function to be fit.
   Output
   RealData      -- Weights as actual std. dev.s and/or covariances.

   odr_error
   odr_stop

Usage information
=================

Introduction
------------

Why Orthogonal Distance Regression (ODR)?  Sometimes one has
measurement errors in the explanatory (a.k.a., "independent")
variable(s), not just the response (a.k.a., "dependent") variable(s).
Ordinary Least Squares (OLS) fitting procedures treat the data for
explanatory variables as fixed, i.e., not subject to error of any kind.
Furthermore, OLS procedures require that the response variables be an
explicit function of the explanatory variables; sometimes making the
equation explicit is impractical and/or introduces errors.  ODR can
handle both of these cases with ease, and can even reduce to the OLS
case if that is sufficient for the problem.

ODRPACK is a FORTRAN-77 library for performing ODR with possibly
non-linear fitting functions.  It uses a modified trust-region
Levenberg-Marquardt-type algorithm [1]_ to estimate the function
parameters.  The fitting functions are provided by Python functions
operating on NumPy arrays.  The required derivatives may be provided
by Python functions as well, or may be estimated numerically.  ODRPACK
can do explicit or implicit ODR fits, or it can do OLS.  Input and
output variables may be multi-dimensional.  Weights can be provided to
account for different variances of the observations, and even
covariances between dimensions of the variables.

odr provides two interfaces: a single function, and a set of
high-level classes that wrap that function; please refer to their
docstrings for more information.  While the docstring of the function
odr does not have a full explanation of its arguments, the classes do,
and arguments of the same name usually have the same requirements.
Furthermore, the user is urged to at least skim the `ODRPACK User's
Guide <http://docs.scipy.org/doc/external/odrpack_guide.pdf>`_ -
"Know Thy Algorithm."

Use
---

See the docstrings of `odr.odrpack` and the functions and classes for
usage instructions.  The ODRPACK User's Guide (linked above) is also
quite helpful.

References
----------
.. [1] P. T. Boggs and J. E. Rogers, "Orthogonal Distance Regression,"
   in "Statistical analysis of measurement error models and
   applications: proceedings of the AMS-IMS-SIAM joint summer research
   conference held June 10-16, 1989," Contemporary Mathematics,
   vol. 112, pg. 186, 1990.

"""
# version: 0.7
# author: Robert Kern <robert.kern@gmail.com>
# date: 2006-09-21

from odrpack import *
from models import *

__all__ = filter(lambda s: not s.startswith('_'), dir())

from numpy.testing import Tester
test = Tester().test
