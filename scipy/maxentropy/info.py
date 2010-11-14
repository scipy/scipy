"""
Routines for fitting maximum entropy models
===========================================

Contains two classes for fitting maximum entropy models (also known
as "exponential family" models) subject to linear constraints on the
expectations of arbitrary feature statistics.  One class, "model", is
for small discrete sample spaces, using explicit summation. The other,
"bigmodel", is for sample spaces that are either continuous (and
perhaps high-dimensional) or discrete but too large to sum over, and
uses importance sampling.  conditional Monte Carlo methods.

The maximum entropy model has exponential form

..
   p(x) = exp(theta^T f(x)) / Z(theta)

.. math::
   p\\left(x\\right)=\\exp\\left(\\frac{\\theta^{T}f\\left(x\\right)}
                                  {Z\\left(\\theta\\right)}\\right)

with a real parameter vector theta of the same length as the feature
statistic f(x), For more background, see, for example, Cover and
Thomas (1991), *Elements of Information Theory*.

See the file bergerexample.py for a walk-through of how to use these
routines when the sample space is small enough to be enumerated.

See bergerexamplesimulated.py for a a similar walk-through using
simulation.

Copyright: Ed Schofield, 2003-2006
License: BSD-style (see LICENSE.txt in main source directory)

Modules
-------

.. autosummary::
   :toctree: generated/

   maxentropy -
   maxentutils -

Classes
-------

.. autosummary::
   :toctree: generated/

   DivergenceError -
   basemodel -
   bigmodel -
   conditionalmodel -
   model -

Functions
---------

.. autosummary::
   :toctree: generated/

   arrayexp -
   arrayexpcomplex -
   columnmeans -
   columnvariances -
   densefeaturematrix -
   densefeatures -
   dotprod -
   flatten -
   innerprod -
   innerprodtranspose -
   logsumexp -
   logsumexp_naive -
   robustlog -
   rowmeans -
   sample_wr -
   sparsefeaturematrix -
   sparsefeatures -

Objects
-------

.. autosummary::
   :toctree: generated/

   division -

"""

postpone_import = 1
depends = ['optimize']
