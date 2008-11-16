================================================
Maximum entropy models (:mod:`scipy.maxentropy`)
================================================

.. automodule:: scipy.maxentropy


Models
======

.. autoclass:: model

.. autosummary::
   :toctree: generated/

   model.beginlogging
   model.endlogging
   model.clearcache
   model.crossentropy
   model.dual
   model.fit
   model.grad
   model.log
   model.logparams
   model.normconst
   model.reset
   model.setcallback
   model.setparams
   model.setsmooth
   model.expectations
   model.lognormconst
   model.logpmf
   model.pmf_function
   model.setfeaturesandsamplespace

.. autoclass:: bigmodel

.. autosummary::
   :toctree: generated/

   bigmodel.estimate
   bigmodel.logpdf
   bigmodel.pdf
   bigmodel.pdf_function
   bigmodel.resample
   bigmodel.setsampleFgen
   bigmodel.settestsamples
   bigmodel.stochapprox
   bigmodel.test
   
.. autoclass:: conditionalmodel

.. autosummary::
   :toctree: generated/

   conditionalmodel.dual
   conditionalmodel.expectations
   conditionalmodel.fit
   conditionalmodel.lognormconst
   conditionalmodel.logpmf

Utilities
=========

.. autosummary::
   :toctree: generated/

   arrayexp
   arrayexpcomplex
   columnmeans
   columnvariances
   densefeaturematrix
   densefeatures
   dotprod
   flatten
   innerprod
   innerprodtranspose
   logsumexp
   logsumexp_naive
   robustlog
   rowmeans
   sample_wr
   sparsefeaturematrix
   sparsefeatures
