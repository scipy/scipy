================================================
Maximum entropy models (:mod:`scipy.maxentropy`)
================================================

.. automodule:: scipy.maxentropy

Models
======
.. autoclass:: scipy.maxentropy.basemodel

.. autosummary::
   :toctree: generated/

   basemodel.beginlogging
   basemodel.endlogging
   basemodel.clearcache
   basemodel.crossentropy
   basemodel.dual
   basemodel.fit
   basemodel.grad
   basemodel.log
   basemodel.logparams
   basemodel.normconst
   basemodel.reset
   basemodel.setcallback
   basemodel.setparams
   basemodel.setsmooth

.. autoclass:: scipy.maxentropy.model

.. autosummary::
   :toctree: generated/

   model.expectations
   model.lognormconst
   model.logpmf
   model.pmf_function
   model.setfeaturesandsamplespace

.. autoclass:: scipy.maxentropy.bigmodel

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

.. autoclass:: scipy.maxentropy.conditionalmodel

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
