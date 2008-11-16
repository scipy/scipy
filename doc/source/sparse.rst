=====================================
Sparse matrices (:mod:`scipy.sparse`)
=====================================

.. automodule:: scipy.sparse


Sparse matrix classes
=====================

.. autosummary::
   :toctree: generated/

   csc_matrix
   csr_matrix
   bsr_matrix
   lil_matrix
   dok_matrix
   coo_matrix
   dia_matrix


Functions
=========

Building sparse matrices:

.. autosummary::
   :toctree: generated/

   eye
   identity
   kron
   kronsum
   lil_eye
   lil_diags
   spdiags
   tril
   triu
   bmat
   hstack
   vstack

Identifying sparse matrices:

.. autosummary::
   :toctree: generated/

   issparse
   isspmatrix
   isspmatrix_csc
   isspmatrix_csr
   isspmatrix_bsr
   isspmatrix_lil
   isspmatrix_dok
   isspmatrix_coo
   isspmatrix_dia

Exceptions
==========

.. autoexception:: SparseEfficiencyWarning

.. autoexception:: SparseWarning
