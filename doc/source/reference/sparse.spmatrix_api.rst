.. _spmatrix_api:

====================================================
Sparse Matrix Legacy Interface (:mod:`scipy.sparse`)
====================================================

.. currentmodule:: scipy.sparse

SciPy 2-D sparse matrix interface for numeric data.

.. warning::

   The SciPy sparse package is switching to an array interface, compatible
   with NumPy arrays, from the older numpy.matrix interface.  We recommend
   that you use the array objects (`csr_array`, `coo_array`, etc.).

   When using the array interface, please note that:

   - ``x * y`` performs element-wise multiplication (just like NumPy arrays).
     To make code work with both arrays and matrices, use ``x @ y`` for
     matrix multiplication.
   - Operations, like ``sum``, produce np.matrix for spmatrix while they produce
     np.array for sparray. Multiplication of np.matrix and np.array
     differ in the same way as between spmatrix and sparray.
   - Sparse arrays use array style *slicing* operations, returning scalars,
     or 1D or 2D sparse arrays according to the numpy standards for indexing.
     If you need always 2D results, use an appropriate index with sparse arrays.
     For example, ``A[:, i, None]`` or ``A[:, [i]]``.

   The construction utilities (`eye`, `kron`, `random`, `diags`, etc.)
   have appropriate replacements (see :ref:`sparse-construction-functions`).

   For more information see
   :ref:`Migration from spmatrix to sparray <migration_to_sparray>`.


Sparse matrix classes
=====================

.. autosummary::
   :toctree: generated/

   bsr_matrix - Block Sparse Row matrix
   coo_matrix - A sparse matrix in COOrdinate format
   csc_matrix - Compressed Sparse Column matrix
   csr_matrix - Compressed Sparse Row matrix
   dia_matrix - Sparse matrix with DIAgonal storage
   dok_matrix - Dictionary Of Keys based sparse matrix
   lil_matrix - Row-based list of lists sparse matrix
   spmatrix - Sparse matrix base class

Building sparse matrices
------------------------

.. autosummary::
   :toctree: generated/

   eye - Sparse MxN matrix whose k-th diagonal is all ones
   identity - Identity matrix in sparse matrix format
   diags - Return a sparse matrix from diagonals
   spdiags - Return a sparse matrix from diagonals
   bmat - Build a sparse matrix from sparse sub-blocks
   random - Random values in a given shape matrix
   rand - Random values in a given shape matrix (old interface)

**To combine matrices use the same functions as for** :ref:`combining-arrays`.

Identifying sparse matrices
---------------------------

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


