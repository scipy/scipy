.. _migration_to_sparray:

Migration from spmatrix to sparray
==================================

This document provides guidance for converting code from sparse *matrices*
to sparse *arrays* in ``scipy.sparse``.

The change from sparse matrices to sparse arrays mirrors conversion from
``np.matrix`` to ``np.ndarray``. Essentially we must move from an all-2D
matrix-multiplication-centric ``matrix`` object to a 1D or 2D “array”
object that supports the matrix multiplication operator and elementwise
computation.

Notation: For this guide we denote the sparse array classes generally as
``sparray`` and the sparse matrix classes ``spmatrix``. Dense numpy
arrays are denoted ``np.ndarray`` and dense matrix classes are
``np.matrix``. Supported sparse formats are denoted BSR, COO, CSC, CSR,
DIA, DOK, LIL and all formats are supported by both sparray and
spmatrix. The term ``sparse`` refers to either ``sparray`` or
``spmatrix``, while ``dense`` refers to either ``np.ndarray`` or
``np.matrix``.

Overview and big picture
------------------------

-  The constructor names ``*_matrix``, e.g. ``csr_matrix``, are changed
   to ``*_array``.
-  spmatrix ``M`` is always 2D (rows x columns) even e.g. ``M.min(axis=0)``.
   sparray ``A`` can be 1D or 2D.
   Numpy scalars are returned for full (0D) reductions, i.e. ``M.min()``.
-  Iterating over a sparray gives 1D sparrays. Iterating spmatrix gives 2D row spmatrices
-  Operators that change behavior are: ``*, @, *=, @=, **``

   -  Scalar multiplication, e.g. ``5 * A``, uses ``*``, and ``5 @ A`` is not
      implemented.
   -  sparrays use ``*`` for elementwise multiplication and ``@`` for
      matrix multiplication while spmatrices use either operator
      ``*`` or ``@`` for matrix multiplication. Either can use
      ``A.multiply(B)`` for elementwise multiplication.
   -  Scalar exponents, e.g. ``A**2``, use elementwise power for sparray and
      matrix power for spmatrix. To get matrix power for sparrays use
      ``scipy.sparse.linalg.matrix_power(A, n)``.

-  Checking the sparse type and format:

   -  ``issparse(A)`` returns ``True`` for any sparse array/matrix.
   -  ``isspmatrix(M)`` returns ``True`` for any sparse matrix.
   -  ``isspmatrix_csr(M)`` checks for a sparse matrix with specific format.
      It should be replaced with an array compatible version such as:
   -  ``issparse(A) and A.format == 'csr'`` looks for a CSR sparse
      array/matrix.

Recommended steps for migration
-------------------------------

-  First pass (leaving spmatrix in the code):

   -  In your spmatrix code, change ``*`` to ``@`` for matrix
      multiplication. Note that scalar multiplication with sparse should
      use ``*``.
   -  Matrix powers, e.g. ``M**3``, should be converted to
      ``scipy.sparse.linalg.matrix_power(A, 3)``
   -  Implement alternatives to unsupported functions/methods like
      ``A.getnnz()`` -> ``A.nnz`` (see :ref:`sparse-migration-removed-methods`
      below).
   -  Change any logic regarding ``issparse()`` and ``isspmatrix()`` as
      needed. Usually, this means replacing ``isspmatrix`` with ``issparse``,
      and ``isspmatrix_csr(G)`` with ``issparse(G) and G.format == "csr"``.
      Moreover ``isspmatrix_csr(G) or isspmatrix_csc(G)`` becomes
      ``issparse(G) and G.format in ['csr', 'csc']``.
   -  Convert all ``spdiags`` calls to ``dia_matrix``.
      See docs in :func:`spdiags<scipy.sparse.spdiags>`.
   -  Run all your tests on the resulting code. You are still using
      spmatrix, not sparray. But your code and tests are prepared for
      the change.

-  Second pass (switching to sparray):

   -  Convert construction functions like ``diags`` and ``triu`` to the
      array version (see :ref:`sparse-migration-construction` below).
   -  Check all functions/methods for which migration causes 1D return
      values. These are mostly indexing and the reduction functions
      (see :ref:`sparse-migration-shapes-reductions` below).
   -  Check all places where you iterate over spmatrices and change them
      to account for the sparrays yielding 1D sparrays rather than 2D spmatrices.
   -  Find and change places where your code makes use of ``np.matrix``
      features. Convert those to ``np.ndarray`` features.
   -  Rename all ``*_matrix`` constructor calls to ``*_array``.
   -  Test your code. And **read** your code. You have migrated to
      sparray.

.. _sparse-migration-construction:

Details: construction functions
-------------------------------

These four functions are new and only handle sparrays:
:func:`~scipy.sparse.block_array`, :func:`~scipy.sparse.diags_array`,
:func:`~scipy.sparse.eye_array`, and :func:`~scipy.sparse.random_array`.
Their signatures are::

   def block_array(blocks, format=None, dtype=None):
   def diags_array(diagonals, /, *, offsets=0, shape=None, format=None, dtype=None):
   def eye_array(m, n=None, *, k=0, dtype=float, format=None):
   def random_array(m, n, density=0.01, format='coo', dtype=None, rng=None, data_random_state=None):

Existing functions that need careful migration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These functions return sparray or spmatrix, depending on the input types they
receive: :func:`~scipy.sparse.kron`, :func:`~scipy.sparse.kronsum`,
:func:`~scipy.sparse.hstack`, :func:`~scipy.sparse.vstack`,
:func:`~scipy.sparse.block_diag`, :func:`~scipy.sparse.tril`, and
:func:`~scipy.sparse.triu`. Their signatures are::

   def kron(A, B, format=None):
   def kronsum(A, B, format=None):
   def hstack(blocks, format=None, dtype=None):
   def vstack(blocks, format=None, dtype=None):
   def block_diag(mats, format=None, dtype=None):
   def tril(A, k=0, format=None):
   def triu(A, k=0, format=None):

Use of these functions should be examined and inputs adjusted to ensure return
values are sparrays. And in turn the outputs should be treated as sparrays.
To return sparrays, at least one input must be an sparray. If you use
list-of-lists or numpy arrays as input you should convert one of them
to a sparse array to get sparse arrays out.

Functions that changed names for the migration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   =========  =============
   Function    New function
   =========  =============
   eye         eye_array
   identity    eye_array
   diags       diags_array
   spdiags     dia_array
   bmat        block
   rand        random_array
   random      random_array
   =========  =============

.. _sparse-migration-shapes-reductions:

Details: shape changes and reductions
-------------------------------------

-  Construction using 1d-list of values:

   -  ``csr_matrix([1, 2, 3]).shape == (1, 3)`` creates a 2D matrix.
   -  ``csr_array([1, 2, 3]).shape == (3,)`` creates a 1D array

-  Indexing:

   -  Indexing of sparray allows 1D objects which can be made 2D using
      ``np.newaxis`` or ``None``. E.g., ``A[3, None, :]`` gives a 2D
      row. Indexing of 2D sparray with implicit (not given) column index
      gives a 1D result, e.g. ``A[3]`` (note: best not to do this - write it as
      ``A[3, :]`` instead). If you need a 2D result, use ``np.newaxis``, or
      ``None`` in your index, or wrap the integer index as a list for which
      fancy indexing gives 2D, e.g. ``A[[3], :]``.
   -  Iteration over sparse object: ``next(M)`` yields a sparse 2D row matrix,
      ``next(A)`` yields a sparse 1D array.

-  Reduction operations along an axis reduce the shape:

   -  ``M.sum(axis=1)`` returns a 2D row matrix by summing along axis 1.
   -  ``A.sum(axis=1)`` returns a 1D ``coo_array`` summing along axis 1.
      Some reductions return dense arrays/matrices instead of sparse ones:

      ============  =========
      Method        Result
      ============  =========
      sum(axis)     dense
      mean(axis)    dense
      argmin(axis)  dense
      argmax(axis)  dense
      min(axis)     sparse
      max(axis)     sparse
      nanmin(axis)  sparse
      nanmax(axis)  sparse
      ============  =========

   Generally, 2D sparray inputs lead to 1D results. 2D spmatrix
   inputs lead to 2D results.

-  Some reductions return a scalar. Those should behave as they did
   before and shouldn’t need to be considered during migration. E.g.
   ``A.sum()``

.. _sparse-migration-removed-methods:

Removed methods and attributes
------------------------------

-  The methods ``get_shape``, ``getrow``, ``getcol``, ``asfptype``, ``getnnz``,
   ``getH`` and the attributes ``.A`` and ``.H`` are only present on spmatrices,
   not sparrays. It is recommended that you replace usage of them with
   alternatives before starting the shift to sparray.

       ===============  ====================
       Function         Alternative
       ===============  ====================
       M.get_shape()    A.shape
       M.getformat()    A.format
       M.asfptype(…)    A.astype(…)
       M.getmaxprint()  A.maxprint
       M.getnnz()       A.nnz
       M.getnnz(axis)   A.count_nonzero(axis)
       M.getH()         A.conj().T
       M.getrow(i)      A[i, :]
       M.getcol(j)      A[:, j]
       M.A              A.toarray()
       M.H              A.conj().T
       ===============  ====================

-  Shape assignment (``M.shape = (2, 6)``) is not permitted for sparray.
   Instead you should use ``A.reshape``.

-  ``M.getnnz()`` returns the number of stored values – not the number
   of non-zeros. ``A.nnz`` does the same. To get the number of
   non-zeros, use ``A.count_nonzero()``. This is not new to the
   migration, but can be confusing.

   To migrate from the ``axis`` parameter of ``M.getnnz(axis=...)``,
   you can use ``A.count_nonzero(axis=...)``
   but it is not an exact replacement because it counts nonzero
   values instead of stored values. The difference is the number
   of explicitly stored zero values. If you really want the number
   of stored values by axis you will need to use some numpy tools.

   The numpy tools approach works for COO, CSR, CSC formats, so convert
   to one of them. For CSR and CSC, the major axis is compressed and
   ``np.diff(A.indptr)`` returns a dense 1D array with the number of
   stored values for each major axis value (row for CSR and column
   for CSC). The minor axes can be computed using
   ``np.bincount(A.indices, minlength=N)`` where ``N`` is the length
   of the minor axis (e.g. ``A.shape[1]`` for CSR). The ``bincount``
   function works for any axis of COO format using ``A.coords[axis]``
   in place of ``A.indices``.

Other
-----

-  If you provide compressed data to a constructor,
   e.g. ``csr_array((data, indices, indptr))`` sparse arrays set the index
   dtype by only checking the index arrays dtype, while sparse matrices check the
   index values too and may downcast to int32
   (see `gh-18509 <https://github.com/scipy/scipy/pull/18509>`__ for more details).

-  Binary operations with operators ``+, -, *, /, @, !=, >`` and sparse and/or
   dense operands:

   -  If all inputs are sparse, the output is usually sparse as well. The
      exception being ``/`` which returns dense (dividing by the default
      value ``0`` is ``nan``).

   -  If inputs are mixed sparse and dense, the result is usually dense
      (i.e., ``np.ndarray``). Exceptions are ``*`` which is sparse, and ``/``
      which is not implemented for ``dense / sparse``, and returns sparse for
      ``sparse / dense``.

-  Binary operations with operators ``+, -, *, /, @, !=, >`` and array
   and/or matrix operands:

   -  If all inputs are arrays, the outputs are arrays and the same is true for
      matrices.

   -  When mixing sparse arrays with sparse matrices, the leading operand
      provides the type for the output, e.g. ``sparray + spmatrix`` gives a
      sparse array while reversing the order gives a sparse matrix.

   -  When mixing dense matrices with sparse arrays, the results are usually
      arrays with exceptions for comparisons, e.g. ``>`` which return dense
      matrices.

   -  When mixing dense arrays with sparse matrices, the results are usually
      matrices with an exception for ``array @ sparse matrix`` which returns a
      dense array.
