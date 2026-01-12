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
      matrix power for spmatrix. Matrix power for sparrays uses
      ``scipy.sparse.linalg.matrix_power(A, n)``.
-  When index arrays are provided to the constructor functions, spmatrix
   selects a dtype based on dtype and values of the incoming arrays, while
   sparray only considers the dtype of the incoming arrays. For example,
   ``M=csr_matrix((data, indices, indptr))`` results in ``int32`` dtype for
   ``M.indices`` so long as the values in ``indices`` and ``indptr`` are small,
   even if the ``dtype`` of the incoming arrays are ``int64``. In contrast,
   ``A=csr_array((data, indices, indptr))`` results in ``int64`` dtype for
   ``A.indices`` when the input arrays are ``int64``. This provides more
   predictable, often larger, index dtypes in sparrays and less casting
   to match dtypes.

-  Checking the sparse type and format:

   -  ``issparse(A)`` returns ``True`` for any sparse array/matrix.
   -  ``isspmatrix(M)`` returns ``True`` for any sparse matrix.
   -  ``isspmatrix_csr(M)`` checks for a sparse matrix with specific format.
      It should be replaced with an array compatible version such as:
   -  ``issparse(A) and A.format == 'csr'`` which checks for a CSR sparse
      array/matrix.

-  Handling your software package API with sparse input/output:

   -  Inputs are fairly easy to make work with either spmatrix or sparray. So
      long as you use ``A.multiply(B)`` for elementwise and ``A @ B`` for matrix
      multiplication, and you use ``sparse.linalg.matrix_power`` for matrix
      power, you should be fine after you complete the "first pass" of the
      migration steps described in the next section. Your code will handle
      both types of inputs interchangeably.
   -  Migrating sparse outputs from your functions requires a little more thought.
      Make a list of all your public functions that return spmatrix objects.
      Check whether you feel OK returning sparrays instead. That depends on
      your library and its users. If you want to allow these functions to
      continue to return spmatrix or sparray objects, you can often do that
      using a sparse input that also serves as a signal for what type of output
      should be returned. Design your function to return the type that was input.
      That approach can be extended to dense inputs. If the input is an np.matrix
      or a masked array with np.matrix as its ``._baseclass`` attribute, then
      return spmatrix. Otherwise return an sparray. Without those inputs, two
      other approaches are to create a keyword argument to signal which to return,
      or create a new function (like we have done with, e.g. ``eye_array``) that
      has the same basic syntax, but returns sparray. Which method you choose
      should depend on your library and your users and your preferences.

Recommended steps for migration
-------------------------------

-  First pass (leaving spmatrix in the code):

   -  In your spmatrix code, change ``*`` to ``@`` for matrix
      multiplication. Note that scalar multiplication with sparse should
      use ``*``. (See helper-code :ref:`sparse-migration-star-vs-at` below)
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
      The git search idiom ``git grep 'isspm[a-z_]*('`` can help find these.
   -  Convert all ``spdiags`` calls to ``dia_matrix``.
      See docs in :func:`spdiags<scipy.sparse.spdiags>`.
      A search for ``spdiags`` is all you need here.
   -  Run all your tests on the resulting code. You are still using
      spmatrix, not sparray. But your code and tests are prepared for
      the change and you should be able to take sparrays as input to your
      code and have them mostly "just work".

-  Second pass (switching to sparray):

   -  Convert construction functions like ``diags`` and ``triu`` to the
      array version (see :ref:`sparse-migration-construction` below).
   -  Rename all ``*_matrix`` constructor calls to ``*_array``.
   -  Check all functions/methods for which migration causes 1D return
      values. These are mostly indexing and the reduction functions
      (see :ref:`sparse-migration-shapes-reductions` below).
   -  Check all places where you iterate over spmatrices and change them
      to account for the sparrays yielding 1D sparrays rather than 2D spmatrices.
   -  Find and change places where your code makes use of ``np.matrix``
      features. Convert those to ``np.ndarray`` features.
   -  If your code reads sparse from files with ``mmread``, ``hb_read``
      or ``loadmat``, use the new keyword argument ``spmatrix=False``
      in those functions to read to sparray.
   -  If you use sparse libraries that only accept ``int32`` index arrays
      for sparse representations, we suggest using just-in-time conversion.
      Convert to ``int32`` just before you call the code that requires ``int32``.
   -  ``sparray`` selects index dtype based on the dtype of the input array instead
      of the values in the array. So if you want your index arrays to be ``int32``,
      you will need to ensure an ``int32`` dtype for each index array like ``indptr``
      that you pass to ``csr_array``. With ``spmatrix`` it is tempting to use the
      default int64 dtype for ``numpy`` arrays and rely on the sparse constructor
      to downcast if the values were small. But this downcasting leads to extra
      recasting when working with other matrices, slices or arithmetic expressions.
      For ``sparray`` you can still rely on the constructors to choose dtypes. But
      you are also given the power to choose your index dtype via the dtype of the
      incoming index arrays rather than their values. So, if you want ``int32``,
      set the dtype, e.g. ``indices = np.array([1,3,6], dtype=np.int32)`` or
      ``indptr = np.arange(9, dtype=np.int32)``, when creating the index arrays.
      See :ref:`sparse-migration-index-array-dtypes` below for more info.
      In many settings, the index array dtype isn't crucial and you can just let
      the constructors choose the dtype for both sparray and spmatrix.
   -  Test your code. And **read** your code. You have migrated to sparray.

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
   def random_array(shape, density=0.01, format='coo', dtype=None, rng=None, data_sampler=None):

The ``random_array`` function has a ``shape`` (2-tuple) arg rather than
two integers. And the ``rng`` arg defaults to NumPy's new ``default_rng()``.
This differs from the spmatrix ``rand`` and ``random`` which default to
the global RandomState instance. If you don't care much about these things,
leaving it as the default should work fine.  If you care about seeding your
random numbers, you should probably add a ``rng=...`` keyword argument
to this call when you switch functions. In summary, to migrate to ``random_array``
change the function name, switch the shape argument to a single tuple argument,
leave any other parameters as before, and think about what
sort of ``rng=`` argument should be used, if any.

The `diags_array` function uses keyword-only rules for arguments. So you have
to type the `offsets=` in front of the offsets arguments. That seems like a pain
during migration from using `diags`, but it helps avoid confusion and eases reading.
A single shape parameter replaces two integers for this migration as well.

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

   =========  =============  ==================================
   Function    New function   Comments
   =========  =============  ==================================
   eye         eye_array
   identity    eye_array
   diags       diags_array    keyword-only input
   spdiags     dia_array      shape as 2-tuple
   bmat        block
   rand        random_array   shape as 2-tuple and default_rng
   random      random_array   shape as 2-tuple and default_rng
   =========  =============  ==================================

.. _sparse-migration-shapes-reductions:

Details: shape changes and reductions
-------------------------------------

-  Construction using 1d-list of values:

   -  ``csr_array([1, 2, 3]).shape == (3,)`` 1D input makes a 1D array.
   -  ``csr_matrix([1, 2, 3]).shape == (1, 3)`` 1D input makes a 2D matrix.

-  Indexing and iteration:

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

   -  ``M.min(axis=1)`` returns a 2D row matrix of the min along axis 1.
   -  ``A.min(axis=1)`` returns a 1D ``coo_array`` of the min along axis 1.
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
   ``A.min()``

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

.. _sparse-migration-star-vs-at:

Use tests to find * and ** spots
--------------------------------

-  It can be tricky to distinguish scalar multiplication ``*`` from
   matrix multiplciation ``*`` as you migrate your code. Python solved
   this, in theory, by introducing the matrix multiplication operator
   ``@``. ``*`` is used for scalar multiplication while ``@`` for matrix
   multiplication. But converting expressions that use ``*`` for both
   can be tricky and cause eye strain. Luckily, if your code has a
   test suite that covers the expressions you need to convert, you
   can use it to find places where ``*`` is being used for matrix
   multiplication involving sparse matrices. Change those to ``@``.

   The approach monkey-patches the spmatrix class dunder methods
   to raise an exception when ``*`` is used for matrix multiplication
   (and not raise for scalar multiplication). The test suite will
   flag a failure at these locations. And a test failure is a success
   here because it shows where to make changes. Change the offending
   ``*`` to ``@``, look nearby for other similar changes, and run the
   tests again. Similarly, this approach helps find where ``**`` is
   used for matrix power. SciPy raises an exception when ``@`` is
   used with for scalar multiplication, so that will catch places where
   you change when you shouldn't have. So the test suite with this
   monkey-patch checks the corrections too.

   Add the following code to your ``conftest.py`` file.
   Then run your tests locally. If there are many matrix expressions,
   you might want to test one section of your codebase at a time.
   A quick read of the code shows that it raises a ``ValueError`` whenever
   ``*`` is used between two matrix-like objects (sparse or dense),
   and whenever ``**`` is used for matrix power. It also produces a warning
   whenever sum/mean/min/max/argmin/argmax are used with an axis so the
   output will be 2D with spmatrix and 1D with sparray. That means you
   check that the code will handle either 1D or 2D output via
   ``flatten``/``ravel``, ``np.atleast_2d`` or indexing.

   .. code-block:: python

        #================== Added to check spmatrix usage ========================
        import scipy
        from warnings import warn

        def flag_this_call(*args, **kwds):
            raise ValueError("Old spmatrix function names for rand/spdiags called")

        scipy.sparse._construct.rand = flag_this_call
        scipy.sparse._construct.spdiags = flag_this_call

        class _strict_mul_mixin:
            def __mul__(self, other):
                if not scipy.sparse._sputils.isscalarlike(other):
                    raise ValueError('Operator * used here! Change to @?')
                return super().__mul__(other)

            def __rmul__(self, other):
                if not scipy.sparse._sputils.isscalarlike(other):
                    raise ValueError('Operator * used here! Change to @?')
                return super().__rmul__(other)

            def __imul__(self, other):
                if not scipy.sparse._sputils.isscalarlike(other):
                    raise ValueError('Operator * used here! Change to @?')
                return super().__imul__(other)

            def __pow__(self, *args, **kwargs):
                raise ValueError('spmatrix ** found! Use linalg.matrix_power?')

            @property
            def A(self):
                raise TypeError('spmatrix A property found! Use .toarray()')

            @property
            def H(self):
                raise TypeError('spmatrix H property found! Use .conjugate().T')

            def asfptype(self):
                raise TypeError('spmatrix asfptype found! rewrite needed')

            def get_shape(self):
                raise TypeError('spmatrix get_shape found! Use .shape')

            def getformat(self):
                raise TypeError('spmatrix getformat found! Use .format')

            def getmaxprint(self):
                raise TypeError('spmatrix getmaxprint found! Use .shape')

            def getnnz(self):
                raise TypeError('spmatrix getnnz found! Use .nnz')

            def getH(self):
                raise TypeError('spmatrix getH found! Use .conjugate().T')

            def getrow(self):
                raise TypeError('spmatrix getrow found! Use .row')

            def getcol(self):
                raise TypeError('spmatrix getcol found! Use .col')

            def sum(self, *args, **kwds):
                axis = args[0] if len(args)==1 else args if args else kwds.get("axis", None)
                if axis is not None:
                    warn(f"\nMIGRATION WARNING: spmatrix sum found using axis={axis}. "
                         "\nsparray with a single axis will produce 1D output. "
                         "\nCheck nearby to ensure 1D output is handled OK in this spot.\n")
                print(f"{args=} {axis=} {kwds=}")
                return super().sum(*args, **kwds)

            def mean(self, *args, **kwds):
                axis = args[0] if len(args)==1 else args if args else kwds.get("axis", None)
                if axis is not None:
                    warn(f"\nMIGRATION WARNING: spmatrix mean found using axis={axis}."
                         "\nsparray with a single axis will produce 1D output.\n"
                         "Check nearby to ensure 1D output is handled OK in this spot.\n")
                return super().mean(*args, **kwds)

            def min(self, *args, **kwds):
                axis = args[0] if len(args)==1 else args if args else kwds.get("axis", None)
                if axis is not None:
                    warn(f"\nMIGRATION WARNING: spmatrix min found using axis={axis}."
                         "\nsparray with a single axis will produce 1D output. "
                         "Check nearby to ensure 1D output is handled OK in this spot.\n")
                return super().min(*args, **kwds)

            def max(self, *args, **kwds):
                axis = args[0] if len(args)==1 else args if args else kwds.get("axis", None)
                if axis is not None:
                    warn(f"\nMIGRATION WARNING: spmatrix max found using axis={axis}."
                         "\nsparray with a single axis will produce 1D output. "
                         "Check nearby to ensure 1D output is handled OK in this spot.\n")
                return super().max(*args, **kwds)

            def argmin(self, *args, **kwds):
                axis = args[0] if len(args)==1 else args if args else kwds.get("axis", None)
                if axis is not None:
                    warn(f"\nMIGRATION WARNING: spmatrix argmin found using axis={axis}."
                         "\nsparray with a single axis will produce 1D output. "
                         "Check nearby to ensure 1D output is handled OK in this spot.\n")
                return super().argmin(*args, **kwds)

            def argmax(self, *args, **kwds):
                axis = args[0] if len(args)==1 else args if args else kwds.get("axis", None)
                if axis is not None:
                    warn(f"\nMIGRATION WARNING: spmatrix argmax found using axis={axis}."
                         "\nsparray with a single axis will produce 1D output. "
                         "Check nearby to ensure 1D output is handled OK in this spot.\n")
                return super().argmax(*args, **kwds)


        class coo_matrix_strict(_strict_mul_mixin, scipy.sparse.coo_matrix):
            pass

        class bsr_matrix_strict(_strict_mul_mixin, scipy.sparse.bsr_matrix):
            pass

        class csr_matrix_strict(_strict_mul_mixin, scipy.sparse.csr_matrix):
            pass

        class csc_matrix_strict(_strict_mul_mixin, scipy.sparse.csc_matrix):
            pass

        class dok_matrix_strict(_strict_mul_mixin, scipy.sparse.dok_matrix):
            pass

        class lil_matrix_strict(_strict_mul_mixin, scipy.sparse.lil_matrix):
            pass

        class dia_matrix_strict(_strict_mul_mixin, scipy.sparse.dia_matrix):
            pass

        scipy.sparse.coo_matrix = scipy.sparse._coo.coo_matrix = coo_matrix_strict
        scipy.sparse.bsr_matrix = scipy.sparse._bsr.bsr_matrix = bsr_matrix_strict
        scipy.sparse.csr_matrix = scipy.sparse._csr.csr_matrix = csr_matrix_strict
        scipy.sparse.csc_matrix = scipy.sparse._csc.csc_matrix = csc_matrix_strict
        scipy.sparse.dok_matrix = scipy.sparse._dok.dok_matrix = dok_matrix_strict
        scipy.sparse.lil_matrix = scipy.sparse._lil.lil_matrix = lil_matrix_strict
        scipy.sparse.dia_matrix = scipy.sparse._dia.dia_matrix = dia_matrix_strict

        scipy.sparse._compressed.csr_matrix = csr_matrix_strict

        scipy.sparse._construct.bsr_matrix = bsr_matrix_strict
        scipy.sparse._construct.coo_matrix = coo_matrix_strict
        scipy.sparse._construct.csc_matrix = csc_matrix_strict
        scipy.sparse._construct.csr_matrix = csr_matrix_strict
        scipy.sparse._construct.dia_matrix = dia_matrix_strict

        scipy.sparse._extract.coo_matrix = coo_matrix_strict

        scipy.sparse._matrix.bsr_matrix = bsr_matrix_strict
        scipy.sparse._matrix.coo_matrix = coo_matrix_strict
        scipy.sparse._matrix.csc_matrix = csc_matrix_strict
        scipy.sparse._matrix.csr_matrix = csr_matrix_strict
        scipy.sparse._matrix.dia_matrix = dia_matrix_strict
        scipy.sparse._matrix.dok_matrix = dok_matrix_strict
        scipy.sparse._matrix.lil_matrix = lil_matrix_strict

        del coo_matrix_strict
        del bsr_matrix_strict
        del csr_matrix_strict
        del csc_matrix_strict
        del dok_matrix_strict
        del lil_matrix_strict
        del dia_matrix_strict
        #==========================================

.. _sparse-migration-index-array-dtypes:

Index Array DTypes
------------------

If you provide compressed indices to a constructor,
e.g. ``csr_array((data, indices, indptr))`` sparse arrays set the index
dtype by only checking the index arrays dtype, while sparse matrices check the
index values too and may downcast to int32
(see `gh-18509 <https://github.com/scipy/scipy/pull/18509>`__ for more details).
This means you may get int64 indexing when you used to get int32.
You can control this by setting the ``dtype`` before instantiating, or by
recasting after construction.

Two sparse utility functions can help with handling the index dtype.
Use ``get_index_dtype(arrays, maxval, check_contents)`` while creating indices
to find an appropriate dtype (int32 or int64) to use for your compressed indices.

Use ``safely_cast_index_arrays(A, idx_dtype)`` for recasting after construction,
while making sure you con't create overflows during downcasting.
This function doesn't actually change the input array. The cast arrays are returned.
And copies are only made when needed. So you can check if casting was done using
``if indices is not A.indices:``

The function signatures are::

    def get_index_dtype(arrays=(), maxval=None, check_contents=False):
    def safely_cast_index_arrays(A, idx_dtype=np.int32, msg=""):

Example idioms include the following for ``get_index_dtype``::

   .. code-block:: python

       # select index dtype before construction based on shape
       shape = (3, 3)
       idx_dtype = scipy.sparse.get_index_dtype(maxval=max(shape))
       indices = np.array([0, 1, 0], dtype=idx_dtype)
       indptr = np.arange(3, dtype=idx_dtype)
       A = csr_array((data, indices, indptr), shape=shape)

and for ``safely_cast_index_arrays``::

   .. code-block:: python

       # rescast after construction, raising exception if shape too big
       indices, indptr = scipy.sparse.safely_cast_index_arrays(B, np.int32)
       B.indices, B.indptr = indices, indptr

Other
-----

-  Binary operators ``+, -, *, /, @, !=, >`` act on sparse and/or dense operands:

   -  If all inputs are sparse, the output is usually sparse as well. The
      exception being ``/`` which returns dense (dividing by the default
      value ``0`` is ``nan``).

   -  If inputs are mixed sparse and dense, the result is usually dense
      (i.e., ``np.ndarray``). Exceptions are ``*`` which is sparse, and ``/``
      which is not implemented for ``dense / sparse``, and returns sparse for
      ``sparse / dense``.

-  Binary operators ``+, -, *, /, @, !=, >`` with array and/or matrix operands:

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
