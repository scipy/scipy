.. _dev-arrayapi:

Support for the array API standard
==================================

.. note:: Array API standard support is still experimental and hidden behind an
          environment variable. Only a small part of the public API is covered
          right now.

This guide describes how to **use** and **add support for** the
`Python array API standard <https://data-apis.org/array-api/latest/index.html>`_.
This standard allows users to use any array API compatible array library
with parts of SciPy out of the box.

The `RFC`_ defines how SciPy implements support for the standard, with the main
principle being *"array type in equals array type out"*. In addition, the
implementation does more strict validation of allowed array-like inputs, e.g.
rejecting numpy matrix and masked array instances, and arrays with object
dtype.

In the following, an array API compatible namespace is noted as ``xp``.


Using array API standard support
--------------------------------

To enable the array API standard support, an environment variable must be set
before importing SciPy:

.. code:: bash

   export SCIPY_ARRAY_API=1

This both enables array API standard support and the more strict input
validation for array-like arguments. *Note that this environment variable is
meant to be temporary, as a way to make incremental changes and merge them into
``main`` without affecting backwards compatibility immediately. We do not
intend to keep this environment variable around long-term.*

This clustering example shows usage with PyTorch tensors as inputs and return
values:

.. code:: python

    >>> import torch
    >>> from scipy.cluster.vq import vq
    >>> code_book = torch.tensor([[1., 1., 1.],
    ...                           [2., 2., 2.]])
    >>> features  = torch.tensor([[1.9, 2.3, 1.7],
    ...                           [1.5, 2.5, 2.2],
    ...                           [0.8, 0.6, 1.7]])
    >>> code, dist = vq(features, code_book)
    >>> code
    tensor([1, 1, 0], dtype=torch.int32)
    >>> dist
    tensor([0.4359, 0.7348, 0.8307])

Note that the above example works for PyTorch CPU tensors. For GPU tensors or
CuPy arrays, the expected result for ``vq`` is a ``TypeError``, because ``vq``
uses compiled code in its implementation, which won't work on GPU.

More strict array input validation will reject ``np.matrix`` and
``np.ma.MaskedArray`` instances, as well as arrays with ``object`` dtype:

.. code:: python

    >>> import numpy as np
    >>> from scipy.cluster.vq import vq
    >>> code_book = np.array([[1., 1., 1.],
    ...                       [2., 2., 2.]])
    >>> features  = np.array([[1.9, 2.3, 1.7],
    ...                       [1.5, 2.5, 2.2],
    ...                       [0.8, 0.6, 1.7]])
    >>> vq(features, code_book)
    (array([1, 1, 0], dtype=int32), array([0.43588989, 0.73484692, 0.83066239]))

    >>> # The above uses numpy arrays; trying to use np.matrix instances or object
    >>> # arrays instead will yield an exception with `SCIPY_ARRAY_API=1`:
    >>> vq(np.asmatrix(features), code_book)
    ...
    TypeError: 'numpy.matrix' are not supported

    >>> vq(np.ma.asarray(features), code_book)
    ...
    TypeError: 'numpy.ma.MaskedArray' are not supported

    >>> vq(features.astype(np.object_), code_book)
    ...
    TypeError: object arrays are not supported


Currently supported functionality
`````````````````````````````````

The following modules provide array API standard support when the environment
variable is set:

- `scipy.cluster`
- `scipy.constants`
- `scipy.datasets`
- `scipy.fft`
- `scipy.io`
- `scipy.ndimage`
- `scipy.special`
- `scipy.stats`

Individual functions in the above modules provide a capability table in the
documentation like the one below. If the table is absent, the function does not
yet support backends other than NumPy.


Example capabilities table
--------------------------

=========  =========  =========
Library    CPU        GPU
=========  =========  =========
NumPy      ✅         n/a
CuPy       n/a        ✅
PyTorch    ✅         ✅
JAX        ⚠️ no JIT  ⛔
Dask       ⛔         n/a
=========  =========  =========

In the example above, the feature has some support for NumPy, CuPy, PyTorch, and JAX
arrays, but no support for Dask arrays. Some backends, like JAX and PyTorch, natively
support multiple devices (CPU and GPU), but SciPy support for such arrays may be
limited; for instance, this SciPy feature is only expected to work with JAX arrays
located on the CPU. Additionally, some backends can have major caveats; in the example
the function will fail when running inside ``jax.jit``.
Additional caveats may be listed in the docstring of the function.

While the elements of the table marked with "n/a" are inherently out of scope, we are
continually working on filling in the rest.
Dask wrapping around backends other than NumPy (notably, CuPy) is currently out of scope
but it may change in the future.

Please see `the tracker issue`_ for updates.


Implementation notes
--------------------

A key part of the support for the array API standard and specific compatibility
functions for Numpy, CuPy and PyTorch is provided through
`array-api-compat <https://github.com/data-apis/array-api-compat>`_.
This package is included in the SciPy codebase via a git submodule (under
``scipy/_lib``), so no new dependencies are introduced.

``array-api-compat`` provides generic utility functions and adds aliases such
as ``xp.concat`` (which, for numpy, mapped to ``np.concatenate`` before NumPy added
``np.concat`` in NumPy 2.0). This allows using a uniform API across NumPy, PyTorch,
CuPy and JAX (with other libraries, such as Dask, being worked on).

When the environment variable isn't set and hence array API standard support in
SciPy is disabled, we still use the wrapped version of the NumPy namespace,
which is ``array_api_compat.numpy``. That should not change behavior of SciPy
functions, as it's effectively the existing ``numpy`` namespace with a number of
aliases added and a handful of functions amended/added for array API standard
support. When support is enabled, ``xp = array_namespace(input)`` will
be the standard-compatible namespace matching the input array type to a
function (e.g., if the input to `cluster.vq.kmeans` is a PyTorch tensor, then
``xp`` is ``array_api_compat.torch``).


Adding array API standard support to a SciPy function
-----------------------------------------------------

As much as possible, new code added to SciPy should try to follow as closely as
possible the array API standard (these functions typically are best-practice
idioms for NumPy usage as well). By following the standard, effectively adding
support for the array API standard is typically straightforward, and we ideally
don't need to maintain any customization.

Various helper functions are available in ``scipy._lib._array_api`` - please see
the ``__all__`` in that module for a list of current helpers, and their docstrings
for more information.

To add support to a SciPy function which is defined in a ``.py`` file, what you
have to change is:

1. Input array validation,
2. Using ``xp`` rather ``np`` functions,
3. When calling into compiled code, convert the array to a NumPy array before
   and convert it back to the input array type after.

Input array validation uses the following pattern::

   xp = array_namespace(arr) # where arr is the input array
   # alternatively, if there are multiple array inputs, include them all:
   xp = array_namespace(arr1, arr2)

   # replace np.asarray with xp.asarray
   arr = xp.asarray(arr)
   # uses of non-standard parameters of np.asarray can be replaced with _asarray
   arr = _asarray(arr, order='C', dtype=xp.float64, xp=xp)

Note that if one input is a non-NumPy array type, all array-like inputs have to
be of that type; trying to mix non-NumPy arrays with lists, Python scalars or
other arbitrary Python objects will raise an exception. For NumPy arrays, those
types will continue to be accepted for backwards compatibility reasons.

If a function calls into a compiled code just once, use the following pattern::

   x = np.asarray(x)  # convert to numpy right before compiled call(s)
   y = _call_compiled_code(x)
   y = xp.asarray(y)  # convert back to original array type

If there are multiple calls to compiled code, ensure doing the conversion just
once to avoid too much overhead.

Here is an example for a hypothetical public SciPy function ``toto``::

  def toto(a, b):
      a = np.asarray(a)
      b = np.asarray(b, copy=True)

      c = np.sum(a) - np.prod(b)

      # this is some C or Cython call
      d = cdist(c)

      return d

You would convert this like so::

  def toto(a, b):
      xp = array_namespace(a, b)
      a = xp.asarray(a)
      b = xp_copy(b, xp=xp)  # our custom helper is needed for copy

      c = xp.sum(a) - xp.prod(b)

      # this is some C or Cython call
      c = np.asarray(c)
      d = cdist(c)
      d = xp.asarray(d)

      return d

Going through compiled code requires going back to a NumPy array, because
SciPy's extension modules only work with NumPy arrays (or memoryviews in the
case of Cython). For arrays on CPU, the
conversions should be zero-copy, while on GPU and other devices the attempt at
conversion will raise an exception. The reason for that is that silent data
transfer between devices is considered bad practice, as it is likely to be a
large and hard-to-detect performance bottleneck.


Adding tests
------------

To run a test on multiple array backends, you should add the ``xp`` fixture to it,
which is valued to the currently tested array namespace. 

The following pytest markers are available:

* ``skip_xp_backends(backend=None, reason=None, np_only=False, cpu_only=False, eager_only=False, exceptions=None)``:
  skip certain backends or categories of backends.
  See docstring of ``scipy.conftest.skip_or_xfail_xp_backends`` for information on how
  to use this marker to skip tests.
* ``xfail_xp_backends(backend=None, reason=None, np_only=False, cpu_only=False, eager_only=False, exceptions=None)``:
  xfail certain backends or categories of backends.
  See docstring of ``scipy.conftest.skip_or_xfail_xp_backends`` for information on how
  to use this marker to xfail tests.
* ``skip_xp_invalid_arg`` is used to skip tests that use arguments which
  are invalid when ``SCIPY_ARRAY_API`` is enabled. For instance, some tests of
  `scipy.stats` functions pass masked arrays to the function being tested, but
  masked arrays are incompatible with the array API. Use of the
  ``skip_xp_invalid_arg`` decorator allows these tests to protect against
  regressions when ``SCIPY_ARRAY_API`` is not used without resulting in failures
  when ``SCIPY_ARRAY_API`` is used. In time, we will want these functions to emit
  deprecation warnings when they receive array API invalid input, and this
  decorator will check that the deprecation warning is emitted without it
  causing the test to fail. When ``SCIPY_ARRAY_API=1`` behavior becomes the
  default and only behavior, these tests (and the decorator itself) will be
  removed.
* ``array_api_backends``: this marker is automatically added by the ``xp`` fixture to
  all tests that use it. This is useful e.g. to select all and only such tests::

    python dev.py test -b all -m array_api_backends

``scipy._lib._array_api`` contains array-agnostic assertions such as ``xp_assert_close``
which can be used to replace assertions from `numpy.testing`.

When these assertions are executed within a test that uses the ``xp`` fixture, they
enforce that the namespaces of both the actual and desired arrays match the namespace
which was set by the fixture. Tests without the ``xp`` fixture infer the namespace from
the desired array. This machinery can be overridden by explicitly passing the ``xp=``
parameter to the assertion functions.

The following examples demonstrate how to use the markers::

  from scipy.conftest import skip_xp_invalid_arg
  from scipy._lib._array_api import xp_assert_close
  ...
  @pytest.mark.skip_xp_backends(np_only=True, reason='skip reason')
  def test_toto1(self, xp):
      a = xp.asarray([1, 2, 3])
      b = xp.asarray([0, 2, 5])
      xp_assert_close(toto(a, b), a)
  ...
  @pytest.mark.skip_xp_backends('array_api_strict', reason='skip reason 1')
  @pytest.mark.skip_xp_backends('cupy', reason='skip reason 2')
  def test_toto2(self, xp):
      ...
  ...
  # Do not run when SCIPY_ARRAY_API is used
  @skip_xp_invalid_arg
  def test_toto_masked_array(self):
      ...

Passing names of backends into ``exceptions`` means that they will not be skipped
by ``cpu_only=True`` or ``eager_only=True``. This is useful when delegation
is implemented for some, but not all, non-CPU backends, and the CPU code path 
requires conversion to NumPy for compiled code::

  # array-api-strict and CuPy will always be skipped, for the given reasons.
  # All libraries using a non-CPU device will also be skipped, apart from
  # JAX, for which delegation is implemented (hence non-CPU execution is supported).
  @pytest.mark.skip_xp_backends(cpu_only=True, exceptions=['jax.numpy'])
  @pytest.mark.skip_xp_backends('array_api_strict', reason='skip reason 1')
  @pytest.mark.skip_xp_backends('cupy', reason='skip reason 2')
  def test_toto(self, xp):
      ...

After applying these markers, ``dev.py test`` can be used with the new option
``-b`` or ``--array-api-backend``::

  python dev.py test -b numpy -b torch -s cluster

This automatically sets ``SCIPY_ARRAY_API`` appropriately. To test a library
that has multiple devices with a non-default device, a second environment
variable (``SCIPY_DEVICE``, only used in the test suite) can be set. Valid
values depend on the array library under test, e.g. for PyTorch, valid values are
``"cpu", "cuda", "mps"``. To run the test suite with the PyTorch MPS
backend, use: ``SCIPY_DEVICE=mps python dev.py test -b torch``.

Note that there is a GitHub Actions workflow which tests with array-api-strict,
PyTorch, and JAX on CPU.


Testing the JAX JIT compiler
----------------------------
The `JAX JIT compiler <https://jax.readthedocs.io/en/latest/jit-compilation.html>`_
introduces special restrictions to all code wrapped by `@jax.jit`, which are not
present when running JAX in eager mode. Notably, boolean masks in `__getitem__`
and `.at` aren't supported, and you can't materialize the arrays by applying
`bool()`, `float()`, `np.asarray()` etc. to them.

To properly test scipy with JAX, you need to wrap the tested scipy functions
with `@jax.jit` before they are called by the unit tests.
To achieve this, you should tag them as follows in your test module::

  from scipy._lib._lazy_testing import lazy_xp_function
  from scipy.mymodule import toto

  lazy_xp_function(toto)

  def test_toto(xp):
      a = xp.asarray([1, 2, 3])
      b = xp.asarray([0, 2, 5])
      # When xp==jax.numpy, toto is wrapped with @jax.jit
      xp_assert_close(toto(a, b), a)

See full documentation in `scipy/_lib/_lazy_testing.py`.


Additional information
----------------------

Here are some additional resources which motivated some design decisions and
helped during the development phase:

* Initial `PR <https://github.com/tupui/scipy/pull/24>`__ with some discussions
* Quick started from this `PR <https://github.com/scipy/scipy/pull/15395>`__ and
  some inspiration taken from
  `scikit-learn <https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/utils/_array_api.py>`__.
* `PR <https://github.com/scikit-learn/scikit-learn/issues/22352>`__ adding Array
  API support to scikit-learn
* Some other relevant scikit-learn PRs:
  `#22554 <https://github.com/scikit-learn/scikit-learn/pull/22554>`__ and
  `#25956 <https://github.com/scikit-learn/scikit-learn/pull/25956>`__

.. _RFC: https://github.com/scipy/scipy/issues/18286
.. _the tracker issue: https://github.com/scipy/scipy/issues/18867
