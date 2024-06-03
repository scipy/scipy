.. _dev-arrayapi:

Support for the array API standard
==================================

.. note:: Array API standard support is still experimental and hidden behind an
          environment variable. Only a small part of the public API is covered
          right now.

This guide describes how to **use** and **add support for** the
`Python array API standard <https://data-apis.org/array-api/latest/index.html>`_.
This standard allows users to use any array API compatible array library
with SciPy out of the box.

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
is not a pure Python function and hence won't work on GPU.

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

- `scipy.cluster.hierarchy`
- `scipy.cluster.vq`
- `scipy.constants`
- `scipy.fft`

Support is provided in `scipy.special` for the following functions:
`scipy.special.log_ndtr`, `scipy.special.ndtr`, `scipy.special.ndtri`,
`scipy.special.erf`, `scipy.special.erfc`, `scipy.special.i0`,
`scipy.special.i0e`, `scipy.special.i1`, `scipy.special.i1e`,
`scipy.special.gammaln`, `scipy.special.gammainc`, `scipy.special.gammaincc`,
`scipy.special.logit`, and `scipy.special.expit`.

Support is provided in `scipy.stats` for the following functions:
`scipy.stats.pearsonr` and `scipy.stats.moment`.


Implementation notes
--------------------

A key part of the support for the array API standard and specific compatibility
functions for Numpy, CuPy and PyTorch is provided through
`array-api-compat <https://github.com/data-apis/array-api-compat>`_.
This package is included in the SciPy code base via a git submodule (under
``scipy/_lib``), so no new dependencies are introduced.

``array-api-compat`` provides generic utility functions and adds aliases such
as ``xp.concat`` (which, for numpy, maps to ``np.concatenate``). This allows
using a uniform API across NumPy, PyTorch, CuPy and JAX (with other libraries,
such as Dask, coming in the future).

When the environment variable isn't set and hence array API standard support in
SciPy is disabled, we still use the "augmented" version of the NumPy namespace,
which is ``array_api_compat.numpy``. That should not change behavior of SciPy
functions, it's effectively the existing ``numpy`` namespace with a number of
aliases added and a handful of functions amended/added for array API standard
support. When support is enabled, depending on the type of arrays, ``xp`` will
return the standard-compatible namespace matching the input array type to a
function (e.g., if the input to `cluster.vq.kmeans` is a PyTorch array, then
``xp`` is ``array_api_compat.torch``).


Adding array API standard support to a SciPy function
-----------------------------------------------------

As much as possible, new code added to SciPy should try to follow as closely as
possible the array API standard (these functions typically are best-practice
idioms for NumPy usage as well). By following the standard, effectively adding
support for the array API standard is typically straightforward, and we ideally
don't need to maintain any customization.

Three helper functions are available:

* ``array_namespace``: return the namespace based on input arrays and do some
  input validation (like refusing to work with masked arrays, please see the
  `RFC`_.)
* ``_asarray``: a drop-in replacement for ``asarray`` with the additional
  parameters ``check_finite`` and ``order``. As stated above, try to limit
  the use of non-standard features. In the end we would want to upstream our
  needs to the compatibility library. Passing ``xp=xp`` avoids duplicate calls
  of ``array_namespace`` internally.
* ``copy``: an alias for ``_asarray(x, copy=True)``.
  The ``copy`` parameter was only introduced to ``np.asarray`` in NumPy 2.0,
  so use of the helper is needed to support ``<2.0``. Passing ``xp=xp`` avoids
  duplicate calls of ``array_namespace`` internally.

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

   # uses of non-standard parameters of np.asarray can be replaced with _asarray
   arr = _asarray(arr, order='C', dtype=xp.float64, xp=xp)

Note that if one input is a non-numpy array type, all array-like inputs have to
be of that type; trying to mix non-numpy arrays with lists, Python scalars or
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
      b = copy(b, xp=xp)  # our custom helper is needed for copy

      c = xp.sum(a) - xp.prod(b)

      # this is some C or Cython call
      c = np.asarray(c)
      d = cdist(c)
      d = xp.asarray(d)

      return d

Going through compiled code requires going back to a NumPy array, because
SciPy's extension modules only work with NumPy arrays (or memoryviews in the
case of Cython), but not with other array types. For arrays on CPU, the
conversions should be zero-copy, while on GPU and other devices the attempt at
conversion will raise an exception. The reason for that is that silent data
transfer between devices is considered bad practice, as it is likely to be a
large and hard-to-detect performance bottleneck.


Adding tests
------------

The following pytest markers are available:

* ``array_api_compatible -> xp``: use a parametrisation to run a test on
  multiple array backends.
* ``skip_xp_backends(*backends, reasons=None, np_only=False, cpu_only=False)``:
  skip certain backends and/or devices. ``np_only`` skips tests for all backends
  other than the default NumPy backend.
  ``@pytest.mark.usefixtures("skip_xp_backends")`` must be used alongside this
  marker for the skipping to apply.
* ``skip_xp_invalid_arg`` is used to skip tests that use arguments which
  are invalid when ``SCIPY_ARRAY_API`` is used. For instance, some tests of
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

The following is an example using the markers::

  from scipy.conftest import array_api_compatible, skip_xp_invalid_arg
  ...
  @pytest.mark.skip_xp_backends(np_only=True,
                                 reasons=['skip reason'])
  @pytest.mark.usefixtures("skip_xp_backends")
  @array_api_compatible
  def test_toto1(self, xp):
      a = xp.asarray([1, 2, 3])
      b = xp.asarray([0, 2, 5])
      toto(a, b)
  ...
  @pytest.mark.skip_xp_backends('array_api_strict', 'cupy',
                                 reasons=['skip reason 1',
                                          'skip reason 2',])
  @pytest.mark.usefixtures("skip_xp_backends")
  @array_api_compatible
  def test_toto2(self, xp):
      a = xp.asarray([1, 2, 3])
      b = xp.asarray([0, 2, 5])
      toto(a, b)
  ...
  # Do not run when SCIPY_ARRAY_API is used
  @skip_xp_invalid_arg
  def test_toto_masked_array(self):
      a = np.ma.asarray([1, 2, 3])
      b = np.ma.asarray([0, 2, 5])
      toto(a, b)

Passing a custom reason to ``reasons`` when ``cpu_only=True`` is unsupported
since ``cpu_only=True`` can be used alongside passing ``backends``. Also,
the reason for using ``cpu_only`` is likely just that compiled code is used
in the function(s) being tested.

When every test function in a file has been updated for array API
compatibility, one can reduce verbosity by telling ``pytest`` to apply the
markers to every test function using ``pytestmark``::

    from scipy.conftest import array_api_compatible

    pytestmark = [array_api_compatible, pytest.mark.usefixtures("skip_xp_backends")]
    skip_xp_backends = pytest.mark.skip_xp_backends
    ...
    @skip_xp_backends(np_only=True, reasons=['skip reason'])
    def test_toto1(self, xp):
        a = xp.asarray([1, 2, 3])
        b = xp.asarray([0, 2, 5])
        toto(a, b)

After applying these markers, ``dev.py test`` can be used with the new option
``-b`` or ``--array-api-backend``::

  python dev.py test -b numpy -b pytorch -s cluster

This automatically sets ``SCIPY_ARRAY_API`` appropriately. To test a library
that has multiple devices with a non-default device, a second environment
variable (``SCIPY_DEVICE``, only used in the test suite) can be set. Valid
values depend on the array library under test, e.g. for PyTorch (currently the
only library with multi-device support that is known to work) valid values are
``"cpu", "cuda", "mps"``. So to run the test suite with the PyTorch MPS
backend, use: ``SCIPY_DEVICE=mps python dev.py test -b pytorch``.

Note that there is a GitHub Actions workflow which runs ``pytorch-cpu``.


Additional information
----------------------

Here are some additional resources which motivated some design decisions and
helped during the development phase:

* Initial `PR <https://github.com/tupui/scipy/pull/24>`__ with some discussions
* Quick started from this `PR <https://github.com/scipy/scipy/pull/15395>`__ and
  some inspiration taken from
  `scikit-learn <https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/utils/_array_api.py>`__.
* `PR <https://github.com/scikit-learn/scikit-learn/issues/22352>`__ adding Array
  API surpport to scikit-learn
* Some other relevant scikit-learn PRs:
  `#22554 <https://github.com/scikit-learn/scikit-learn/pull/22554>`__ and
  `#25956 <https://github.com/scikit-learn/scikit-learn/pull/25956>`__

.. _RFC: https://github.com/scipy/scipy/issues/18286
