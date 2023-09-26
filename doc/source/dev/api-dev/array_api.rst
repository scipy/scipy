Support for Array API
=====================

.. note:: Array API support is still experimental and hidden behind an
          environment variable. Only a small part of the public API is covered
          right now.

This guide describes how to **use** and **add support for** the
`Python Array API standard <https://data-apis.org/array-api/latest/index.html>`_.
This standard allows users to use any Array API compatible array library
with SciPy out of the box.

The `RFC`_ defines how SciPy implements support for the standard, with the main
principle being *"array type in equals array type out"*. In addition, the
implementation does more strict validation of allowed array-like inputs, e.g.
rejecting numpy matrix and masked array instances, and arrays with object
dtype.

In the following, an Array API compatible namespace is noted as ``xp``.


Using Array API support
-----------------------

To enable the Array API standard support, an environment variable must be set
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

The following modules provide Array API standard support when the environment
variable is set:

- ``scipy.cluster.hierarchy``
- ``scipy.cluster.vq``
- ``scipy.fft``


Implementation notes
--------------------

A key part of the support for the array API standard and specific compatibility
functions for Numpy, CuPy and PyTorch is provided through
`array-api-compat <https://github.com/data-apis/array-api-compat>`_.
This package is included in the SciPy code base via a git submodule (under
``scipy/_lib``), so no new dependencies are introduced.

``array-api_compat`` provides generic utility functions and adds aliases such
as ``xp.concat`` (which, for numpy, maps to ``np.concatenate``). This allows
using a uniform API across NumPy, PyTorch and CuPy (as of right now; support
for other libraries like JAX is expected to be added in the future).

When the environment variable isn't set and hence Array API support in SciPy is
disabled, we still use the "augmented" version of the NumPy namespace, which is
``array_api_compat.numpy``. That should not change behavior of SciPy functions,
it's effectively the existing ``numpy`` namespace with a number of aliases
added and a handful of functions amended/added for array API standard support.
When support is enabled, depending on the type of arrays, ``xp`` will return the
standard-compatible namespace matching the input array type to a function (e.g.,
if the input to ``cluster.vq.kmeans`` is a PyTorch array, then ``xp`` is
``array_api_compat.torch``).


Adding Array API support to a SciPy function
--------------------------------------------

As much as possible, new code added to SciPy should try to follow as closely as
possible the Array API standard (these functions typically are best-practice
idioms for NumPy usage as well). By following the standard, effectively adding
support for Array API is typically straightforward, and we ideally don't need
to maintain any customization.

Two helper functions are available:

* ``array_namespace``: detect the namespace based on input arrays and do some
  input validation (like refusing to work with masked arrays, please see the
  `RFC`_.)
* ``as_xparray``: a drop-in replacement for ``np.asarray`` with additional
  features like ``copy, check_finite``. As stated above, try to limit the use
  of non standard features. In the end we would want to upstream our needs to
  the compatibility library.

To add support to a SciPy function which is defined in a ``.py`` file, what you
have to change is:

1. Input array validation,
2. Using ``xp`` rather ``np`` functions,
3. When calling into compiled code, convert the array to a NumPy array before
   and convert it back to the input array type after.

Input array validation uses the following pattern::

   xp = array_namespace(arr)  # where `arr` is the first input array
   # Do this for each input array, it applies all the validation steps (reject
   # matrix, etc.) as well as the conversion to a numpy array if it's a
   # sequence, or preserve the non-numpy array type:
   arr = as_xparray(arr, xp=xp)

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
      b = as_xparray(b, copy=True)  # our custom helper is needed for copy

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
* ``skip_if_array_api``: don't run a test if ``SCIPY_ARRAY_API`` is on.
* ``skip_if_array_api_gpu``: don't run a test if GPU is involved (also applies
  to PyTorch's MPS mode).
* ``skip_if_array_api_backend(backend)``: don't run a test for a specific
  backend

The following is an example using the main decorator responsible of the
namespace parametrization::

  from scipy.conftest import array_api_compatible
  ...
  @array_api_compatible
  def test_toto(self, xp):
      a = xp.asarray([1, 2, 3])
      b = xp.asarray([0, 2, 5])
      toto(a, b)

Then ``dev.py`` can be used with he new option ``-b`` or
``--array-api-backend``::

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
