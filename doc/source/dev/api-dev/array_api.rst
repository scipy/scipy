Support for Array API
=====================

.. note:: Array API support is still experimental and only a fraction of the
          public API is covered

This guide describes how to **use** and **add support** for the
`Array API <https://data-apis.org/array-api/latest/index.html>`_ standard.
This standard allow users to use any Array API compatible array without having
to do anything.

The `RFC`_ defines how SciPy implements the standard.

In the following, a Array API compatible namespace is noted as ``xp``.

Use Array API
-------------

There are two global variables controlling the behaviour of functions which
are supporting the Array API standard.

* ``SCIPY_ARRAY_API``: either a bool, ``"all"``, or a list of supported
  backend. For now these are ``"cupy", "pytorch", "numpy", "numpy.array_api"``.
  The global variable can be adjusted dynamically. This is an optional flag.
  By default the backend is ``numpy`` and it only affects functions modified
  with the helpers described in the next section.
* ``SCIPY_DEVICE``: can be ``"cpu", "gpu", "mps"`` (something supported by the
  backend)

A new runtime dependency is introduced via a submodule:
`array-api-compat <https://github.com/data-apis/array-api-compat>`_.

This library allows to convert arrays from various libraries (NumPy, PyTorch,
Cupy) into Array API compatible arrays. It add aliases such as ``xp.concat``.

When the global flag is off, we still use the "augmented" version of NumPy
arrays, which is ``array_api_compat.numpy``. When the flag is on, depending
on the type of arrays, it will return the corresponding compatible namespace.

Add support to a function
-------------------------

As much as possible, new code should try to follow as closely as possible the
Array API. By following the standard, effectively adding support for Array API
is trivial and we ideally don't need to maintain any customization.

Two helper functions are available:

* ``array_namespace``: detect the namespace based on input arrays and do some
  input validation (like refusing to work with masked arrays, please see the
  `RFC`_.)
* ``as_xparray``: a drop-in replacement for ``np.asarray`` with additional
  features like ``copy, check_finite``. As stated above, try to limit the use
  of non standard features. In the end we would want to upstream our needs to
  the compatibility library.

Here are some practical guidelines using an example::

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


The key is that, going to compiled code requires to go back to a NumPy array.
Meaning that it will raise if you are on a GPU as the conversion would be
silent otherwise-which is not good. For CPU interactions, it should be a
zero-copy operation.

Adding tests
------------

Some testing fixture:

* ``array_api_compatible -> xp``: use a parametrisation to run a test on
  multiple Array backend.
* ``skip_if_array_api``: don't run a test if ``SCIPY_ARRAY_API`` is on.
* ``skip_if_array_api_gpu``: don't run a test if GPU is involved (also applies
  to PyTorch's MPS mode.)
* ``skip_if_array_api_backend(backend)``: don't run a test for a specific
  backend

Following is an example using the main decorator responsible of the namespace
parametrization::

  @array_api_compatible
  def test_toto(self, xp):
      a = xp.asarray([1, 2, 3])
      b = xp.asarray([0, 2, 5])
      toto(a, b)


Then ``dev.py`` can be used with he new option ``-b`` or
``--array-api-backend``::

  python dev.py test -b numpy -b pytorch -s cluster

This automatically set ``SCIPY_ARRAY_API`` appropriately. And if wanted, to
test a different device, ``SCIPY_DEVICE`` can be manually set. e.g.
``SCIPY_DEVICE=mps python dev.py ...``.

Finally, there is a GitHub action workflow which runs ``pytorch-cpu``.


Additional information
----------------------

Here are some additional resources which motivated and helped during the
development phase:

* Initial `PR <https://github.com/tupui/scipy/pull/24>`__ with some discussions
* Quick started from this `PR <https://github.com/scipy/scipy/pull/15395>`__ and
  take some inspiration from
  `scikit-learn <https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/utils/_array_api.py>`__.
* `PR <https://github.com/scikit-learn/scikit-learn/issues/22352>`__ adding Array
  API surpport to scikit-learn
* Some other PRs from scikit-learn
  `22554 <https://github.com/scikit-learn/scikit-learn/pull/22554>`__ and
  `25956 <https://github.com/scikit-learn/scikit-learn/pull/25956>`__

.. _RFC: https://github.com/scipy/scipy/issues/18286
