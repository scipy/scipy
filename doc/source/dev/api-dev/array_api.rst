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
uses compiled code in its implementation, which won't work on GPU, and there
are currently no GPU specific implementations to delegate to.

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
but this may change in the future.

Please see `the tracker issue`_ for updates.


.. _dev-arrayapi_implementation_notes:

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

.. _dev-arrayapi_adding_support:

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
      b = xp.asarray(b, copy=True)

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

In some cases, compiled code can be supported on alternative backends
through delegation to native implementations. Such delegation has currently been set
up in `~scipy.fft`, `~scipy.ndimage`, `~scipy.signal`
(see :ref:`array_api_support_signal_caveats`),
and `~scipy.special`, though there is not yet a
standard approach, and support in each module has mostly evolved separately.

There is also some effort being put into expanding access to native
implementations, such as the `xsf project <https://github.com/scipy/xsf/issues/1>`_
to establish a library of mathematical special function implementations which support
both CPU and GPU.

.. _dev-arrayapi_jax_support:

A note on JAX support
`````````````````````

JAX was designed with deliberate restrictions to make code easier to reason about and exploits
this to better support features like JIT-compilation and autodifferentiation. The most
relevant restrictions for SciPy developers are:

* JAX arrays are immutable. Rather than performing in-place updates of arrays, one
  can use the `at <https://docs.jax.dev/en/latest/_autosummary/jax.numpy.ndarray.at.html>`_
  property to transform an array in an equivalent way. Inside a JIT compiled function,
  an expression like ``x = x.at[idx].set(y)`` will be applied in-place under the hood.

* Functions using the JAX JIT must be functionally pure. They cannot have side
  effects, cannot mutate data, and their outputs must be determined completely by their
  inputs. Raising a Python exception is a side-effect that is not permitted within a JITed
  function.

* Within the JIT, value based control flow with Python ``if`` statements is not permitted.
  Only static properties of arrays such as their ``shape`` and ``dtype`` are permitted to be
  used with ``if``. `xp.where <https://data-apis.org/array-api/2024.12/API_specification/generated/array_api.where.html#where>`_
  and `array_api_extra.apply_where <https://data-apis.org/array-api-extra/generated/array_api_extra.apply_where.html>`_ are
  provide some basic control flow that works with the JIT.

* Within the JIT, the shapes of output arrays cannot depend dynamically on the *values* in input arrays.

See `Common Gotchas in JAX <https://docs.jax.dev/en/latest/notebooks/Common_Gotchas_in_JAX.html>`_ for more information.

**Recommendations for developers:**

* To work around the mutability restriction, developers adding JAX support
  to SciPy functions which make in-place updates should use
  `array_api_extra.at <https://data-apis.org/array-api-extra/generated/array_api_extra.at.html>`_
  which works for all array API compatible backends, delegating to JAX's ``at`` for
  JAX arrays and performing regular in-place operations with ``__setitem__`` for other kinds
  of arrays.

* The restriction that functions be functionally pure to support the JAX JIT necessitates
  that input-validation that raises with bad input must be skipped when ``xp`` is JAX.

* Compiled functions ``f`` which cannot be supported on JAX through delegation to a native
  implementation can potentially still be supported through the use of
  `array_api_extra.lazy_apply <https://data-apis.org/array-api-extra/generated/array_api_extra.lazy_apply.html>`_, which uses JAX's `pure_callback <https://docs.jax.dev/en/latest/_autosummary/jax.pure_callback.html>`_
  mechanism to enable calling Python functions within JIT-ed JAX functions.

  Using ``lazy_apply``, the example function ``toto`` might be made compatible
  with the JAX JIT in the following way::

    def toto(a, b):
        xp = array_namespace(a, b)
        a = xp.asarray(a)
        b = xp.asarray(b, copy=True)

        c = xp.sum(a) - xp.prod(b)

        # this is some C or Cython call
        # as_numpy=True tells lazy_apply to convert to and from NumPy.
        d = xpx.lazy_apply(cdist, c, as_numpy=True)

        return d

  ``lazy_apply`` can be used so long as ``f`` is a pure function for which the output
  shape can be determined knowing only the input shapes.

* `xp.where <https://data-apis.org/array-api/2024.12/API_specification/generated/array_api.where.html#where>`_
  and `array_api_extra.apply_where <https://data-apis.org/array-api-extra/generated/array_api_extra.apply_where.html>`_ provide a level of basic control flow that works with the JIT
  and in some cases these can be used to replace the value dependent use of ``if``. In
  some cases its also possible to wrap code using ``if`` within a pure function and use
  ``lazy_apply``.

**JAX Eager:**

It is also possible to run JAX in eager-mode without the JIT (in fact this is the
default behavior when ``@jax.jit`` is not applied). Eager-mode comes with serious
performance limitations and is typically only used to debug functions which are
ultimately intended to be run with the JIT. Do not be tempted to attempt to distinguish
whether JAX is being used with the JIT in order to bypass some restrictions while running
with eager-mode. There is no way to make this distinction using JAX's public API, and any means
of determining if JAX is running with the JIT would necessarily involve implementation
details that SciPy should not rely upon. In general, support for eager-mode is not a high-value
target, and it is not considered a good use of developer time to put significant effort
into enabling eager-only support.

.. _dev-arrayapi_default_dtype:

Default Datatypes
`````````````````

The Array API Standard allows conforming libraries to have
`default datatypes <https://data-apis.org/array-api/latest/API_specification/data_types.html#default-data-types>`_
for integers and real and complex floating point numbers which differ from the
``int64``, ``float64``, ``complex128`` defaults used by NumPy. Our aim is to
have array API supporting SciPy functions with array inputs have behavior which
is independent of the default dtype to the extent that this is practical. This means
that any when using array creation functions from the ``xp`` namespace such as ``xp.zeros``
or ``xp.arange``, one should take care to explicitly set a dtype with the ``dtype``
kwarg; otherwise, the result will depend on the default dtype.

Note that SciPy is currently only tested in CI on platforms/backends/backend-configurations
where a ``float64`` dtype is available. At the moment, ``float32``-only support
varies from function to function and is not well-documented. Some examples of
``float32``-only backends are JAX in its
`default configuration <https://docs.jax.dev/en/latest/default_dtypes.html>`_,
and PyTorch with the
`Metal performance shader backend <https://pytorch.org/blog/introducing-accelerated-pytorch-training-on-mac/>`_
on ARM Mac. We are open to expanding and improving ``float32``-only
support in cases where this is feasible and there is sufficient user interest.


Array creation functions without array inputs
`````````````````````````````````````````````

For array creation functions without array inputs,
adding array API standard support can be accomplished by adding keyword
only arguments ``xp`` and ``device`` to specify the desired backend and
device respectively. See for instance `~scipy.signal.buttap` which
constructs the analog prototype of an Nth-order Butterworth filter. It may also
be desirable to add a ``dtype`` kwarg to control the output ``dtype`` for
such functions.

Note that none of these keyword arguments are necessary when
there are array arguments from which the backend, device, and desired output dtype
can be inferred (ideally, output dtype should follow input dtypes through the expected
`type promotion rules <https://data-apis.org/array-api/latest/API_specification/type_promotion.html>`_).
For the sake of API simplicity and consistency and in the
spirit of "There should be one-- and preferably only one --obvious way to do it.",
it is recommended to avoid the use of these kwargs in functions which take
at least one array argument.

It is still under debate how array creation functions without array inputs should
behave with respect to :ref:`default dtypes <dev-arrayapi_default_dtype>`.
Should they respect default dtype or should the output dtype be fixed across
backends and defaults? Should there be a ``dtype`` kwarg for controlling the output
dtype or is being able to apply ``xp.astype`` on the output sufficient?
Since there is not yet a consistent pattern to follow, for now its
important to clearly document how such functions behave with respect to the
default dtype in the :ref:`extra_note <dev-arrayapi_extra_note>` described below.


Documenting array API standard support
--------------------------------------

Support for alternative array API standard backends can be registered and
documented using the ``xp_capabilities`` decorator which has the following
signature::

  def xp_capabilities(
      *,
      # Alternative capabilities table.
      # Used only for testing this decorator.
      capabilities_table=None,
      # Generate pytest.mark.skip/xfail_xp_backends.
      # See documentation in conftest.py.
      # lists of tuples [(module name, reason), ...]
      skip_backends=(), xfail_backends=(),
      cpu_only=False, np_only=False, reason=None,
      out_of_scope=False, exceptions=(),
      # lists of tuples [(module name, reason), ...]
      warnings=(),
      # xpx.testing.lazy_xp_function kwargs.
      # Refer to array-api-extra documentation.
      allow_dask_compute=False, jax_jit=True,
      # Extra note to inject into the docstring
      extra_note=None,
  ):

This is available in ``scipy._lib._array_api`` and can be applied to functions,
methods, and classes to declare the current extent of their array API standard
support. For the sake of brevity, in the remainder of this document, we write
as if ``xp_capabilities`` only applies to functions.

The ``xp_capabilities`` decorator is what inserts the capabilities table into
docstrings. It also allows developers to tag tests with the
``@make_xp_test_case`` decorator or apply ``make_xp_pytest_param`` to pytest
parameters to automatically generate backend specific SKIP/XFAIL markers, and
setting up testing that functions work with the JAX JIT or work in Dask lazily
(i.e. without materializing arrays with ``dask.compute`` or otherwise triggering
computation with ``dask.persist``).

.. warning::

   The modification of docstrings by ``xp_capabilities`` can potentially break
   intersphinx references because it currently has the side effect of replacing
   implicit roles with ``:func:``. This can be avoided by explicitly
   setting the role for references to classes and methods that are
   outside of SciPy. The following snippet is taken from the docstring for
   `~scipy.signal.detrend` where the role ``:meth:`` for a class method is
   needed.

   .. code-block:: rst

       See Also
       --------
       :meth:`numpy.polynomial.polynomial.Polynomial.fit` : Create least squares fit polynomial

Note that in some modules a systematic process for delegation to native
implementations is set up, where functions are replaced with wrappers
that perform delegation. In this case, ``xp_capabilities`` is not always
applied as a decorator with ``@`` syntax, but may instead be applied
programatically on the wrappers. When working on array API standard
support within a module, its important to be aware of how such delegation
is set up, if any, and how ``xp_capabilities`` is being applied. A common
practice currently is to have a file, ``_support_alternative_backends.py``
within a module that sets up such delegation. See for instance
`scipy/signal/_support_alternative_backends.py <https://github.com/scipy/scipy/blob/main/scipy/signal/_support_alternative_backends.py>`_.


Basic behavior
``````````````

Using ``xp_capabilities``  with no arguments, like this::

  @xp_capabilities()
  def my_function(x):
      ...

declares that a function works on all supported backends, on JAX with the JIT
and lazily in Dask. This is most likely to occur if a function is written
entirely in terms of the array API standard as described earlier in this
document. Such functions are commonly referred to as array-agnostic. For
functions which are written mostly in terms of the array API standard, but
include calls to compiled code sandwiched between conversions to and from NumPy,
``xp_capabilities`` should be given the ``cpu_only=True`` option. Backends
which are supported on GPU by such a function ``f`` through delegation to a
native implementation can be specified with the ``exceptions`` kwarg, which
in this case takes a list of strings specifying GPU-capable backends. The
currently supported string values when using ``cpu_only=True`` are
``"cupy"``, ``"torch"``, and ``"jax.numpy"``.

It is recommended to reserve ``cpu_only=False`` (the default) for
array-agnostic functions which are expected to work on all array API compliant
backends, including ones not tested in SciPy and ones that do not yet exist.
If a function is supported on GPU on all tested backends through delegation to
respective native implementations, one should use ``cpu_only=True`` while listing
each backend in the list of ``exceptions`` like so::

  @xp_capabilities(cpu_only=True, exceptions=["cupy"])

When setting ``cpu_only=True``, one may list a reason by passing a string with
the ``reason`` kwarg. This can be helpful for documenting why something is not
supported for other developers. The reason will appear in the ``pytest`` output
when the SciPy test suite is run with ``pytest`` in verbose mode
(with the ``-v`` flag).


JAX JIT
```````

One may declare a function as not supporting the JAX JIT with the option
``jax_jit=False``. See the earlier note on :ref:`JAX support <dev-arrayapi_jax_support>`
for more information.


Unsupported functions
`````````````````````

Functions which do not support the array API standard through the means
described earlier in this document should either receive the ``np_only=True``
option or the ``out_of_scope=True`` option. The former should be used for
functions which are not currently supported but which are considered in-scope
for future support.

Functions for which array API support has *not* been added following the
procedures described earlier in this document, but for which delegation to a
native implementation has been set up for one or more array API backends, should
still use ``np_only=True`` in their ``xp_capabilities`` entries. Just as for
``cpu_only=True``, exceptions can be passed with the ``exceptions`` kwarg (and
also just as for ``cpu_only=True`` one can pass a reason with the ``reason``
kwarg)::

  @xp_capabilities(
      np_only=True, exceptions=["cupy"],
      reason="not converted yet but has CuPy delegation."
  )

Valid strings to pass in the exceptions list are ``"array_api_strict"``,
``"cupy"``, ``"dask.array"``, ``"jax.numpy"``, and ``"torch"``.

If ``np_only=True`` and ``"torch"`` or ``"jax.numpy"`` is added to
the lists of exceptions, it will be declared as supported on both CPU and
GPU.

.. _dev-arrayapi_jax_jit_no_gpu:

.. dropdown:: Functions with JAX JIT support but no GPU support
  :icon: alert
  :color: warning
	      
  It's possible for a function to be natively available in JAX,
  support ``jax.jit``, but not be supported on GPU. Thus, it may be
  possible for JAX delegation to be set up in a function which has
  not yet received the array API standard compatibility treatment,
  and for the JIT to be supported but not the GPU.

  Because ``exceptions`` does double
  duty declaring exceptions to ``cpu_only=True`` and ``np_only=True``, it is not
  possible to express this situation using ``xp_capabilities`` in the way
  described above. This is not too serious of an issue because the intention is
  that ``np_only=True`` is only a temporary state. Through the means described
  above in the section on :ref:`adding array API support <dev-arrayapi_adding_support>`,
  it is a reasonable goal for all functions in SciPy's public API to at least reach
  the state ``cpu_only=True``. For functions still waiting in the ``np_only=True`` state,
  ``xp_capabilities``'s ``skip_backends`` kwarg can be used as an escape hatch to
  allow more fine grained declaration of capabilities. See the section on
  :ref:`skip_backends and xfail_backends <dev-arrayapi_skip_xfail_backends>`.

``out_of_scope=True`` signals that there is no intention to ever provide array
API support for a given function. There is not yet a formal policy for which
functions should be considered out-of-scope. Some general rules of thumb that are
being followed are to exclude:

* functions which do not operate on arrays such as :doc:`scipy.constants.value <../../reference/generated/scipy.constants.value>`
* functions which are too implementation specific such as those in `scipy.linalg.blas` which give direct wrappers to low-level BLAS routines.
* functions which would inherently be very difficult or even impossible to compute efficiently on accelerated computing devices.

As an example, the contents of `scipy.odr` are considered out-of-scope for a
combination of reasons 2 and 3 above. `scipy.odr` essentially provides a direct
wrapper of the monolithic ODRPACK Fortran library, and its API is tied
to the structure of this monolithic library. An efficient GPU
accelerated implementation of nonlinear weighted orthogonal distance regression
would benefit from not having to support an API so tightly coupled to ODRPACK
but is also a challenging problem in its own right.

(Since the previous paragraph was written `scipy.odr` has been slated for
deprecation. Things that are deprecated are inherently out-of-scope).

Considerations of what to consider in-scope are evolving, and something which is now
considered out-of-scope may be decided to be in-scope in the future if sufficient user
interest and feasability are demonstrated.

.. _dev-arrayapi_skip_xfail_backends:

``skip_backends`` and ``xfail_backends``
````````````````````````````````````````
One may pass lists of tuples of backend string, reason pairs to ``xp_capabilities``
with the ``skip_backends`` and ``xfail_backends`` kwargs. The valid backend strings
are ``"array_api_strict"``, ``"cupy"``, ``"dask.array"``, ``"jax.numpy"`` and ``"torch"``
(note that one should almost never want to skip tests for
`array_api_strict <https://data-apis.org/array-api-strict/>`_ as failures with this
backend most likely indicate a failure to correctly follow the array API standard).

Any backend passed in such a way with either kwarg is declared as unsupported with both
CPU and GPU. The difference between ``skip_backends`` and ``xfail_backends`` is that for
tests using the ``xp`` fixture, ``skip_backends`` adds ``pytest.skip`` markers for
backends and the corresponding tests are skipped entirely, while with ``xfail_backends``,
``pytest.xfail`` markers are added, and tests are still run but expected to fail.

One example in which it is pertinent to use
``skip_backends`` is for functions which otherwise support the array API standard, but
use features which are not available on a particular backend, such as mutation of
arrays through item assignment, which is not supported in JAX. For instance the following
can be used to signify a function which is otherwise array-agnostic, but uses
in-place item-assignment::

  @xp_capabilities(
      skip_backends=[("jax.numpy", "in-place item assignment")]
  )
  def function_with_internal_mutation(x):
      ...


Another case is when there is a bug in the corresponding array library, in which
case the ``reason`` string should contain a link to the upstream issue.

In the caveat above about functions with JAX JIT support but no JAX GPU support
we discussed the edge-case of a function which has not been given array
API standard support in the usual way, is available on JAX through delegation to
a native implementation which supports ``jax.jit``, but does not work on the GPU.
For now, such situations can in principle be handled by using ``cpu_only=True``
and passing in any backends which are not supported even on CPU to ``skip_backends``::

  @xp_capabilities(
      cpu_only=True,
      skip_backends=[
          ("array_api_strict", "not supported"),
	  ("cupy", "not supported"),
	  ("dask.array", "not supported"),
	  ("torch", "not supported"),
      ]
  )
  def oddball_function(x):
      ...

Such situations are hopefully rare enough that special handling isn't needed.
``xp_capabilities`` has evolved naturally over time to meet developer needs; good
suggestions for ways to improve developer ergonomics are welcome.

Dask Compute
````````````

The default, ``allow_dask_compute=False`` declares that a function works lazily
in Dask and will not materialize any Dask arrays with ``dask.compute`` or
otherwise initiate computation with ``dask.persist``. Use
``allow_dask_compute=True`` to declare that a function supports Dask arrays but
not lazily. Developers can also pass an integer to give a cap for the number of
combined calls to ``dask.compute`` and ``dask.persist`` that are allowed. If a function
is not array-agnostic, then it will typically be the case that
``allow_dask_compute=True`` should be set, unless Dask specific codepaths have been added.

Dask support is currently deprioritized due to structural barriers
that make the development of meaningful Dask support particularly challenging. At present,
developers should feel free to reflexively add ``skip_backends=[("dask.array", "deprioritized")]``
to the ``xp_capabilities`` entry of any function they are working on. Reprioritization may
be considered in the future if a champion emerges and the structural outlook improves.
See `RFC: Should Dask support remain a priority? #24205 <https://github.com/scipy/scipy/issues/24205>`_
for relevant discussion.

.. _dev-arrayapi_extra_note:

``extra_note``
``````````````
Some functions may be supported on an alternative backend, but only in particular
cases, perhaps only for some values of a kwarg, for real arrays but not complex ones,
or only for arrays with fewer than a given number of dimensions. Such caveats should
be supplied with the ``extra_note`` kwarg of ``xp_capabilities``. Note that the
implementation of ``extra_note`` simply inserts a string directly into the Notes section
of the docstring, and one must be careful about whitespace. This is perhaps
best demonstrated with an example::

  uses_choose_conv_extra_note = (
    """CuPy does not support inputs with ``ndim>1`` when ``method="auto"``
    but does support higher dimensional arrays for ``method="direct"``
    and ``method="fft"``.

    """
  )

.. _dev-arrayapi_adding_tests:


Adding tests
------------

To run a test on multiple array backends, you should add the ``xp`` fixture to
it. ``xp`` currently supports testing with the following backends:

* `array_api_strict <https://data-apis.org/array-api-strict/>`_
* `cupy <https://cupy.dev/>`_
* `dask.array <https://www.dask.org/>`_
* `jax.numpy <https://docs.jax.dev/en/latest/>`_,
* `numpy <https://numpy.org/>`_
* `torch <https://pytorch.org/>`_

``xp`` is a
`pytest fixture <https://docs.pytest.org/en/6.2.x/fixture.html>`_
which is parameterized over all currently installed backends among
those listed above. Note that ``xp`` takes values from the set of "raw"
namespaces, not the wrapped namespaces from
:ref:`array_api_compat <dev-arrayapi_implementation_notes>`.

``scipy._lib._array_api`` provides the ``make_xp_test_case``
decorator, and the ``make_xp_pytest_param`` and ``make_xp_pytest_marks``
functions to declare which functions are being tested by a test.  These draw on
the ``xp_capabilities`` entries for a function (or in some cases those for a
list of functions) to insert the relevant backend specific skip and xfail
markers.

**make_xp_test_case:**

In most cases, developers should use ``make_xp_test_case``, which is applied as a
decorator to a test function, test method, or entire test class. Applying it to a
test class is equivalent to applying it to each method separately. The decorator can
be applied at both the class and method level as below::

  @make_xp_test_case(my_function)
  class TestMyFunction:
      def test1(self, xp):
          ...

      @make_xp_test_case(other_function)
      def test_integration_with_other_function(self, xp)
          ...

Applying ``@make_xp_test_case(my_function)`` to ``TestMyFunction`` causes
all skips and xfails from the ``xp_capabilities`` entry for ``my_function``
to be applied to all methods in the class. Additional applications of
``@make_xp_test_case`` to individual methods will add additional skips and
xfails and not override the class level decorator. Below is an equivalent
way to write the above example. This style can become unwieldy when there
are many methods in a class testing the same function.::

  class TestMyFunction:
      @make_xp_test_case(my_function)
      def test1(self, xp):
          ...

      @make_xp_test_case(my_function, other_function)
      def test_integration_with_other_function(self, xp)
          ...

**make_xp_pytest_param:**

``make_xp_pytest_param`` covers the situation where a common test body is
parametrized over a list of functions using ``pytest.mark.parametrize``.
It is not used as frequently as ``make_xp_test_case`` but this kind of
situation is not too uncommon.::

  @pytest.mark.parametrize(
      "func",
      [make_xp_pytest_param(func) for func in tested_funcs]
  )
  def test_foo(func, xp):
      ...

Without access to ``make_xp_pytest_param``, one might instead have to do
something like::

  @make_xp_test_case(*test_funcs)
  @pytest.mark.parametrize(
      "func", tested_funcs
  )
  def test_foo(func, xp):
      ...

But then ``test_foo`` would take on the collective skips and xfails
for all of the functions in ``test_funcs`` taken together, leading to
tests being run with unnecessary skips and xfails.

Unlike ``make_xp_test_case``, only a single function can be passed to any given
call to ``make_xp_pytest_param``. Additional arguments specify additional
parameters for ``pytest.mark.parametrize``, such as in the example
below::

  @pytest.mark.parametrize(
      "func,norm",
      [
          make_xp_pytest_param(func, norm)
	  for func, norm in it.product(tested_functions, [True, False])
      ]
   )
   def test_normed_foo(func, norm, xp):
       ...

**make_xp_pytest_marks:**

``make_xp_pytest_marks`` is rarely used. It directly returns a list of
pytest marks which can be used with the ``pytestmark = ...`` variable
to set marks for all tests in a file.

**Strict checks:**

The ``xp`` fixture should almost always be used along with ``make_xp_test_case``
or one of the similar functions listed above and the ``xp`` fixture has
strict checks to enforce this. If one had accidentally written::

  @pytest.mark.parametrize(
      "func", tested_funcs
  )
  def test_foo(func, xp):
      ...

without using ``make_xp_pytest_param`` then running this test would result
in an error with the the message::

  ERROR scipy/my_module/tests/test_foo.py::test_foo[numpy] - UserWarning: test uses `xp`
  fixture without drawing from `xp_capabilities`  but is not explicitly marked with ``pytest.mark...

Since ``xp_capabilities`` is used to declare alternative backend support for the
purpose of both testing and documentation, this strict check in the ``xp``
fixture ensures that documentation of tested array API capabilities does not
become out-of-date. There may be cases where one intentionally does cannot or
does not want to use ``make_xp_test_case`` or an equivalent, such as for private
functions which do not have associated ``xp_capabilities`` entries. To bypass
the strict checks, one can explicitly mark a test with
``@pytest.mark.uses_xp_capabilities(False)``. An optional ``reason`` string can
be passed to this mark. Tests of private functionality for which there are no
relevant ``xp_capabilities`` entries, one should use ``reason="private"``.::

  pytest.mark.uses_xp_capabilities(False, reason="private")
  def test_private_toto_helper(xp):
      ...

Directly adding pytest markers
``````````````````````````````

Though most of the time it's sufficient to use ``make_xp_test_case`` and
``make_xp_pytest_param``, the following ``pytest`` markers are available and can
be added directly to tests. (``make_xp_test_case`` and its equivalents provide a
declarative means of adding ``skip_xp_backends`` and ``xfail_xp_backends``
markers).

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

    spin test -b all -m array_api_backends
* ``uses_xp_capabilities(status, funcs=None, reason=None)``: discussed above.
  ``make_xp_test_case`` and its equivalents apply the marker
  ``uses_xp_capabilities(True)`` and direct use of ``uses_xp_capabilities(False)``
  can be used to declare a test intentionally does not use ``make_xp_test_case``
  or one of its equivalents.

Test specific skips and xfails
``````````````````````````````

For a public function ``f``, ``skip_xp_backends`` and ``xfail_xp_backends`` should
only be used directly for backend related skips and xfails which are needed for
the specific test but which do not reflect the general capabilities of
``f``. Reasons to directly use ``skip_xp_backends`` include when:

1. the test body itself contains unsupported functionality (though one should
   try to avoid this whenever possible, see the subsection on testing
   practice below).
2. ``f`` is only partially supported on a backend and the test relies on
   cases which are not supported, e.g. tests involving complex values for
   functions which only support real values on a given backend, tests involving
   higher dimensional arrays for functions which only support arrays of size 2d
   or less on a given backend.
3. the test exposes a bug in ``f`` on a given backend which crashes test
   execution.

For tests exposing bugs on alternative backends that do not crash test
execution, such as bugs that lead to numerical errors, it is preferable to use
``xfail_xp_backends`` so we can be notified with an ``XPASS`` when the
upstream bug is fixed.

``xfail_xp_backends`` should not be used for test failures for an alternative
backend which are actually unrelated to ``f`` but are instead due to bugs
outside ``f`` exposed by other parts of the test body. To avoid such situations,
we recommend as a general practice to attempt to isolate use of the alternative
backend only to the function ``f`` being tested with a caveat that there are
situations where or it is necessary or desired to do otherwise: see the section
on :ref:`backend isolation <dev-arrayapi_backend_isolation>` below for more
information.

Tests which are inherently NumPy only should avoid the ``xp`` fixture
altogether rather than using it with an ``np_only=True`` skip marker.

Note that, in one case, ``xp_capabilities`` offers more granularity than
``skip_xp_backends`` and ``xfail_xp_backends``. ``xp_capabilities`` allows
developers to separately declare support for the JAX JIT and support for lazy
computation with Dask with the respective ``jax_jit`` and ``allow_dask_compute``
kwargs.  ``skip_xp_backends`` (``xfail_xp_backends``) offers only an
``eager_only`` kwarg which can only add skips (xfails) for both the JAX JIT and
lazy Dask together. The current state is that one cannot add test specific skips
(xfails) for the JAX JIT without also adding them for lazy Dask and vice
versa. This is a known limitation that arose through the historical process
through which ``xp_capabilities``, ``skip_xp_backends``, and
``xfail_xp_backends`` were developed, and it may be addressed in the future if
there is sufficient developer need.


Array-agnostic assertions
`````````````````````````

``scipy._lib._array_api`` contains array-agnostic assertions such as ``xp_assert_close``
which can be used to replace assertions from `numpy.testing`.

When these assertions are executed within a test that uses the ``xp`` fixture, they
enforce that the namespaces of both the actual and desired arrays match the namespace
which was set by the fixture. Tests without the ``xp`` fixture infer the namespace from
the desired array. This machinery can be overridden by explicitly passing the ``xp=``
parameter to the assertion functions.


Examples
````````

The following examples demonstrate how to use direct markers together with
``make_xp_test_case``::

  from scipy.conftest import skip_xp_invalid_arg
  from scipy._lib._array_api import xp_assert_close, make_xp_test_case

  @make_xp_test_case(toto)
  class TestToto:
      def test_toto_list_input(self):
      # This test is inherently NumPy only so avoids the xp fixture altogether.
          a = [1., 2., 3.]
          b = [0., 2., 5.]
          xp_assert_close(toto(a, b), np.array(a))
  ...
      @pytest.mark.skip_xp_backends(
          'cupy',
	  reason='cupy does not support inputs with ndim>2'
      )
      def test_toto2(self, xp):
          ...
  ...
      # Do not run when SCIPY_ARRAY_API=1 is used since calling toto on masked
      # arrays will raise in this case.
      @skip_xp_invalid_arg
      def test_toto_masked_array(self):
          ...


Running tests
`````````````

After applying these markers, either through ``make_xp_test_case`` or one of its
equvilents, or directly,
``spin test`` can be used with the option ``-b`` or ``--array-api-backend``::

  spin test -b numpy -b torch -s cluster

This automatically sets ``SCIPY_ARRAY_API`` appropriately and will cause
tests with the ``xp`` fixture to run only for the selected backends to be
collected. Valid backends are ``numpy``, ``array_api_strict``, ``cupy``,
``dask.array``, ``jax.numpy``, and ``torch``. One may also use the
``-m array_api_backends`` option to restrict collection to only tests using
the ``xp`` fixture. For instance the following command causes pytest to only
collect tests using the ``xp`` fixture with the CuPy backend::

  spin test -b cupy -m array_api_backends

To test a library
that has multiple devices with a non-default device, a second environment
variable (``SCIPY_DEVICE``, only used in the test suite) can be set. Valid
values depend on the array library under test, e.g. for PyTorch, valid values are
``"cpu", "cuda", "mps"``. To run the test suite with the PyTorch MPS
backend, use: ``SCIPY_DEVICE=mps spin test -b torch``.

Note that in SciPy's GitHub Actions workflows, there are regular tests
with array-api-strict, Dask, PyTorch, and JAX on CPU, and tests with
CuPy, PyTorch, and JAX on GPU.

A third environment variable (``SCIPY_DEFAULT_DTYPE``, again only used in the
test suite) can be used to control the :ref:`default dtype <dev-arrayapi_default_dtype>`
used by ``torch`` in tests. Valid values are ``"float64"`` and ``"float32"``.
If ``SCIPY_DEFAULT_DTYPE`` is unset, then ``torch``'s default dtype will be ``float64``.

The intention behind testing with different default dtypes is primarily to catch
subtle bugs that can arise with the ``torch`` backend due to internal array creation
that does not explicitly specify a dtype. The intention is not to implicitly test
that functions are numerically accurate with both ``float32`` and ``float64`` inputs,
or that input dtype controls output dtype. These tasks, if done, should instead be accomplished
mindfully by explicitly setting dtypes within tests. For the sake of consistency,
tests intended to be run with the ``torch`` backend should not use array creation
functions without explicitly setting the dtype. At the time of writing, there are many
tests in the test suite which do not follow this practice, and this could be a good source
of first issues for new contributors.

.. _dev-arrayapi_backend_isolation:

Backend isolation in tests
``````````````````````````

In most cases, it's important that for any supported function ``f``, there exist
tests using the ``xp`` fixture that restrict use of alternative backends to only
the function ``f`` being tested. Other functions evaluated within a test, for
the purpose of producing reference values, inputs, round-trip calculations,
etc. should instead use the NumPy backend. This helps ensure that any failures
that occur on a backend actually relate to the function of interest, and avoids
the need to skip backends due to lack of support for functions other than
``f``. Property based integration tests which check that some invariant holds
using the same alternative backend across different functions can also have
value, giving a window into the general health of backend support for a module,
but in order to ensure the test suite actually reflects the state of backend
support for each function, it's usually vital to have tests which isolate use
of the alternative backend only to the function being tested.

To help facilitate such backend isolation, there is a function
``_xp_copy_to_numpy`` in ``scipy._lib._array_api`` which can copy an arbitrary
``xp`` array to a NumPy array, bypassing any device transfer guards, while
preserving dtypes. It is essential that this function is only used in
tests. Attempts to copy a device array to NumPy outside of tests should fail,
because otherwise it is opaque as to whether a function is working on GPU or
not. Creation of input arrays and reference output arrays, and computations that
verify that the output of the function being tested satisfies an invariant (such
as round trip tests that a function composed with its inverse gives the identity
function), should all be done with NumPy (using the ``_xp_copy_to_numpy``
function if necessary).

Such backend isolation should not be applied blindly. Consider for example a
vectorized root finding function like `scipy.optimize.elementwise.find_root`.
When testing such a function on alternative backends, isolating use of the
alternative backend only to ``find_root`` by using an input callable ``f`` (the
function for which roots are sought) that converts to and from NumPy would not
be desirable since since ``find_root`` and ``f`` are so tightly coupled in this
case. In other cases, a function ``h`` used in the tests of a function ``g`` may
be known to be so simple and rock solid that there is no point in going through
the trouble of backend isolation. Maintainers are free to use their discretion to
decide whether backend isolation is necessary or desirable.

Testing the JAX JIT compiler (and lazy evaluation with Dask)
------------------------------------------------------------
The `JAX JIT compiler <https://jax.readthedocs.io/en/latest/jit-compilation.html>`_
introduces special restrictions to all code wrapped by ``@jax.jit``, which are not
present when running JAX in eager mode. Notably, boolean masks in ``__getitem__``
and ``.at`` aren't supported, and you can't materialize the arrays by applying
``bool()``, ``float()``, ``np.asarray()`` etc. to them.

To properly test scipy with JAX, the tested scipy functions must be wrapped
with ``@jax.jit`` before they are called by the unit tests. This is done
automatically when using ``make_xp_test_case`` and its friends when the
associated ``xp_capabilities`` entry (or entries) have ``jax_jit=True``::

  from scipy._lib._array_api import make_xp_test_case, xp_assert_close
  from scipy.mymodule import toto

  @make_xp_test_case(toto)
  def test_toto(xp):
     a = xp.asarray([3, 10, 5, 16, 8, 4, 2, 1, ])
     b = xp.asarray([3, 5, 8, 4, 2, 1])
     # When xp==jax.numpy, toto is wrapped with @jax.jit
     # so long as the xp_capabilities entry for toto has
     # jax_jit=True.
     xp_assert_close(toto(a), b)

To achieve this for private functions without ``xp_capabilities`` entries,
you should tag them as follows in your test module::

  from scipy._lib._array_api import xp_assert_close
  from scipy._lib.array_api_extra.testing import lazy_xp_function
  from scipy.mymodule import _private_toto_helper

  lazy_xp_function(_private_toto_helper)

  @pytest.mark.uses_xp_capabilities(False, reason="private")
  def test_private_toto_helper(xp):
      a = xp.asarray([1, 2, 3])
      b = xp.asarray([0, 2, 5])
      # When xp==jax.numpy, _private_toto_helper is wrapped with @jax.jit
      xp_assert_close(_private_toto_helper(a, b), a)

.. warning::

   If instead of importing the functions from ``scipy.mymodule``, the above example
   imported ``mymodule`` and called ``toto`` through the qualified name
   ``mymodule.toto``, ``@jax.jit`` would not be applied to ``toto``.  This is due to an
   implementation specific quirk which limits the application of ``@jax.jit`` only
   to functions in the globals of the module that defines the current test.
   If one wishes to use a pattern like ``mymodule.toto`` in a test, one must define a
   variable ``lazy_xp_modules`` at the top of the test file to specify additional places
   the testing framework should look for functions tagged with ``lazy_xp_function``::

     import scipy.mymodule as mymodule
     from scipy._lib._array_api import make_xp_test_case, xp_assert_close

     lazy_xp_modules = [mymodule]

     @make_xp_test_case(mymodule.toto)
     def test_toto(xp):
         a = xp.asarray([3, 10, 5, 16, 8, 4, 2, 1, ])
         b = xp.asarray([3, 5, 8, 4, 2, 1])
         # When xp==jax.numpy, toto is wrapped with @jax.jit
         # so long as the xp_capabilities entry for toto has
         # jax_jit=True.
         xp_assert_close(toto(a), b)

   This can be slightly annoying to remember at first, but in practice isn't too bad
   once one gets in the habit of checking for this. The essential complexity of
   ``lazy_xp_function`` is actually quite high, and the current design trades off on
   developer ergonomics to allow for a simpler implementation.

Testing lazy evaluation with Dask works similarly, except ``lazy_xp_function`` wraps
functions with a decorator that disables ``compute()`` and ``persist()`` and ensures
that exceptions and warnings are raised eagerly. Similarly as for the JAX JIT,
``make_xp_test_case`` and friends will automatically do this when the associated
``xp_capabilities`` entry has ``allow_dask_compute=False``. The same warning about
requiring ``lazy_xp_modules`` applies for tests Dask works with lazy evaluation just
as it does for tests of the JAX JIT.

See full documentation `here <https://data-apis.org/array-api-extra/generated/array_api_extra.testing.lazy_xp_function.html>`_.

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

.. _dev-arrayapi_coverage:

API Coverage
------------
The below tables show the current state of alternative backend support across
SciPy's modules. Currently only public functions and function-like callable
objects are included in the tables, but it is planned to eventually also include
relevant public classes. Functions which are deemed out-of-scope are excluded
from consideration. If a module or submodule contains no in-scope functions, it
is excluded from the tables. For example, `scipy.spatial.transform` is currently
excluded because it's API contains no functions, but may be included in the future
when the scope expands to include classes. `scipy.odr` and `scipy.datasets` are excluded
because their contents are considered out-of-scope.

.. toctree::
   :hidden:

   array_api_modules_tables/cluster_vq
   array_api_modules_tables/cluster_hierarchy
   array_api_modules_tables/constants
   array_api_modules_tables/differentiate
   array_api_modules_tables/fft
   array_api_modules_tables/integrate
   array_api_modules_tables/interpolate
   array_api_modules_tables/io
   array_api_modules_tables/linalg
   array_api_modules_tables/linalg_interpolative
   array_api_modules_tables/ndimage
   array_api_modules_tables/optimize
   array_api_modules_tables/optimize_elementwise
   array_api_modules_tables/signal
   array_api_modules_tables/signal_windows
   array_api_modules_tables/sparse
   array_api_modules_tables/sparse_linalg
   array_api_modules_tables/sparse_csgraph
   array_api_modules_tables/spatial
   array_api_modules_tables/spatial_distance
   array_api_modules_tables/special
   array_api_modules_tables/stats
   array_api_modules_tables/stats_contingency
   array_api_modules_tables/stats_qmc

Support on CPU
``````````````

.. array-api-support-per-module::
   :backend_type: cpu
   :cluster.vq: array_api_support_cluster_vq_cpu
   :cluster.hierarchy: array_api_support_cluster_hierarchy_cpu
   :constants: array_api_support_constants_cpu
   :differentiate: array_api_support_differentiate_cpu
   :fft: array_api_support_fft_cpu
   :integrate: array_api_support_integrate_cpu
   :interpolate: array_api_support_interpolate_cpu
   :io: array_api_support_io_cpu
   :linalg: array_api_support_linalg_cpu
   :linalg.interpolative: array_api_support_linalg_interpolative_cpu
   :ndimage: array_api_support_ndimage_cpu
   :optimize: array_api_support_optimize_cpu
   :optimize.elementwise: array_api_support_optimize_elementwise_cpu
   :signal: array_api_support_signal_cpu
   :signal.windows: array_api_support_signal_windows_cpu
   :sparse: array_api_support_sparse_cpu
   :sparse.linalg: array_api_support_sparse_linalg_cpu
   :sparse.csgraph: array_api_support_sparse_csgraph_cpu
   :spatial: array_api_support_spatial_cpu
   :spatial.distance: array_api_support_spatial_distance_cpu
   :special: array_api_support_special_cpu
   :stats: array_api_support_stats_cpu
   :stats.contingency: array_api_support_stats_contingency_cpu
   :stats.qmc: array_api_support_stats_qmc_cpu

Support on GPU
``````````````

.. array-api-support-per-module::
   :backend_type: gpu
   :cluster.vq: array_api_support_cluster_vq_gpu
   :cluster.hierarchy: array_api_support_cluster_hierarchy_gpu
   :constants: array_api_support_constants_gpu
   :differentiate: array_api_support_differentiate_gpu
   :fft: array_api_support_fft_gpu
   :integrate: array_api_support_integrate_gpu
   :interpolate: array_api_support_interpolate_gpu
   :io: array_api_support_io_gpu
   :linalg: array_api_support_linalg_gpu
   :linalg.interpolative: array_api_support_linalg_interpolative_gpu
   :ndimage: array_api_support_ndimage_gpu
   :optimize: array_api_support_optimize_gpu
   :optimize.elementwise: array_api_support_optimize_elementwise_gpu
   :signal: array_api_support_signal_gpu
   :signal.windows: array_api_support_signal_windows_gpu
   :sparse: array_api_support_sparse_gpu
   :sparse.linalg: array_api_support_sparse_linalg_gpu
   :sparse.csgraph: array_api_support_sparse_csgraph_gpu
   :spatial: array_api_support_spatial_gpu
   :spatial.distance: array_api_support_spatial_distance_gpu
   :special: array_api_support_special_gpu
   :stats: array_api_support_stats_gpu
   :stats.contingency: array_api_support_stats_contingency_gpu
   :stats.qmc: array_api_support_stats_qmc_gpu

Support with JIT
````````````````

.. array-api-support-per-module::
   :backend_type: jit
   :cluster.vq: array_api_support_cluster_vq_jit
   :cluster.hierarchy: array_api_support_cluster_hierarchy_jit
   :constants: array_api_support_constants_jit
   :differentiate: array_api_support_differentiate_jit
   :fft: array_api_support_fft_jit
   :integrate: array_api_support_integrate_jit
   :interpolate: array_api_support_interpolate_jit
   :io: array_api_support_io_jit
   :linalg: array_api_support_linalg_jit
   :linalg.interpolative: array_api_support_linalg_interpolative_jit
   :ndimage: array_api_support_ndimage_jit
   :optimize: array_api_support_optimize_jit
   :optimize.elementwise: array_api_support_optimize_elementwise_jit
   :signal: array_api_support_signal_jit
   :signal.windows: array_api_support_signal_windows_jit
   :sparse: array_api_support_sparse_jit
   :sparse.linalg: array_api_support_sparse_linalg_jit
   :sparse.csgraph: array_api_support_sparse_csgraph_jit
   :spatial: array_api_support_spatial_jit
   :spatial.distance: array_api_support_spatial_distance_jit
   :special: array_api_support_special_jit
   :stats: array_api_support_stats_jit
   :stats.contingency: array_api_support_stats_contingency_jit
   :stats.qmc: array_api_support_stats_qmc_jit
