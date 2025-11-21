.. _array_api_support_signal:

Array API Standard Support: ``signal``
======================================
This page explains some caveats of the `~scipy.signal` module and provides (currently
incomplete) tables about the
:ref:`CPU <array_api_support_signal_cpu>`,
:ref:`GPU <array_api_support_signal_gpu>` and
:ref:`JIT <array_api_support_signal_jit>` support.


.. _array_api_support_signal_caveats:

Caveats
-------
`JAX <https://docs.jax.dev/en/latest/jax.scipy.html>`__ and `CuPy
<https://docs.cupy.dev/en/stable/reference/scipy_signal.html>`__ provide alternative
implementations for some `~scipy.signal` functions. When such a function is called, a
decorator decides which implementation to use by inspecting the `xp` parameter.

Hence, there can be, especially during CI testing, discrepancies in behavior between
the default NumPy-based implementation and the JAX and CuPy backends. Skipping the
incompatible backends in unit tests, as described in the
:ref:`dev-arrayapi_adding_tests` section, is the currently recommended workaround.

The functions are decorated inline with ``_dispatchable`` from
``scipy/signal/_support_alternative_backends.py``:

.. literalinclude:: ../../../../../scipy/signal/_support_alternative_backends.py
    :lineno-match:

which is implemented via [``spatch``](https://scientific-python.github.io/spatch/api/for_libraries.html).

If ``SCIPY_ARRAY_API`` is set, this will set up dispatching based on the
signature and whether ``cupy=`` and/or ``jax=`` are passed in the decorator.
The decorator consists of:

* The signature description, either listing all array parameters or a callable
  which returns a tuple of array parameters.
* ``cupy=`` argument which defaults to ``True`` and signals if CuPy supports
  the function. It can be set to a string if CuPy uses a different name.
* ``jax=`` argument which defaults to ``False`` and has to be set
  to ``True`` or a string if the function is supported on JAX.

Other arguments are forwarded to ``spatch``.

.. _array_api_support_signal_cpu:

Support on CPU
--------------

.. array-api-support-per-function::
    :module: signal
    :backend_type: cpu

.. _array_api_support_signal_gpu:

Support on GPU
--------------

.. array-api-support-per-function::
    :module: signal
    :backend_type: gpu

.. _array_api_support_signal_jit:

Support with JIT
----------------

.. array-api-support-per-function::
    :module: signal
    :backend_type: jit
