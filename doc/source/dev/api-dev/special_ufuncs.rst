Adding vectorized ufuncs in ``scipy.special``
=============================================

.. highlight:: none

Many of the functions in ``special`` are vectorized versions of scalar
functions. The scalar functions are written by hand and the necessary
loops for vectorization are generated automatically. This section
discusses the steps necessary to add a new vectorized special
function.

The first step in adding a new vectorized function is writing the
corresponding scalar function. This is done in the
`xsf repository <https://github.com/scipy/xsf>`_, see
`here <https://github.com/scipy/xsf/blob/main/CONTRIBUTING.md>`_ for further details.

After implementing the scalar function, create the ufunc in
``_special_ufuncs.cpp`` using ``xsf::numpy::ufunc``. Add the ufunc's
docstring to ``_special_ufuncs_docs.cpp``, and add its name to the
``special_ufuncs`` list in ``_generate_pyx.py``. For a public function,
also add wrappers to ``cython_special_wrappers.cpp`` and
``cython_special_wrappers.h``, and add the function to
``cython_special.pyx``. Look in these files for examples.

When writing the parameters section of the documentation for ufuncs,
the type of an argument should be ``array_like``. Discussion of
whether an argument can be e.g. real or complex-valued should be saved
for the description. So for example, if we were to document the
parameters for the Gamma function then it should look like this::

  Parameters
  ----------
  z : array_like
      Real or complex valued argument

When documenting the returns section, the type of the returned value
should be ``scalar or ndarray`` since ufuncs return scalars when given
scalars as arguments. Also keep in mind that providing a ``name`` for
the return value is optional, and indeed is often not helpful for
special functions. So for the Gamma function we might have something
like this::

  Returns
  -------
  scalar or ndarray
      Values of the Gamma function
