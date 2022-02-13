Adding vectorized ufuncs in ``scipy.special``
=============================================

.. highlight:: none

Many of the functions in ``special`` are vectorized versions of scalar
functions. The scalar functions are written by hand and the necessary
loops for vectorization are generated automatically. This section
discusses the steps necessary to add a new vectorized special
function.

The first step in adding a new vectorized function is writing the
corresponding scalar function. This can be done in Cython, C, C++, or
Fortran. If starting from scratch then Cython should be preferred
because the code is easier to maintain for developers only familiar
with Python. If the primary code is in Fortran then it is necessary to
write a C wrapper around the code; for examples of such wrappers see
``specfun_wrappers.c``.

After implementing the scalar function, register the new function by
adding an entry to ``functions.json``. The docstring in
``generate_ufuncs.py`` explains the format. Also add documentation for
the new function by adding an entry to ``add_newdocs.py``; look in the
file for examples.

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
