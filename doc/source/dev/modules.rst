Module-Specific Instructions
============================

Some SciPy modules have specific development workflows that it is
useful to be aware of while contributing.

``scipy.special``
-----------------

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
adding a line to the ``FUNC`` string in ``generate_ufuncs.py``. The
docstring for that file explains the format. Also add documentation
for the new function by adding an entry to ``add_newdocs.py``; look in
the file for examples.
