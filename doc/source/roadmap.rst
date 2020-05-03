.. _scipy-roadmap:

SciPy Roadmap
=============

This roadmap page contains only the most important ideas and needs for SciPy
going forward.  For a more detailed roadmap, including per-subpackage status,
many more ideas, API stability and more, see :ref:`scipy-roadmap-detailed`.


Implement sparse arrays in addition to sparse matrices
------------------------------------------------------

The sparse matrix formats are mostly feature-complete, however the main issue
is that they act like ``numpy.matrix`` (which will be deprecated in NumPy at
some point).  What we want is sparse *arrays* that act like ``numpy.ndarray``.
This is being worked on in https://github.com/pydata/sparse, which is quite far
along.  The tentative plan is:

- Start depending on ``pydata/sparse`` once it's feature-complete enough (it
  still needs a CSC/CSR equivalent) and okay performance-wise.
- Add support for ``pydata/sparse`` to ``scipy.sparse.linalg`` (and perhaps to
  ``scipy.sparse.csgraph`` after that).
- Indicate in the documentation that for new code users should prefer
  ``pydata/sparse`` over sparse matrices.
- When NumPy deprecates ``numpy.matrix``, vendor that or maintain it as a
  stand-alone package.


Support for distributed arrays and GPU arrays
---------------------------------------------

NumPy is splitting its API from its execution engine with
``__array_function__`` and ``__array_ufunc__``.  This will enable parts of SciPy
to accept distributed arrays (e.g. ``dask.array.Array``) and GPU arrays (e.g.
``cupy.ndarray``) that implement the ``ndarray`` interface.  At the moment it is
not yet clear which algorithms will work out of the box, and if there are
significant performance gains when they do.  We want to create a map of which
parts of the SciPy API work, and improve support over time.

In addition to making use of NumPy protocols like ``__array_function__``, we can
make use of these protocols in SciPy as well.  That will make it possible to
(re)implement SciPy functions like, e.g., those in ``scipy.signal`` for Dask
or GPU arrays (see
`NEP 18 - use outside of NumPy <http://www.numpy.org/neps/nep-0018-array-function-protocol.html#use-outside-of-numpy>`__).


Statistics enhancements
-----------------------

The `scipy.stats` enhancements listed in the :ref:`scipy-roadmap-detailed` are of
particularly high importance to the project.

- Improve the options for fitting a probability distribution to data.
- Expand the set of hypothesis tests.  In particular, include all the basic
  variations of analysis of variance.
- Add confidence intervals for all statistical tests.
