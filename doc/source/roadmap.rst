.. _scipy-roadmap:

SciPy Roadmap
=============

This roadmap page contains only the most important ideas and needs for SciPy
going forward.  For a more detailed roadmap, including per-subpackage status,
many more ideas, API stability and more, see :ref:`scipy-roadmap-detailed`.


Evolve BLAS and LAPACK support
------------------------------

The Python and Cython interfaces to BLAS and LAPACK in ``scipy.linalg`` are one
of the most important things that SciPy provides. In general ``scipy.linalg``
is in good shape, however we can make a number of improvements:

1. Library support. Our released wheels now ship with OpenBLAS, which is
currently the only feasible performant option (ATLAS is too slow, MKL cannot be
the default due to licensing issues, Accelerate support is dropped because
Apple doesn't update Accelerate anymore). OpenBLAS isn't very stable though,
sometimes its releases break things and it has issues with threading (currently
the only issue for using SciPy with PyPy3).  We need at the very least better
support for debugging OpenBLAS issues, and better documentation on how to build
SciPy with it.  An option is to use BLIS for a BLAS interface (see `numpy
gh-7372 <https://github.com/numpy/numpy/issues/7372>`__).

2. Support for newer LAPACK features.  In SciPy 1.2.0 we increased the minimum
supported version of LAPACK to 3.4.0.  Now that we dropped Python 2.7, we can
increase that version further (MKL + Python 2.7 was the blocker for >3.4.0
previously) and start adding support for new features in LAPACK.


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


Fourier transform enhancements
------------------------------

The new ``scipy.fft`` subpackage should be extended to add a backend system with
support for PyFFTW and mkl-fft.


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


Improve source builds on Windows
--------------------------------

SciPy critically relies on Fortran code. This is still problematic on Windows.
There are currently only two options: using Intel Fortran, or using
MSVC + gfortran.  The former is expensive, while the latter works (it's what we
use for releases) but is quite hard to do correctly.  For allowing contributors
and end users to reliably build SciPy on Windows, using the Flang compiler
looks like the best way forward long-term.  Until Flang support materializes,
we need to streamline and better document the MSVC + gfortran build.


Statistics enhancements
-----------------------

The `scipy.stats` enhancements listed in the :ref:`scipy-roadmap-detailed` are of
particularly high importance to the project.

- Improve the options for fitting a probability distribution to data.
- Expand the set of hypothesis tests.  In particular, include all the basic
  variations of analysis of variance.
- Add confidence intervals for all statistical tests.
