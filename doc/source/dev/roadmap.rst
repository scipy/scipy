.. _scipy-roadmap:

SciPy Roadmap
=============

This roadmap page contains only the most important ideas and needs for SciPy
going forward.  For a more detailed roadmap, including per-subpackage status,
many more ideas, API stability and more, see :ref:`scipy-roadmap-detailed`.


Support for distributed arrays and GPU arrays
---------------------------------------------

NumPy has split its API from its execution engine with
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
`NEP 18 - use outside of NumPy <http://www.numpy.org/neps/nep-0018-array-function-protocol.html#use-outside-of-numpy>`__).  NumPy's features in this areas are still evolving,
see e.g. `NEP 37 - A dispatch protocol for NumPy-like modules <https://numpy.org/neps/nep-0037-array-module.html>`__,
and SciPy is an important "client" for those features.


Performance improvements
------------------------

Speed improvements, lower memory usage and the ability to parallelize
algorithms are beneficial to most science domains and use cases.  We have
established an API design pattern for multiprocessing - using the ``workers``
keyword - that can be adopted in many more functions.

Making it easier for users to use Numba's ``@njit`` in their code that relies
on SciPy functionality would unlock a lot of performance gain.  That needs a
strategy though, all solutions are still maturing (see for example
`this overview <https://fluiddyn.netlify.app/transonic-vision.html>`__).

Finally, many individual functions can be optimized for performance.
``scipy.optimize`` and ``scipy.interpolate`` functions are particularly often
requested in this respect.


Support for more hardware platforms
-----------------------------------

SciPy now has continuous integration for ARM64 (or ``aarch64``) and POWER8/9
(or ``ppc64le``), and binaries are available via
`Miniforge <https://github.com/conda-forge/miniforge>`__.  Wheels on PyPI for
these platforms are now also possible (with the ``manylinux2014`` standard),
and requests for those are becoming more frequent.

Additionally, having IBM Z (or ``s390x``) in CI is now possible with TravisCI
but not yet done - and ``manylinux2014`` wheels for that platform are also
possible then.  Finally, resolving open AIX build issues would help users.


Implement sparse arrays in addition to sparse matrices
------------------------------------------------------

SciPy sparse matrices are being replaced by sparse arrays.
The sparse matrix formats are mostly feature-complete, however their main issue
is that they act like ``numpy.matrix`` (which will be deprecated in NumPy at
some point). What we want is sparse *arrays* that act like ``numpy.ndarray``
(See discussion at `gh-18915 <https://github.com/scipy/scipy/issues/18915>`_).
Sparse arrays support all features of sparse matrices as of 1.15.
In addition to 2D arrays, 1D sparse arrays are supported in DOK, COO, CSR formats.
Further functionality e.g. nD array support and broadcasting for some operations
is being developed.  The future plan is:

- Extend sparse array API to nD arrays:
    - COO, CSR and DOK formats. COO format already partially in place.
    - The nD formats use 2D CSR code to do nD things like
      indexing/min-max/arithmetic.
- Sparse array binary operations will support broadcasting in some settings.
  Broadcasting is tricky for sparse arrays because it leans heavily on the strided
  memory model of dense arrays, and so does not always fit sparse data formats.
  Our optimistic goal is to support broadcasting for all operations where that
  makes sense for sparse data structures. We start with binary operations like `A + B`.
- Help other libraries convert to sparse arrays from sparse matrices.
  Create transition guide and helpful scripts to flag code that needs changing.
- Deprecate and then remove "sparse matrix" in favor of "sparse array".
- Work with NumPy on deprecation/removal of ``numpy.matrix``.
