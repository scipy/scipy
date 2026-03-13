.. _scipy-roadmap-detailed:

Detailed SciPy Roadmap
======================

Most of this roadmap is intended to provide a high-level view on what is
most needed per SciPy submodule in terms of new functionality, bug fixes, etc.
Besides important "business as usual" changes, it contains ideas for major new
features - those are marked as such, and are expected to take significant
dedicated effort.  Things not mentioned in this roadmap are
not necessarily unimportant or out of scope, however we (the SciPy developers)
want to provide to our users and contributors a clear picture of where SciPy is
going and where help is needed most.

.. note:: This is the detailed roadmap.  A very high-level overview with only
   the most important ideas is :ref:`scipy-roadmap`.


General
-------
This roadmap will be evolving together with SciPy. Updates can be submitted as
pull requests. For large or disruptive changes you may want to discuss
those first on the scipy-dev forum.


API changes
```````````
In general, we want to evolve the API to remove known warts as much as possible,
*however as much as possible without breaking backwards compatibility*.


Test coverage
`````````````
Test coverage of code added in the last few years is quite good, and we aim for
a high coverage for all new code that is added.  However, there is still a
significant amount of old code for which coverage is poor.  Bringing that up to
the current standard is probably not realistic, but we should plug the biggest
holes.

Besides coverage there is also the issue of correctness - older code may have a
few tests that provide decent statement coverage, but that doesn't necessarily
say much about whether the code does what it says on the box.  Therefore code
review of some parts of the code (``stats``, ``signal`` and ``ndimage`` in
particular) is necessary.


Documentation
`````````````
The documentation is in good shape. Expanding of current docstrings - adding
examples, references, and better explanations - should continue.  Most modules
also have a tutorial in the reference guide that is a good introduction,
however there are a few missing or incomplete tutorials - this should be fixed.


Benchmarks
``````````
The ``asv``-based benchmark system is in reasonable shape.  It is quite easy to
add new benchmarks, however running the benchmarks is not very intuitive.
Making this easier is a priority.


Use of Cython
`````````````
Cython's old syntax for using NumPy arrays should be removed and replaced with
Cython memoryviews.

Binary sizes of extensions built from Cython code are large, and compile times
are long. We should aim to combine extension modules where possible (e.g.,
``stats._boost`` contains many extension modules now), and limit the use of
Cython to places where it's the best choice. Note that conversion of Cython
to C++ is ongoing in ``scipy.special``.


Use of Pythran
``````````````
Pythran is still an optional build dependency, and can be disabled with
``-Duse-pythran=false``. The aim is to make it a hard dependency - for that to
happen it must be clear that the maintenance burden is low enough.


Use of Fortran libraries
````````````````````````
SciPy owes a lot of its success to relying on wrapping well established
Fortran libraries (QUADPACK, FITPACK, ODRPACK, ODEPACK etc). The Fortran 77
that these libraries are written in is quite hard to maintain, and the use
of Fortran is problematic for many reasons; e.g., it makes our wheel builds
much harder to maintain, it has repeatedly been problematic for supporting
new platforms like macOS arm64 and Windows on Arm, and it is highly problematic
for Pyodide's SciPy support. Our goal is to remove all Fortran code from SciPy
by replacing the functionality with code written in other languages. Progress
towards this goal is tracked in `gh-18566 <https://github.com/scipy/scipy/issues/18566>`__.


Continuous integration
``````````````````````
Continuous integration currently covers 32/64-bit Windows, macOS on x86-64/arm,
32/64-bit Linux on x86, and Linux on aarch64 - as well as a range of versions
of our dependencies and building release quality wheels. Reliability of CI has
not been good recently (H1 2023), due to the large amount of configurations to
support and some CI jobs needing an overhaul. We aim to reduce build times by
removing the remaining distutils-based jobs when we drop that build system
and make the set of configurations in CI jobs more orthogonal.


Size of binaries
````````````````
SciPy binaries are quite large (e.g. an unzipped manylinux wheel for 1.7.3 is
39 MB on PyPI and 122 MB after installation), and this can be problematic - for
example for use in AWS Lambda, which has a 250 MB size limit. We aim to keep
binary size as low as possible; when adding new compiled extensions, this needs
checking. Stripping of debug symbols in ``multibuild`` can perhaps be improved
(see `this issue <https://github.com/multi-build/multibuild/issues/162>`__).
An effort should be made to slim down where possible, and not add new large
files. In the future, things that are being considered (very tentatively) and
may help are separating out the bundled ``libopenblas`` and removing support
for ``long double``.


Modules
-------

cluster
```````
``dendrogram`` needs a rewrite, it has a number of hard to fix open issues and
feature requests.


constants
`````````
This module needs only periodic updates to the numerical values.

differentiate
`````````````
This module was added with the understanding that its scope would be limited.
The goal is to provide basic differentiation of black-box functions, not
replicate all features of existing differentiation tools. With that in mind,
items for future work include:

- Improve support for callables that accept additional arguments (e.g. add
  ``*args`` to ``jacobian`` and ``hessian``). Note that this is not trivial
  due to the way elements of the arrays are eliminated when their corresponding
  calculations have converged.
- Improve implementation of `scipy.differentiate.hessian`: rather than chaining
  first-order differentiation, use a second-order differentiation stencil.
- Consider the addition of an option to use relative step sizes.
- Consider generalizing the approach to use "polynomial extrapolation"; i.e.,
  rather than estimating derivatives of a given order from the minimal number
  of function evaluations, use a least-squares approach to improve robustness
  to numerical noise.

fft
````
This module is in good shape.


integrate
`````````
- Complete the feature set of the new ``cubature`` function, and add an interface
  tailored to one-dimensional integrals.
- Migrate functions for generating quadrature rule points and weights from `special`,
  improve their reliability, and add support for other important rules.
- Complete the feature set of ``solve_ivp``, including all functionality of the
  ``odeint`` interface.
- Gracefully sunset superseded functions and classes.


interpolate
```````````


*FITPACK replacements*

We need to finish providing functional equivalents of the FITPACK Fortran library.
The aim is to reimplement the algorithmic content of the FITPACK monolith in a
modular fashion, to allow better user control and extensibility.
To this end, we need

- 1D splines: add the periodic spline fitting to `make_splrep` function.
- 2D splines: provide functional replacements of the `BivariateSpline` family
  of classes for scattered and gridded data.  

Once we have a feature-complete set of replacements, we can gracefully sunset
our use of FITPACK.

*Ideas for new features*


- add the ability to construct smoothing splines with variable number of knots
  to `make_smoothing_spline`. Currently, `make_smoothing_spline` always uses all
  data points for the knot vector.
- experiment with user-selectable smoothing criteria: add the P-splines penalty
  form to complement the 2nd derivative penalty of `make_smoothing_spline` and
  jumps of the 3rd derivative of `make_splrep`
- investigate ways of vectorizing the spline construction for batches of the
  independent variable. This class of features has been requested by
  `scipy.stats` developers
- allow specifying the boundary conditions for spline fitting.
- investigate constrained spline fitting problems which the FITPACK library was
  providing in Fortran, but were never exposed to the python interfaces
- NURBS support can be added if there is sufficient user interest

*Scalability and performance*

We need to check performance on large data sets and improve where lacking. This
is relevant for both regular grid N-D interpolators and QHull-based scattered
N-D interpolators. For the latter ones, we additionally need to investigate if
we can improve on the 32-bit integer limitations of the QHull library.


io
``
wavfile:

- PCM float will be supported, for anything else use ``audiolab`` or other
  specialized libraries.
- Raise errors instead of warnings if data not understood.

Other sub-modules (matlab, netcdf, idl, harwell-boeing, arff, matrix market)
are in good shape.


linalg
``````
``scipy.linalg`` is in good shape.

Needed:

- Reduce duplication of functions with ``numpy.linalg``, make APIs consistent.
- ``get_lapack_funcs`` should always use ``flapack``
- Wrap more LAPACK functions
- One too many funcs for LU decomposition, remove one

Ideas for new features:

- Add type-generic wrappers in the Cython BLAS and LAPACK
- Make many of the linear algebra routines into gufuncs
- Complete support for batched operations (see
  `gh-21466 <https://github.com/scipy/scipy/issues/21466>`__)

**BLAS and LAPACK**

The Python and Cython interfaces to BLAS and LAPACK in ``scipy.linalg`` are one
of the most important things that SciPy provides. In general ``scipy.linalg``
is in good shape, however we can make a number of improvements:

1. Add support for ILP64 (64-bit) BLAS and LAPACK (see
   `gh-21889 <https://github.com/scipy/scipy/issues/21889>`__)
2. Unify the two sets of low-level BLAS/LAPACK wrappers, probably dropping the
   ``f2py``-based ones (see
   `gh-20682 <https://github.com/scipy/scipy/issues/20682>`__)
3. Improve and document the various ways we link to BLAS and LAPACK from C
   and C++ code internally in SciPy (see
   `gh-20002 <https://github.com/scipy/scipy/issues/20002>`__ and
   `gh-21130 <https://github.com/scipy/scipy/pull/21130>`__)


misc
````
All features have been removed from ``scipy.misc``, and the namespace itself
will eventually be removed.


ndimage
```````
Underlying ``ndimage`` is a powerful interpolation engine.  Users come
with an expectation of one of two models: a pixel model with ``(1,
1)`` elements having centers ``(0.5, 0.5)``, or a data point model,
where values are defined at points on a grid.  Over time, we've become
convinced that the data point model is better defined and easier to
implement, but this should be clearly communicated in the documentation.

More importantly, still, SciPy implements one *variant* of this data
point model, where datapoints at any two extremes of an axis share a
spatial location under *periodic wrapping* mode.  E.g., in a 1D array,
you would have ``x[0]`` and ``x[-1]`` co-located.  A very common
use-case, however, is for signals to be periodic, with equal spacing
between the first and last element along an axis (instead of zero
spacing).  Wrapping modes for this use-case were added in
`gh-8537 <https://github.com/scipy/scipy/pull/8537>`__, next the
interpolation routines should be updated to use those modes.
This should address several issues, including gh-1323, gh-1903, gh-2045
and gh-2640.

The morphology interface needs to be standardized:

- binary dilation/erosion/opening/closing take a "structure" argument,
  whereas their grey equivalent take size (has to be a tuple, not a scalar),
  footprint, or structure.
- a scalar should be acceptable for size, equivalent to providing that same
  value for each axis.
- for binary dilation/erosion/opening/closing, the structuring element is
  optional, whereas it's mandatory for grey.  Grey morphology operations
  should get the same default.
- other filters should also take that default value where possible.


odr
```
This module is in reasonable shape, although it could use a bit more
maintenance.  No major plans or wishes here.


optimize
````````
We aim to continuously improve the set of optimizers provided by this module.
For large scale problems, the state of the art continues to advance and we aim
to keep up by leveraging implementations from domain-specific libraries like
HiGHS and PRIMA. Other areas for future work include the following.

- Improve the interfaces of existing optimizers (e.g. ``shgo``).
- Improve usability of the benchmark system, and add features for comparing
  results more easily (e.g. summary plots).

signal
``````
*Convolution and correlation*: (Relevant functions are convolve, correlate,
fftconvolve, convolve2d, correlate2d, and sepfir2d.) Eliminate the overlap with
`ndimage` (and elsewhere).  From ``numpy``, ``scipy.signal`` and ``scipy.ndimage``
(and anywhere else we find them), pick the "best of class" for 1-D, 2-D and n-d
convolution and correlation, put the implementation somewhere, and use that
consistently throughout SciPy.

*B-splines*: (Relevant functions are gauss_spline,
cspline1d, qspline1d, cspline2d, qspline2d, cspline1d_eval, and spline_filter.)
Move the good stuff to `interpolate` (with appropriate API changes to match how
things are done in `interpolate`), and eliminate any duplication.

*Filter design*: merge `firwin` and `firwin2` so `firwin2` can be removed.

*Continuous-Time Linear Systems*: Further improve the performance of ``ltisys``
(fewer internal transformations between different representations). Fill gaps in lti
system conversion functions.

*Second Order Sections*: Make SOS filtering equally capable as existing
methods. This includes ltisys objects, an `lfiltic` equivalent, and numerically
stable conversions to and from other filter representations. SOS filters could
be considered as the default filtering method for ltisys objects, for their
numerical stability.


sparse
``````
The sparse matrix formats are mostly feature-complete, however the main issue
is that they act like ``numpy.matrix`` (which will be deprecated in NumPy at
some point).

What we want is sparse arrays that act like ``numpy.ndarray``. Initial
support for a new set of classes (``csr_array`` et al.) was added in SciPy
``1.8.0`` and stabilized in ``1.12.0`` with construction functions for
arrays, ``1.14.0`` with 1D array support and ``1.15.0`` with 1D indexing.
The sparse array codebase now supports all sparse matrix features and in
addition supports 1D arrays and the first steps toward nD arrays.
There is a transition guide to help users and libraries convert their code
to sparse arrays.

Next steps toward sparse array conversion:

- Extend sparse array API to nD arrays.
    - Support for COO, CSR and DOK formats.
    - Some COO features exist in 1.15.
- Introduce support for broadcasting in operations where sparse formats
  can effectively do that.
- Help other libraries convert to sparse arrays from sparse matrices.
  NetworkX, dipy, scikit-image, pyamg, cvxpy and scikit-learn are in
  progress or have completed conversion to sparse arrays.
- Add deprecation warnings for sparse matrix.
- Work with NumPy on deprecation/removal of ``numpy.matrix``.
- Deprecate and then remove sparse matrix in favor of sparse array.
- Start API shift of construction function names (``diags``, ``block``, etc.)
    - Note: as a whole, the construction functions undergo two name shifts.
      Once to move from matrix creation to new functions for array creation
      (i.e. ``eye`` -> ``eye_array``). Then a second move to change names to match
      the ``array_api`` name (i.e. ``eye_array`` to ``eye``) after sparse matrices are
      removed. We will keep the ``*_array`` versions for a long time as
      (maybe hidden) aliases.
- Add construction function names matching ``array_api`` names.
- Deprecate the transition construction function names.

An alternative (more ambitious, and unclear if it will materialize)
plan is being worked on in https://github.com/pydata/sparse.
To support that effort we aim to support PyData Sparse in all functions that
accept sparse arrays.  Support for ``pydata/sparse`` in ``scipy.sparse.linalg``
and ``scipy.sparse.csgraph`` is mostly complete.

Regarding the different sparse matrix formats: there are a lot of them.  These
should be kept, but improvements/optimizations should go into CSR/CSC, which
are the preferred formats. LIL may be the exception, it's inherently
inefficient. It could be dropped if DOK is extended to support all the
operations LIL currently provides.


sparse.csgraph
``````````````
This module is in good shape.


sparse.linalg
`````````````
There are a significant number of open issues for ``_arpack`` and ``lobpcg``.

``_isolve``:

- callback keyword is inconsistent
- tol keyword is broken, should be relative tol
- Fortran code not re-entrant (but we don't solve, maybe reuse from
  PyKrilov)

``_dsolve``:

- add license-compatible sparse Cholesky or incomplete Cholesky
- add license-compatible sparse QR
- improve interface to SuiteSparse UMFPACK
- add interfaces to SuiteSparse CHOLMOD and SPQR


spatial
```````
QHull wrappers are in good shape, as is ``KDTree``.

A rewrite of ``spatial.distance`` metrics in C++ is in progress - this should
improve performance, make behaviour (e.g., for various non-float64 input
dtypes) more consistent, and fix a few remaining issues with definitions of the
math implement by a few of the metrics.


special
```````
Though there are still a lot of functions that need improvements in precision,
probably the only show-stoppers are hypergeometric functions, parabolic cylinder
functions, and spheroidal wave functions. Three possible ways to handle this:

1. Get good double-precision implementations. This is doable for parabolic
   cylinder functions (in progress). I think it's possible for hypergeometric
   functions, though maybe not in time. For spheroidal wavefunctions this is
   not possible with current theory.

2. Port Boost's arbitrary precision library and use it under the hood to get
   double precision accuracy. This might be necessary as a stopgap measure
   for hypergeometric functions; the idea of using arbitrary precision has
   been suggested before by @nmayorov and in
   `gh-5349 <https://github.com/scipy/scipy/issues/5349>`__.  Likely
   necessary for spheroidal wave functions, this could be reused:
   https://github.com/radelman/scattering.

3. Add clear warnings to the documentation about the limits of the existing
   implementations.


stats
`````

The ``scipy.stats`` subpackage aims to provide fundamental statistical
methods as might be covered in standard statistics texts such as Johnson's
"Miller & Freund's Probability and Statistics for Engineers", Sokal & Rohlf's
"Biometry", or Zar's "Biostatistical Analysis".  It does not seek to duplicate
the advanced functionality of downstream packages (e.g. StatsModels,
LinearModels, PyMC3); instead, it can provide a solid foundation on which
they can build.  (Note that these are rough guidelines, not strict rules.
"Advanced" is an ill-defined and subjective term, and "advanced" methods
may also be included in SciPy, especially if no other widely used and
well-supported package covers the topic.  Also note that *some* duplication
with downstream projects is inevitable and not necessarily a bad thing.)

The following improvements will help SciPy better serve this role.

- Improve statistical tests: include methods for generating confidence
  intervals, and implement exact and randomized p-value calculations -
  considering the possibility of ties - where computationally feasible.
- Extend the new univariate distribution infrastructure, adding support
  for discrete distributions and circular continuous distributions.
  Add a selection of the most widely used statistical distributions
  under the new infrastructure, performing rigorous accuracy testing
  and documenting the results. Enable users to create custom distributions
  that leverage the new infrastructure.
- Improve the multivariate distribution infrastructure to ensure a
  consistent API, thorough testing, and complete documentation.
- Continue to make the APIs of functions more consistent, with standard
  support for batched calculations, NaN policies, and dtype preservation.
