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
This roadmap will be evolving together with SciPy.  Updates can be submitted as
pull requests.  For large or disruptive changes you may want to discuss
those first on the scipy-dev mailing list.


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


Use of venerable Fortran libraries
``````````````````````````````````
SciPy owes a lot of its success to relying on wrapping well established
Fortran libraries (QUADPACK, FITPACK, ODRPACK, ODEPACK etc). Some of these
libraries are aging well, others less so. We should audit our use of these
libraries with respect to the maintenance effort, the functionality, and the
existence of (possibly partial) alternatives, *including those inside SciPy*.


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
may help are separating out the bundled` ``libopenblas`` and removing support
for ``long double``.


Modules
-------

cluster
```````
``dendrogram`` needs a rewrite, it has a number of hard to fix open issues and
feature requests.


constants
`````````
This module is basically done, low-maintenance and without open issues.


fft
````
This module is in good shape.


integrate
`````````
Needed for ODE solvers:

- Documentation is pretty bad, needs fixing
- A new ODE solver interface  (``solve_ivp``) was added in SciPy 1.0.0.
  In the future we can consider (soft-)deprecating the older API.

The numerical integration functions are in good shape.  Support for integrating
complex-valued functions and integrating multiple intervals (see `gh-3325
<https://github.com/scipy/scipy/issues/3325>`__) could be added.


interpolate
```````````

*Spline fitting*: we need spline fitting routines with better user control. This
includes 

    - user-selectable alternatives for the smoothing criteria (manual,
      cross-validation etc); gh-16653 makes a start in this direction;

    - several strategies for knot placement, both manual and automatic (using
      algorithms by Dierckx, de Boor, possibly other). 

Once we have a reasonably feature complete set, we can start taking a long look
at the future of the venerable FITPACK Fortran library, which currently is the
only way of constructing smoothing splines in SciPy.

*Tensor-product splines*: `RegularGridInterpolator` provides a minimal
implementation. We want to evolve it both for new features (e.g. derivatives),
performance and API (possibly provide a transparent N-dimensional tensor-product
B-spline object).

*Scalability and performance*: For the FITPACK-based functionality, the data
size is limited by 32-bit Fortran integer size (for non-ILP64 builds).
For N-D scattered interpolators (which are QHull based) and N-D regular grid
interpolators we need to check performance on large data sets and improve
where lacking (gh-16483 makes progress in this direction).

*Ideas for new features*: NURBS support could be added.


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

**BLAS and LAPACK**

The Python and Cython interfaces to BLAS and LAPACK in ``scipy.linalg`` are one
of the most important things that SciPy provides. In general ``scipy.linalg``
is in good shape, however we can make a number of improvements:

1. Library support. Our released wheels now ship with OpenBLAS, which is
   currently the only feasible performant option (ATLAS is too slow, MKL cannot
   be the default due to licensing issues, Accelerate support is dropped
   because Apple doesn't update Accelerate anymore). OpenBLAS isn't very stable
   though, sometimes its releases break things and it has issues with threading
   (currently the only issue for using SciPy with PyPy3).  We need at the very
   least better support for debugging OpenBLAS issues, and better documentation
   on how to build SciPy with it.  An option is to use BLIS for a BLAS
   interface (see `numpy gh-7372 <https://github.com/numpy/numpy/issues/7372>`__).

2. Support for newer LAPACK features.  In SciPy 1.2.0 we increased the minimum
   supported version of LAPACK to 3.4.0.  Now that we dropped Python 2.7, we
   can increase that version further (MKL + Python 2.7 was the blocker for
   >3.4.0 previously) and start adding support for new features in LAPACK.


misc
````
``scipy.misc`` will be removed as a public module.  Most functions in it have
been moved to another submodule or deprecated.  The few that are left:

- ``derivative``, ``central_diff_weight`` : remove, possibly replacing them
  with more extensive functionality for numerical differentiation.
- ``ascent``, ``face``, ``electrocardiogram`` : remove or move to the
  appropriate subpackages (e.g. ``scipy.ndimage``, ``scipy.signal``).


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
Overall this module is in good shape. Two good global optimizers were added in
1.2.0; large-scale optimizers is still a gap that could be filled.  Other
things that are needed:

- Many ideas for additional functionality (e.g. integer constraints) in
  ``linprog``, see `gh-9269 <https://github.com/scipy/scipy/issues/9269>`__.
- Add functionality to the benchmark suite to compare results more easily
  (e.g. with summary plots).
- deprecate the ``fmin_*`` functions in the documentation, ``minimize`` is
  preferred.
- ``scipy.optimize`` has an extensive set of benchmarks for accuracy and speed of
  the global optimizers. That has allowed adding new optimizers (``shgo`` and
  ``dual_annealing``) with significantly better performance than the existing
  ones.  The ``optimize`` benchmark system itself is slow and hard to use
  however; we need to make it faster and make it easier to compare performance of
  optimizers via plotting performance profiles.


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

*Wavelets*: what's there now doesn't make much sense.  Continuous wavelets
only at the moment - decide whether to completely rewrite or remove them.
Discrete wavelet transforms are out of scope (PyWavelets does a good job
for those).


sparse
``````
The sparse matrix formats are mostly feature-complete, however the main issue
is that they act like ``numpy.matrix`` (which will be deprecated in NumPy at
some point).

What we want is sparse arrays, that act like ``numpy.ndarray``. In SciPy
``1.8.0`` a new set of classes (``csr_array`` et al.) has been added - these
need testing in the real world, as well as a few extra features like 1-D array
support.
An alternative (more ambitious, and unclear if it will materialize at this
point) plan is being worked on in https://github.com/pydata/sparse.  The
tentative plan for that was/is:

- Start depending on ``pydata/sparse`` once it's feature-complete enough (it
  still needs a CSC/CSR equivalent) and okay performance-wise.
- Add support for ``pydata/sparse`` to ``scipy.sparse.linalg`` (and perhaps to
  ``scipy.sparse.csgraph`` after that).
- Indicate in the documentation that for new code users should prefer
  ``pydata/sparse`` over sparse matrices.
- When NumPy deprecates ``numpy.matrix``, vendor that or maintain it as a
  stand-alone package.

Regarding the different sparse matrix formats: there are a lot of them.  These
should be kept, but improvements/optimizations should go into CSR/CSC, which
are the preferred formats.  LIL may be the exception, it's inherently
inefficient.  It could be dropped if DOK is extended to support all the
operations LIL currently provides.


sparse.csgraph
``````````````
This module is in good shape.


sparse.linalg
`````````````
There are a significant number of open issues for ``_arpack`` and ``lobpcg``.
``_propack`` is new in 1.8.0, TBD how robust it will turn out to be.

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

In addition to the items described in the :ref:`scipy-roadmap`, the following
improvements will help SciPy better serve this role.

- Add fundamental and widely used hypothesis tests, such as:

  - post hoc tests (e.g. Dunnett's test)
  - the various types of analysis of variance (ANOVA):

    - two-way ANOVA (single replicate, uniform number of replicates, variable
      number of replicates)
    - multiway ANOVA (i.e. generalize two-way ANOVA)
    - nested ANOVA
    - analysis of covariance (ANCOVA)

  Also, provide an infrastructure for implementing hypothesis tests.
- Add additional tools for meta-analysis
- Add tools for survival analysis
- Speed up random variate sampling (method ``rvs``) of distributions, 
  leveraging ``scipy.stats.sampling`` where appropriate
- Expand QMC capabilities and performance
- Enhance the `fit` method of the continuous probability distributions:

  - Expand the options for fitting to include:

    - maximal product spacings
    - method of L-moments / probability weighted moments

  - Include measures of goodness-of-fit in the results
  - Handle censored data (e.g. merge `gh-13699 <https://github.com/scipy/scipy/pull/13699>`__)

- Implement additional widely used continuous and discrete probability
  distributions, e.g. mixture distributions.

- Improve the core calculations provided by SciPy's probability distributions
  so they can robustly handle wide ranges of parameter values.  Specifically,
  replace many of the PDF and CDF methods from the Fortran library CDFLIB
  used in scipy.special with Boost implementations as in
  `gh-13328 <https://github.com/scipy/scipy/pull/13328>`__.

In addition, we should:

- Continue work on making the function signatures of ``stats`` and
  ``stats.mstats`` more consistent, and add tests to ensure that that
  remains the case.
- Improve statistical tests: return confidence intervals for the test
  statistic, and implement exact p-value calculations - considering the
  possibility of ties - where computationally feasible.
