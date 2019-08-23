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

Also, it should be made (even) more clear what is public and what is private in
SciPy.  Everything private should be named starting with an underscore as much
as possible.


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
The documentation is in good shape.  Expanding of current docstrings and
putting them in the standard NumPy format should continue, so the number of
reST errors and glitches in the html docs decreases.  Most modules also have a
tutorial in the reference guide that is a good introduction, however there are
a few missing or incomplete tutorials - this should be fixed.


Benchmarks
``````````
The ``asv``-based benchmark system is in reasonable shape.  It is quite easy to
add new benchmarks, however running the benchmarks is not very intuitive.
Making this easier is a priority.  In addition, we should run them in our CI
(gh-8779 is an ongoing attempt at this).


Other
`````

Regarding Cython code:

- It's not clear how much functionality can be Cythonized without making the
  .so files too large.  This needs measuring.
- Cython's old syntax for using NumPy arrays should be removed and replaced
  with Cython memoryviews.

Regarding build environments:

- SciPy builds from source on Windows now with a MSVC + MinGW-w64 gfortran
  toolchain, which we're using for official releases.
  MSVC + Intel Fortran + MKL works as well, and is easier for users (as long
  as they have access to ifort and MKL of course).  This mainly needs better
  documentation at the moment.
- We're aiming to gradually increase the minimum version of LAPACK that is
  required, so we can use newer features.  Support for Accelerate on macOS
  has been dropped.  We do rely quite heavily on OpenBLAS, and its stability
  is a worry (often only one of the recent releases works without test
  failures) - improvements in testing and build documentation at least are
  needed.

Continuous integration is in good shape, it covers Windows, macOS and Linux, as well
as a range of versions of our dependencies and building release quality wheels.


Modules
-------

cluster
```````
This module is in good shape.


constants
`````````
This module is basically done, low-maintenance and without open issues.


fft
````

Ideas for new features:

- Add a backend/plugin system.  At the moment pyFFTW is monkeypatching SciPy,
  and ``mkl_fft`` provides ``fftpack``-compatible functions as well.  We should
  provide a method to support such packages.

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

Ideas for new features:

- Spline fitting routines with better user control.
- Transparent tensor-product splines.
- NURBS support.
- Mesh refinement and coarsening of B-splines and corresponding tensor products.

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
``scipy.linalg`` is in good shape.  We have started requiring more recent
LAPACK versions (minimum version increases from 3.1.0 to 3.4.0 in SciPy 1.2.0);
we want to add support for newer features in LAPACK.

Needed:

- Reduce duplication of functions with ``numpy.linalg``, make APIs consistent.
- ``get_lapack_funcs`` should always use ``flapack``
- Wrap more LAPACK functions
- One too many funcs for LU decomposition, remove one

Ideas for new features:

- Add type-generic wrappers in the Cython BLAS and LAPACK
- Make many of the linear algebra routines into gufuncs


misc
````
``scipy.misc`` will be removed as a public module.  Most functions in it have
been moved to another submodule or deprecated.  The few that are left:

- ``info``, ``who`` : these are NumPy functions
- ``derivative``, ``central_diff_weight`` : remove, possibly replacing them
  with more extensive functionality for numerical differentiation.


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

- Many ideas for additional functionality (e.g. integer constraints, sparse
  matrix support, performance improvements) in ``linprog``, see
  `gh-9269 <https://github.com/scipy/scipy/issues/9269>`__.
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

*B-splines*: (Relevant functions are bspline, cubic, quadratic, gauss_spline,
cspline1d, qspline1d, cspline2d, qspline2d, cspline1d_eval, and spline_filter.)
Move the good stuff to `interpolate` (with appropriate API changes to match how
things are done in `interpolate`), and eliminate any duplication.

*Filter design*: merge `firwin` and `firwin2` so `firwin2` can be removed.

*Continuous-Time Linear Systems*: remove `lsim2`, `impulse2`, `step2`.  The
`lsim`, `impulse` and `step` functions now "just work" for any input system.
Further improve the performance of ``ltisys`` (fewer internal transformations
between different representations). Fill gaps in lti system conversion functions.

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
some point).  What we want is sparse arrays, that act like ``numpy.ndarray``.
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
Arpack is in good shape.

isolve:

- callback keyword is inconsistent
- tol keyword is broken, should be relative tol
- Fortran code not re-entrant (but we don't solve, maybe re-use from
  PyKrilov)

dsolve:

- add sparse Cholesky or incomplete Cholesky
- look at CHOLMOD


Ideas for new features:

- Wrappers for PROPACK for faster sparse SVD computation.


spatial
```````
QHull wrappers are in good shape, as is ``cKDTree``.

Needed:

- ``KDTree`` will be removed, and ``cKDTree`` will be renamed to ``KDTree``
  in a backwards-compatible way.
- ``distance_wrap.c`` needs to be cleaned up (maybe rewrite in Cython).


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

- Add fundamental and widely used hypothesis tests:

  - Alexander-Govern test
  - Somers' D
  - Kendall's tau-c
  - Page's L-test
  - Tukey-Kramer test
  - Dunnett's test
  - the various types of analysis of variance (ANOVA):

    - two-way ANOVA (single replicate, uniform number of replicates, variable
      number of replicates)
    - multiway ANOVA (i.e. generalize two-way ANOVA)
    - nested ANOVA
    - analysis of covariance (ANCOVA)

- Where appropriate, include confidence intervals for the statistic in the
  results of any statistical test.
- Add additional tools for meta-analysis; currently we have just `combine_pvalues`.
- Enhance the `fit` method of the continuous probability distributions:

  - Expand the options for fitting to include:

    - method of moments
    - maximal product spacings
    - method of L-moments / probability weighted moments

  - Include measures of goodness-of-fit in the results
  - Handle censored data

- Implement additional widely used continuous and discrete probability
  distributions:

  - noncentral hypergeometric distribution (both Fisher's and Wallenius')
  - negative hypergeometric distribution
  - multivariate hypergeometric distribution
  - multivariate t distribution
  - mixture distributions

- Improve the core calculations provided by SciPy's probability distributions 
  so they can robustly handle wide ranges of parameter values.  Specifically,
  replace many of the PDF and CDF methods from the Fortran library CDFLIB
  used in scipy.special with better code, perhaps ported from the Boost C++
  library.
  
In addition, we should:

- Continue work on making the function signatures of ``stats`` and
  ``stats.mstats`` more consistent, and add tests to ensure that that
  remains the case.
- Return ``Bunch`` objects from functions that now return many values, and for
  functions for which extra return values are desired (see
  `gh-3665 <https://github.com/scipy/scipy/issues/3665>`__).
- Improve statistical tests (p-value calculation, alternative hypothesis), for
  example implement an exact two-sided KS test (see
  `gh-8341 <https://github.com/scipy/scipy/issues/8341>`__) or a one-sided
  Wilcoxon test (see `gh-9046 <https://github.com/scipy/scipy/issues/9046>`__).
- Address the various issues regarding ``stats.mannwhitneyu``, and pick up the 
  stalled PR in `gh-4933 <https://github.com/scipy/scipy/pull/4933>`__.
