.. _building-blas-and-lapack:

BLAS and LAPACK
===============

Selecting BLAS and LAPACK libraries
-----------------------------------

Using pkg-config to detect libraries in a nonstandard location
--------------------------------------------------------------


How do I deal with Fortran ABI mismatch?
----------------------------------------

Some linear algebra libraries are built with g77 ABI and others with
GFortran ABI, and these two ABIs are incompatible. Therefore, if you
build SciPy with ``gfortran`` and link to a linear algebra library, like
MKL, which is built with g77 ABI, then there'll be an exception or a
segfault. SciPy fixes this by using the CBLAS API for the few
functions in the BLAS API that suffers from this issue.

Note that SciPy needs to know at build time, what needs to be done and
the build system will automatically check whether linear algebra
library is MKL and if so, use the CBLAS API instead of the BLAS API.
If autodetection fails or if the user wants to override this
autodetection mechanism, use the following:

Use the ``-Duse-g77-abi=true`` build option. E.g.,::

    $ meson setup build -Duse-g77-abi=true

A more complete example, also configuring the BLAS/LAPACK libraries and picking
a better Python install behavior (this is what conda-forge could be using for
example)::

    $ meson setup builddir -Duse-g77-abi=true -Dblas=blas -Dlapack=lapack -Dpython.install_env=auto
    $ meson install -C builddir

Work-in-progress
----------------

These options are planned to be fully supported, but currently not usable out
of the box:

- ILP64 (64-bit integer size) builds: large parts of SciPy support using ILP64
  BLAS/LAPACK. Note that support is still incomplete, so SciPy *also* requires
  LP64 (32-bit integer size) BLAS/LAPACK.
- Automatically selecting from multiple possible BLAS and LAPACK options, with
  a user-provided order of precedence

