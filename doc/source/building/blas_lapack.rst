.. _building-blas-and-lapack:

BLAS and LAPACK
===============

.. _blas-lapack-selection:

Selecting BLAS and LAPACK libraries
-----------------------------------

BLAS and LAPACK library selection, other than the OpenBLAS default, is
implemented via Meson `build options
<https://mesonbuild.com/Build-options.html#build-options>`__. For example, to
select plain ``libblas`` and ``liblapack`` (this is typically Netlib
BLAS/LAPACK on Linux distros, and can be dynamically switched between
implementations on conda-forge), use::

    $ # for a development build
    $ spin build -S-Dblas=blas -S-Dlapack=lapack

    $ # to build and install a wheel
    $ python -m build -Csetup-args=-Dblas=blas -Csetup-args=-Dlapack=lapack
    $ pip install dist/scipy*.whl

    $ # Or, with pip>=23.1, this works too:
    $ python -m pip -Csetup-args=-Dblas=blas -Csetup-args=-Dlapack=lapack

Other options that should work (as long as they're installed with
``pkg-config`` or CMake support) include:

- Accelerate: ``-Dblas=accelerate``,
- MKL: use the name of the relevant pkg-config file, e.g.
  ``-Dblas=mkl-dynamic-lp64-gomp``,
- BLIS: ``-Dblas=blis`` (is BLAS-only, so will need a ``-Dlapack=`` as well),
- ATLAS: ``-Dblas=blas-atlas`` (name of the pkg-config file may vary)

Note that both Accelerate and ``scipy-openblas`` have flags in ``spin``
that are easier to remember, since they're commonly used for development::

    $ spin build --with-accelerate
    $ spin build --with-scipy-openblas=32

The ``-Dlapack`` flag isn't needed for Accelerate, MKL or ``scipy-openblas``,
since we can be sure that BLAS and LAPACK are the same for those options.
E.g., to create a wheel with Accelerate, use::

    $ python -m build -Csetup-args=-Dblas=accelerate


Using pkg-config to detect libraries in a nonstandard location
--------------------------------------------------------------

The way BLAS and LAPACK detection works under the hood is that Meson tries
to discover the specified libraries first with ``pkg-config``, and then
with CMake. If all you have is a standalone shared library file (e.g.,
``armpl_lp64.so`` in ``/a/random/path/lib/`` and a corresponding header
file in ``/a/random/path/include/``), then what you have to do is craft
your own pkg-config file. It should have a matching name (so in this
example, ``armpl_lp64.pc``) and may be located anywhere. The
``PKG_CONFIG_PATH`` environment variable should be set to point to the
location of the ``.pc`` file. The contents of that file should be::

    libdir=/path/to/library-dir      # e.g., /a/random/path/lib
    includedir=/path/to/include-dir  # e.g., /a/random/path/include
    version=1.2.3                    # set to actual version
    extralib=-lm -lpthread -lgfortran   # if needed, the flags to link in dependencies
    Name: armpl_lp64
    Description: ArmPL - Arm Performance Libraries
    Version: ${version}
    Libs: -L${libdir} -larmpl_lp64      # linker flags
    Libs.private: ${extralib}
    Cflags: -I${includedir}

To check that this works as expected, you should be able to run::

    $ pkg-config --libs armpl_lp64
    -L/path/to/library-dir -larmpl_lp64
    $ pkg-config --cflags armpl_lp64
    -I/path/to/include-dir


Specifying the Fortran ABI to use
---------------------------------

Some linear algebra libraries are built with the ``g77`` ABI (also known as
"the ``f2c`` calling convention") and others with GFortran ABI, and these two
ABIs are incompatible. Therefore, if you build SciPy with ``gfortran`` and link
to a linear algebra library like MKL, which is built with a ``g77`` ABI,
there'll be an exception or a segfault. SciPy fixes this by using ABI wrappers
which rely on the CBLAS API for the few functions in the BLAS API that suffer
from this issue.

Note that SciPy needs to know at build time what needs to be done.
The build system will automatically check whether linear algebra
library is MKL or Accelerate (which both always use the ``g77`` ABI) and if so,
use the CBLAS API instead of the BLAS API. If autodetection fails or if the
user wants to override this autodetection mechanism for building against plain
``libblas``/``liblapack`` (this is what conda-forge does for example), use the
``-Duse-g77-abi=true`` build option. E.g.,::

    $ python -m build -C-Duse-g77-abi=true -Csetup-args=-Dblas=blas -Csetup-args=-Dlapack=lapack


64-bit integer (ILP64) BLAS/LAPACK
----------------------------------

Support for ILP64 BLAS and LAPACK, as of version 1.18.0, requires that also LP64
symbols are available in the same library. Hence, ILP64 support is only available
MKL and Accelerate. *Note: after the deprecated scipy.odr (the last Fortran
module) is removed, this restriction will be removed and OpenBLAS and other
libraries will also be supported.*


To build SciPy from source with ILP64 support, use ``-Duse-ilp64=true``. E.g.,
to use Accelerate on macOS::

    $ python -m build --wheel -Csetup-args=-Dblas=accelerate -C-Duse-ilp64=true -Dcython-blas-abi=lp64

And to use MKL on x86-64::

    # Note the "lp64" in the `-Dblas=` argument is not a mistake; this is
    # necessary as long as the cython_blas ABI is set to "lp64" (see next section)
    $ python -m build --wheel -Csetup-args=-Dblas=mkl-dynamic-lp64-seq -C-Duse-ilp64=true -Dcython-blas-abi=lp64

.. note::

   These instructions will change for scipy>=1.19.0, when the
   ``scipy.linalg.cython_blas`` ABI doesn't need to be forced to LP64 anymore.
   Please see the next section for more details.

Building with ``-Duse-ilp64=true`` by default flips the low-level Python and
Cython APIs in ``scipy.linalg`` to ILP64 as well. This may require downstream usage
to be made compatible.

From Python, low-level BLAS and LAPACK functions are available from ``scipy.linalg.blas``
and ``scipy.linalg.lapack`` namespaces::

    >>> from scipy.linalg.blas import dgemm  # this may be an LP64 or an ILP64 function

To choose the variant of a low-level routine, use ``get_blas_funcs`` and
``get_lapack_funcs`` functions::

    >>> from scipy.linalg.blas import get_blas_funcs
    >>> daxpy = get_blas_funcs('axpy', (np.ones(3),), ilp64='preferred')
    >>> daxpy.int_dtype
    dtype('int64')       # depends on the build option

High-level linear algebra functions (``norm``, ``solve`` and so on) should use this
mechanism under the hood, so they work transparently with either LP64 or ILP64
routines.


Cython BLAS/LAPACK integer ABI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Cython BLAS and LAPACK API (``scipy.linalg.cython_blas`` and
``scipy.linalg.cython_lapack``) uses the ``blas_int`` type for all integer
parameters. By default, ``blas_int`` follows the ``use-ilp64`` setting: it is
``int`` (32-bit) for LP64 builds and ``int64_t`` (64-bit) for ILP64 builds.

This can be overridden with the ``-Dcython-blas-abi`` build option, which
accepts three values:

- ``auto`` (default): follows the ``use-ilp64`` setting
- ``lp64``: always use 32-bit integers, even when ``use-ilp64=true``
- ``ilp64``: always use 64-bit integers

For example, to build with ILP64 BLAS/LAPACK but keep the Cython API at LP64
(for downstream compatibility)::

    $ spin build -S-Dblas=accelerate -S-Duse-ilp64=true -S-Dcython-blas-abi=lp64

This is most useful for Accelerate, which has good support for using both
LP64 and ILP64 at the same time.

MKL also supports this (it has ILP64 symbols with ``_64`` symbol suffix
available in ``mkl-dynamic-lp64-*`` shared libraries). MKL's support will work
fine for a Python-only stack, as long as you ensure to not link to different
shared libraries - avoid mixing ``mkl_rt.so``, ``mkl-dynamic-lp64-*.so`` and
``mkl-dynamic-ilp64-*.so``.

It may be useful when starting to introduce ILP64 usage to have SciPy itself
use ILP64 but keeping the Cython API at LP64, because downstream packages may
not yet support 64-bit integers in their Cython BLAS/LAPACK calls.

The build configuration can be checked at runtime via
``scipy.show_config()`` — look for the ``'blas cython ilp64'`` entry.

Some LAPACK functions use boolean variables (Fortran ``logical``). For booleans,
``cython_lapack`` uses the ``blas_bint`` type: when ``blas_int`` resolves to C ``int``,
``blas_bint`` resolves to Cython's ``bint`` type; when ``blas_int`` resolves to
``int64_t``, so does ``blas_bint``.


Using Cython BLAS/LAPACK ABI in downstream packages
"""""""""""""""""""""""""""""""""""""""""""""""""""

Downstream packages which consume ``cython_blas`` or ``cython_lapack`` interfaces
should ideally directly use the ``blas_int`` integer type at all call sites.

Some packages, however, might prefer to continue using ``int`` types, and manually
map between ``int`` and ``blas_int`` types (It is convenient to localize this mapping
in a single internal wrapper which converts ``int`` inputs to ``blas_int`` before
calling a LAPACK function, and then converts its ``blas_int`` outputs back to ``int``).
We stress that doing this limits the array sizes to ``< INT_MAX`` even if the
LAPACK itself is ILP64 enabled.

Consult `a worked example`_ which illustrates both approaches.

.. _a worked example: https://github.com/scipy/scipy/tree/main/scipy/linalg/tests/_cython_examples/ilp64_test_package


Work-in-progress
----------------

These options are planned to be fully supported, but currently not usable out
of the box:

- ILP64 (64-bit integer size) builds: making ``-Duse-ilp64=true`` work without
  using another option to avoid issues with the deprecated ``scipy.odr``.
  Those options are: ``-Dcython-blas-abi=lp64``, or
  ``-D_without-fortran=true``. The latter is a temporary option that doesn't
  use a Fortran compiler at all and doesn't install ``scipy.odr``.
- Automatically selecting from multiple possible BLAS and LAPACK options, with
  a user-provided order of precedence
