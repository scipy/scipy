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
    $ python dev.py build -C-Dblas=blas -C-Dlapack=lapack

    $ # to build and install a wheel
    $ python -m build -Csetup-args=-Dblas=blas -Csetup-args=-Dlapack=lapack
    $ pip install dist/scipy*.whl

    $ # Or, with pip>=23.1, this works too:
    $ python -m pip -Csetup-args=-Dblas=blas -Csetup-args=-Dlapack=lapack

Other options that should work (as long as they're installed with
``pkg-config`` or CMake support) include ``mkl``, ``atlas``, ``blis`` and
``accelerate``.

Note that both Accelerate and ``scipy-openblas`` have flags in ``dev.py``
that are easier to remember, since they're commonly used for development::

    $ python dev.py build --with-accelerate
    $ python dev.py build --with-scipy-openblas

The ``-Dlapack`` flag isn't needed for Accelerate, MKL or ``scipy-openblas``,
since we can be sure that BLAS and LAPACK are the same for those options.
E.g., to create a wheel with Accelerate (on macOS >=13.3 only), use::

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

Note that SciPy needs to know at build time, what needs to be done and
the build system will automatically check whether linear algebra
library is MKL or Accelerate (which both always use the ``g77`` ABI) and if so,
use the CBLAS API instead of the BLAS API. If autodetection fails or if the
user wants to override this autodetection mechanism for building against plain
``libblas``/``liblapack`` (this is what conda-forge does for example), use the
``-Duse-g77-abi=true`` build option. E.g.,::

    $ python -m build -C-Duse-g77-abi=true -Csetup-args=-Dblas=blas -Csetup-args=-Dlapack=lapack 


Work-in-progress
----------------

These options are planned to be fully supported, but currently not usable out
of the box:

- ILP64 (64-bit integer size) builds: large parts of SciPy support using ILP64
  BLAS/LAPACK. Note that support is still incomplete, so SciPy *also* requires
  LP64 (32-bit integer size) BLAS/LAPACK.
- Automatically selecting from multiple possible BLAS and LAPACK options, with
  a user-provided order of precedence

