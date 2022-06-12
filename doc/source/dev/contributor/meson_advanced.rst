.. _meson-advanced:

===========================
Advanced Meson build topics
===========================

.. _blas-lapack-selection:

Select a different BLAS or LAPACK library
=========================================

BLAS and LAPACK library selection, other than the OpenBLAS default, is
implemented via Meson `build options
<https://mesonbuild.com/Build-options.html#build-options>`__. For example, to
select plain ``libblas`` and ``liblapack`` (this is typically Netlib
BLAS/LAPACK on Linux distros, and can be dynamically switched between
implementations on conda-forge), use::

    $ # for a development build
    $ meson setup build -Dblas=blas -Dlapack=lapack --prefix=$PWD/build-install
    $ python dev.py

    $ # to build and install a wheel
    $ python -m build -C-Dblas=blas -C-Dlapack=lapack
    $ pip install dist/scipy*.whl

Other options that should work (as long as they're installed with
``pkg-config`` support) include ``mkl`` and ``blis``.


Use different build types with Meson
====================================

Meson provides different build types while configuring the project. You can see
the available options for build types in
`the "core options" section of the Meson documentation <https://mesonbuild.com/Builtin-options.html#core-options>`__.

Assuming that you are building from scratch (do ``git clean -xdf`` if needed),
you can configure the build as following to use the ``debug`` build type::

    meson setup build --buildtype debug  --prefix=$PWD/build-install

Now, you can use the ``dev.py`` interface for further building, installing and
testing SciPy::

    python dev.py -s linalg

This will work because after initial configuration, Meson will remember the
config options.


Use GCC and Clang builds in parallel
====================================

It may be useful to have several builds of SciPy in the same repo, for example
to compare the differences between two compilers for diagnosing an issue. As
discussed, Meson is fully out-of-place, so different builds will not interfere
with each other. We assume in the rest of this section that GCC is the default.
For example, let us build using GCC and Clang.

1. Build with GCC::

    python dev.py --build-only

Using the above command, meson will build with the (default) GCC compilers in
the ``build`` directory.  It will then install SciPy into
``$PWD/build-install/lib/python3.x/site-packages/``.

2. Build with Clang::

    CC=clang CXX=clang++ FC=gfortran python dev.py --build-only --build-dir=build-clang

Using the above commands, Meson will build with the Clang, Clang++ and Gfortran
compilers in the ``build-clang`` directory.  It will then install SciPy into
``$PWD/build-clang-install/lib/python3.x/site-packages/``.

Meson will remember the compiler selection for the ``build-clang`` directory and
it cannot be changed, so each future invocation of
``python dev.py --build-dir=build-clang`` it will automatically use Clang.
Tip: use an alias to make this easier to use, e.g.,
``alias dev-clang="python dev.py --build-dir=build-clang"``.

A common reason to have two builds is to compare between them. For example,
to run the ``scipy.linalg`` tests for builds with both compilers, do::

    python dev.py -s linalg  # will run the tests for the GCC build
    python dev.py --build-dir build-clang -s linalg  # will run the tests for the Clang build

