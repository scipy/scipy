.. _meson-advanced:

Use different build types with Meson
====================================

Meson provides different build types while configuring the project. You can see
the available options for build types at
`meson documentation <https://mesonbuild.com/Builtin-options.html#core-options>`.

Assuming that you are building from scratch(do ``git clean -xdf`` if needed), you can
configure the build as following to use the ``debug`` build type::

    meson setup build --buildtype debug  --prefix=$PWD/build-install

Now, you can use the ``dev.py`` interface for further building, installing and testing SciPy::

    python dev.py -s linalg

This will work because after initial configuration, Meson will remember the config options.


Use GCC and Clang build in parallel
===================================

As discussed, Meson is fully out-of-place, so different builds will not interfere
with each other. We assume in the rest of this section that GCC is the default.
For example, let us build using GCC and Clang.

1. Build with GCC::

    python dev.py --build-only

Using the above command, meson will build with gcc compiler in ``build`` directory.
It will then install SciPy into ``$PWD/build-install/lib/python3.x/site-packages/``.

2. Build with Clang::

    CC=clang CXX=clang++ FC=gfortran python dev.py --build-only --build-dir=build-clang

Using the above commands, meson will build with clang compiler in ``build-clang`` directory.
It will then install SciPy into ``$PWD/build-clang-install/lib/python3.x/site-packages/``.
Meson will remember the compiler selection for the ``build-clang`` directory and
it cannot be changed, so each future invocation of
``python dev.py --build-dir=build-clang`` will automatically use Clang.

A common reason to have two builds is to compare between them. For example,
to run the ``scipy.linalg`` tests for builds with both compilers, do::

    python dev.py -s linalg                   # will run the tests for the GCC build
    python dev.py --build-dir build-clang -s linalg    # will run the tests for the Clang build

Note: python3.x represents the current python3 version in your machine.
