Compiler selection and customizing a build
==========================================

Selecting a specific compiler
-----------------------------

Meson supports the standard environment variables ``CC``, ``CXX`` and ``FC`` to
select specific C, C++ and/or Fortran compilers. These environment variables are
documented in `the reference tables in the Meson docs
<https://mesonbuild.com/Reference-tables.html#compiler-and-linker-flag-environment-variables>`__.

Note that environment variables only get applied from a clean build, because
they affect the configure stage (i.e., ``meson setup``). An incremental rebuild
does not react to changes in environment variables - you have to run ``git
clean -xdf`` and do a full rebuild, or run ``meson setup --reconfigure``.


Adding a custom compiler or linker flag
---------------------------------------

Meson by design prefers builds being configured through command-line options
passed to ``meson setup``. It provides many built-in options:

- For enabling a debug build and the optimization level, see the next section
  on "build types",
- Enabling ``-Werror`` in a portable manner is done via ``-Dwerror=true``,
- Enabling warning levels is done via ``-Dwarning_level=<val>``, with ``<val>``
  one of ``{0, 1, 2, 3, everything}``,
- There are many other builtin options, from activating Visual Studio
  (``-Dvsenv=true``) and building with link time optimization (``-Db_lto``) to
  changing the default C++ language level (``-Dcpp_std='c++17'``) or linker
  flags (``-Dcpp_link_args='-Wl,-z,defs'``).

For a comprehensive overview of options, see `Meson's builtin options docs page
<https://mesonbuild.com/Builtin-options.html>`__.

Meson also supports the standard environment variables ``CFLAGS``,
``CXXFLAGS``, ``FFLAGS`` and ``LDFLAGS`` to inject extra flags - with the same
caveat as in the previous section about those environment variables being
picked up only for a clean build and not an incremental build.


Using different build types with Meson
--------------------------------------

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


Controlling build parallelism
-----------------------------

By default, ``ninja`` will launch ``2*n_cpu + 2``, with ``n_cpu`` the number of
physical CPU cores, parallel build jobs. This is fine in the vast majority of
cases, and results in close to optimal build times. In some cases, on machines
with a small amount of RAM relative to the number of CPU cores, this leads to a
job running out of memory. In case that happens, lower the number of jobs ``N``
such that you have at least 2 GB RAM per job. For example, to launch 6 jobs::

    python -m pip install . -Ccompile-args="-j6"

or::

    python dev.py build -j6


Use GCC and Clang builds in parallel
------------------------------------

It may be useful to have several builds of SciPy in the same repo, for example
to compare the differences between two compilers for diagnosing an issue. As
discussed, Meson is fully out-of-place, so different builds will not interfere
with each other. We assume in the rest of this section that GCC is the default.
For example, let us build using GCC and Clang.

1. Build with GCC::

    python dev.py build

Using the above command, meson will build with the (default) GCC compilers in
the ``build`` directory, and install to the ``build-install`` directory.

2. Build with Clang::

    CC=clang CXX=clang++ FC=gfortran python dev.py --build-dir=build-clang build

Using the above commands, Meson will build with the Clang, Clang++ and Gfortran
compilers in the ``build-clang`` directory, and then install SciPy into
``build-clang-install``.

Meson will remember the compiler selection for the ``build-clang`` directory and
it cannot be changed, so each future invocation of
``python dev.py --build-dir=build-clang <command>`` it will automatically use Clang.

Tip: use an alias to make this easier to use, e.g.,
``alias dev-clang="python dev.py --build-dir=build-clang"`` and then
``dev-clang build``.

A common reason to have two builds is to compare between them. For example,
to run the ``scipy.linalg`` tests for builds with both compilers, do::

    python dev.py -s linalg                          # run tests for the GCC build
    python dev.py --build-dir build-clang -s linalg  # run tests for the Clang build


