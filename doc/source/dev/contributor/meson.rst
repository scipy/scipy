.. _meson:

=============================
How to build SciPy with Meson
=============================

.. warning::

   The change over from a `numpy.distutils` to a Meson based build happened
   recently, so there may still be rough edges. All platforms are supported,
   however niche configurations may not yet be tested.  These instructions
   should work reliably on Linux and macOS with a conda environment and
   OpenBLAS as the BLAS/LAPACK library.
   Building SciPy on Windows has always been a little tricky - it does work
   provided that all build dependencies are installed (most importantly,
   OpenBLAS and a compatible set of compilers: either Mingw-w64 or MSVC + Intel
   Fortran works).  If you have issues on Windows, please look at
   `this CI job <https://github.com/scipy/scipy/blob/main/.github/workflows/windows.yml>`_
   for details.


Quickstart from scratch
=======================

Clone the repo if you haven't done so yet, and initialize the git submodules::

  git clone git@github.com:scipy/scipy.git
  git submodule update --init

Create a conda development environment, build SciPy with Meson and run the test
suite::

  conda env create -f environment.yml
  conda activate scipy-dev
  python dev.py test


Full details and explanation
============================

To build SciPy, let's start with a clean repo::

  git clone git@github.com:scipy/scipy.git
  git submodule update --init

We will use conda here, because it's the easiest way to get a fully
reproducible environment. If you do not have a conda environment yet, the
recommended installer is
`Mambaforge <https://github.com/conda-forge/miniforge#mambaforge>`__
(``mamba`` is basically a much faster ``conda``).

To create a development environment::

  conda env create -f environment.yml  # `mamba` works too for this command
  conda activate scipy-dev

Support for Cython in Meson is very new, and we also need some recent bug
fixes and new features in Meson - hence we need a ``>=0.60.x`` release
(automatically installed via use of ``environment.yml`` above).

Meson uses a configure and a build stage. To configure it for putting the build
artifacts in ``build/`` and a local install under ``build-install/`` and then
build::

  meson setup build --prefix=$PWD/build-install
  ninja -C build

In the command above, ``-C`` is followed by the name of the build directory.
You can have multiple builds at the same time. Meson is fully out-of-place, so
those builds will not interfere with each other. You can for example have a GCC
build, a Clang build and a debug build in different directories.

To then install SciPy into the prefix (``build-install/`` here, but note that
that's just an arbitrary name we picked here)::

  meson install -C build

It will then install to ``build-install/lib/python3.9/site-packages/scipy``, which
is not on your Python path, so to add it do (*note, having to use ``PYTHONPATH``
is temporary, this will be changed once we merge support for building wheels*)::

  export PYTHONPATH=$PWD/build-install/lib/python3.9/site-packages/

Now we should be able to import ``scipy`` and run the tests. Remembering that
we need to move out of the root of the repo to ensure we pick up the package
and not the local ``scipy/`` source directory::

  cd doc
  python -c "from scipy import constants as s; s.test()"

The above runs the tests for a single module, ``constants``. Other ways of
running the tests should also work, for example::

  pytest --pyargs scipy

The full test suite should pass, without any build warnings on Linux (with GCC
9 at least) and a moderate amount on the other platforms.

.. _the-dev-py-interface:

The ``dev.py`` interface
========================

The above configure-build-install-test docs are useful to understand how the
Meson build works, and for working on build improvements.
If you want the "all-in-one" command for all of the above, run::

  python dev.py test

This interface has many options, allowing you to perform all regular
development-related tasks (building, running tests, building docs, running
benchmarks, etc.). Here we document a few of the most commonly used options;
run ``python dev.py --help`` or ``--help`` on each of the subcommands for more
details.

Use the following command to build and install SciPy::

  python dev.py build

To run the tests use::

  python dev.py test

To run the tests for a particular submodule(let's say ``optimize``), you can use::

  python dev.py test -s optimize


To learn more about Meson
=========================

It's worth pointing out that Meson has `very good documentation <https://mesonbuild.com/>`__;
it pays off to read it, and is often the best source of answers for "how to do X". Furthermore, an extensive pdf book on Meson can be obtained for free at https://nibblestew.blogspot.com/2021/12/this-year-receive-gift-of-free-meson.html

To learn more about the design principles Meson uses, the recent talks linked
from `mesonbuild.com/Videos <https://mesonbuild.com/Videos.html>`__ are also a
good resource.

For running the Linux Meson CI job locally, one can use the ``act`` tool, see
:ref:`using-act`.

.. _meson-faq:

Meson frequently asked questions
================================

**Q: What are the changes in dependencies when switching to Meson?**

There are a couple of new dependencies:

- ``meson``: the Meson build system, installable as a pure Python package from
  PyPI or conda-forge
- ``ninja``: the build tool invoked by Meson to do the actual building (e.g.
  invoking compilers). Installable also from PyPI (on all common platforms) or
  conda-forge.
- ``pkg-config``: the tool used for discovering dependencies (in particular
  BLAS/LAPACK). Available on conda-forge (and Homebrew, Chocolatey, and Linux
  package managers), but not packaged on PyPI.

In case your ``pkg-config`` is not on the ``PATH`` and you don't want to add
it, you can set an environment variable to let Meson find it. For example for
Homebrew:
``export PKG_CONFIG_PATH="/opt/homebrew/opt/openblas/lib/pkgconfig"``.

Note that we are also losing dependencies, namely ``numpy.distutils`` and
``setuptools``. Overall we are (a) switching build systems, and (b) adding
``pkg-config`` for more reliable dependency discovery than the hardcoded paths
that ``numpy.distutils`` used.

**Q: I currently use in-place builds, how is my workflow changing?**

Meson by design does not support in-place builds. This has advantages (e.g.,
one can use multiple parallel builds, caching becomes easier, etc.) - however
it does mean that one current workflow is no longer supported.

The recommended workflow is to use ``python dev.py``. This works exactly the
same way as ``python runtests.py`` worked before. What it does is rebuild if
needed, and then install SciPy to a private directory (default is
``build-install/`` in-tree) before running tests or other development tasks. This
way modifications to pure Python code get picked up.

If you use an IDE with, e.g., a "Run" button for scripts which were pointing to
an in-place build, and you would really like to continue using that same
workflow instead of ``python dev.py``, then you have a few options:

- After modifying pure Python code in the SciPy repo, install it on the command
  line with ``python dev.py build``, or with ``meson install -C build``
  before running your script.
- If your IDE supports it, customize what the "Run" button does before running
  the script, to do the install each time (this is expected to take 2-3 sec)
  For this to work, the install directory
  (``build-install/lib/python3.X/site-packages`` inside the repo by default)
  should be added to ``PYTHONPATH``. *Note that the Spyder IDE does not yet
  support this; its developers are looking at implementing support before the
  SciPy 1.9.0 release).*

**Q: I'm seeing a warning "Broken python installation detected. ..."**

Please ignore these warnings, they are innocuous. They indicate that the
install path is outside of a ``site-packages`` directory (which we prefer as
the default for ``python dev.py``).

Upgrading Meson to 0.62.0 will make the warning go away.

**Q: How do the current build/install commands change?**

*Old workflows (numpy.distutils based):*

1. ``python runtests.py``
2. ``python setup.py build_ext -i`` + ``export
   PYTHONPATH=/home/username/path/to/scipy/reporoot`` (and then edit pure
   Python code in SciPy and run it with ``python some_script.py``).
3. ``python setup.py develop`` - this is similar to (2), except in-place build
   is made permanently visible in env.
4. ``python setup.py bdist_wheel`` + ``pip install dist/scipy*.whl`` - build
   wheel in current env (i.e. uses installed numpy, etc.) and install it.
5. ``pip install .`` - build wheel in an isolated build env against deps in
   ``pyproject.toml`` and install it. *Note: be careful, this is usually not
   the correct command for development installs - typically you want to use (4)
   or* ``pip install . -v --no-build-isolation``.

*New workflows (Meson based):*

Note that currently (29 Dec 2021) only (1) is implemented. The rest is to be
added/documented in follow-up PRs over the next few weeks to months.

1. ``python dev.py``
2. *no direct equivalent for in-place builds (but see FAQ entry on in-place
   builds)*
3. *same as (2)*
4. ``python -m build --no-isolation`` + ``pip install dist/scipy*.whl`` - see
   `pypa/build <https://pypa-build.readthedocs.io/en/latest/>`_; it's also
   possible Meson will gain the capability to build wheels directly, but
   ``python -m build`` is going to become the standard way of doing this.
5. ``pip install .`` - this will work unchanged after switching the default in
   ``pyproject.toml`` to Meson.

