.. _meson:

=============================
How to build SciPy with Meson
=============================

*note: these instructions are for Linux and macOS. Windows is not officially
supported yet (if you want to use it on Windows anyway, please look at
`this CI job <https://github.com/rgommers/scipy/blob/meson/.github/workflows/windows.yml>`__
for details)! Other Unix-like OSes may work, but are untested (please open an
issue if you have tested and something seems broken)*


Quickstart from scratch
=======================

Clone the repo if you haven't done so yet, and initialize the git submodules::

  git clone git@github.com:scipy/scipy.git
  git submodule update --init

Create a conda development environment, build SciPy with Meson and run the test
suite::

  conda env create -f environment_meson.yml
  conda activate scipy-meson
  python dev.py


Full details and explanation
============================

To build SciPy, we need the SciPy ``master`` branch. Note that further work
on Meson integration is being done in the ``meson`` branch from ``@rgommers``'s
fork. We stay with SciPy master here::

  git clone git@github.com:scipy/scipy.git
  git submodule update --init

We will use conda here, because it's the easiest way to get a fully
reproducible environment. If you do not have a conda environment yet, the
recommended installer is
`Mambaforge <https://github.com/conda-forge/miniforge#mambaforge>`__
(``mamba`` is basically a much faster ``conda``).

To create a development environment::

  conda env create -f environment_meson.yml  # `mamba` works too for this command
  conda activate scipy-meson

Support for Cython in Meson is very new, and we also need some recent bug
fixes and new features in Meson - hence we need the ``0.60.0`` release
(automatically installed via use of ``environment_meson.yml`` above).

Meson uses a configure and a build stage. To configure it for putting the build
artifacts in ``build/`` and a local install under ``installdir/`` and then
build::

  meson setup build --prefix=$PWD/installdir
  ninja -C build

In the command above, ``-C`` is followed by the name of the build directory.
You can have multiple builds at the same time. Meson is fully out-of-place, so
those builds will not interfere with each other. You can for example have a GCC
build, a Clang build and a debug build in different directories.

To then install SciPy into the prefix (``installdir/`` here, but note that
that's just an arbitrary name we picked here)::

  meson install -C build

It will then install to ``installdir/lib/python3.9/site-packages/scipy``, which
is not on your Python path, so to add it do (*note, having to use ``PYTHONPATH``
is temporary, this will be changed once we merge support for building wheels*)::

  export PYTHONPATH=$PWD/installdir/lib/python3.9/site-packages/

Now we should be able to import ``scipy`` and run the tests. Remembering that
we need to move out of the root of the repo to ensure we pick up the package
and not the local ``scipy/`` source directory::

  cd doc
  python -c "from scipy import constants as s; s.test()"

The above runs the tests for a single module, ``constants``. Other ways of
running the tests should also work, for example::

  pytest --pyargs scipy

Current status (24 Dec '21) is that the full test suite passes on Linux, macOS
and Windows with OpenBLAS, without any build warnings on Linux (with GCC 9 at
least) and a moderate amount on the other platforms. There is CI (one job in
SciPy master, and more on ``@rgommers``'s fork) to keep it that way.
The current status is already good enough to work on both build related issues
(e.g. build warnings, debugging some C/C++ extension) and on general SciPy
development. It is already a much smoother/faster experience than
working with the default ``distutils``-based build one gets with
``python setup.py develop`` - especially when working on compiled code.


The ``dev.py`` interface
========================

The above configure-build-install-test docs are useful to understand how the
Meson build works, and for working on build improvements.
If you want the "all-in-one" command for all of the above, run::

  python dev.py

This interface has many options, allowing you to perform all regular
development-related tasks (building, running tests, building docs, running
benchmarks, etc.). Here we document a few of the most commonly used options;
run ``python dev.py --help`` for more details.

Use the following command to build and install SciPy::

  python dev.py --build-only

To run the tests use (``-n`` is short for ``--no-build``)::

  python dev.py -n

To run the tests for a particular submodule(let's say ``optimize``), you can use::

  python dev.py -n -s optimize


To learn more about Meson
=========================

It's worth pointing out that Meson has `very good documentation <https://mesonbuild.com/>`__;
it pays off to read it, and is often the best source of answers for "how to do X".

To learn more about the design principles Meson uses, the recent talks linked
from `mesonbuild.com/Videos <https://mesonbuild.com/Videos.html>`__ are a good
resource.

For running the Linux Meson CI job locally, one can use the ``act`` tool, see
:ref:`using-act`.
